#include "instruction.h"
#include "prog.h"
#include "platform.h"
#include "tools/perfect_hash.h"
#include "tools/json_struct.h"
#include "tools/softmc_utils.h"
#include "tools/ProgressBar.hpp"

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <list>
#include <cassert>
#include <bitset>
#include <chrono>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
using namespace boost::program_options;
using namespace boost::filesystem;

#include <array>
#include <algorithm>
#include <numeric>
#include <regex>

using namespace std;

#define CASR 0
#define BASR 1
#define RASR 2

#define NUM_SOFTMC_REGS 16
#define FPGA_PERIOD 1.5015f // ns

#define RED_TXT "\033[31m"
#define GREEN_TXT "\033[32m"
#define YELLOW_TXT "\033[33m"
#define BLUE_TXT "\033[34m"
#define MAGENTA_TXT "\033[35m"
#define NORMAL_TXT "\033[0m"

/*** DRAM Organization Parameters - UPDATE here if the organization of your DRAM module differs ***/
int NUM_BANKS = 16; // this is the total number of banks in the chip
int NUM_BANK_GROUPS = 4;
int NUM_ROWS  = 32768;
int ROW_SIZE = 8192;
int NUM_COLS_PER_ROW = 128;
/******/

/*** DRAM Timing Parameters - UPDATE the timing parameters to match the timings of your module ***/
float DEFAULT_TRCD = 13.5f; // ns
float DEFAULT_TRAS = 35.0f; // ns
float DEFAULT_TRP = 13.5f; // ns
float DEFAULT_TWR = 15.0f; // ns
float DEFAULT_TRFC = 260.0f; // ns
float DEFAULT_TRRDS = 5.3f; // ns (ACT-ACT to different bank groups)
float DEFAULT_TRRDL = 6.4f; // ns (ACT-ACT to same bank group)
/******/


int trcd_cycles = (int) ceil(DEFAULT_TRCD/FPGA_PERIOD);
int tras_cycles = (int) ceil(DEFAULT_TRAS/FPGA_PERIOD);
int trp_cycles = (int) ceil(DEFAULT_TRP/FPGA_PERIOD);
int twr_cycles = (int) ceil(DEFAULT_TWR/FPGA_PERIOD);
int trfc_cycles = (int) ceil(DEFAULT_TRFC/FPGA_PERIOD);
int trrds_cycles = (int) ceil(DEFAULT_TRRDS/FPGA_PERIOD);
int trrdl_cycles = (int) ceil(DEFAULT_TRRDL/FPGA_PERIOD);

// Retention Profiler Parameters
uint RETPROF_NUM_ITS = 100; // When a candidate row group is found, the profiler repeats the retention time test on the on the row num_test_iterations number of times to make sure the row is reliably weak
float RETPROF_RETTIME_STEP = 1.0f; // defines by how much to increase the target retention time if sufficient row groups not found in the previous iteration
float RETPROF_RETTIME_MULT_H = 1.2f;

vector<uint32_t> reserved_regs{CASR, BASR, RASR};

typedef struct RowData {
    bitset<512> input_data_pattern;
    uint pattern_id;
    string label;
} RowData;

typedef struct WeakRow {
    uint row_id;
    std::vector<uint> bitflip_locs;
    WeakRow(uint _row_id, vector<uint> _bitflip_locs) : 
            row_id(_row_id), bitflip_locs(_bitflip_locs) {}
} WeakRow;

JS_OBJECT_EXTERNAL(WeakRow,
                JS_MEMBER(row_id),
                JS_MEMBER(bitflip_locs));
                // JS_MEMBER(data_pattern_type));

typedef struct WeakRowSet {
    std::vector<WeakRow> row_group;
    uint bank_id;
    uint ret_ms;
    uint data_pattern_type;
    uint rowdata_ind;
    WeakRowSet(std::vector<WeakRow> _weak_rows, uint _bank_id, uint _ret_ms, uint _data_pattern_type, uint _rowdata_ind) : 
        row_group(_weak_rows), bank_id(_bank_id), ret_ms(_ret_ms), data_pattern_type(_data_pattern_type), rowdata_ind(_rowdata_ind)  {}
    bool contains_rows_in(const WeakRowSet& other) {
        for(auto& wr : row_group) {
            auto it = std::find_if(other.row_group.begin(), other.row_group.end(), [&](const WeakRow& w) {return w.row_id == wr.row_id;});

            if(it != other.row_group.end())
                return true;
        }

        return false;
    }
    std::string rows_as_str() {
        std::string ret = "(";

        for(auto& wr : row_group) {
            ret = ret + to_string(wr.row_id) + ", ";
        }

        // remove ', ' at the end
        if(ret.size() > 1) {
            ret = ret.substr(0, ret.size() - 2);
        }

        ret = ret + ")";

        return ret;
    }
} WeakRowSet;

JS_OBJECT_EXTERNAL(WeakRowSet,
                JS_MEMBER(row_group),
                JS_MEMBER(bank_id),
                JS_MEMBER(ret_ms),
                JS_MEMBER(data_pattern_type));


// returns a vector of bit positions that experienced bitflips
void collect_bitflips(vector<uint>& bitflips, const char* read_data, const RowData& rh_row) {

    bitset<512> read_data_bitset;

    uint32_t* iread_data = (uint32_t*) read_data;

    // check for bitflips in each cache line
    for(int cl = 0; cl < ROW_SIZE/64; cl++) {

        read_data_bitset.reset();
        for(int i = 0; i < 512/32; i++) {
            bitset<512> tmp_bitset = iread_data[cl*(512/32) + i];

            read_data_bitset |= (tmp_bitset << i*32);
        }

        // compare and print errors
        bitset<512> error_mask = read_data_bitset ^ rh_row.input_data_pattern;

        if(error_mask.any()) {
            // there is at least one bitflip in this cache line

            for(uint i = 0; i < error_mask.size(); i++){
                if(error_mask.test(i)){
                    bitflips.push_back(cl*CACHE_LINE_BITS + i);
                }
            }
        }
    }
}

bool in_same_bg(int bank1, int bank2);

int get_total_ACT_cycs(const vector<int>& target_banks) {
    int total_cycs = 4*target_banks.size(); // for SMC_LI that sets the bank id register

    for (int ind = 0; ind < (target_banks.size() - 1); ind++) {
        if(in_same_bg(target_banks[ind], target_banks[ind + 1]))
            total_cycs += 1 + trrdl_cycles;
        else
            total_cycs += 1 + trrds_cycles;
    }

    return total_cycs;
}

// update this function when changing the hammering part of the SoftMC
// program
int iters_per_round(int ref_interval_ns, int num_aggr_rows, const vector<int>& target_banks) {

    int iters = 0;

    int total_hammer_cycles = (target_banks.size() + get_total_ACT_cycs(target_banks) + tras_cycles)*num_aggr_rows + trp_cycles*(num_aggr_rows - 1) + 1 + 24*2;
    int total_hammer_ns = total_hammer_cycles * FPGA_PERIOD;
    iters = floor(ref_interval_ns/total_hammer_ns);

    cout << "RowHammer iterations to perform within " << ref_interval_ns << "ns refresh interval: " << iters << endl;

    return iters;
}

template<typename T>
bool vec_contains(const vector<T>& vec, T val) {

    return (find(vec.cbegin(), vec.cend(), val) != vec.cend());

}

// figures out whether two banks are in the same bank group
bool in_same_bg(int bank1, int bank2) {
    int bg1 = bank1/NUM_BANK_GROUPS;
    int bg2 = bank2/NUM_BANK_GROUPS;

    return bg1 == bg2;
}

void activateBanks(Program& program, const vector<int>& target_banks, int BANK_ADDR_REG, int ROW_ADDR_REG,
        bool is_fake_hammering = false) {

    // activate rows in the target banks
    int remaining_cycs = 0;
    int total_act_cycles = 0;
    for (int ind = 0; ind < target_banks.size(); ind++) {
      int bank_id = target_banks[ind];
      program.add_inst(SMC_LI(bank_id, BANK_ADDR_REG));
      total_act_cycles += 4;

      int cur_rrd = 0;
      if ((ind + 1) < target_banks.size()){
          if(in_same_bg(bank_id, target_banks[ind + 1]))
              cur_rrd = trrdl_cycles;
          else
              cur_rrd = trrds_cycles;
      }

      cur_rrd = min(0, cur_rrd-4); // -4 because we use SMC_LI to update the BANK_ADDR_REG
    
      remaining_cycs = add_op_with_delay(program, is_fake_hammering ? SMC_NOP() : SMC_ACT(BANK_ADDR_REG, 0, ROW_ADDR_REG, 0), remaining_cycs, cur_rrd);
      total_act_cycles += 4;
    }

    if(trcd_cycles > total_act_cycles)
      remaining_cycs = trcd_cycles - total_act_cycles;

    if(remaining_cycs > 0)
      add_op_with_delay(program, SMC_NOP(), remaining_cycs, 0);
}

void writeToDRAM(Program& program, const uint target_bank, const uint start_row, 
        const uint row_batch_size, const vector<RowData>& rows_data) {

    const int REG_TMP_WRDATA = 15;
    const int REG_BANK_ADDR = 12;
    const int REG_ROW_ADDR = 13;
    const int REG_COL_ADDR = 14;
    const int REG_NUM_COLS = 11;

    const int REG_BATCH_IT = 6;
    const int REG_BATCH_SIZE = 5;


    bitset<512> bitset_int_mask(0xFFFFFFFF);
    int remaining_cycs = 0;

    // ===== BEGIN SoftMC Program =====
  
    program.add_inst(SMC_LI(start_row, REG_ROW_ADDR));
    program.add_inst(SMC_LI(target_bank, REG_BANK_ADDR));

    add_op_with_delay(program, SMC_PRE(REG_BANK_ADDR, 0, 1), 0, 0); // precharge all banks
    
    program.add_inst(SMC_LI(NUM_COLS_PER_ROW*8, REG_NUM_COLS));

    program.add_inst(SMC_LI(8, CASR)); // Load 8 into CASR since each READ reads 8 columns
    program.add_inst(SMC_LI(1, BASR)); // Load 1 into BASR
    program.add_inst(SMC_LI(1, RASR)); // Load 1 into RASR


    /* ==== Initialize data of rows in the batch ==== */
    program.add_inst(SMC_LI(0, REG_BATCH_IT));
    program.add_inst(SMC_LI(row_batch_size, REG_BATCH_SIZE));

    assert(row_batch_size % rows_data.size() == 0 && "Data patterns to initialize consecutive rows with must be multiple of the batch of row to initialize at once.");

    program.add_label("INIT_BATCH");

        int row_it = 0;
        for(auto& row_data : rows_data) {
            // set up the input data in the wide register
    	    for (int pos = 0; pos < 16; pos++) {
              program.add_inst(SMC_LI((((row_data.input_data_pattern >> 32*pos) & bitset_int_mask).to_ulong() & 0xFFFFFFFF), REG_TMP_WRDATA));
      	      program.add_inst(SMC_LDWD(REG_TMP_WRDATA, pos));
      	    }

            remaining_cycs -= (16*4*2 + 4);
            assert(remaining_cycs <= 0 && "I should add some delay here");

            // activate the next row and increment the row address register
            add_op_with_delay(program, SMC_ACT(REG_BANK_ADDR, 0, REG_ROW_ADDR, 1), 0, trcd_cycles - 1);
            
            // write data to the row and precharge
      	    program.add_inst(SMC_LI(0, REG_COL_ADDR));

            string new_lbl = "INIT_ROW" + to_string(row_it++);
      	    program.add_label(new_lbl);
            add_op_with_delay(program, SMC_WRITE(REG_BANK_ADDR, 0, REG_COL_ADDR, 1, 0, 0), 0, 0);
            remaining_cycs = 0;
      	    program.add_branch(program.BR_TYPE::BL, REG_COL_ADDR, REG_NUM_COLS, new_lbl);

            // Wait for t(write-precharge)
            // & precharge the open bank
            remaining_cycs = add_op_with_delay(program, SMC_PRE(REG_BANK_ADDR, 0, 0), 0, trp_cycles);
        }

    program.add_inst(SMC_ADDI(REG_BATCH_IT, rows_data.size(), REG_BATCH_IT));
    program.add_branch(program.BR_TYPE::BL, REG_BATCH_IT, REG_BATCH_SIZE, "INIT_BATCH");

    program.add_inst(SMC_END());
}


void writeToDRAM(Program& program, SoftMCRegAllocator& reg_alloc, const uint target_bank, const WeakRowSet& wrs, const vector<RowData>& rows_data) {

    SMC_REG REG_TMP_WRDATA = reg_alloc.allocate_SMC_REG();
    SMC_REG REG_BANK_ADDR = reg_alloc.allocate_SMC_REG();
    SMC_REG REG_ROW_ADDR = reg_alloc.allocate_SMC_REG();
    SMC_REG REG_COL_ADDR = reg_alloc.allocate_SMC_REG();
    SMC_REG REG_NUM_COLS = reg_alloc.allocate_SMC_REG();

    bitset<512> bitset_int_mask(0xFFFFFFFF);

    // ===== BEGIN SoftMC Program =====
  
    program.add_inst(SMC_LI(target_bank, REG_BANK_ADDR));
    add_op_with_delay(program, SMC_PRE(REG_BANK_ADDR, 0, 1), 0, 0); // precharge all banks
    program.add_inst(SMC_LI(NUM_COLS_PER_ROW*8, REG_NUM_COLS));

    program.add_inst(SMC_LI(8, CASR)); // Load 8 into CASR since each READ reads 8 columns
    program.add_inst(SMC_LI(1, BASR)); // Load 1 into BASR
    program.add_inst(SMC_LI(1, RASR)); // Load 1 into RASR

    /* ==== Initialize data of rows in the row group ==== */

    // initialize the entire range that corresponds to the psysical row ids according to the row_group_pattern
    PhysicalRowID first_phys_row_id = to_physical_row_id(wrs.row_group.front().row_id);
    PhysicalRowID last_phys_row_id = to_physical_row_id(wrs.row_group.back().row_id);
    assert(last_phys_row_id >= first_phys_row_id);

    std::vector<LogicalRowID> rows_to_init;
    rows_to_init.reserve(last_phys_row_id - first_phys_row_id + 1);

    for(uint i = first_phys_row_id; i <= last_phys_row_id; i++) {
        rows_to_init.push_back(to_logical_row_id(i));
    }


    // set up the input data in the wide register
    for (int pos = 0; pos < 16; pos++) {
        program.add_inst(SMC_LI((((rows_data[wrs.rowdata_ind].input_data_pattern >> 32*pos) & bitset_int_mask).to_ulong() & 0xFFFFFFFF), REG_TMP_WRDATA));
        program.add_inst(SMC_LDWD(REG_TMP_WRDATA, pos));
    }

    for(uint row_id : rows_to_init){
        program.add_inst(SMC_LI(row_id, REG_ROW_ADDR));
        // activate the next row and increment the row address register
        add_op_with_delay(program, SMC_ACT(REG_BANK_ADDR, 0, REG_ROW_ADDR, 1), 0, trcd_cycles - 1);
        
        // write data to the row and precharge
        program.add_inst(SMC_LI(0, REG_COL_ADDR));

        string new_lbl = createSMCLabel("INIT_ROW_DATA");
        program.add_label(new_lbl);
            add_op_with_delay(program, SMC_WRITE(REG_BANK_ADDR, 0, REG_COL_ADDR, 1, 0, 0), 0, 0);
            add_op_with_delay(program, SMC_WRITE(REG_BANK_ADDR, 0, REG_COL_ADDR, 1, 0, 0), 0, 0);
            add_op_with_delay(program, SMC_WRITE(REG_BANK_ADDR, 0, REG_COL_ADDR, 1, 0, 0), 0, 0);
            add_op_with_delay(program, SMC_WRITE(REG_BANK_ADDR, 0, REG_COL_ADDR, 1, 0, 0), 0, 0);
        program.add_branch(program.BR_TYPE::BL, REG_COL_ADDR, REG_NUM_COLS, new_lbl);

        // Wait for t(write-precharge)
        // & precharge the open bank
        add_op_with_delay(program, SMC_PRE(REG_BANK_ADDR, 0, 0), 0, trp_cycles - 1);
    }
    

    program.add_inst(SMC_END());

    reg_alloc.free_SMC_REG(REG_TMP_WRDATA);
    reg_alloc.free_SMC_REG(REG_BANK_ADDR);
    reg_alloc.free_SMC_REG(REG_ROW_ADDR);
    reg_alloc.free_SMC_REG(REG_COL_ADDR);
    reg_alloc.free_SMC_REG(REG_NUM_COLS);
}

void readFromDRAM(Program& program, const uint target_bank, const uint start_row, const uint row_batch_size) {

    const int REG_BANK_ADDR = 12;
    const int REG_ROW_ADDR = 13;
    const int REG_COL_ADDR = 14;
    const int REG_NUM_COLS = 11;

    const int REG_BATCH_IT = 6;
    const int REG_BATCH_SIZE = 5;

    int remaining_cycs = 0;

    // ===== BEGIN SoftMC Program =====
  
    program.add_inst(SMC_LI(start_row, REG_ROW_ADDR));
    program.add_inst(SMC_LI(target_bank, REG_BANK_ADDR));

    add_op_with_delay(program, SMC_PRE(REG_BANK_ADDR, 0, 1), 0, 0); // precharge all banks
    
    //program.add_inst(SMC_LI(NUM_ROWS, REG_NUM_ROWS));
    program.add_inst(SMC_LI(NUM_COLS_PER_ROW*8, REG_NUM_COLS));

    program.add_inst(SMC_LI(8, CASR)); // Load 8 into CASR since each READ reads 8 columns
    program.add_inst(SMC_LI(1, BASR)); // Load 1 into BASR
    program.add_inst(SMC_LI(1, RASR)); // Load 1 into RASR

    /* ==== Read the data of rows in the batch ==== */
    program.add_inst(SMC_LI(0, REG_BATCH_IT));
    program.add_inst(SMC_LI(row_batch_size, REG_BATCH_SIZE));

    program.add_label("READ_BATCH");

    // activate the next row and increment the row address register
    add_op_with_delay(program, SMC_ACT(REG_BANK_ADDR, 0, REG_ROW_ADDR, 1), 0, trcd_cycles - 1);
    
    // issue read cmds to read out the entire row and precharge
    program.add_inst(SMC_LI(0, REG_COL_ADDR));

    string new_lbl = "READ_ROW";
    program.add_label(new_lbl);
    add_op_with_delay(program, SMC_READ(REG_BANK_ADDR, 0, REG_COL_ADDR, 1, 0, 0), remaining_cycs, 4);
    remaining_cycs = 0;
    program.add_branch(program.BR_TYPE::BL, REG_COL_ADDR, REG_NUM_COLS, new_lbl);

    // Wait for t(write-precharge)
    // & precharge the open bank
    remaining_cycs = add_op_with_delay(program, SMC_PRE(REG_BANK_ADDR, 0, 0), 0, trp_cycles);

    program.add_inst(SMC_ADDI(REG_BATCH_IT, 1, REG_BATCH_IT));
    program.add_branch(program.BR_TYPE::BL, REG_BATCH_IT, REG_BATCH_SIZE, "READ_BATCH");

    program.add_inst(SMC_END());
}

void readFromDRAM(Program& program, const uint target_bank, const WeakRowSet& wrs) {

    const int REG_BANK_ADDR = 12;
    const int REG_ROW_ADDR = 13;
    const int REG_COL_ADDR = 14;
    const int REG_NUM_COLS = 11;

    int remaining_cycs = 0;

    // ===== BEGIN SoftMC Program =====
  
    program.add_inst(SMC_LI(target_bank, REG_BANK_ADDR));

    add_op_with_delay(program, SMC_PRE(REG_BANK_ADDR, 0, 1), 0, 0); // precharge all banks
    
    //program.add_inst(SMC_LI(NUM_ROWS, REG_NUM_ROWS));
    program.add_inst(SMC_LI(NUM_COLS_PER_ROW*8, REG_NUM_COLS));

    program.add_inst(SMC_LI(8, CASR)); // Load 8 into CASR since each READ reads 8 columns
    program.add_inst(SMC_LI(1, BASR)); // Load 1 into BASR
    program.add_inst(SMC_LI(1, RASR)); // Load 1 into RASR

    /* ==== Read the data of rows in the row group ==== */

    for (auto& wr : wrs.row_group) {
        program.add_inst(SMC_LI(wr.row_id, REG_ROW_ADDR));
    
        // activate the next row and increment the row address register
        add_op_with_delay(program, SMC_ACT(REG_BANK_ADDR, 0, REG_ROW_ADDR, 0), 0, trcd_cycles - 1);
        
        // issue read cmds to read out the entire row and precharge
        program.add_inst(SMC_LI(0, REG_COL_ADDR));

        string new_lbl = createSMCLabel("READ_ROW");
        program.add_label(new_lbl);
        add_op_with_delay(program, SMC_READ(REG_BANK_ADDR, 0, REG_COL_ADDR, 1, 0, 0), remaining_cycs, 4);
        remaining_cycs = 0;
        program.add_branch(program.BR_TYPE::BL, REG_COL_ADDR, REG_NUM_COLS, new_lbl);

        // Wait for t(write-precharge)
        // & precharge the open bank
        remaining_cycs = add_op_with_delay(program, SMC_PRE(REG_BANK_ADDR, 0, 0), 0, trp_cycles);
    }

    program.add_inst(SMC_END());
}

uint determineRowBatchSize(const uint retention_ms, const uint num_data_patterns) {

    uint pcie_cycles = ceil(5000/FPGA_PERIOD); // assuming 5us pcie transfer latency
    uint setup_cycles = 36;
    uint pattern_loop_cycles = 64/*write reg init*/ + 1 /*ACT*/ + trcd_cycles + 4 + 
        (4 + 24)*NUM_COLS_PER_ROW /*row write*/ + 1 + trp_cycles;


    // cycles(retention_ms) = pcie_cycles + setup_cycles +
    // (X/NUM_PATTERNS)*(NUM_PATTERNS*pattern_loop_cycles + 28)
    //
    // X: batch size
    //
    // X = ((retention_cycles - pcie_cycles -
    // setup_cycles)/(NUM_PATTERNS*pattern_loop_cycles + 28))*NUM_PATTERNS
    ulong retention_cycles = floor((retention_ms*1000000)/FPGA_PERIOD);

    uint batch_size = ((retention_cycles - pcie_cycles - setup_cycles)/(num_data_patterns*pattern_loop_cycles + 28))*num_data_patterns;

    // cout << "Calculated initial batch size as " << batch_size << " for " << retention_ms << " ms" << endl;
    // cout << "Rounding batch_size to the previous power-of-two number" << endl;

    assert(NUM_ROWS % num_data_patterns == 0 && "Number of specified data patterns must be a divisor of NUM_ROWS, i.e., power of two");

    // rounding
    batch_size = min(1 << (uint)(log2(batch_size)), NUM_ROWS);

    // cout << "The final batch_size: " << batch_size << endl;

    return batch_size;
}

void checkForLeftoverPCIeData(SoftMCPlatform& platform) {
    // checking if there is more data to receive
    uint additional_bytes = 0;

    char* buf = new char[4];

    while (additional_bytes += platform.receiveData(buf, 4)) { // try to receive 4 bytes
        cout << "Received total of " << additional_bytes << " additional bytes" << endl;
        cout << "Data: " << *(int*)buf << endl;
    }

    delete[] buf;
}

// pairs of row_id and number of bitflips
static std::vector<std::pair<int, std::vector<uint>>> bitflip_history;
static std::vector<uint> locs_to_check;
void init_row_pattern_fitter(const std::string row_group_pattern) {
    // initialize with row id -1 and empty bitflip vectors
    bitflip_history = std::vector<std::pair<int, vector<uint>>>(row_group_pattern.size(), std::pair<int, std::vector<uint>>(-1, std::vector<uint>()));

    for(uint i = 0; i < row_group_pattern.size(); i++) {
        char c = row_group_pattern[i];
        if (c == 'R')
            locs_to_check.push_back(i);
    }
    
}

void clear_bitflip_history() {
    std::fill(bitflip_history.begin(), bitflip_history.end(), std::pair<int, std::vector<uint>>(-1, std::vector<uint>()));
}

bool fits_into_row_pattern(const vector<uint>& bitflips, const uint row_id) {
    // check if receiving row_id in order
    if(bitflip_history.back().first != -1 && bitflip_history.back().first != (row_id - 1)){
        std::cerr << RED_TXT << "ERROR: Did not profile rows in order. Got row id " << row_id << " after row id " << bitflip_history.back().first << NORMAL_TXT << std::endl;
        exit(-1);
    }

    // remove the oldest entry in bitflip_history
    bitflip_history.erase(bitflip_history.begin());

    // insert a new entry
    bitflip_history.push_back(std::pair<uint, vector<uint>>(row_id, bitflips));

    // check if the current bitflip_history matches the pattern
    for(auto loc : locs_to_check) {
        if(bitflip_history[loc].second.size() == 0)
            return false;
    }
    return true;
}

void build_WeakRowSet(vector<WeakRowSet>& wrs, const std::string& row_group_pattern, const vector<RowData>& rows_data, 
                        const uint target_bank, const uint batch_ind, const uint batch_first_row, const uint retention_ms) {

    vector<WeakRow> row_group;
    uint vec_size = std::count(row_group_pattern.begin(), row_group_pattern.end(), 'R');
    row_group.reserve(vec_size);
    
    for(auto loc : locs_to_check) {
        auto& bfh = bitflip_history[loc];
        row_group.emplace_back(to_logical_row_id(bfh.first), bfh.second);
    }

    assert(rows_data.size() == 1); // remove this if you are trying to enable support for different input data patterns for different rows
    wrs.emplace_back(row_group, target_bank, retention_ms, rows_data[0].pattern_id, 0);

    // to prevent a row being part of multiple WeakRowSets
    clear_bitflip_history();
}

void test_retention(SoftMCPlatform& platform, const uint retention_ms, const uint target_bank, const uint first_row_id, 
                    const uint row_batch_size, const vector<RowData>& rows_data, const std::string& row_group_pattern, char* buf, vector<WeakRowSet>& row_group) {
    
    Program writeProg;
    writeToDRAM(writeProg, target_bank, first_row_id, row_batch_size, rows_data);

    // execute the program
    auto t_start_issue_prog = chrono::high_resolution_clock::now();
    platform.execute(writeProg);
    auto t_end_issue_prog = chrono::high_resolution_clock::now();

    chrono::duration<double, milli> prog_issue_duration(t_end_issue_prog - t_start_issue_prog);

    // cout << "Issuing the DRAM write program took: " << prog_issue_duration.count() << " ms." << endl; 
    waitMS(retention_ms - prog_issue_duration.count());

    // at this point we expect writing data to DRAM to be finished
    auto t_end_ret_wait = chrono::high_resolution_clock::now();

    chrono::duration<double, milli> write_ret_duration(t_end_ret_wait - t_start_issue_prog);
    // cout << "Writing to DRAM + waiting for the target " << retention_ms << " ms retention time took: " << write_ret_duration.count() << endl;

    
    // READ DATA BACK AND CHECK ERRORS 
    auto t_prog_started = chrono::high_resolution_clock::now();
    Program readProg;
    readFromDRAM(readProg, target_bank, first_row_id, row_batch_size);
    platform.execute(readProg);
    //checkForLeftoverPCIeData(platform);
    platform.receiveData(buf, ROW_SIZE*row_batch_size); // reading all RH_NUM_ROWS at once
    //t_two_rows_recvd = chrono::high_resolution_clock::now();
    //elapsed = t_two_rows_recvd - t_prog_started;
    //cout << "Time for reading two rows: " << elapsed.count()*1000 << "ms" << endl;

    bool check_time = false;
    if(check_time) {
        check_time = false;
        auto t_two_rows_recvd = chrono::high_resolution_clock::now();

        auto elapsed = t_two_rows_recvd - t_prog_started;
        cout << "Time interval for reading back " << row_batch_size << "rows: " << elapsed.count()*1000.0f << "ms" << endl;
    }

    //t_prog_started = chrono::high_resolution_clock::now();
    
    // go over physical row IDs in order
    for (int i = 0; i < row_batch_size; i++) {
        PhysicalRowID phys_row_id = first_row_id + i;
        LogicalRowID log_row_id = to_logical_row_id(phys_row_id);

        // std::cout << "first_row_id: " << first_row_id << std::endl;
        // std::cout << "row_batch_size: " << row_batch_size << std::endl;
        // std::cout << "log_row_id: " << log_row_id << std::endl;
        assert(log_row_id < (first_row_id + row_batch_size) && log_row_id >= first_row_id &&
                "ERROR: The used Logical to Physical row address mapping results in logical address out of bounds of the row_batch size. Consider revising the code.");

        vector<uint> bitflips; 
        collect_bitflips(bitflips, buf + (log_row_id - first_row_id)*ROW_SIZE, rows_data[(log_row_id - first_row_id) % rows_data.size()]);

        if (fits_into_row_pattern(bitflips, phys_row_id)) {
            build_WeakRowSet(row_group, row_group_pattern, rows_data, target_bank, i, first_row_id, retention_ms);
        }
    }

    //t_two_rows_recvd = chrono::high_resolution_clock::now();

    //elapsed = t_two_rows_recvd - t_prog_started;
    //cout << "Error checking time: " << elapsed.count()*1000.0f << "ms" << endl;

    //cout << "Finished testing rows " << start_row << "-" << start_row + row_batch_size - 1 << endl;
}

// return true if the same bit locations in WeakRowSet wrs experience bitflips
bool check_retention_failute_repeatability(SoftMCPlatform& platform, const uint retention_ms, const uint target_bank, WeakRowSet& wrs, 
                    const vector<RowData>& rows_data, char* buf, bool filter_out_failures = false) {
    
    Program writeProg;
    SoftMCRegAllocator reg_alloc(NUM_SOFTMC_REGS, reserved_regs);
    writeToDRAM(writeProg, reg_alloc, target_bank, wrs, rows_data);

    // execute the program
    auto t_start_issue_prog = chrono::high_resolution_clock::now();
    platform.execute(writeProg);
    auto t_end_issue_prog = chrono::high_resolution_clock::now();

    chrono::duration<double, milli> prog_issue_duration(t_end_issue_prog - t_start_issue_prog);

    // cout << "Issuing the DRAM write program took: " << prog_issue_duration.count() << " ms." << endl; 
    waitMS(retention_ms - prog_issue_duration.count());

    // at this point we expect writing data to DRAM to be finished
    auto t_end_ret_wait = chrono::high_resolution_clock::now();

    chrono::duration<double, milli> write_ret_duration(t_end_ret_wait - t_start_issue_prog);
    // cout << "Writing to DRAM + waiting for the target " << retention_ms << " ms retention time took: " << write_ret_duration.count() << endl;

    
    // READ DATA BACK AND CHECK ERRORS 
    auto t_prog_started = chrono::high_resolution_clock::now();
    Program readProg;
    readFromDRAM(readProg, target_bank, wrs);
    platform.execute(readProg);
    //checkForLeftoverPCIeData(platform);
    platform.receiveData(buf, ROW_SIZE*wrs.row_group.size()); // reading all RH_NUM_ROWS at once
    //t_two_rows_recvd = chrono::high_resolution_clock::now();
    //elapsed = t_two_rows_recvd - t_prog_started;
    //cout << "Time for reading two rows: " << elapsed.count()*1000 << "ms" << endl;

    bool check_time = false;
    if(check_time) {
        check_time = false;
        auto t_two_rows_recvd = chrono::high_resolution_clock::now();

        auto elapsed = t_two_rows_recvd - t_prog_started;
        cout << "Time interval for reading back " << wrs.row_group.size() << "rows: " << elapsed.count()*1000.0f << "ms" << endl;
    }

    //t_prog_started = chrono::high_resolution_clock::now();
    
    for (int i = 0; i < wrs.row_group.size(); i++) {
        vector<uint> bitflips;
        collect_bitflips(bitflips, buf + i*ROW_SIZE, rows_data[wrs.rowdata_ind]);

        // filter out bit locations that flipped
        if(filter_out_failures) {
            for (uint bf_loc : bitflips) {
                auto it = std::find(wrs.row_group[i].bitflip_locs.begin(), wrs.row_group[i].bitflip_locs.end(), bf_loc);
                if(it != wrs.row_group[i].bitflip_locs.end())
                    wrs.row_group[i].bitflip_locs.erase(it);
            }
        } else { // filter out bit locations that did not flip
            for (auto it = wrs.row_group[i].bitflip_locs.begin(); it != wrs.row_group[i].bitflip_locs.end(); it++) {
                auto it_find = std::find(bitflips.begin(), bitflips.end(), *it);
                if(it_find == bitflips.end())
                    wrs.row_group[i].bitflip_locs.erase(it--);
            }
        }

        // return false if all bitflip locations in a weak row were filtered out
        if(wrs.row_group[i].bitflip_locs.size() == 0)
            return false;
    }

    return true;

    //t_two_rows_recvd = chrono::high_resolution_clock::now();

    //elapsed = t_two_rows_recvd - t_prog_started;
    //cout << "Error checking time: " << elapsed.count()*1000.0f << "ms" << endl;

    //cout << "Finished testing rows " << start_row << "-" << start_row + row_batch_size - 1 << endl;
}

// check if the candicate row groups have repeatable retention bitflips according to the RETPROF configuration parameters
// clears candidate_weaks
void analyze_weaks(SoftMCPlatform& platform, const vector<RowData>& rows_data, 
            vector<WeakRowSet>& candidate_weaks, vector<WeakRowSet>& row_group, const uint weak_rows_needed) {

    // vector<WeakRow> multi_it_weaks(RETPROF_NUM_ITS);

    for(auto& wr : candidate_weaks) {
        std::cout << BLUE_TXT << "Checking retention time consistency of row(s) " << wr.rows_as_str() << NORMAL_TXT << std::endl;
        char buf[ROW_SIZE*wr.row_group.size()];
        
        // Setting up a progress bar
        progresscpp::ProgressBar progress_bar(RETPROF_NUM_ITS, 70, '#', '-');
        progress_bar.display();
        
        bool success = true;
        for(uint i = 0; i < RETPROF_NUM_ITS; i++) {

            // std::cout << "Iteration: " << i + 1 << "/" << RETPROF_NUM_ITS << endl;

            // test whether the row experiences bitflips with RETPROF_RETTIME_MULT_H higher retention time
            if(!check_retention_failute_repeatability(platform, (int)wr.ret_ms*RETPROF_RETTIME_MULT_H, wr.bank_id, wr, rows_data, buf)) {
                progress_bar.done();
                std::cout << RED_TXT << "HIGH RETENTION CHECK FAILED" << NORMAL_TXT << std::endl;
                success = false;
                break;
            }

            // std::cout << YELLOW_TXT << "HIGH RETENTION CHECK SUCCEEDED" << NORMAL_TXT << std::endl;

            // test whether the row never experiences bitflips with RETPROF_RETTIME_MULT_L lower retention time
            if(!check_retention_failute_repeatability(platform, (int)wr.ret_ms*RETPROF_RETTIME_MULT_H*0.5f, wr.bank_id, wr, rows_data, buf, true)){
                progress_bar.done();
                std::cout << RED_TXT << "LOW RETENTION CHECK FAILED" << NORMAL_TXT << std::endl;
                success = false;
                break;
            }

            // std::cout << YELLOW_TXT << "LOW RETENTION CHECK SUCCEEDED" << NORMAL_TXT << std::endl;

            ++progress_bar;
            progress_bar.display();
        }

        if(!success)
            continue;

        progress_bar.done();

        std::cout << MAGENTA_TXT << "PASSED" << NORMAL_TXT << std::endl;
        row_group.push_back(std::move(wr));

        if(row_group.size() == weak_rows_needed)
            break;
    }

    candidate_weaks.clear();

    // Sort the rows in each wrs based on the physical row IDs
    for (auto& wr : row_group){
        std::sort(wr.row_group.begin(), wr.row_group.end(), 
            [] (const WeakRow& lhs, const WeakRow& rhs) {
                return to_physical_row_id(lhs.row_id) < to_physical_row_id(rhs.row_id);
            }
        );
    }
}

std::string wrs_to_string(const WeakRowSet& wrs) { 
    return JS::serializeStruct(wrs);
}

int main(int argc, char** argv)
{

    string out_filename = "./out.txt";
    int test_mode = 0;
    int target_bank = 1;
    int target_row = -1;
    int starting_ret_time = 64;
    int num_row_groups = 1;
    string row_group_pattern = "R-R"; // to search for rows that have specific distances among each other. 
    // For example, "R-R" (default) makes RowScout search for two rows that 1) are one row address apart and 2) have similar retention times.
    // Similarly, "RR" makes RowScout search for two rows that 1) have consecutive row addresses and 2) have similar retention times.
    // "R" makes RowScout search for any row that would experience a retention failure.

    int input_data_pattern = 1;
    vector<int> row_range{-1, -1};

    uint arg_log_phys_conv_scheme = 0;

    bool append_output = false;

    // try{
    options_description desc("RowScout Options");
    desc.add_options()
        ("help,h", "Prints this usage statement.")
        ("out,o", value(&out_filename)->default_value(out_filename), "Specifies a path for the output file.")
        ("bank,b", value(&target_bank)->default_value(target_bank), "Specifies the address of the bank to be profiled.")
        ("range", value<vector<int>>(&row_range)->multitoken(), "Specifies a range of row addresses (start and end values are both inclusive) to be profiled. By default, the range spans an entire bank.")
        ("init_ret_time,r", value(&starting_ret_time)->default_value(starting_ret_time), "Specifies the initial retention time (in milliseconds) to test the rows specified by --bank and --range. When RowScout cannot find a set of rows that satisfy the requirements specified by other options, RowScout increases the retention time used in profiling and repeats the profiling procedure.")
        ("row_group_pattern", value(&row_group_pattern)->default_value(row_group_pattern), "Specifies the distances among rows in a row group that RowScout must find. Must include only 'R' and '-'. Example values: R-R (two one-row-address-apart rows with similar retention times) , RR (two consecutively-addressed rows with similar retention times).")
        ("num_row_groups,w", value(&num_row_groups)->default_value(num_row_groups), "Specifies the number of row groups that RowScout must find.")
        ("log_phys_scheme", value(&arg_log_phys_conv_scheme)->default_value(arg_log_phys_conv_scheme), "Specifies how to convert logical row IDs to physical row ids and the other way around. Pass 0 (default) for sequential mapping, 1 for the mapping scheme typically used in Samsung chips.")
        ("input_data,i", value(&input_data_pattern)->default_value(input_data_pattern), "Specifies the data pattern to initialize rows with for profiling. Defined value are 0: random, 1: all ones, 2: all zeros, 3: colstripe (0101), 4: inverse colstripe (1010), 5: checkered (0101, 1010), 6: inverse checkered (1010, 0101)")
        ("append", bool_switch(&append_output), "When specified, the output is appended to the --out file (if it exists). Otherwise the --out file is cleared.")
        ;

    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);

    if (vm.count("help")) {
        cout << desc << endl;
        return 0;
    }

    assert(row_range.size() == 2 && "--row_range must specify exactly two tokens. E.g., <--row_range  0 5>");
    if(row_range[0] > row_range[1])
        swap(row_range[0], row_range[1]);

    if (row_range[0] < 0)
        row_range[0] = 0;

    if (row_range[1] < 0)
        row_range[1] = NUM_ROWS - 1;

    if (row_range[1] > (NUM_ROWS - 1)){
        cout << "Specified row range exceeds the number of rows in the DRAM module. Adjusting the range accordingly." << endl;
        cout << "Specified: " << row_range[1] << ", actual num rows: " << NUM_ROWS << endl;
        row_range[1] = NUM_ROWS - 1;
    }

    // make sure row_group_pattern contains only R(r) or -
    if(row_group_pattern.find_first_not_of("Rr-") != std::string::npos) {
        cerr << RED_TXT << "ERROR: --row_group_pattern should contain only R or -" << NORMAL_TXT << std::endl;
        exit(-1);
    }

    // make sure row_group_pattern is uppercase
    for (auto& c : row_group_pattern) c = toupper(c);

    // the pattern should start with R and end with R
    row_group_pattern = row_group_pattern.substr(row_group_pattern.find_first_of('R'), row_group_pattern.find_last_of('R') - row_group_pattern.find_first_of('R') + 1);

    init_row_pattern_fitter(row_group_pattern);


    path out_dir(out_filename);
    out_dir = out_dir.parent_path();

    if (!(exists(out_dir))) {
        if (!create_directory(out_dir)) {
            cerr << "Cannot create directory: " << out_dir << ". Exiting..." << endl;
            return -1;
        }
    }

    boost::filesystem::ofstream out_file;
    if(append_output)
        out_file.open(out_filename, boost::filesystem::ofstream::app);
    else
        out_file.open(out_filename);

    vector<RowData> rows_data;
    
    SoftMCPlatform platform;
    int err;

    if((err = platform.init()) != SOFTMC_SUCCESS){
        cerr << "Could not initialize SoftMC Platform: " << err << endl;
        return err;
    }

    platform.reset_fpga();
    // sleep(1);
  
    // disable refresh
    platform.set_aref(false);

    // init random data generator
    srand(0);


    assert(arg_log_phys_conv_scheme < uint(LogPhysRowIDScheme::MAX));
    logical_physical_conversion_scheme = (LogPhysRowIDScheme) arg_log_phys_conv_scheme;

  
    bitset<512> bitset_int_mask(0xFFFFFFFF);

    auto t_prog_started = chrono::high_resolution_clock::now();
    auto t_two_rows_recvd = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed;
    bool check_time;

    target_row = 0;

    const uint default_data_patterns[] = {0x0, 0xFFFFFFFF, 0x00000000, 0x55555555, 0xAAAAAAAA, 0xAAAAAAAA, 0x55555555};

    // this is ugly but I am leaving it like this to make things easier in case we decide to use different input data patterns for different rows.
    vector<int> input_data_patterns = {input_data_pattern};
    for (int inp_pat : input_data_patterns) {
        RowData rd;
        bitset<512> rdata;

        switch(inp_pat) {
            case 0: {//random
                // GENERATING RANDOM TEST DATA
                uint32_t rand_int;

                for (int pos = 0; pos < 16; pos ++) {
                    rdata <<= 32;
                        rand_int = (rand() << 16) | (0x0000FFFF & rand());
                        //cout << "generated random 32-bit: " << hex << rand_int << dec << endl;
                    rdata |= rand_int;
                }
                break; }
            case 1:
            case 2: // 1's for victim rows, 0's for aggressor rows
            case 3: 
            case 4: 
            case 5: 
            case 6: { 
                for (int pos = 0; pos < 16; pos ++) {
                    rdata <<= 32;
                    rdata |= default_data_patterns[inp_pat];
                }

                break; }
            default: {
                cerr << "Undefined input data pattern mode: " << inp_pat << endl;
                return -1;
                break; }
        }

        rd.input_data_pattern = rdata;
        rd.pattern_id = inp_pat;
        rows_data.push_back(rd);
    }

    int retention_ms = starting_ret_time;
    uint64_t buf_size = 0;
    char* buf = nullptr;
    vector<WeakRowSet> candidate_weaks;
    vector<WeakRowSet> row_group;

    uint num_wrs_written_out = 0;

    uint last_num_weak_rows = 0;

    // write out profiler configuration to the output file
    // out_file << "RETPROF_NUM_ITS: " << RETPROF_NUM_ITS << std::endl;
    // // out_file << "RETPROF_SUCCESSFUL_ITS_THRESH: " << RETPROF_SUCCESSFUL_ITS_THRESH << std::endl;
    // out_file << "RETPROF_RETTIME_MULT_H: " << RETPROF_RETTIME_MULT_H << std::endl;
    // out_file << "RETPROF_RETTIME_STEP: " << RETPROF_RETTIME_STEP << std::endl;
    // out_file << "============" << std::endl;

    while(true) {

        std::cout << "Profiling with " << retention_ms << " ms retention time" << std::endl;

        uint max_row_batch_size = determineRowBatchSize(retention_ms, rows_data.size());
        uint target_region_size = row_range[1] - row_range[0] + 1;
        uint row_batch_size = min(max_row_batch_size, target_region_size);

        // check the size of the buffer to read the data to and increase its size if needed
        if(buf_size < row_batch_size*ROW_SIZE) {
            if(buf != nullptr)
                delete[] buf;

            buf = new char[row_batch_size*ROW_SIZE];
            buf_size = row_batch_size*ROW_SIZE;
        }

        // apply the retention time to the corresponding row region
        uint num_profiled_rows = 0;
        while(num_profiled_rows < target_region_size) {
            clear_bitflip_history();
            test_retention(platform, retention_ms, target_bank, row_range[0], row_batch_size, rows_data, row_group_pattern, buf, candidate_weaks);

            if(candidate_weaks.size() > 0) {
                // remove rows already identified as weak from candidate_weaks
                for (auto& wr : row_group) {
                    for (auto it = candidate_weaks.begin(); it != candidate_weaks.end(); it++) {
                        if(wr.contains_rows_in(*it))
                            candidate_weaks.erase(it--);
                    }
                }
            }

            // analyze the candidate row groups to ensure the bitflips are repeatable and the retention time is determined accurately (i.e., we do not want a cell to fail for periods smaller than the determined retention time)
            if(candidate_weaks.size() > 0) {
                std::cout << RED_TXT << "Found " << candidate_weaks.size() << " new candidate row groups." << NORMAL_TXT << std::endl;
                analyze_weaks(platform, rows_data, candidate_weaks, row_group, num_row_groups);
            }

            while (num_wrs_written_out < row_group.size()) {
                out_file << wrs_to_string(row_group[num_wrs_written_out++]) << std::endl;
            }

            if(row_group.size() >= num_row_groups)
                break;

            num_profiled_rows += row_batch_size;
        }

        auto cur_time = chrono::high_resolution_clock::now();
        elapsed = cur_time - t_prog_started;
        //cout << "Time for reading two rows: " << elapsed.count()*1000 << "ms" << endl;
        std::cout << GREEN_TXT << "[" << (int) elapsed.count() << " s] Found " << row_group.size() - last_num_weak_rows << " new (" << 
            row_group.size() << " total) row groups" << NORMAL_TXT << std::endl;
        last_num_weak_rows = row_group.size();

        if(row_group.size() >= num_row_groups)
            break;

        retention_ms += (int)(starting_ret_time*RETPROF_RETTIME_STEP); 
    }

    // checkForLeftoverPCIeData(platform);
    out_file.close();

    std::cout << "The test has finished!" << endl;

    

    return 0;
}
