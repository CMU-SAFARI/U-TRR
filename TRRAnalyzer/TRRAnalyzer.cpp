#include "instruction.h"
#include "prog.h"
#include "platform.h"
#include "tools/perfect_hash.h"
#include "tools/json_struct.h"
#include "tools/softmc_utils.h"
#include "tools/ProgressBar.hpp"

#include <string>
#include <fstream>
#include <iostream>
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

// #define PRINT_SOFTMC_PROGS

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
float DEFAULT_TREFI = 7800.0f;
/******/

int trcd_cycles = (int) ceil(DEFAULT_TRCD/FPGA_PERIOD);
int tras_cycles = (int) ceil(DEFAULT_TRAS/FPGA_PERIOD);
int trp_cycles = (int) ceil(DEFAULT_TRP/FPGA_PERIOD);
int twr_cycles = (int) ceil(DEFAULT_TWR/FPGA_PERIOD);
int trfc_cycles = (int) ceil(DEFAULT_TRFC/FPGA_PERIOD);
int trrds_cycles = (int) ceil(DEFAULT_TRRDS/FPGA_PERIOD);
int trrdl_cycles = (int) ceil(DEFAULT_TRRDL/FPGA_PERIOD);
int trefi_cycles = (int) ceil(DEFAULT_TREFI/FPGA_PERIOD);

// TRR Analyzer Parameters
const float TRR_RETTIME_MULT          = 1.2f;
const uint  TRR_CHECK_HAMMERS         = 500000;
const uint  TRR_DUMMY_ROW_DIST        = 2; // the minimum row distance between dummy rows
const uint  TRR_WEAK_DUMMY_DIST        = 5000; // the minimum row distance between weak and dummy aggressor rows
const uint  TRR_ALLOWED_RET_TIME_DIFF = 64;

const uint default_data_patterns[] = {0x0, 0xFFFFFFFF, 0x00000000, 0x55555555, 0xAAAAAAAA, 0xAAAAAAAA, 0x55555555};

vector<uint32_t> reserved_regs{CASR, BASR, RASR};

typedef struct RowData {
    bitset<512> input_data_pattern;
    uint pattern_id;
    string label;
} RowData;


typedef struct WeakRow {
    uint row_id;
    std::vector<uint> bitflip_locs;
    
    WeakRow(){}
    WeakRow(uint _row_id, vector<uint> _bitflip_locs) : 
            row_id(_row_id), bitflip_locs(_bitflip_locs) {}
} WeakRow;

typedef struct HammerableRowSet {
    std::vector<uint> victim_ids;
    std::vector<std::vector<uint>> vict_bitflip_locs;
    std::vector<std::vector<uint>> uni_bitflip_locs;
    std::vector<uint> aggr_ids;
    std::vector<uint> uni_ids;
    bitset<512> data_pattern;
    uint bank_id;
    uint ret_ms;
} HammerableRowSet;

JS_OBJECT_EXTERNAL(WeakRow,
                JS_MEMBER(row_id),
                JS_MEMBER(bitflip_locs));

typedef struct WeakRowSet {
    std::vector<WeakRow> row_group;
    uint bank_id;
    uint ret_ms;
    uint index_in_file;
    uint data_pattern_type;
    uint rowdata_ind;
    WeakRowSet(){}
    WeakRowSet(std::vector<WeakRow> _weak_rows, uint _bank_id, uint _ret_ms, uint _data_pattern_type, uint _rowdata_ind) : 
            row_group(_weak_rows), bank_id(_bank_id), ret_ms(_ret_ms), data_pattern_type(_data_pattern_type), rowdata_ind(_rowdata_ind) {}
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
void collect_bitflips(vector<uint>& bitflips, const char* read_data, const bitset<512>& input_data_pattern, const vector<uint> bitflips_loc) {

    bitset<512> read_data_bitset;
    uint bit_loc;

    uint32_t* iread_data = (uint32_t*) read_data;

    // check for bitflips in each cache line
    if(bitflips_loc.size() == 0){
        // check for bitflips in each cache line
        for(int cl = 0; cl < ROW_SIZE/64; cl++) {

            read_data_bitset.reset();
            for(int i = 0; i < 512/32; i++) {
                bitset<512> tmp_bitset = iread_data[cl*(512/32) + i];

                read_data_bitset |= (tmp_bitset << i*32);
            }

            // compare and print errors
            bitset<512> error_mask = read_data_bitset ^ input_data_pattern;

            if(error_mask.any()) {
                // there is at least one bitflip in this cache line
                for(uint i = 0; i < error_mask.size(); i++){
                    if(error_mask.test(i)){
                        bit_loc = cl*CACHE_LINE_BITS + i;
                        bitflips.push_back(bit_loc);
                    }
                }
            }
        }
    }else{
        for(auto bitflip: bitflips_loc){
            uint cl = floor(bitflip/CACHE_LINE_BITS);
            uint offset = bitflip%CACHE_LINE_BITS;

            read_data_bitset.reset();
            for(int i = 0; i < 512/32; i++) {
                bitset<512> tmp_bitset = iread_data[cl*(512/32) + i];
                read_data_bitset |= (tmp_bitset << i*32);
            }

            bitset<512> error_mask = read_data_bitset ^ input_data_pattern;

            if(error_mask.test(offset)){
                bitflips.push_back(bitflip);
            }
        }
    }
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

void waitForRetentionTime(const uint ret_time_ms) {
    
    static constexpr chrono::duration<double, milli> min_sleep_duration(1);
    auto start = chrono::high_resolution_clock::now();
    while (chrono::duration<double, milli>(chrono::high_resolution_clock::now() - start).count() < ret_time_ms) {
        this_thread::sleep_for(min_sleep_duration);
    }
}

void readFromDRAM(Program& program, const uint target_bank, const uint start_row, const uint row_batch_size) {

    // const int REG_TMP_WRDATA = 15;
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

std::string wrs_to_string(const WeakRowSet& wrs) { 
    return JS::serializeStruct(wrs);
}

bool getNextJSONObj(boost::filesystem::ifstream& f_row_groups, string& s_weak) {
    
    std::string cur_line;
    std::getline(f_row_groups, cur_line);

    s_weak = "";
    while (cur_line != "}") {
        if(f_row_groups.eof())
            return false;

        s_weak += cur_line;
        std::getline(f_row_groups, cur_line);
    }
    s_weak += cur_line;

    return true;
}

void parse_all_weaks(boost::filesystem::ifstream& f_row_groups, vector<WeakRowSet>& row_groups) {
    string s_weak;

    uint i = 0;
    while(getNextJSONObj(f_row_groups, s_weak)) {
        JS::ParseContext context(s_weak);
        WeakRowSet wr;
        context.parseTo(wr);
        wr.index_in_file = i++;
        row_groups.push_back(wr);
    }
}

void pick_weaks(boost::filesystem::ifstream& f_row_groups, vector<WeakRowSet>& all_weak_rows, vector<WeakRowSet>& picked_weak_rows, const uint num_row_groups) {

    for(uint i = picked_weak_rows.size(); i < num_row_groups; i++){

        auto it = all_weak_rows.begin();

        if(it == all_weak_rows.end()) {
            std::cerr << RED_TXT << "ERROR: The weaks rows file does not contain a sufficient number of hammerable weak rows" << NORMAL_TXT << std::endl;
            std::cerr << RED_TXT << "Needed: " << num_row_groups << ", found: " << picked_weak_rows.size() << NORMAL_TXT << std::endl;
            exit(-1);
        }

        if(picked_weak_rows.size() > 0) {
            // remove picked rows that have different retention times
            for(auto it_picked = picked_weak_rows.begin(); it_picked != picked_weak_rows.end(); it_picked++) {
                if(std::abs((int)it_picked->ret_ms - (int)it->ret_ms) > TRR_ALLOWED_RET_TIME_DIFF) {
                    picked_weak_rows.erase(it_picked--);
                    i--;
                }
            }
        }

        picked_weak_rows.push_back(*it);
        all_weak_rows.erase(it);
    }
}

void init_row_data(Program& prog, SoftMCRegAllocator& reg_alloc, const SMC_REG reg_bank_addr, const SMC_REG reg_num_cols, 
    const vector<uint>& rows_to_init, const vector<bitset<512>>& data_patts) {

    uint initial_free_regs = reg_alloc.num_free_regs();

    bitset<512> bitset_int_mask(0xFFFFFFFF);

    prog.add_inst(SMC_LI(8, CASR)); // Load 8 into CASR since each WRITE writes 8 columns
    SMC_REG reg_row_addr = reg_alloc.allocate_SMC_REG();
    SMC_REG reg_col_addr = reg_alloc.allocate_SMC_REG();

    assert(rows_to_init.size() == data_patts.size());
    for(uint i = 0; i < rows_to_init.size(); i++){
        uint target_row = rows_to_init[i];

        prog.add_inst(SMC_LI(target_row, reg_row_addr));

        if(i == 0 || data_patts[i-1] != data_patts[i]) {
            // set up the input data in the wide register
            SMC_REG reg_wrdata = reg_alloc.allocate_SMC_REG();
            for (int pos = 0; pos < 16; pos++) {
                prog.add_inst(SMC_LI((((data_patts[i] >> 32*pos) & bitset_int_mask).to_ulong() & 0xFFFFFFFF), reg_wrdata));
                prog.add_inst(SMC_LDWD(reg_wrdata, pos));
            }
            reg_alloc.free_SMC_REG(reg_wrdata);
        }

        // activate the target row
        uint remaining = add_op_with_delay(prog, SMC_ACT(reg_bank_addr, 0, reg_row_addr, 0), 0, trcd_cycles - 5);

        // write data to the row and precharge
        add_op_with_delay(prog, SMC_LI(0, reg_col_addr), remaining, 0);

        string new_lbl = createSMCLabel("INIT_ROW");
        prog.add_label(new_lbl);
        add_op_with_delay(prog, SMC_WRITE(reg_bank_addr, 0, reg_col_addr, 1, 0, 0), 0, 0);
        prog.add_branch(prog.BR_TYPE::BL, reg_col_addr, reg_num_cols, new_lbl);
        
        // precharge the open bank
        add_op_with_delay(prog, SMC_PRE(reg_bank_addr, 0, 0), 0, trp_cycles);
    }

    reg_alloc.free_SMC_REG(reg_row_addr);
    reg_alloc.free_SMC_REG(reg_col_addr);

    assert(reg_alloc.num_free_regs() == initial_free_regs);
}

void init_data_row_range(Program& prog, SoftMCRegAllocator& reg_alloc, const SMC_REG reg_bank_addr, const SMC_REG reg_num_cols, 
    const uint first_row_id, const uint last_row_id, const bitset<512>& data_patt) {

    uint initial_free_regs = reg_alloc.num_free_regs();

    bitset<512> bitset_int_mask(0xFFFFFFFF);

    prog.add_inst(SMC_LI(8, CASR)); // Load 8 into CASR since each WRITE writes 8 columns
    SMC_REG reg_row_addr = reg_alloc.allocate_SMC_REG();
    SMC_REG reg_col_addr = reg_alloc.allocate_SMC_REG();

    SMC_REG reg_last_row_id = reg_alloc.allocate_SMC_REG();
    prog.add_inst(SMC_LI(first_row_id, reg_row_addr));
    prog.add_inst(SMC_LI(last_row_id + 1, reg_last_row_id)); // +1 because we don't have BL_TYPE::BLE

    // initialize the wide register that contains data to write to DRAM
    SMC_REG reg_wrdata = reg_alloc.allocate_SMC_REG();
    for (int pos = 0; pos < 16; pos++) {
        prog.add_inst(SMC_LI((((data_patt >> 32*pos) & bitset_int_mask).to_ulong() & 0xFFFFFFFF), reg_wrdata));
        prog.add_inst(SMC_LDWD(reg_wrdata, pos));
    }

    std::string lbl_init_row_it = createSMCLabel("INIT_ROW_IT");
    prog.add_label(lbl_init_row_it);
        // activate the target row
        uint remaining = add_op_with_delay(prog, SMC_ACT(reg_bank_addr, 0, reg_row_addr, 1), 0, trcd_cycles - 5);

        // write data to the row and precharge
        add_op_with_delay(prog, SMC_LI(0, reg_col_addr), remaining, 0);

        string new_lbl = createSMCLabel("INIT_ROW");
        prog.add_label(new_lbl);
        add_op_with_delay(prog, SMC_WRITE(reg_bank_addr, 0, reg_col_addr, 1, 0, 0), 0, 0);
        prog.add_branch(prog.BR_TYPE::BL, reg_col_addr, reg_num_cols, new_lbl);
        
        // precharge the open bank
        add_op_with_delay(prog, SMC_PRE(reg_bank_addr, 0, 0), 0, 0); // no need for tRP since branch latency is 24 cycles
    prog.add_branch(prog.BR_TYPE::BL, reg_row_addr, reg_last_row_id, lbl_init_row_it);
        
    reg_alloc.free_SMC_REG(reg_row_addr);
    reg_alloc.free_SMC_REG(reg_col_addr);
    reg_alloc.free_SMC_REG(reg_last_row_id);
    reg_alloc.free_SMC_REG(reg_wrdata);

    assert(reg_alloc.num_free_regs() == initial_free_regs);
}

void init_row_data(Program& prog, SoftMCRegAllocator& reg_alloc, const SMC_REG reg_bank_addr, const SMC_REG reg_num_cols, 
    const uint target_row, const bitset<512>& data_pattern) {

    vector<uint> rows_to_init{target_row};
    vector<bitset<512>> data_patts{data_pattern};
    init_row_data(prog, reg_alloc, reg_bank_addr, reg_num_cols, rows_to_init, data_patts);
}

void hammer_aggressors(Program& prog, SoftMCRegAllocator& reg_alloc, const SMC_REG reg_bank_addr, const vector<uint>& rows_to_hammer,
                        const std::vector<uint>& num_hammers, const bool cascaded_hammer, const uint hammer_duration);

void init_hammerable_row_set(Program& prog, SoftMCRegAllocator& reg_alloc, const SMC_REG reg_bank_addr, const SMC_REG reg_num_cols, 
        const HammerableRowSet& hr, const bool init_aggrs_first, const bool ignore_aggrs, const bool init_only_victims){
    vector<uint> rows_to_init;
    vector<bitset<512>> data_patts;

    uint total_rows = hr.victim_ids.size() + hr.aggr_ids.size() + hr.uni_ids.size();
    rows_to_init.reserve(total_rows);
    data_patts.reserve(total_rows);

    for(auto& vict : hr.victim_ids) {
        rows_to_init.push_back(vict);
        data_patts.push_back(hr.data_pattern);
    }

    for(auto& uni : hr.uni_ids) {
        rows_to_init.push_back(uni);
        data_patts.push_back(hr.data_pattern);
    }

    bitset<512> aggr_data_patt = hr.data_pattern;
    aggr_data_patt.flip();
    assert(aggr_data_patt != hr.data_pattern);

    if(!ignore_aggrs && !init_only_victims) {
        for(auto& aggr : hr.aggr_ids) {
            if(init_aggrs_first) {
                rows_to_init.insert(rows_to_init.begin(), aggr);
                data_patts.insert(data_patts.begin(), aggr_data_patt);
            } else {
                rows_to_init.push_back(aggr);
                data_patts.push_back(aggr_data_patt);
            }
        }
    }

    init_row_data(prog, reg_alloc, reg_bank_addr, reg_num_cols, rows_to_init, data_patts);
}


void init_HRS_data(Program& prog, SoftMCRegAllocator& reg_alloc, const SMC_REG reg_bank_addr, const SMC_REG reg_num_cols, const HammerableRowSet& hr,
                    const bool init_aggrs_first, const bool ignore_aggrs, const bool init_only_victims) {

    // // initialize data of the two rows on the sides as well (if they are not out of bounds)
    // // we do this because these rows can be aggressor rows in the TRRAnalyzer experiments and
    // // they may affect the retention of the first and last rows. Since we collect the first retention failures
    // // using a big batch of rows (where these side rows are likely to be initialized), we should do the same here

    // // As a slight change to the note above, let's not initialize the aggressors on the sides (unless the rh_type forces that)
    // // to make it easier to find which ACTs are sampling in Hynix modules 

    // init_data_row_range(prog, reg_alloc, reg_bank_addr, reg_num_cols, first_row_id, last_row_id, hr.data_pattern);

    // UPDATE 09.12.2020 - Instead of initializing a range of rows, now we initilize the victims first and then the aggressors one by one.
    // This is to make aggressors the last row activated prior to a refresh when no rows are hammering during the hammering phase. 
    // This change is useful for analyzing the sampling method of Hynix modules
    init_hammerable_row_set(prog, reg_alloc, reg_bank_addr, reg_num_cols, hr, init_aggrs_first, ignore_aggrs, init_only_victims);

    
}


void hammer_aggressors(Program& prog, SoftMCRegAllocator& reg_alloc, const SMC_REG reg_bank_addr, const vector<uint>& rows_to_hammer,
                        const std::vector<uint>& num_hammers, const bool cascaded_hammer, const uint hammer_duration) {

    if(rows_to_hammer.size() < 1)
        return; // nothing to hammer

    uint initial_free_regs = reg_alloc.num_free_regs();
    uint remaining_cycs = 0;

    SMC_REG reg_row_addr = reg_alloc.allocate_SMC_REG();
    SMC_REG reg_cur_hammers = reg_alloc.allocate_SMC_REG();
    SMC_REG reg_num_hammers = reg_alloc.allocate_SMC_REG();

    if(!cascaded_hammer){

        // it is complicated to efficiently hammer rows different number of times while activating them one after another
        // We implement the following algorithm:
        // 1. If there is a non-zero element in hammer_per_ref, hammer the rows corresponding to those elements using the smallest non-zero hammer_per_ref value. If all hammer_per_ref elements are zero, exit
        // 2. decrement all non-zero elements of hammer_per_ref vector by the smallest value
        // 3. go back to 1

        auto hammers_per_round = num_hammers;

        while (1) {

            auto min_non_zero = std::min_element(hammers_per_round.begin(), hammers_per_round.end(), 
                    [](const uint& a, const uint& b) {return ((a > 0) && (a < b)) || (b == 0);}
                );

            if (min_non_zero == hammers_per_round.end() || *min_non_zero == 0) {
                break;
            }

            uint min_elem = *min_non_zero;

            // perform hammering
            prog.add_inst(SMC_LI(min_elem, reg_num_hammers));
            prog.add_inst(SMC_LI(0, reg_cur_hammers));
            std::string lbl_rh = createSMCLabel("ROWHAMMERING");
            prog.add_label(lbl_rh);
            for (int ind_row = 0; ind_row < rows_to_hammer.size(); ind_row++) {
                if(hammers_per_round[ind_row] == 0) // do not anymore hammer a row that has 0 remaining hammers
                    continue;

                int row_id = rows_to_hammer[ind_row];
                prog.add_inst(SMC_LI(row_id, reg_row_addr));

                if(hammer_duration < 20)
                    remaining_cycs = add_op_with_delay(prog, SMC_ACT(reg_bank_addr, 0, reg_row_addr, 0), 0, tras_cycles + hammer_duration - 1);
                else {
                    remaining_cycs = add_op_with_delay(prog, SMC_ACT(reg_bank_addr, 0, reg_row_addr, 0), 0, hammer_duration % 4);
                    remaining_cycs = add_op_with_delay(prog, SMC_SLEEP(std::floor(hammer_duration/4.0f)), remaining_cycs, tras_cycles - 1);
                }
                    
                remaining_cycs = add_op_with_delay(prog, SMC_PRE(reg_bank_addr, 0, 0), 0, trp_cycles - 5);
            }

            prog.add_inst(SMC_ADDI(reg_cur_hammers, 1, reg_cur_hammers));
            prog.add_branch(Program::BR_TYPE::BL, reg_cur_hammers, reg_num_hammers, lbl_rh);


            // this subtracts min_elem from every non-zero element
            std::for_each(hammers_per_round.begin(), hammers_per_round.end(), [&](uint& a) {if (a > 0) a -= min_elem;});
        }

        
    } else { // cascaded_hammer == true
        for (int ind_row = 0; ind_row < rows_to_hammer.size(); ind_row++) {
            int row_id = rows_to_hammer[ind_row];

            if(num_hammers[ind_row] == 0) // do not hammer rows with 0 hammer count
                continue;

            prog.add_inst(SMC_LI(row_id, reg_row_addr));
            prog.add_inst(SMC_LI(num_hammers[ind_row], reg_num_hammers));
            prog.add_inst(SMC_LI(0, reg_cur_hammers));

            string lbl_rh = createSMCLabel("ROWHAMMERING");
            prog.add_label(lbl_rh);

            if(hammer_duration < 20)
                remaining_cycs = add_op_with_delay(prog, SMC_ACT(reg_bank_addr, 0, reg_row_addr, 0), 0, tras_cycles + hammer_duration - 1);
            else {
                remaining_cycs = add_op_with_delay(prog, SMC_ACT(reg_bank_addr, 0, reg_row_addr, 0), 0, 0);
                remaining_cycs = add_op_with_delay(prog, SMC_SLEEP(std::ceil(hammer_duration/4.0f)), 0, tras_cycles - 5);
            }

            remaining_cycs = add_op_with_delay(prog, SMC_PRE(reg_bank_addr, 0, 0), 0, 0);
            remaining_cycs = 0;
            prog.add_inst(SMC_ADDI(reg_cur_hammers, 1, reg_cur_hammers));
            prog.add_branch(Program::BR_TYPE::BL, reg_cur_hammers, reg_num_hammers, lbl_rh);
        }
    }

    reg_alloc.free_SMC_REG(reg_row_addr);
    reg_alloc.free_SMC_REG(reg_cur_hammers);
    reg_alloc.free_SMC_REG(reg_num_hammers);

    assert(reg_alloc.num_free_regs() == initial_free_regs);
}


void init_HRS_data(SoftMCPlatform& platform, const std::vector<HammerableRowSet>& vec_hr,
                    const bool init_aggrs_first, const bool ignore_aggrs, const bool init_only_victims,
                    const uint num_pre_init_bank0_hammers, const uint pre_init_nops,
                    Program* prog = nullptr, SoftMCRegAllocator* reg_alloc = nullptr) {

    bool exec_prog_and_clean = false;
    if (prog == nullptr) {
        prog = new Program();
        reg_alloc = new SoftMCRegAllocator(NUM_SOFTMC_REGS, reserved_regs);
        exec_prog_and_clean = true;
    }
    
    if(exec_prog_and_clean)
        add_op_with_delay(*prog, SMC_PRE(0, 0, 1), 0, 0); // precharge all banks

    SMC_REG reg_bank_addr = reg_alloc->allocate_SMC_REG();
    SMC_REG reg_num_cols = reg_alloc->allocate_SMC_REG();

    if(num_pre_init_bank0_hammers > 0) {
        prog->add_inst(SMC_LI(0, reg_bank_addr));
        std::vector<uint> rows_to_hammer = std::vector<uint>{0};
        std::vector<uint> bank0_hammers_per_ref = std::vector<uint>{num_pre_init_bank0_hammers};
        hammer_aggressors(*prog, *reg_alloc, reg_bank_addr, rows_to_hammer, bank0_hammers_per_ref, true, 0);
    }

    if(pre_init_nops > 0) {
        if(pre_init_nops < 3) {
            for(uint i = 0; i < pre_init_nops; i++)
                prog->add_inst(__pack_mininsts(SMC_NOP(), SMC_NOP(), SMC_NOP(), SMC_NOP()));
        } else {
            prog->add_inst(SMC_SLEEP(pre_init_nops));
        }
    }

    prog->add_inst(SMC_LI(NUM_COLS_PER_ROW*8, reg_num_cols));

    for(const auto& hrs : vec_hr) {
        prog->add_inst(SMC_LI(hrs.bank_id, reg_bank_addr));
        init_HRS_data(*prog, *reg_alloc, reg_bank_addr, reg_num_cols, hrs, init_aggrs_first, ignore_aggrs, init_only_victims);
    }

    reg_alloc->free_SMC_REG(reg_bank_addr);
    reg_alloc->free_SMC_REG(reg_num_cols);

    if(exec_prog_and_clean) {
        prog->add_inst(SMC_END());
        platform.execute(*prog);
        #ifdef PRINT_SOFTMC_PROGS
        std::cout << "--- SoftMCProg: Initializing Victims and Aggressors ---" << std::endl;
        prog->pretty_print();
        #endif

        delete prog;
        delete reg_alloc;
    }
}

void read_row_data(Program& prog, SoftMCRegAllocator& reg_alloc, const SMC_REG reg_bank_addr, const SMC_REG reg_num_cols, 
                    const vector<uint>& rows_to_read) {

    uint initial_free_regs = reg_alloc.num_free_regs();

    prog.add_inst(SMC_LI(8, CASR)); // Load 8 into CASR since each READ reads 8 columns
    SMC_REG reg_row_addr = reg_alloc.allocate_SMC_REG();

    for(auto target_row : rows_to_read) {
        prog.add_inst(SMC_LI(target_row, reg_row_addr));

        // activate the victim row
        add_op_with_delay(prog, SMC_ACT(reg_bank_addr, 0, reg_row_addr, 0), 0, trcd_cycles - 1);
        
        // read data from the row and precharge
        SMC_REG reg_col_addr = reg_alloc.allocate_SMC_REG();
        prog.add_inst(SMC_LI(0, reg_col_addr));

        string new_lbl = createSMCLabel("READ_ROW");
        prog.add_label(new_lbl);
        add_op_with_delay(prog, SMC_READ(reg_bank_addr, 0, reg_col_addr, 1, 0, 0), 0, 0);
        prog.add_branch(prog.BR_TYPE::BL, reg_col_addr, reg_num_cols, new_lbl);
        reg_alloc.free_SMC_REG(reg_col_addr);

        // precharge the open bank
        add_op_with_delay(prog, SMC_PRE(reg_bank_addr, 0, 0), 0, trp_cycles);
    }

    reg_alloc.free_SMC_REG(reg_row_addr);
    assert(reg_alloc.num_free_regs() == initial_free_regs);
}

void read_row_data(Program& prog, SoftMCRegAllocator& reg_alloc, const SMC_REG reg_bank_addr, const SMC_REG reg_num_cols, 
                    const WeakRowSet wrs) {

    vector<uint> rows_to_read;
    rows_to_read.reserve(wrs.row_group.size());
    for (auto& wr : wrs.row_group)
        rows_to_read.push_back(wr.row_id);
    read_row_data(prog, reg_alloc, reg_bank_addr, reg_num_cols, rows_to_read);
}

bitset<512> setup_data_pattern(const uint data_pattern_type) {
    
    bitset<512> data_pattern;
    switch(data_pattern_type) {
        case 0: {//random
            // GENERATING RANDOM TEST DATA
            uint32_t rand_int;

            for (int pos = 0; pos < 16; pos ++) {
                data_pattern <<= 32;
                    rand_int = (rand() << 16) | (0x0000FFFF & rand());
                    //cout << "generated random 32-bit: " << hex << rand_int << dec << endl;
                data_pattern |= rand_int;
            }
            break; }
        case 1:
        case 2: // 1's for victim rows, 0's for aggressor rows
        case 3: 
        case 4: 
        case 5: 
        case 6: { 
            for (int pos = 0; pos < 16; pos ++) {
                data_pattern <<= 32;
                data_pattern |= default_data_patterns[data_pattern_type];
            }

            break; }
        default: {
            std:: cerr << RED_TXT << "ERROR: Undefined input data pattern mode: " << data_pattern_type << NORMAL_TXT << std::endl;
            exit(-1);
            }
    }
    
    return data_pattern;
}

HammerableRowSet toHammerableRowSet(const WeakRowSet& wrs, const std::string row_layout) {
    HammerableRowSet hr;

    hr.bank_id = wrs.bank_id;
    hr.ret_ms = wrs.ret_ms;

    hr.victim_ids.reserve(wrs.row_group.size());

    uint wrs_ind = 0;
    uint vict_to_aggr_dist = 1;
    for(const char row_type : row_layout) {
        switch (row_type)
        {
        case 'r':
        case 'R':
            hr.victim_ids.push_back(wrs.row_group[wrs_ind].row_id);
            hr.vict_bitflip_locs.push_back(wrs.row_group[wrs_ind++].bitflip_locs);
            vict_to_aggr_dist = 1;
            break;
        case 'a':
        case 'A':
            hr.aggr_ids.push_back(to_logical_row_id(to_physical_row_id(hr.victim_ids.back()) + vict_to_aggr_dist));

            if(hr.aggr_ids.back() == wrs.row_group[wrs_ind].row_id) // advance wrs_ind if the corresponding row is profiled as a retention weak row but we would like to use it as an aggressor row
                wrs_ind++;

            vict_to_aggr_dist = 1;
            break;
        case 'u':
        case 'U':
            hr.uni_ids.push_back(wrs.row_group[wrs_ind].row_id);
            hr.uni_bitflip_locs.push_back(wrs.row_group[wrs_ind++].bitflip_locs);
            vict_to_aggr_dist = 1;
            break;
        case '-':
            vict_to_aggr_dist++;
            break;
        
        default:
            std::cerr << RED_TXT << "ERROR: Unexpected character in row_layout. row_layout: " << row_layout << ", unexpected char: " << row_type << NORMAL_TXT << std::endl;
            break;
        }
    }

    hr.data_pattern = setup_data_pattern(wrs.data_pattern_type);
    // for aggressor rows we use the bit inverse of victim's data pattern

    return hr;
}

// runs a small test that checks whether the weak rows in wrs can be hammered using aggressor rows determined based on the rh_type.
// If the aggressor rows are physically close to the weak rows, then we should observe RowHammer bitflips.
bool is_hammerable(SoftMCPlatform& platform, const WeakRowSet& wrs, const std::string row_layout, const bool cascaded_hammer) {
    
    Program p_testRH;
    uint target_bank = wrs.bank_id;

    /****************************/
    /*** Initialize DRAM data ***/
    /****************************/

    int remaining_cycs = 0;

    SoftMCRegAllocator reg_alloc(NUM_SOFTMC_REGS, reserved_regs);

    SMC_REG reg_bank_addr = reg_alloc.allocate_SMC_REG();
    p_testRH.add_inst(SMC_LI(target_bank, reg_bank_addr));

    add_op_with_delay(p_testRH, SMC_PRE(reg_bank_addr, 0, 1), 0, 0); // precharge all banks
    
    SMC_REG reg_num_cols = reg_alloc.allocate_SMC_REG();
    p_testRH.add_inst(SMC_LI(NUM_COLS_PER_ROW*8, reg_num_cols));

    p_testRH.add_inst(SMC_LI(8, CASR)); // Load 8 into CASR since each READ reads 8 columns
    p_testRH.add_inst(SMC_LI(1, BASR)); // Load 1 into BASR
    p_testRH.add_inst(SMC_LI(1, RASR)); // Load 1 into RASR

    HammerableRowSet hr = toHammerableRowSet(wrs, row_layout);

    // no need to check if the WRS is hammerable if there are no aggressor rows
    // A WRS won't have aggressor rows when using the 'R' row layout and BETWEEN_VICTIMS RHType
    if((hr.aggr_ids.size() + hr.uni_ids.size()) == 0)
        return true;

    // write to the victim and aggressor row(s)
    // first all victims, and then the aggressors are initialized
    init_hammerable_row_set(p_testRH, reg_alloc, reg_bank_addr, reg_num_cols, hr, false, false, false);


    /*************************/
    /*** Perform Hammering ***/
    /*************************/

    vector<uint> hammers(hr.aggr_ids.size(), TRR_CHECK_HAMMERS);
    hammer_aggressors(p_testRH, reg_alloc, reg_bank_addr, hr.aggr_ids, hammers, cascaded_hammer, 0);


    /********************************/
    /*** issue DRAM READ commands ***/
    /********************************/
    read_row_data(p_testRH, reg_alloc, reg_bank_addr, reg_num_cols, hr.victim_ids);
    p_testRH.add_inst(SMC_END());

    platform.execute(p_testRH);
    #ifdef PRINT_SOFTMC_PROGS
    std::cout << "--- SoftMCProg: Checking if the victims are hammerable ---" << std::endl;
    p_testRH.pretty_print();
    #endif

    /*********************************************/
    /*** read PCIe data and check for bitflips ***/
    /*********************************************/
    char buf[ROW_SIZE*hr.victim_ids.size()];
    platform.receiveData(buf, ROW_SIZE*hr.victim_ids.size());
    vector<uint> bitflips;

    // we expect all victim rows to be hammerable
    bool all_victims_have_bitflips = true;
    for(uint i = 0; i < hr.victim_ids.size(); i++) {
        bitflips.clear();
        collect_bitflips(bitflips, buf + ROW_SIZE*i, hr.data_pattern, vector<uint> {});

        if(bitflips.size() == 0) {
            all_victims_have_bitflips = false;
            std::cout << RED_TXT << "No RH bitflips found in row " << hr.victim_ids[i] << NORMAL_TXT << std::endl;
        }
    }

    return all_victims_have_bitflips;
}

void pick_dummy_aggressors(vector<uint>& dummy_aggrs, const uint dummy_aggrs_bank, const uint num_dummies, const vector<WeakRowSet>& weak_row_sets,
                            const uint dummy_ids_offset) {

    uint cur_dummy = TRR_DUMMY_ROW_DIST % NUM_ROWS;
    if(weak_row_sets.size() != 0)
        cur_dummy = (weak_row_sets[0].row_group[0].row_id + TRR_WEAK_DUMMY_DIST + dummy_ids_offset) % NUM_ROWS;

    const uint MAX_TRIES = 1000000;
    uint cur_try = 0;

    while(dummy_aggrs.size() != num_dummies) {
        if(cur_try == MAX_TRIES) {
            std::cerr << RED_TXT << "ERROR: Failed to pick dummy aggressor rows after " << MAX_TRIES << " tries" << NORMAL_TXT << std::endl;
            std::cerr << "Consider reducing the number of weak or dummy rows" << std::endl;
            exit(-1);
        }

        // check if there is a victim close to cur_dummy
        if(dummy_aggrs_bank == weak_row_sets[0].bank_id){
            for(auto& wrs : weak_row_sets){
                for(auto& wr : wrs.row_group){
                    if(std::abs((int)cur_dummy - (int)wr.row_id) < TRR_WEAK_DUMMY_DIST) {
                        cur_dummy = (cur_dummy + TRR_WEAK_DUMMY_DIST) % NUM_ROWS;
                        cur_try++;
                        continue;
                    }
                }
            }
        }

        dummy_aggrs.push_back(cur_dummy);
        cur_dummy = (cur_dummy + TRR_DUMMY_ROW_DIST) % NUM_ROWS;
    }
}

void perform_refresh(Program& prog, SoftMCRegAllocator& reg_alloc, const uint num_refs_per_round,
                        const uint pre_ref_delay) {

    SMC_REG reg_num_refs_per_cycle = reg_alloc.allocate_SMC_REG();
    SMC_REG reg_it_refs_per_cycle = reg_alloc.allocate_SMC_REG();

    prog.add_inst(SMC_LI(num_refs_per_round, reg_num_refs_per_cycle));
    prog.add_inst(SMC_LI(0, reg_it_refs_per_cycle));

    if (pre_ref_delay >= 8) {
        prog.add_inst(SMC_SLEEP(std::ceil(pre_ref_delay/4.0f)));
    }

    std::string lbl_issue_per_cycle_refs = createSMCLabel("PER_CYCLE_REFS");
    prog.add_label(lbl_issue_per_cycle_refs);    
        add_op_with_delay(prog, SMC_REF(), 0, 0);
        add_op_with_delay(prog, SMC_SLEEP(ceil((trfc_cycles - 1 - 24 - 4)/4.0f)), 0, 0);
        add_op_with_delay(prog, SMC_ADDI(reg_it_refs_per_cycle, 1, reg_it_refs_per_cycle), 0, 0);
    prog.add_branch(prog.BR_TYPE::BL, reg_it_refs_per_cycle, reg_num_refs_per_cycle, lbl_issue_per_cycle_refs);

    reg_alloc.free_SMC_REG(reg_num_refs_per_cycle);
    reg_alloc.free_SMC_REG(reg_it_refs_per_cycle);
}

void issue_REFs(SoftMCPlatform& platform, const uint num_refs,
                Program* prog = nullptr, SoftMCRegAllocator* reg_alloc = nullptr) {

    bool exec_prog_and_clean = false;
    if (prog == nullptr) {
        prog = new Program();
        reg_alloc = new SoftMCRegAllocator(NUM_SOFTMC_REGS, reserved_regs);
        exec_prog_and_clean = true;
    }

    SMC_REG reg_num_refs = reg_alloc->allocate_SMC_REG();
    SMC_REG reg_issued_refs = reg_alloc->allocate_SMC_REG();
  
    prog->add_inst(SMC_LI(num_refs, reg_num_refs));
    prog->add_inst(SMC_LI(0, reg_issued_refs));

    int remaining_cycs = 0;
    if(exec_prog_and_clean)
        remaining_cycs = add_op_with_delay(*prog, SMC_PRE(0, 0, 1), 0, trp_cycles); // precharge all banks

    if(remaining_cycs > 0) { 
        remaining_cycs = add_op_with_delay(*prog, SMC_NOP(), remaining_cycs, 0);
    }

    std::string lbl_issue_refs = createSMCLabel("ISSUE_REFS");
    prog->add_label(lbl_issue_refs);

    add_op_with_delay(*prog, SMC_REF(), 0, 0);

    // the SLEEP function waits for 1 SoftMC frontend logic cycle (4 * FPGA_PERIOD). Therefore, we divide by 4
    add_op_with_delay(*prog, SMC_SLEEP(ceil((trefi_cycles - 4 - 24 - 3)/4.0f)), 0, 0);

    prog->add_inst(SMC_ADDI(reg_issued_refs, 1, reg_issued_refs));
    prog->add_branch(prog->BR_TYPE::BL, reg_issued_refs, reg_num_refs, lbl_issue_refs);

    if(exec_prog_and_clean) {
        // perform a dummy read to signal the end of a program

        SMC_REG reg_bank_id = reg_alloc->allocate_SMC_REG();
        SMC_REG reg_row_id = reg_alloc->allocate_SMC_REG();
        SMC_REG reg_col_id = reg_alloc->allocate_SMC_REG();
        add_op_with_delay(*prog, SMC_LI(0, reg_bank_id), 0, 0);
        add_op_with_delay(*prog, SMC_LI(0, reg_row_id), 0, 0);
        add_op_with_delay(*prog, SMC_LI(0, reg_col_id), 0, 0);

        add_op_with_delay(*prog, SMC_ACT(reg_bank_id, 0, reg_row_id, 0), 0, trcd_cycles);
        add_op_with_delay(*prog, SMC_READ(reg_bank_id, 0, reg_col_id, 0, 0, 0), 0, tras_cycles - trcd_cycles);
        add_op_with_delay(*prog, SMC_PRE(reg_bank_id, 0, 0), 0, trp_cycles);

        prog->add_inst(SMC_END());
        platform.execute(*prog);
        #ifdef PRINT_SOFTMC_PROGS
        std::cout << "--- SoftMCProg: Issuing REFs ---" << std::endl;
        prog->pretty_print();
        #endif

        reg_alloc->free_SMC_REG(reg_bank_id);
        reg_alloc->free_SMC_REG(reg_row_id);
        reg_alloc->free_SMC_REG(reg_col_id);

        char cl[64];
        platform.receiveData(cl, 64);
    }    

    reg_alloc->free_SMC_REG(reg_num_refs);
    reg_alloc->free_SMC_REG(reg_issued_refs);

    if(exec_prog_and_clean) {
        delete prog;
        delete reg_alloc;
    }
}

// performing REF at nominal rate, i.e., a REF cmd is issued once every 7.8us
// dummy rows are hammered between the REF cmds
void hammer_dummies(SoftMCPlatform& platform, const uint bank_id, const vector<uint>& dummy_aggrs, const uint num_refs,
                    Program* prog = nullptr, SoftMCRegAllocator* reg_alloc = nullptr) {

    bool exec_prog_and_clean = false;
    if (prog == nullptr) {
        prog = new Program();
        reg_alloc = new SoftMCRegAllocator(NUM_SOFTMC_REGS, reserved_regs);
        exec_prog_and_clean = true;
    }

    SMC_REG reg_num_refs = reg_alloc->allocate_SMC_REG();
    SMC_REG reg_issued_refs = reg_alloc->allocate_SMC_REG();

    SMC_REG reg_row_id = reg_alloc->allocate_SMC_REG();
    SMC_REG reg_bank_id = reg_alloc->allocate_SMC_REG();
    prog->add_inst(SMC_LI(bank_id, reg_bank_id));
  
    prog->add_inst(SMC_LI(num_refs, reg_num_refs));
    prog->add_inst(SMC_LI(0, reg_issued_refs));

    int remaining_cycs = 0;
    if(exec_prog_and_clean)
        remaining_cycs = add_op_with_delay(*prog, SMC_PRE(0, 0, 1), 0, 0); // precharge all banks

    uint num_dummies = dummy_aggrs.size();
    uint cycs_hammer_dummies_once = (tras_cycles + trp_cycles)*num_dummies + 28;
    uint hammers_per_round = floor((trefi_cycles - trfc_cycles)/cycs_hammer_dummies_once);

    if(remaining_cycs > 0) { 
        remaining_cycs = add_op_with_delay(*prog, SMC_NOP(), remaining_cycs, 0);
    }

    std::string lbl_issue_refs = createSMCLabel("ISSUE_REFS");
    prog->add_label(lbl_issue_refs);
    
    add_op_with_delay(*prog, SMC_REF(), 0, trfc_cycles - 1);

        // hammer the dummy rows
        SMC_REG reg_hammer_it = reg_alloc->allocate_SMC_REG();
        SMC_REG reg_hammers_per_ref = reg_alloc->allocate_SMC_REG();
        prog->add_inst(SMC_LI(0, reg_hammer_it));
        prog->add_inst(SMC_LI(hammers_per_round, reg_hammers_per_ref));

        std::string lbl_hammer = createSMCLabel("AFTER_INIT_DUMMY_HAMMERING");
        prog->add_label(lbl_hammer);
            remaining_cycs = 0;
            for (uint dummy_row_id : dummy_aggrs) {
                prog->add_inst(SMC_LI(dummy_row_id, reg_row_id));

                remaining_cycs = add_op_with_delay(*prog, SMC_ACT(reg_bank_id, 0, reg_row_id, 0), remaining_cycs, tras_cycles - 1);
                remaining_cycs = add_op_with_delay(*prog, SMC_PRE(reg_bank_id, 0, 0), remaining_cycs, trp_cycles - 1);
            }

            add_op_with_delay(*prog, SMC_ADDI(reg_hammer_it, 1, reg_hammer_it), 0, 0);
        prog->add_branch(prog->BR_TYPE::BL, reg_hammer_it, reg_hammers_per_ref, lbl_hammer);


    prog->add_inst(SMC_ADDI(reg_issued_refs, 1, reg_issued_refs));
    prog->add_branch(prog->BR_TYPE::BL, reg_issued_refs, reg_num_refs, lbl_issue_refs);

    if(exec_prog_and_clean) {
        // perform a dummy read to signal the end of a program
    
        SMC_REG reg_col_id = reg_alloc->allocate_SMC_REG();
        add_op_with_delay(*prog, SMC_LI(0, reg_bank_id), 0, 0);
        add_op_with_delay(*prog, SMC_LI(0, reg_row_id), 0, 0);
        add_op_with_delay(*prog, SMC_LI(0, reg_col_id), 0, 0);

        add_op_with_delay(*prog, SMC_ACT(reg_bank_id, 0, reg_row_id, 0), 0, trcd_cycles);
        add_op_with_delay(*prog, SMC_READ(reg_bank_id, 0, reg_col_id, 0, 0, 0), 0, tras_cycles - trcd_cycles);
        add_op_with_delay(*prog, SMC_PRE(reg_bank_id, 0, 0), 0, trp_cycles);

        prog->add_inst(SMC_END());
        platform.execute(*prog);
        #ifdef PRINT_SOFTMC_PROGS
        std::cout << "--- SoftMCProg: Hammering Dummy Aggressors ---" << std::endl;
        prog->pretty_print();
        #endif

        reg_alloc->free_SMC_REG(reg_col_id);

        char cl[64];
        platform.receiveData(cl, 64);
    }

    reg_alloc->free_SMC_REG(reg_num_refs);
    reg_alloc->free_SMC_REG(reg_issued_refs);
    reg_alloc->free_SMC_REG(reg_bank_id);
    reg_alloc->free_SMC_REG(reg_row_id);
    reg_alloc->free_SMC_REG(reg_hammer_it);
    reg_alloc->free_SMC_REG(reg_hammers_per_ref);    

    if(exec_prog_and_clean) {
        delete prog;
        delete reg_alloc;
    }
}


void hammer_hrs(SoftMCPlatform& platform,
                const vector<HammerableRowSet>& hammerable_rows, const std::vector<uint>& hammers_per_round, const bool cascaded_hammer,
                const uint num_rounds, const bool skip_hammering_aggr, const bool ignore_dummy_hammers, 
                const uint hammer_duration, const uint num_refs_per_round, const uint pre_ref_delay,
                const vector<uint>& dummy_aggrs, const uint dummy_aggrs_bank, const bool hammer_dummies_first, const bool hammer_dummies_independently,
                const uint num_bank0_hammers = 0, Program* prog = nullptr, SoftMCRegAllocator* reg_alloc = nullptr) {

    bool exec_prog_and_clean = false;
    if (prog == nullptr) {
        prog = new Program();
        reg_alloc = new SoftMCRegAllocator(NUM_SOFTMC_REGS, reserved_regs);
        exec_prog_and_clean = true;
    }

    SMC_REG reg_bank_addr = reg_alloc->allocate_SMC_REG();
    prog->add_inst(SMC_LI(hammerable_rows[0].bank_id, reg_bank_addr));

    if(exec_prog_and_clean)
        add_op_with_delay(*prog, SMC_PRE(reg_bank_addr, 0, 1), 0, 0); // precharge all banks


    SMC_REG reg_cur_its = reg_alloc->allocate_SMC_REG();
    SMC_REG reg_refresh_cycles = reg_alloc->allocate_SMC_REG();
    prog->add_inst(SMC_LI(0, reg_cur_its));
    prog->add_inst(SMC_LI(num_rounds, reg_refresh_cycles));

    string lbl_hammer_loop = createSMCLabel("HAMMER_MAIN_LOOP");
    prog->add_label(lbl_hammer_loop);

    // hammering all aggressor and dummy rows
    std::vector<uint> all_rows_to_hammer;
    std::vector<uint> weak_rows_to_hammer;

    if(hammer_dummies_first && !hammer_dummies_independently && !ignore_dummy_hammers)
        all_rows_to_hammer.insert(all_rows_to_hammer.end(), dummy_aggrs.begin(), dummy_aggrs.end());

    
    for(auto& hrs : hammerable_rows) {
        if(!skip_hammering_aggr) {
            all_rows_to_hammer.insert(all_rows_to_hammer.end(), hrs.aggr_ids.begin(), hrs.aggr_ids.end());
            all_rows_to_hammer.insert(all_rows_to_hammer.end(), hrs.uni_ids.begin(), hrs.uni_ids.end());
        }
        weak_rows_to_hammer.insert(weak_rows_to_hammer.end(), hrs.aggr_ids.begin(), hrs.aggr_ids.end());
        weak_rows_to_hammer.insert(weak_rows_to_hammer.end(), hrs.uni_ids.begin(), hrs.uni_ids.end());
    }

    if(!hammer_dummies_first & !hammer_dummies_independently & !ignore_dummy_hammers)
        all_rows_to_hammer.insert(all_rows_to_hammer.end(), dummy_aggrs.begin(), dummy_aggrs.end());

    
    uint total_hammers_per_ref = 0;
    
    std::vector<uint> t_dummy_hammers_per_ref;
    std::vector<uint> t_aggr_hammers_per_ref;
    if(hammer_dummies_first) {
        t_dummy_hammers_per_ref = std::vector<uint>(hammers_per_round.begin(), hammers_per_round.begin() + dummy_aggrs.size());
        t_aggr_hammers_per_ref = std::vector<uint>(hammers_per_round.begin() + dummy_aggrs.size(), hammers_per_round.end());
    } else {
        t_dummy_hammers_per_ref = std::vector<uint>(hammers_per_round.begin() + weak_rows_to_hammer.size(), hammers_per_round.end());
        t_aggr_hammers_per_ref = std::vector<uint>(hammers_per_round.begin(), hammers_per_round.begin() + weak_rows_to_hammer.size());
    }

    std::vector<uint> new_hammers_per_ref;

    if(hammer_dummies_first && !ignore_dummy_hammers && !hammer_dummies_independently)
        new_hammers_per_ref.insert(new_hammers_per_ref.end(), t_dummy_hammers_per_ref.begin(), t_dummy_hammers_per_ref.end());

    if(!skip_hammering_aggr)
        new_hammers_per_ref.insert(new_hammers_per_ref.end(), t_aggr_hammers_per_ref.begin(), t_aggr_hammers_per_ref.end());

    if(!hammer_dummies_first && !ignore_dummy_hammers && !hammer_dummies_independently)
        new_hammers_per_ref.insert(new_hammers_per_ref.end(), t_dummy_hammers_per_ref.begin(), t_dummy_hammers_per_ref.end());

    if(!hammer_dummies_independently)
        assert(new_hammers_per_ref.size() == all_rows_to_hammer.size());
        

    for (uint i = 0; i < new_hammers_per_ref.size(); i++) {
        uint hammer_count = new_hammers_per_ref[i];
        total_hammers_per_ref += hammer_count;
    }

    if (total_hammers_per_ref == 0 && num_refs_per_round == 0)
        return; // no need to launch SoftMC program since there is nothing to do (i.e., to rows to hammer and no refs to do)


    if(hammer_dummies_independently && hammer_dummies_first && !ignore_dummy_hammers){
        prog->add_inst(SMC_LI(dummy_aggrs_bank, reg_bank_addr));

        std::vector<uint32_t> dummy_hammers_per_round;
        
        if(hammer_dummies_first)
            dummy_hammers_per_round = std::vector<uint32_t>(hammers_per_round.begin(), hammers_per_round.begin() + dummy_aggrs.size());
        else
            dummy_hammers_per_round = std::vector<uint32_t>(hammers_per_round.begin() + weak_rows_to_hammer.size(), hammers_per_round.end());
        
        hammer_aggressors(*prog, *reg_alloc, reg_bank_addr, dummy_aggrs, dummy_hammers_per_round, cascaded_hammer, 0);
    }

    if(total_hammers_per_ref > 0){

        if(dummy_aggrs_bank != hammerable_rows[0].bank_id) {
            assert(cascaded_hammer && "ERROR: Dummy aggressors can be in a different bank only when using cascaded hammering. This feature is not yet implemented for sequential hammering.");

            if(hammer_dummies_first && !hammer_dummies_independently && !ignore_dummy_hammers){
                prog->add_inst(SMC_LI(dummy_aggrs_bank, reg_bank_addr));
                auto dummy_hammers_per_round = std::vector<uint32_t>(hammers_per_round.begin(), hammers_per_round.begin() + dummy_aggrs.size());
                hammer_aggressors(*prog, *reg_alloc, reg_bank_addr, dummy_aggrs, dummy_hammers_per_round, cascaded_hammer, 0);
            }

            if(!skip_hammering_aggr) {
                prog->add_inst(SMC_LI(hammerable_rows[0].bank_id, reg_bank_addr));
                
                auto weak_rows_hammers_per_ref = std::vector<uint32_t>(
                    (hammer_dummies_first & !hammer_dummies_independently) ? hammers_per_round.begin() + dummy_aggrs.size() : hammers_per_round.begin(),
                    (hammer_dummies_first & !hammer_dummies_independently) ? hammers_per_round.end() : hammers_per_round.begin() + weak_rows_to_hammer.size());
                hammer_aggressors(*prog, *reg_alloc, reg_bank_addr, weak_rows_to_hammer, weak_rows_hammers_per_ref, cascaded_hammer, hammer_duration);
            }

            if(!hammer_dummies_first && !hammer_dummies_independently && !ignore_dummy_hammers){
                prog->add_inst(SMC_LI(dummy_aggrs_bank, reg_bank_addr));
                auto dummy_hammers_per_round = std::vector<uint32_t>(hammers_per_round.begin() + weak_rows_to_hammer.size(), hammers_per_round.end());
                hammer_aggressors(*prog, *reg_alloc, reg_bank_addr, dummy_aggrs, dummy_hammers_per_round, cascaded_hammer, 0);
            }

        } else {
            hammer_aggressors(*prog, *reg_alloc, reg_bank_addr, all_rows_to_hammer, new_hammers_per_ref, cascaded_hammer, hammer_duration);
        }
    }

    if(hammer_dummies_independently && !hammer_dummies_first && !ignore_dummy_hammers){
        prog->add_inst(SMC_LI(dummy_aggrs_bank, reg_bank_addr));

        std::vector<uint32_t> dummy_hammers_per_round;
        
        if(hammer_dummies_first)
            dummy_hammers_per_round = std::vector<uint32_t>(hammers_per_round.begin(), hammers_per_round.begin() + dummy_aggrs.size());
        else
            dummy_hammers_per_round = std::vector<uint32_t>(hammers_per_round.begin() + weak_rows_to_hammer.size(), hammers_per_round.end());
        
        hammer_aggressors(*prog, *reg_alloc, reg_bank_addr, dummy_aggrs, dummy_hammers_per_round, cascaded_hammer, 0);
    }

    
    if (num_bank0_hammers > 0) {
        prog->add_inst(SMC_LI(0, reg_bank_addr));
        all_rows_to_hammer = std::vector<uint>{0};
        std::vector<uint> bank0_hammers_per_ref = std::vector<uint>{num_bank0_hammers};
        hammer_aggressors(*prog, *reg_alloc, reg_bank_addr, all_rows_to_hammer, bank0_hammers_per_ref, cascaded_hammer, hammer_duration);
    }

    // 4) issue a REF
    if(num_rounds > 0 && num_refs_per_round > 0) {
        perform_refresh(*prog, *reg_alloc, num_refs_per_round, pre_ref_delay);

        prog->add_inst(SMC_ADDI(reg_cur_its, 1, reg_cur_its));
        prog->add_branch(Program::BR_TYPE::BL, reg_cur_its, reg_refresh_cycles, lbl_hammer_loop); // 5) repeat hammering + REF num_rounds times
    }

    if(exec_prog_and_clean) {
        prog->add_inst(SMC_END());
        platform.execute(*prog);
        #ifdef PRINT_SOFTMC_PROGS
        std::cout << "--- SoftMCProg: Hammering the Aggressor Rows ---" << std::endl;
        prog->pretty_print();
        #endif
    }

    reg_alloc->free_SMC_REG(reg_cur_its);
    reg_alloc->free_SMC_REG(reg_refresh_cycles);
    reg_alloc->free_SMC_REG(reg_bank_addr);

    if(exec_prog_and_clean) {
        delete prog;
        delete reg_alloc;
    }
}

void perform_dummy_read(SoftMCPlatform& platform) {

    Program prog;
    SoftMCRegAllocator reg_alloc = SoftMCRegAllocator(NUM_SOFTMC_REGS, reserved_regs);

    
    SMC_REG reg_bank_id = reg_alloc.allocate_SMC_REG();
    SMC_REG reg_row_id = reg_alloc.allocate_SMC_REG();
    SMC_REG reg_col_id = reg_alloc.allocate_SMC_REG();

    add_op_with_delay(prog, SMC_LI(0, reg_bank_id), 0, 0);
    add_op_with_delay(prog, SMC_LI(0, reg_row_id), 0, 0);
    add_op_with_delay(prog, SMC_LI(0, reg_col_id), 0, 0);

    add_op_with_delay(prog, SMC_ACT(reg_bank_id, 0, reg_row_id, 0), 0, trcd_cycles);
    add_op_with_delay(prog, SMC_READ(reg_bank_id, 0, reg_col_id, 0, 0, 0), 0, tras_cycles - trcd_cycles);
    add_op_with_delay(prog, SMC_PRE(reg_bank_id, 0, 0), 0, trp_cycles);

    prog.add_inst(SMC_END());
    platform.execute(prog);
    #ifdef PRINT_SOFTMC_PROGS
    std::cout << "--- SoftMCProg: Performing a dummy read ---" << std::endl;
    prog.pretty_print();
    #endif


    // we put 64 bytes to the PCie bus. We need to clear this data
    char cl[64];
    platform.receiveData(cl, 64);
}

void waitMS_softmc(const uint ret_time_ms, Program* prog) {
    // convert milliseconds to SoftMC cycles

    ulong cycs = std::ceil((ret_time_ms*1000000)/FPGA_PERIOD);

    // cycs is the number of DDR cycles now, convert it to FPGA cycles by dividing it by 4
    cycs = std::ceil(cycs/4.0f);

    prog->add_inst(SMC_SLEEP(cycs));
}

vector<vector<uint>> analyzeTRR(SoftMCPlatform& platform, const vector<HammerableRowSet>& hammerable_rows, const vector<uint>& dummy_aggrs, 
                    const uint dummy_aggrs_bank, const uint dummy_hammers_per_round,
                    const bool hammer_dummies_first, const bool hammer_dummies_independently, const bool cascaded_hammer, const std::vector<uint>& hammers_per_round,
                    const float hammer_cycle_time, const uint hammer_duration, const uint num_rounds, const bool skip_hammering_aggr, const uint refs_after_init, const vector<uint>& after_init_dummies,
                    const bool init_aggrs_first, const bool ignore_aggrs, const bool init_only_victims, const bool ignore_dummy_hammers, const bool first_it_aggr_init_and_hammer,
                    const bool refs_after_init_no_dummy_hammer, const uint num_refs_per_round, const uint pre_ref_delay, const std::vector<uint>& hammers_before_wait,
                    const float init_to_hammerbw_delay, const uint num_bank0_hammers, const uint num_pre_init_bank0_hammers,
                    const uint pre_init_nops,
                    const bool use_single_softmc_prog, const uint num_iterations, const bool verbose) {


    Program single_prog;
    SoftMCRegAllocator single_prog_reg_alloc = SoftMCRegAllocator(NUM_SOFTMC_REGS, reserved_regs);

    // create the iteration loop
    SMC_REG reg_iter_counter = single_prog_reg_alloc.allocate_SMC_REG();
    SMC_REG reg_num_iters = single_prog_reg_alloc.allocate_SMC_REG();

    add_op_with_delay(single_prog, SMC_PRE(reg_iter_counter, 0, 1), 0, 0); // precharge all banks
    single_prog.add_inst(SMC_LI(0, reg_iter_counter));
    single_prog.add_inst(SMC_LI(num_iterations, reg_num_iters));

    std::string lbl_iter_loop = createSMCLabel("MAIN_ITERATION_LOOP");
    single_prog.add_label(lbl_iter_loop);


    // // 1) initialize the data of the entire row range from the smallest row id to the largest row id in each HammerableRowSet
    std::string lbl_init_end = createSMCLabel("INIT_ROWS_END");

    auto t_start_init_data = chrono::high_resolution_clock::now();    
    if(!skip_hammering_aggr) {
        if (!use_single_softmc_prog)
            init_HRS_data(platform, hammerable_rows, init_aggrs_first, ignore_aggrs, init_only_victims, num_pre_init_bank0_hammers, pre_init_nops);
        else {
            std::string lbl_init_all = createSMCLabel("INIT_ALL_ROWS");

            if (first_it_aggr_init_and_hammer) {
                SMC_REG reg_zero = single_prog_reg_alloc.allocate_SMC_REG();
                single_prog.add_inst(SMC_LI(0, reg_zero));
                single_prog.add_branch(Program::BR_TYPE::BEQ, reg_iter_counter, reg_zero, lbl_init_all);
                single_prog_reg_alloc.free_SMC_REG(reg_zero);

                init_HRS_data(platform, hammerable_rows, init_aggrs_first, ignore_aggrs, true/*init_only_victims*/, num_pre_init_bank0_hammers, pre_init_nops, &single_prog, &single_prog_reg_alloc);
                
                single_prog.add_branch(Program::BR_TYPE::JUMP, reg_iter_counter, reg_zero, lbl_init_end);
            }

            single_prog.add_label(lbl_init_all);
            init_HRS_data(platform, hammerable_rows, init_aggrs_first, ignore_aggrs, init_only_victims, num_pre_init_bank0_hammers, pre_init_nops, &single_prog, &single_prog_reg_alloc);
        }
    }

    if (use_single_softmc_prog) {
        single_prog.add_label(lbl_init_end);
    }

    auto t_end_issue_prog = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> prog_issue_duration(t_end_issue_prog - t_start_init_data);

    /*** OPTIONAL - issue REF commands after initializing data ***/
    // this is an attempt to reset any REF related state
    if(refs_after_init > 0) {
        // auto t_start_issue_refs = chrono::high_resolution_clock::now();
        if(after_init_dummies.size() == 0){
            if (!use_single_softmc_prog)
                issue_REFs(platform, refs_after_init);
            else
                issue_REFs(platform, refs_after_init, &single_prog, &single_prog_reg_alloc);
        } else {
            if (!use_single_softmc_prog)
                hammer_dummies(platform, dummy_aggrs_bank, after_init_dummies, refs_after_init);
            else
                hammer_dummies(platform, dummy_aggrs_bank, after_init_dummies, refs_after_init, &single_prog, &single_prog_reg_alloc);
        }

        if(refs_after_init_no_dummy_hammer){
            if (!use_single_softmc_prog)
                issue_REFs(platform, refs_after_init);
            else
                issue_REFs(platform, refs_after_init, &single_prog, &single_prog_reg_alloc);
        }
            
        // auto t_end_issue_refs = chrono::high_resolution_clock::now();
        // chrono::duration<double, milli> ref_issue_duration(t_end_issue_refs - t_start_issue_refs);
        // std::cout << YELLOW_TXT << "Completed issuing 8192 REF commands in (ms): " << ref_issue_duration.count() << NORMAL_TXT << std::endl;

        t_end_issue_prog = chrono::high_resolution_clock::now();
    }


    // 2) wait until half of the target retention time, i.e., hammering_start_time = (ret_ms*H_MODIFIER - (num_rounds*refresh_cycle_time))/2
    /* SMC_WAIT() allows sleeping up to only ~6 seconds
       need to use the system timer like in the RetentionProfiler to sleep longer */

    // calculating the time for activating all aggressor rows once

    uint num_all_aggrs = 0;
    for (auto& hr : hammerable_rows){
        num_all_aggrs += hr.aggr_ids.size();
        num_all_aggrs += hr.uni_ids.size();
    }

    uint act_pre_cycles = std::ceil(std::max(DEFAULT_TRAS + DEFAULT_TRP + hammer_duration*FPGA_PERIOD, hammer_cycle_time)/FPGA_PERIOD);
    uint act_all_rows_cycles = act_pre_cycles*(num_all_aggrs + dummy_aggrs.size()) + 24 /*for the branch inst*/;

    // calculating overall time that will be spend on hammering and refreshing
    // NOTE: this is a bit of overestimation as we assume each aggressor and dummy will be activated max. amount of times
    uint max_hammer_acts = *(std::max_element(hammers_per_round.begin(), hammers_per_round.end()));
    max_hammer_acts = std::max(max_hammer_acts, dummy_hammers_per_round);
    uint total_hammer_cycles = max_hammer_acts*act_all_rows_cycles + trfc_cycles*num_refs_per_round + pre_ref_delay;

    float total_hammer_ms = std::ceil(total_hammer_cycles*FPGA_PERIOD)/1000000.0f;

    uint wait_interval_ms = (hammerable_rows[0].ret_ms*TRR_RETTIME_MULT - (num_rounds*total_hammer_ms))/2;

    const uint c_wait_interval_ms = wait_interval_ms;

    if(verbose) {
        std::cout << YELLOW_TXT << "Weak row retention time (ms): " << hammerable_rows[0].ret_ms << NORMAL_TXT << std::endl;
        // std::cout << YELLOW_TXT << "pre_ref_delay (cycles): " << pre_ref_delay << NORMAL_TXT << std::endl;
        // std::cout << YELLOW_TXT << "total_hammer_cycles (cycles): " << total_hammer_cycles << NORMAL_TXT << std::endl;
        std::cout << YELLOW_TXT << "Time to complete hammering phase (ms): " << total_hammer_ms << NORMAL_TXT << std::endl;
        // std::cout << YELLOW_TXT << "TRR_RETTIME_MULT: " << TRR_RETTIME_MULT << NORMAL_TXT << std::endl;
        // std::cout << YELLOW_TXT << "prog_issue_duration (ms): " << prog_issue_duration.count() << NORMAL_TXT << std::endl;
        // std::cout << YELLOW_TXT << "Waiting for (ms): " << wait_interval_ms - prog_issue_duration.count() << NORMAL_TXT << std::endl;
    }


    // hammers the rows here as well if hammers_before_wait contains non zero element
    if(init_to_hammerbw_delay > 0.0f) {
        uint wait_ms = wait_interval_ms*init_to_hammerbw_delay;
        wait_interval_ms -= wait_ms;

        if(!use_single_softmc_prog)
            waitMS(wait_ms);
        else
            waitMS_softmc(wait_ms, &single_prog);
    }

    if(hammers_before_wait.size() > 0){
        if(!use_single_softmc_prog)
            hammer_hrs(platform, hammerable_rows, hammers_before_wait, cascaded_hammer, 1, skip_hammering_aggr | ignore_aggrs, ignore_dummy_hammers,
                    hammer_duration, 0, pre_ref_delay, dummy_aggrs, dummy_aggrs_bank, hammer_dummies_first, hammer_dummies_independently);
        else
            hammer_hrs(platform, hammerable_rows, hammers_before_wait, cascaded_hammer, 1, skip_hammering_aggr | ignore_aggrs, ignore_dummy_hammers,
                    hammer_duration, 0, pre_ref_delay, dummy_aggrs, dummy_aggrs_bank, hammer_dummies_first, hammer_dummies_independently, 0, &single_prog, &single_prog_reg_alloc);
    }
    
    if(!skip_hammering_aggr) {
        if(!use_single_softmc_prog)
            waitMS(wait_interval_ms/* - prog_issue_duration.count()*/);
        else
            waitMS_softmc(wait_interval_ms, &single_prog);
    }

    
    // 3) Perform hammering based on rh_type and hammers_per_round
    auto t_start_hammering = chrono::high_resolution_clock::now();

    if(!use_single_softmc_prog)
        hammer_hrs(platform, hammerable_rows, hammers_per_round, cascaded_hammer, num_rounds, skip_hammering_aggr | ignore_aggrs, ignore_dummy_hammers,
                hammer_duration, num_refs_per_round, pre_ref_delay, dummy_aggrs, dummy_aggrs_bank, hammer_dummies_first, hammer_dummies_independently, num_bank0_hammers);
    else {
        std::string lbl_hammer_all = createSMCLabel("HAMMER_ALL");
        std::string lbl_hammer_end = createSMCLabel("HAMMER_END");
        if (first_it_aggr_init_and_hammer) {
            SMC_REG reg_zero = single_prog_reg_alloc.allocate_SMC_REG();
            single_prog.add_inst(SMC_LI(0, reg_zero));
            
            single_prog.add_branch(Program::BR_TYPE::BEQ, reg_iter_counter, reg_zero, lbl_hammer_all);
            single_prog_reg_alloc.free_SMC_REG(reg_zero);

            hammer_hrs(platform, hammerable_rows, hammers_per_round, cascaded_hammer, num_rounds, true /*skip_hammering_aggr | ignore_aggrs*/, ignore_dummy_hammers,
                hammer_duration, num_refs_per_round, pre_ref_delay, dummy_aggrs, dummy_aggrs_bank, hammer_dummies_first, hammer_dummies_independently, 
                num_bank0_hammers, &single_prog, &single_prog_reg_alloc); 
            
            single_prog.add_branch(Program::BR_TYPE::JUMP, reg_iter_counter, reg_zero, lbl_hammer_end);
        }

        single_prog.add_label(lbl_hammer_all);

        hammer_hrs(platform, hammerable_rows, hammers_per_round, cascaded_hammer, num_rounds, skip_hammering_aggr | ignore_aggrs, ignore_dummy_hammers,
                hammer_duration, num_refs_per_round, pre_ref_delay, dummy_aggrs, dummy_aggrs_bank, hammer_dummies_first, hammer_dummies_independently, 
                num_bank0_hammers, &single_prog, &single_prog_reg_alloc);

        single_prog.add_label(lbl_hammer_end);

    }


    auto t_end_hammering = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> dur_hammering(t_end_hammering - t_start_hammering);
    chrono::duration<double, milli> dur_from_start(t_end_hammering - t_end_issue_prog);

    // 5) wait until the retention time of the weak rows is satisfied (ret_ms * H_MODIFIER)

    if(skip_hammering_aggr) { // no need to wait and read any rows back if we are going to skip operating on the weak rows
        // just returning an empty vector
        return vector<vector<uint>>();
    }


    // std::cout << YELLOW_TXT << "(2nd Wait) Waiting for (ms): " << hammerable_rows[0].ret_ms*TRR_RETTIME_MULT - dur_from_start.count() << NORMAL_TXT << std::endl;
    if(!use_single_softmc_prog)
        waitMS(hammerable_rows[0].ret_ms*TRR_RETTIME_MULT - dur_from_start.count());
    else
        // we cannot use the measured time interval 'dur_from_start' when executing the experiment as a single program
        // Therefore, we use the calculated time here
        waitMS_softmc(c_wait_interval_ms, &single_prog);
    

    // 6) read back the weak rows and check for bitflips
    Program* prog_read = nullptr;
    SoftMCRegAllocator* reg_alloc = nullptr;

    bool exec_prog_and_clean = false;
    if (!use_single_softmc_prog) {
        prog_read = new Program();
        reg_alloc = new SoftMCRegAllocator(NUM_SOFTMC_REGS, reserved_regs);
        exec_prog_and_clean = true;
    } else
    {
        prog_read = &single_prog;
        reg_alloc = &single_prog_reg_alloc;
    }
    

    SMC_REG reg_bank_addr = reg_alloc->allocate_SMC_REG();
    SMC_REG reg_num_cols = reg_alloc->allocate_SMC_REG();
    if(exec_prog_and_clean)
        add_op_with_delay(*prog_read, SMC_PRE(reg_bank_addr, 0, 1), 0, 0); // precharge all banks
    
    prog_read->add_inst(SMC_LI(NUM_COLS_PER_ROW*8, reg_num_cols));

    ulong total_victim_rows = 0;
    for(auto& hrs : hammerable_rows) {
        prog_read->add_inst(SMC_LI(hrs.bank_id, reg_bank_addr));
        auto rows_to_read = hrs.victim_ids;
        rows_to_read.insert(rows_to_read.end(), hrs.uni_ids.begin(), hrs.uni_ids.end());

        read_row_data(*prog_read, *reg_alloc, reg_bank_addr, reg_num_cols, rows_to_read);
        total_victim_rows += rows_to_read.size();
    }

    if(exec_prog_and_clean) {
        prog_read->add_inst(SMC_END());
        platform.execute(*prog_read);
        #ifdef PRINT_SOFTMC_PROGS
        std::cout << "--- SoftMCProg: Reading the Victim rows ---" << std::endl;
        prog_read->pretty_print();
        #endif

        delete prog_read;
        delete reg_alloc;
    } else {
        assert(prog_read == &single_prog);

        // close the iteration loop
        single_prog.add_inst(SMC_ADDI(reg_iter_counter, 1, reg_iter_counter));
        single_prog.add_branch(single_prog.BR_TYPE::BL, reg_iter_counter, reg_num_iters, lbl_iter_loop);

        single_prog.add_inst(SMC_END());
        platform.execute(single_prog);
        #ifdef PRINT_SOFTMC_PROGS
        std::cout << "--- SoftMCProg: Running the experiments as a single SoftMC program ---" << std::endl;
        single_prog.pretty_print();
        #endif

        single_prog_reg_alloc.free_SMC_REG(reg_iter_counter);
        single_prog_reg_alloc.free_SMC_REG(reg_num_iters);
    }

    vector<vector<uint>> loc_bitflips;
    if(!use_single_softmc_prog) {
        // get data from PCIe
        ulong read_data_size = ROW_SIZE*total_victim_rows;
        char buf[read_data_size*2];
        platform.receiveData(buf, read_data_size);
        // std::cout << BLUE_TXT << "Successfully read all the data!" << NORMAL_TXT << std::endl;

        vector<uint> bitflips;
        loc_bitflips.reserve(total_victim_rows);

        uint row_it = 0;
        for (auto& hrs : hammerable_rows) {
            for(uint vict_ind = 0; vict_ind < hrs.victim_ids.size(); vict_ind++) {
                bitflips.clear();
                collect_bitflips(bitflips, buf + row_it*ROW_SIZE, hrs.data_pattern, hrs.vict_bitflip_locs[vict_ind]);
                row_it++;

                loc_bitflips.push_back(bitflips);

                // if(bitflips.size() == 0)
                //     std::cout << RED_TXT << "[Weak Row " << row_groups[i].row_id << "] Did not find any bitflips." << NORMAL_TXT << std::endl;
                // else
                //     std::cout << GREEN_TXT << "[Weak Row " << row_groups[i].row_id << "] Found " << bitflips.size() << " bitflip(s)." << NORMAL_TXT << std::endl;
            }

            for(uint uni_ind = 0; uni_ind < hrs.uni_ids.size(); uni_ind++) {
                bitflips.clear();
                collect_bitflips(bitflips, buf + row_it*ROW_SIZE, hrs.data_pattern, hrs.uni_bitflip_locs[uni_ind]);
                row_it++;

                loc_bitflips.push_back(bitflips);
            }
        }
        assert(row_it == total_victim_rows);
    } else {
        // allocate a large enough buffer to read the entire experiment data at once
        // get data from PCIe

        // UPDATE: we can actually read iteration by iteration since each iteration takes a long time to complete 
        // because of the wait times. So, it is quite unlikely that the PCIe will be a bottleneck
        // Moving this to the main function where we update the progress bar
        
        // ulong read_data_size = ROW_SIZE*total_victim_rows*num_iterations;
        // char* buf = new char[read_data_size];
        // platform.receiveData(buf, read_data_size);

        // vector<uint> bitflips;
        // num_bitflips.reserve(total_victim_rows*num_iterations);

        // for(uint iter_id = 0; iter_id < num_iterations; iter_id++) {
        //     uint row_it = 0;
        //     for (auto& hrs : hammerable_rows) {
        //         for(uint vict : hrs.victim_ids) {
        //             bitflips.clear();
        //             collect_bitflips(bitflips, buf + iter_id*total_victim_rows*ROW_SIZE + row_it*ROW_SIZE, hrs.data_pattern);
        //             row_it++;

        //             num_bitflips.push_back(bitflips.size());

        //             // if(bitflips.size() == 0)
        //             //     std::cout << RED_TXT << "[Weak Row " << row_groups[i].row_id << "] Did not find any bitflips." << NORMAL_TXT << std::endl;
        //             // else
        //             //     std::cout << GREEN_TXT << "[Weak Row " << row_groups[i].row_id << "] Found " << bitflips.size() << " bitflip(s)." << NORMAL_TXT << std::endl;
        //         }
        //     }
        //     assert(row_it == total_victim_rows);
        // }

        // delete[] buf;

        
    }

    return loc_bitflips;
}

// Finds and returns a subset of the retention-profiled rows in wrs that match the provided row_layout
WeakRowSet adjust_wrs(const WeakRowSet& wrs, const std::string& row_layout){

    // calculate row_layout distance vector
    std::vector<uint> wrs_type_dists;
    wrs_type_dists.reserve(row_layout.size());

    int first_r_ind = -1;
    for (uint i = 0; i < row_layout.size(); i++){
        switch (row_layout[i])
        {
            case 'u':
            case 'U':
            case 'r':
            case 'R':{
                if (first_r_ind == -1) {
                    first_r_ind = i;
                    break;
                }

                wrs_type_dists.push_back(i - first_r_ind);
                break;
            }
        }
    }
    assert(first_r_ind != -1 && "ERROR: There must be at least one R or U in the rowlayout (i.e., row_layout)");


    WeakRowSet new_wrs = wrs;
    bool match = false;

    std::vector<uint> wrs_dists;
    wrs_dists.reserve(new_wrs.row_group.size());
    while(true) {
        // finding the row distance between the rows in wrs
        wrs_dists.clear();

        uint first_victim_id = to_physical_row_id(new_wrs.row_group[0].row_id);
        for (uint i = 1; i < new_wrs.row_group.size(); i++)
            wrs_dists.push_back(to_physical_row_id(new_wrs.row_group[i].row_id) - first_victim_id);
        // for RRRRR, wrs_dists would be 1, 2, 3, 4
        // RARAR should match it, which would have dist vector 2 4
        // we should pick weak rows 0, 2, and 4

        if (wrs_dists.size() < wrs_type_dists.size())
            break;

        // the distance vectors are sorted so we can use includes() to determine if wrs_type_dists is a subset of wrs_dists
        match = std::includes(wrs_dists.begin(), wrs_dists.end(), wrs_type_dists.begin(), wrs_type_dists.end());

        if (match)
            break;

        new_wrs.row_group.erase(new_wrs.row_group.begin()); // remove the first row id and check if the remaining match
    }
    
    assert(wrs_dists.size() != 0 || match); // does includes() return true when wrs_type_dists is empty?


    if (match) { // we have a match
        // build a new row id vector that includes only the matching rows

        std::vector<WeakRow> matching_weak_rows;
        matching_weak_rows.reserve(wrs_type_dists.size());
        matching_weak_rows.push_back(new_wrs.row_group[0]);

        uint new_wrs_it = 0;
        for(auto d : wrs_type_dists) { // evict the rows that do not match the dist vector
            for (; new_wrs_it < wrs_dists.size(); new_wrs_it++){
                assert(wrs_dists[new_wrs_it] <= d);
                if (wrs_dists[new_wrs_it] == d){ // found a match
                    matching_weak_rows.push_back(new_wrs.row_group[new_wrs_it + 1]);
                    break;
                }
            }
        }

        new_wrs.row_group = matching_weak_rows;

        return new_wrs;
    }


    std::cerr << RED_TXT << "ERROR: The provided --row_layout does not match the data read from --file_weaks" << NORMAL_TXT << std::endl;
    std::cerr << RED_TXT << "--row_layout: " << row_layout << NORMAL_TXT << std::endl;
    std::cerr << RED_TXT << "Exiting..." << NORMAL_TXT << std::endl;
    exit(-1);

    return new_wrs;
}

void pick_hammerable_row_groups_from_file(SoftMCPlatform& platform, boost::filesystem::ifstream& f_row_groups, vector<WeakRowSet>& row_groups, const uint num_row_groups,
                                        const bool cascaded_hammer, const std::string row_layout) {

    vector<WeakRowSet> all_weaks;
    all_weaks.reserve(100);
    parse_all_weaks(f_row_groups, all_weaks);

    while(row_groups.size() != num_row_groups) {
        // 1) Pick (in order) 'num_weaks' weak rows from 'file_weak_rows' that have the same retention time.
        pick_weaks(f_row_groups, all_weaks, row_groups, num_row_groups);
        
        // cout << "Num picked weak rows: " << row_groups.size() << std::endl;
        // for(auto& wr : row_groups) {
        //     std::cout << "B: " << wr.bank_id << ", R: " << wr.row_id << ", ret: " << wr.ret_ms << " ms" << std::endl;
        // }

        // 2) test whether RowHammer bitflips can be induced on the weak rows
        for (auto it = row_groups.begin(); it != row_groups.end(); it++) {
            if(!is_hammerable(platform, *it, row_layout, cascaded_hammer)) {
                std::cout << RED_TXT << "Candidate victim row set " << it->rows_as_str() << " is not hammerable" << NORMAL_TXT << std::endl;
                row_groups.erase(it--);
                continue;
            }

            std::cout << GREEN_TXT << "Candidate victim row " << it->rows_as_str() << " is hammerable" << NORMAL_TXT << std::endl;

            *it = adjust_wrs(*it, row_layout);
        }
    }
}

void get_row_groups_by_index(boost::filesystem::ifstream& f_row_groups, vector<WeakRowSet>& row_groups, const vector<uint>& ind_weak_rows, const std::string& row_layout) {
    vector<WeakRowSet> all_weaks;
    all_weaks.reserve(100);
    parse_all_weaks(f_row_groups, all_weaks);

    for(uint ind : ind_weak_rows) {
        if(all_weaks.size() <= ind) {
            std::cerr << RED_TXT << "ERROR: The weaks rows file does not contain a sufficient number of hammerable weak rows" << NORMAL_TXT << std::endl;
            std::cerr << RED_TXT << "Needed the weak row at index: " << ind << " but the file contains " << all_weaks.size() << " row_groups" << NORMAL_TXT << std::endl;
            exit(-1);
        }

        WeakRowSet adjusted_wrs = adjust_wrs(all_weaks[ind], row_layout);
        row_groups.push_back(adjusted_wrs);
    }
}

bool check_dummy_vs_rg_collision(const std::vector<uint>& dummy_aggrs, const std::vector<WeakRowSet>& vec_wrs) {
    for(uint dummy : dummy_aggrs) {
        for(const auto& wrs : vec_wrs) {
            for(const auto& weak_row : wrs.row_group) {
                if(dummy == weak_row.row_id)
                    return true;
            }
        }
    }

    return false;
}

void adjust_hammers_per_ref(std::vector<uint>& hammers_per_round, const uint num_aggrs_in_wrs, const bool hammer_rgs_individually, const bool skip_hammering_aggr,
                            const uint num_wrs, const uint total_aggrs, const std::vector<uint>& dummy_aggr_ids, const uint dummy_hammers_per_round,
                            const bool hammer_dummies_first) {
    if(hammers_per_round.size() == 1) { // if a single hammers_per_round is specified, hammer all aggressors in a WRS the same amount
        while(hammers_per_round.size() != num_aggrs_in_wrs)
            hammers_per_round.push_back(hammers_per_round[0]);
    }

    if (!hammer_rgs_individually) { // use the same hammer counts for all WRSs
        uint num_specified_hammers = hammers_per_round.size();

        assert(num_specified_hammers == num_aggrs_in_wrs && 
            "ERROR: --hammers_per_round must specify exactly one hammer count for each aggressor in one WRS unless --hammer_rgs_individually is used to specify hammer counts for all WRS");

        for(uint i = 1; i < num_wrs; i++) {
            for(uint j = 0; j < num_specified_hammers; j++)
                hammers_per_round.push_back(hammers_per_round[j]);
        }
    }

    if(hammers_per_round.size() != total_aggrs) {
        std::cerr << RED_TXT << "ERROR: " << hammers_per_round.size() << " hammers per ref specified for " << 
            total_aggrs << " total aggressors" << NORMAL_TXT << std::endl;
        exit(-1);
    }

    for(auto& dummy_id : dummy_aggr_ids) {
        if(hammer_dummies_first)
            hammers_per_round.insert(hammers_per_round.begin(), dummy_hammers_per_round); // insert at front
        else
            hammers_per_round.push_back(dummy_hammers_per_round);
    }
}

int main(int argc, char** argv)
{
    /* Program options */
    std::string out_filename = "./out.txt";
    uint num_row_groups = 1;
    std::string row_scout_file = "";
    std::string row_layout = "RAR";
    std::vector<uint> hammers_per_round;
    std::vector<uint> hammers_before_wait;
    float init_to_hammerbw_delay = 0.0f;
    bool init_aggrs_first = false;
    bool first_it_aggr_init_and_hammer = false;
    bool first_it_dummy_hammer = false;
    bool hammer_rgs_individually = false;
    uint num_bank0_hammers = 0;
    uint num_pre_init_bank0_hammers = 0;
    uint pre_init_nops = 0;
    uint num_dummy_aggressors = 0;
    int dummy_aggrs_bank = -1;
    uint dummy_hammers_per_round = 1;
    uint dummy_ids_offset = 0;
    vector<uint> arg_dummy_aggr_ids;
    bool hammer_dummies_first = false;
    bool hammer_dummies_independently = false;
    bool cascaded_hammer = false;
    uint num_rounds = 0;
    uint num_refs_per_round = 1;
    uint pre_ref_delay = 0;
    uint num_iterations = 1;
    float hammer_cycle_time = 0.0f; // as nanosec
    uint hammer_duration = 0; // as DDR cycles (1.5ns)
    bool append_output = false;

    vector<uint> row_group_indices;
    bool skip_hammering_aggr = false;

    bool only_pick_rgs = false;

    uint refs_after_init = 0;
    uint num_dummy_after_init = 0;
    bool refs_after_init_no_dummy_hammer = false;

    bool init_only_victims = false;

    bool use_single_softmc_prog = false;
    bool location_out = false;

    uint arg_log_phys_conv_scheme = 0;

    // try{
    options_description desc("TRR Analyzer Options");
    desc.add_options()
        ("help,h", "Prints this usage statement.")
        ("out,o", value(&out_filename)->default_value(out_filename), "Specifies a path for the output file.")
        ("row_scout_file,f", value(&row_scout_file)->required(), "A file containing a list of row groups and their retentions times, i.e., the output of RowScout.")
        ("num_row_groups,w", value(&num_row_groups)->default_value(num_row_groups), "The number of row groups to work with. Row groups are parsed in order from the 'row_scout_file'.")
        ("row_layout", value(&row_layout)->default_value(row_layout), "Specifies how the aggressor rows should be positioned inside a row group. Allowed characters are 'R', 'A', 'U', and '-'. For example, 'RAR' places an aggressor row between two adjacent (victim) rows, as is single-sided RowHammer attacks. 'RARAR' places two aggressor rows to perform double-sided RowHammer attack. '-' specifies a row that is not to be hammered or checked for bit flips. 'U' specifies a (unified) row that will be both hammered and checked for bit flips.")
        ("row_group_indices", value<vector<uint>>(&row_group_indices)->multitoken(), "An optional argument used select which exact row groups in the --row_scout_file to use. When this argument is not provided, TRR Analyzer selects --num_row_groups from the file in order.")
        ("num_rounds", value(&num_rounds)->default_value(num_rounds), "Specifies the number of (hammer + refresh) rounds that the experiment should perform.")
        ("num_iterations", value(&num_iterations)->default_value(num_iterations), "Defines how many times the sequence of {aggr/victim initialization, hammer+ref rounds, reading back and checking for bit flips} should be performed.")

        // aggressor row related args
        ("hammers_per_round", value<vector<uint>>(&hammers_per_round)->multitoken(), "Specifies how many times each of the aggressors in --row_layout will be hammered in a round. You must enter multiple values, one for each aggressor.")
        ("cascaded_hammer", bool_switch(&cascaded_hammer), "When specified, the aggressor and dummy rows are hammered in non-interleaved manner, i.e., one row is hammered --hammers_per_round times and then the next row is hammered. Otherwise, the aggressor and dummy rows get activated one after another --hammers_per_round times.")
        ("hammers_before_wait", value<vector<uint>>(&hammers_before_wait)->multitoken(), "Similar to --hammers_per_round but hammering happens right after data initialization before waiting for half of the retention time.")
        ("hammer_rgs_individually", bool_switch(&hammer_rgs_individually), "When specified, --hammers_per_round specifies hammers for each aggressor row for separately each row group. Otherwise, the same aggressor hammers are applied to all row groups.")
        ("skip_hammering_aggr", bool_switch(&skip_hammering_aggr), "When provided, the aggressor rows are not hammered but just used to pick locations for the dummy rows.")
        ("hammer_duration", value(&hammer_duration)->default_value(hammer_duration), "Specifies the number of additional cycles to wait in row active state while hammering (tRAS + hammer_duration). The default is 0, i.e., tRAS)")
        ("hammer_cycle_time", value(&hammer_cycle_time)->default_value(hammer_cycle_time), "Specifies the time interval between two consecutive activations (the default and the minimum is tRAS + tRP).")
        ("init_aggrs_first", bool_switch(&init_aggrs_first), "When specified, the aggressor rows are initialized with a data pattern before the victim rows.")
        ("first_it_aggr_init_and_hammer", bool_switch(&first_it_aggr_init_and_hammer), "When specified, the aggressor rows are initialized and hammered only during the first iteration.")
        ("init_only_victims", bool_switch(&init_only_victims), "When specified, only the victim rows are initialized at the beginning of an iteration but not the aggressors.")

        // refresh related args
        ("refs_per_round", value(&num_refs_per_round)->default_value(num_refs_per_round), "Specifies how many REF commands to issue at the end of a round, i.e., after hammering.")
        ("refs_after_init", value(&refs_after_init), "Specifies the number of REF commands to issue right after initializing data in DRAM rows.")
                                
        // dummy row related args
        ("num_dummy_aggrs", value(&num_dummy_aggressors)->default_value(num_dummy_aggressors), "Specifies the number of dummy aggressors to hammer in each round. The dummy row addresses are selected such that they are different and in safe distance from the actual aggressor rows.")
        ("dummy_aggrs_bank", value(&dummy_aggrs_bank)->default_value(dummy_aggrs_bank), "Specifies the bank address from which dummy rows should be selected. If not specified, TRR Analyzer picks dummy rows from the same bank as the row groups.")
        ("dummy_aggr_ids", value<vector<uint>>(&arg_dummy_aggr_ids)->multitoken(), "Specifies the exact dummy row addresses to hammer in each round instead of letting TRR Analyzer select the dummy rows.")
        ("dummy_hammers_per_round", value(&dummy_hammers_per_round)->default_value(dummy_hammers_per_round), "Specifies how many times each dummy row to hammer in each round.")
        ("dummy_ids_offset", value(&dummy_ids_offset)->default_value(dummy_ids_offset), "Specifies a value to offset every dummy row address. Useful when there is a need to pick different dummy rows in different runs of TRR Analyzer.")
        ("hammer_dummies_first", bool_switch(&hammer_dummies_first), "When specified, the dummy rows are hammered before hammering the actual aggressor rows.")
        ("hammer_dummies_independently", bool_switch(&hammer_dummies_independently), "When specified, the dummy rows are hammered after the aggressor rows to matter whether --cascaded is used or not. The dummy rows are simply treated as a separate group of rows to hammer after hammering the aggressor rows in interleaved or cascaded way.")
        ("num_dummy_after_init", value(&num_dummy_after_init)->default_value(num_dummy_after_init), "Specifies the number of dummy rows to hammer right after initializing the victim and aggressor rows. These dummy row hammers happen concurrently with --refs_after_init refreshes. Each dummy is hammered as much as possible based on the refresh interval and --refs_after_init.")
        ("refs_after_init_no_dummy_hammer", bool_switch(&refs_after_init_no_dummy_hammer), "When specified, after hammering dummy rows as specified by --num_dummy_after_init, TRR Analyzer also performs another set of refreshes but this time without hammering dummy rows.")
        ("first_it_dummy_hammer", bool_switch(&first_it_dummy_hammer), "When specified, the dummy rows are hammered only during the first iteration.")

        // other. args
        ("init_to_hammerbw_delay", value(&init_to_hammerbw_delay)->default_value(init_to_hammerbw_delay), "A float in range [0,1] that specifies the ratio of time to wait before performing --hammers_before_wait. The default value (0) means all the delay is inserted after performing --hammers_before_wait (if specified)")
        ("num_bank0_hammers", value(&num_bank0_hammers)->default_value(num_bank0_hammers), "Specifies how many times a row from bank 0 should be hammered after hammering the aggressor and dummy rows.")
        ("num_pre_init_bank0_hammers", value(&num_pre_init_bank0_hammers)->default_value(num_pre_init_bank0_hammers), "Specifies how many times a row from bank 0 should be hammered before initializing data in victim and aggressor rows.")
        ("pre_init_nops", value(&pre_init_nops)->default_value(pre_init_nops), "Specifies the number of NOPs (as FPGA cycles, i.e., 4 DRAM cycles) to be inserted before victim/aggressor data initialization.")
        ("pre_ref_delay", value(&pre_ref_delay)->default_value(pre_ref_delay), "Specifies the number of cycles to wait before performing REFs specified by --refs_per_round. Must be 8 or larger if not 0 for this arg to take an effect.")
        ("only_pick_rgs", bool_switch(&only_pick_rgs), "When specified, the test finds hammerable row groups rows in --row_scout_file, but it does not run the TRR analysis.")
        ("log_phys_scheme", value(&arg_log_phys_conv_scheme)->default_value(arg_log_phys_conv_scheme), "Specifies how to convert logical row IDs to physical row ids and the other way around. Pass 0 (default) for sequential mapping, 1 for the mapping scheme typically used in Samsung chips.")
        ("use_single_softmc_prog", bool_switch(&use_single_softmc_prog), "When specified, the entire experiment executes as a single SoftMC program. This is to prevent SoftMC maintenance operations to kick in between multiple SoftMC programs. However, using this option may result in a very large program that may exceed the instruction limit.")
        ("append", bool_switch(&append_output), "When specified, the output of TRR Analyzer is appended to the --out file. Otherwise the --out file is cleared.")
        ("location_out", bool_switch(&location_out), "When specified, the bit flip locations are written to the --out file.")
        ;


    variables_map vm;
    boost::program_options::store(parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
        cout << desc << endl;
        return 0;
    }

    notify(vm);

    if(row_group_indices.size() > 0)
        num_row_groups = row_group_indices.size();

    if(arg_dummy_aggr_ids.size() > 0) {
        num_dummy_aggressors = arg_dummy_aggr_ids.size();
    }

    if(row_layout == "") {
        auto dot_pos = row_scout_file.find_last_of(".");
        if (dot_pos != string::npos)
            row_layout = row_scout_file.substr(dot_pos + 1);
        else {
            std::cerr << RED_TXT << "ERROR: Could not find '.' in the provided --row_scout_file\n" << std::endl;
            exit(-5);
        }
    }

    if(!std::regex_match(row_layout, std::regex("^[RrAaUu-]+$"))) {
        std::cerr << RED_TXT << "ERROR: --row_layout should contain only 'R', 'A', 'U', and '-' characters. Provided: " << row_layout << NORMAL_TXT << std::endl;
        exit(-3);
    }

    if(out_filename != "") {
        path out_dir(out_filename);
        out_dir = out_dir.parent_path();
        if (!(exists(out_dir))) {
            if (!create_directory(out_dir)) {
                cerr << "Cannot create directory: " << out_dir << ". Exiting..." << endl;
                return -1;
            }
        }
    }

    boost::filesystem::ofstream out_file;
    if(out_filename != "") {
        if(append_output)
            out_file.open(out_filename, boost::filesystem::ofstream::app);
        else
            out_file.open(out_filename);
    } else {
        out_file.open("/dev/null");
    }
    
    SoftMCPlatform platform;
    int err;

    if((err = platform.init()) != SOFTMC_SUCCESS){
        cerr << "Could not initialize SoftMC Platform: " << err << endl;
        return err;
    }

    platform.reset_fpga();  
    platform.set_aref(false); // disable refresh

    assert(arg_log_phys_conv_scheme < uint(LogPhysRowIDScheme::MAX));
    logical_physical_conversion_scheme = (LogPhysRowIDScheme) arg_log_phys_conv_scheme;

    // init random data generator
    std::srand(0);
  
    bitset<512> bitset_int_mask(0xFFFFFFFF);

    auto t_prog_started = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed;
    bool check_time;


    vector<WeakRowSet> row_groups;
    vector<uint> picked_weak_indices;
    row_groups.reserve(num_row_groups);
    picked_weak_indices.reserve(num_row_groups);

    boost::filesystem::ifstream f_row_groups;
    boost::filesystem::path p_row_scout_file(row_scout_file);
    if(!boost::filesystem::exists(p_row_scout_file)) {
        std::cerr << RED_TXT << "ERROR: RowScout file not found: " << row_scout_file << NORMAL_TXT << std::endl;
        exit(-1);
    }
    f_row_groups.open(p_row_scout_file);
    
    if(row_group_indices.size() > 0) {
        get_row_groups_by_index(f_row_groups, row_groups, row_group_indices, row_layout);
    }
    else if (num_row_groups > 0) {
        pick_hammerable_row_groups_from_file(platform, f_row_groups, row_groups, num_row_groups, cascaded_hammer, row_layout);
    }
    
    f_row_groups.close();

    if(only_pick_rgs) { // write the picked weak row indices to the output file and exit
        for(auto& rg : row_groups)
            out_file << rg.index_in_file << " ";

        return 0;
    }

    if(dummy_aggrs_bank == -1)
        dummy_aggrs_bank = row_groups[0].bank_id;

    // 3) Pick dummy aggressors rows
    if((num_dummy_aggressors > 0) && arg_dummy_aggr_ids.size() == 0) {
        uint max_dummy_aggrs = num_dummy_aggressors;
        arg_dummy_aggr_ids.reserve(max_dummy_aggrs);
        pick_dummy_aggressors(arg_dummy_aggr_ids, dummy_aggrs_bank, max_dummy_aggrs, row_groups, dummy_ids_offset);

    } else if (arg_dummy_aggr_ids.size() > 0 && (dummy_aggrs_bank == row_groups[0].bank_id)) { // check whether the user provided dummy row ids collide with the aggressor row ids
        if(check_dummy_vs_rg_collision(arg_dummy_aggr_ids, row_groups)) {
            std::cerr << RED_TXT << "ERROR: The user provided dummy aggressor rows collide with victims/aggressor rows. Finishing the test!" << NORMAL_TXT << std::endl;
            return -2;
        }
    }

    std::cout << YELLOW_TXT << "Dummy rows (while hammering): " << std::endl;
    for(auto dummy : arg_dummy_aggr_ids) {
        std::cout << dummy << " ";
    }
    std::cout << std::endl;

    // pick dummy rows that are hammered right after initializing data while performing refresh operations
    std::vector<uint> after_init_dummies;
    if(num_dummy_after_init > 0) {
        std::vector<WeakRowSet> cur_rgs_and_dummies = row_groups;
        WeakRowSet cur_dummies;

        // this is to pick different dummy than those we picked to hammer while hammering the actual aggressor rows
        for (uint dummy_row_id : arg_dummy_aggr_ids)
            cur_dummies.row_group.push_back(WeakRow(dummy_row_id, std::vector<uint>()));
        cur_rgs_and_dummies.push_back(cur_dummies);

        pick_dummy_aggressors(after_init_dummies, dummy_aggrs_bank, num_dummy_after_init, cur_rgs_and_dummies, dummy_ids_offset);
    }

    // std::cout << "Picked the following after init dummies: ";
    // for (auto& dummy_id : after_init_dummies)
    //     std::cout << dummy_id << " ";
    // std::cout << std::endl;


    // 4) Perform TRR analysis for each position of the weak rows among the dummy rows

    // convert row groups to hammerable row set
    std::vector<HammerableRowSet> hrs;
    hrs.reserve(row_groups.size());

    uint total_victims = 0;
    uint total_aggrs = 0;
    for(auto& wrs : row_groups) {
        hrs.push_back(toHammerableRowSet(wrs, row_layout));
        total_victims += hrs.back().victim_ids.size();
        total_aggrs += hrs.back().aggr_ids.size();
    }

    auto aggr_hammers_per_ref = hammers_per_round;

    if(total_aggrs > 0)
        adjust_hammers_per_ref(hammers_per_round, hrs[0].aggr_ids.size(), hammer_rgs_individually, skip_hammering_aggr,
                            row_groups.size(), total_aggrs, arg_dummy_aggr_ids, dummy_hammers_per_round, hammer_dummies_first);

    if(hammers_before_wait.size() > 0)
        adjust_hammers_per_ref(hammers_before_wait, hrs[0].aggr_ids.size(), hammer_rgs_individually, skip_hammering_aggr,
                            row_groups.size(), total_aggrs, arg_dummy_aggr_ids, 0, hammer_dummies_first);

    vector<uint> total_bitflips(total_victims, 0);

    // Setting up a progress bar
    progresscpp::ProgressBar progress_bar(num_iterations, 70, '#', '-');


    std::cout << BLUE_TXT << "Num hammerable row sets: " << hrs.size() << NORMAL_TXT << std::endl;
    uint hr_ind = 0;
    uint hammers_ind = 0;
    for(auto hr : hrs) {
        std::cout << BLUE_TXT << "Hammerable row set " << hr_ind << NORMAL_TXT << std::endl;

        std::cout << BLUE_TXT << "Victims: ";
        for (auto vict_id : hrs[hr_ind].victim_ids) {
            std::cout << vict_id << ", ";
        }
        std::cout << NORMAL_TXT << std::endl;

        std::cout << BLUE_TXT << "Aggressors: ";
        for (auto aggr_id : hrs[hr_ind].aggr_ids) {
            std::cout << aggr_id << " (";
            std::cout << aggr_hammers_per_ref[hammers_ind++] << "), ";
        }
        std::cout << NORMAL_TXT << std::endl;

        std::cout << BLUE_TXT << "Unified Rows: ";
        for (auto uni_id : hrs[hr_ind].uni_ids) {
            std::cout << uni_id << " (";
            std::cout << aggr_hammers_per_ref[hammers_ind++] << "), ";
        }
        std::cout << NORMAL_TXT << std::endl;

        hr_ind++;
    }

    std::cout << BLUE_TXT << "tRAS: " << tras_cycles << " cycles" << NORMAL_TXT << std::endl;

    // printing experiment parameters
    out_file << "row_layout=" << row_layout << std::endl;
    out_file << "--- END OF HEADER ---" << std::endl;


    if(!use_single_softmc_prog) {
        for (uint i = 0; i < num_iterations; i++) {

            bool ignore_aggrs = first_it_aggr_init_and_hammer ? i != 0 : false;
            bool ignore_dummy_hammers = first_it_dummy_hammer ? i != 0 : false;
            bool verbose = (i == 0);
            auto loc_bitflips = analyzeTRR(platform, hrs, arg_dummy_aggr_ids, dummy_aggrs_bank, dummy_hammers_per_round, hammer_dummies_first, hammer_dummies_independently, cascaded_hammer, 
                                                    hammers_per_round, hammer_cycle_time, hammer_duration, num_rounds, skip_hammering_aggr, refs_after_init, after_init_dummies,
                                                    init_aggrs_first, ignore_aggrs, init_only_victims, ignore_dummy_hammers, first_it_aggr_init_and_hammer,
                                                    refs_after_init_no_dummy_hammer, num_refs_per_round, pre_ref_delay, hammers_before_wait, init_to_hammerbw_delay,
                                                    num_bank0_hammers, num_pre_init_bank0_hammers, pre_init_nops, false, 0, verbose);

            ++progress_bar;
            progress_bar.display();
            
            out_file << "Iteration " << i << " bitflips:" << std::endl;

            if(!skip_hammering_aggr) {
                uint bitflips_ind = 0;

                for(auto& hr : hrs) {
                    uint total_rows = hr.victim_ids.size() + hr.uni_ids.size();

                    std::vector<std::string> output_strs_vict, output_strs_uni;

                    for(uint vict : hr.victim_ids) {
                        // out_file << "Victim row " << vict << ": " << num_bitflips[it_vict] << std::endl;
                        string output_str_vict;
                        output_str_vict = "Victim row " + to_string(vict) + ": " + to_string(loc_bitflips[bitflips_ind].size());
                        if(location_out){
                            output_str_vict += ": ";
                            for(auto loc: loc_bitflips[bitflips_ind])
                                output_str_vict += to_string(loc) + ", ";
                        }
                        output_strs_vict.push_back(output_str_vict);
                        total_bitflips[bitflips_ind] += loc_bitflips[bitflips_ind].size();
                        bitflips_ind++;
                    }

                    for(uint uni : hr.uni_ids) {
                        // out_file << "Victim row(U) " << uni << ": " << num_bitflips[it_vict] << std::endl;
                        string output_str_uni;
                        output_str_uni = "Victim row(U) " + to_string(uni) + ": " + to_string(loc_bitflips[bitflips_ind].size());
                        if(location_out){
                            output_str_uni += ": ";
                            for(auto loc: loc_bitflips[bitflips_ind])
                                output_str_uni += to_string(loc) + ", ";
                        }
                        output_strs_uni.push_back(output_str_uni);
                        total_bitflips[bitflips_ind] += loc_bitflips[bitflips_ind].size();
                        bitflips_ind++;
                    }

                    // reordering rows based on their physical row IDs
                    uint uni_ind = 0;
                    for(uint vict_ind = 0; vict_ind < hr.victim_ids.size(); vict_ind++){
                        if(uni_ind != hr.uni_ids.size() && to_physical_row_id(hr.uni_ids[uni_ind]) < to_physical_row_id(hr.victim_ids[vict_ind])){
                            out_file << output_strs_uni[uni_ind++] << std::endl;
                            vict_ind--;
                        } else {
                            out_file << output_strs_vict[vict_ind] << std::endl;
                        }
                    }

                    for(; uni_ind < hr.uni_ids.size(); uni_ind++)
                        out_file << output_strs_uni[uni_ind] << std::endl;
                }


            }
        }
    } else {
        assert(!first_it_dummy_hammer && "ERROR: --first_it_dummy_hammer is not yet supported when running the experiments as a single SoftMC program.");
        // run the experiment as a single SoftMC program
        auto num_bitflips = analyzeTRR(platform, hrs, arg_dummy_aggr_ids, dummy_aggrs_bank, dummy_hammers_per_round, hammer_dummies_first, hammer_dummies_independently, cascaded_hammer, 
                                                    hammers_per_round, hammer_cycle_time, hammer_duration, num_rounds, skip_hammering_aggr, refs_after_init, after_init_dummies,
                                                    init_aggrs_first, false, init_only_victims, false, first_it_aggr_init_and_hammer,
                                                    refs_after_init_no_dummy_hammer, num_refs_per_round, pre_ref_delay, hammers_before_wait, init_to_hammerbw_delay,
                                                    num_bank0_hammers, num_pre_init_bank0_hammers, pre_init_nops, true, num_iterations, true);

        // num_bitflips contains nothing since we have not read data from the PCIe yet
        // receive PCIe data iteration by iteration and keep the out_file format the same

        ulong read_data_size = ROW_SIZE*total_victims;
        char* buf = new char[read_data_size];
        vector<uint> bitflips;

        
        for (uint i = 0; i < num_iterations; i++) {
            if(!skip_hammering_aggr) {
                platform.receiveData(buf, read_data_size);

                out_file << "Iteration " << i << " bitflips:" << std::endl;
                
                uint row_it = 0;
                for (auto& hr : hrs) {
                    for(uint vict_ind = 0; vict_ind < hr.victim_ids.size(); vict_ind++) {
                        bitflips.clear();
                        collect_bitflips(bitflips, buf + row_it*ROW_SIZE, hr.data_pattern, hr.vict_bitflip_locs[vict_ind]);
                        row_it++;


                        out_file << "Victim row " << hr.victim_ids[vict_ind] << ": " << bitflips.size();
                        if(location_out){
                            out_file << ": ";
                            for(auto loc: bitflips)
                                out_file << loc << ", ";
                        }
                        out_file << endl;
                        total_bitflips[row_it] += bitflips.size();

                        // if(bitflips.size() == 0)
                        //     std::cout << RED_TXT << "[Victim Row " << row_groups[i].row_id << "] Did not find any bitflips." << NORMAL_TXT << std::endl;
                        // else
                        //     std::cout << GREEN_TXT << "[Victim Row " << row_groups[i].row_id << "] Found " << bitflips.size() << " bitflip(s)." << NORMAL_TXT << std::endl;
                    }

                    auto aggr_data_pattern = hr.data_pattern;
                    aggr_data_pattern.flip();
                    for(uint uni_ind = 0; uni_ind < hr.uni_ids.size(); uni_ind++) {
                        bitflips.clear();
                        collect_bitflips(bitflips, buf + row_it*ROW_SIZE, aggr_data_pattern, hr.uni_bitflip_locs[uni_ind]);
                        row_it++;


                        out_file << "Victim row(U) " << hr.uni_ids[uni_ind] << ": " << bitflips.size();
                        if(location_out){
                            out_file << ": ";
                            for(auto loc: bitflips)
                                out_file << loc << ", ";
                        }
                        out_file << endl;
                        total_bitflips[row_it] += bitflips.size();
                    }
                }
                assert(row_it == total_victims);
            }

            ++progress_bar;
            progress_bar.display(); // the progress bar is probably not that useful here
        }

        delete[] buf;
    }

    if(!skip_hammering_aggr) {
        out_file << "Total bitflips:" << std::endl;
        uint it_vict = 0;
        for(auto& hr : hrs) {
            for(uint vict : hr.victim_ids) {
                out_file << "Victim row " << vict << ": " << total_bitflips[it_vict++] << std::endl;
            }

            for(uint uni : hr.uni_ids) {
                out_file << "Victim row(U) " << uni << ": " << total_bitflips[it_vict++] << std::endl;
            }
        }
    }

    progress_bar.done();


    std::cout << "The test has finished!" << endl;

    out_file.close();

    return 0;
}