#include "tools/softmc_utils.h"
#include "instruction.h"
#include "prog.h"
#include "platform.h"

#include "tools/perfect_hash.h"
#include "tools/ProgressBar.hpp"

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
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

// #define PRINT_SOFTMC_PROGS
int print_times = 1;

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
int NUM_BANKS = 8; // this is the total number of banks in the chip
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

// TRR Attacker Parameters
const uint DUMMY_ROW_DIST = 2; // the minimum row distance between two dummy rows

const uint default_data_patterns[] = {0x0, 0xFFFFFFFF, 0x00000000, 0x55555555, 0xAAAAAAAA, 0xAAAAAAAA, 0x55555555};

vector<uint32_t> reserved_regs{CASR, BASR, RASR};

// returns a vector of bit positions that experienced bitflips
void collect_bitflips(vector<uint>& bitflips, const char* read_data, const bitset<512>& input_data_pattern) {

    bitflips.clear();

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
        bitset<512> error_mask = read_data_bitset ^ input_data_pattern;

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

// last_row_id is inclusive
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
            uint unroll_factor = 16;
            for(uint i_unroll = 0; i_unroll < unroll_factor; i_unroll++){
                add_op_with_delay(prog, SMC_WRITE(reg_bank_addr, 0, reg_col_addr, 1, 0, 0), 0, 0);
            }
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

std::vector<PhysicalRowID> getRowIDsOfType(const std::string row_layout, const PhysicalRowID row_id_offset, const char row_type) {

    std::vector<PhysicalRowID> row_ids;

    for(uint i = 0; i < row_layout.size(); i++){

        if(row_layout[i] == row_type)
            row_ids.push_back(row_id_offset + i);
    }

    return row_ids;
}

void toLogicalRowIDs(std::vector<uint>& row_ids) {

    for(uint i = 0; i < row_ids.size(); i++){
        row_ids[i] = to_logical_row_id(row_ids[i]);
    }
}

void init_row_data(Program& prog, SoftMCRegAllocator& reg_alloc, const SMC_REG reg_bank_addr, const SMC_REG reg_num_cols, 
    const vector<LogicalRowID>& rows_to_init, const vector<bitset<512>>& data_patts) {

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

void perform_refresh(Program& prog, SoftMCRegAllocator& reg_alloc, const uint num_refs_per_cycle,
                        const uint pre_ref_delay, const bool fake_ref) {

    SMC_REG reg_num_refs_per_cycle = reg_alloc.allocate_SMC_REG();
    SMC_REG reg_it_refs_per_cycle = reg_alloc.allocate_SMC_REG();

    prog.add_inst(SMC_LI(num_refs_per_cycle, reg_num_refs_per_cycle));
    prog.add_inst(SMC_LI(0, reg_it_refs_per_cycle));

    if (pre_ref_delay >= 8) {
        prog.add_inst(SMC_SLEEP(std::ceil(pre_ref_delay/4.0f)));
    }

    std::string lbl_issue_per_cycle_refs = createSMCLabel("PER_CYCLE_REFS");
    prog.add_label(lbl_issue_per_cycle_refs);    
        add_op_with_delay(prog, fake_ref ? SMC_NOP() : SMC_REF(), 0, 0);
        add_op_with_delay(prog, SMC_SLEEP(ceil((trfc_cycles - 1 - 24 - 4)/4.0f)), 0, 0);
        add_op_with_delay(prog, SMC_ADDI(reg_it_refs_per_cycle, 1, reg_it_refs_per_cycle), 0, 0);
    prog.add_branch(prog.BR_TYPE::BL, reg_it_refs_per_cycle, reg_num_refs_per_cycle, lbl_issue_per_cycle_refs);

    reg_alloc.free_SMC_REG(reg_num_refs_per_cycle);
    reg_alloc.free_SMC_REG(reg_it_refs_per_cycle);
}

void issue_REFs(SoftMCPlatform& platform, const uint num_refs) {
    // cout << YELLOW_TXT << "Performing " << num_refs << " refs." << NORMAL_TXT << endl;
    if (num_refs == 0) // do nothing if no REFs are to be performed
        return;

    Program prog;
    SoftMCRegAllocator reg_alloc(NUM_SOFTMC_REGS, reserved_regs);

    add_op_with_delay(prog, SMC_PRE(5, 0, 1), 0, 0); // precharge all banks

    perform_refresh(prog, reg_alloc, num_refs, 0, 0);
    prog.add_inst(SMC_END());

    platform.execute(prog);
    #ifdef PRINT_SOFTMC_PROGS
    prog.pretty_print();
    #endif
}

void perform_hammers(Program& prog, SoftMCRegAllocator& reg_alloc, SMC_REG reg_bank_addr, const std::vector<LogicalRowID> rows_to_hammer, 
                        const std::vector<uint>& num_hammers, const uint num_aggressors, const uint aggr_bank_id,
                        const bool cascaded_hammer, const std::vector<uint> dummy_banks, const bool fake_hammer, const bool fake_dummy_hammer){

    SMC_REG reg_cur_hammers = reg_alloc.allocate_SMC_REG();
    SMC_REG reg_num_hammers = reg_alloc.allocate_SMC_REG();
    SMC_REG reg_row_addr = reg_alloc.allocate_SMC_REG();

    uint remaining_cycs = 0;

    if(!cascaded_hammer){

        // it is complicated to efficiently hammer rows different number of times while activating them one after another
        // We implement the following algorithm:
        // 1. If there is a non-zero element in hammer_per_ref, hammer the rows corresponding to those elements using the smallest non-zero hammer_per_ref value. If all hammer_per_ref elements are zero, exit
        // 2. decrement all non-zero elements of hammer_per_ref vector by the smallest value
        // 3. go back to 1

        auto hammers_per_ref = num_hammers;

        while (1) {

            auto min_non_zero = std::min_element(hammers_per_ref.begin(), hammers_per_ref.end(), 
                    [](const uint& a, const uint& b) {return ((a > 0) && (a < b)) || (b == 0);}
                );

            if (min_non_zero == hammers_per_ref.end() || *min_non_zero == 0) {
                break;
            }

            uint min_elem = *min_non_zero;

            // perform hammering
            prog.add_inst(SMC_LI(min_elem, reg_num_hammers));
            prog.add_inst(SMC_LI(0, reg_cur_hammers));
            std::string lbl_rh = createSMCLabel("ROWHAMMERING");
            prog.add_label(lbl_rh);
            for (int ind_row = 0; ind_row < rows_to_hammer.size(); ind_row++) {
                if(hammers_per_ref[ind_row] == 0) // do not anymore hammer a row that has 0 remaining hammers
                    continue;

                int row_id = rows_to_hammer[ind_row];
                prog.add_inst(SMC_LI(row_id, reg_row_addr));

                // take into account fake_dummy_hammer for the dummy rows. Dummy rows follow the aggressor rows in the rows_to_hammer vector
                bool n_fake_hammer = ind_row >= num_aggressors ? fake_hammer | fake_dummy_hammer : fake_hammer;

                if(ind_row >= num_aggressors) { // hammering a dummy aggressor row
                    // hammer the same dummy row in all dummy_banks
                    const uint RRD_CYCLES = 5; // the spec says 6.4ns when activating rows in the same bank group. 5*1.5ns = 7.5ns
                    for(auto dummy_bank_id : dummy_banks){
                        prog.add_inst(SMC_LI(dummy_bank_id, reg_bank_addr));
                        remaining_cycs = add_op_with_delay(prog, n_fake_hammer ? SMC_NOP() : SMC_ACT(reg_bank_addr, 0, reg_row_addr, 0), remaining_cycs, RRD_CYCLES - 1 - 4);
                    }

                    prog.add_inst(SMC_LI(aggr_bank_id, reg_bank_addr)); // reload the aggressors' bank ID
                    remaining_cycs = add_op_with_delay(prog, n_fake_hammer ? SMC_NOP() : SMC_PRE(reg_bank_addr, 0, 1), tras_cycles - 5, 0); // precharge all banks
                } else { // hammering an aggressor row
                    remaining_cycs = add_op_with_delay(prog, n_fake_hammer ? SMC_NOP() : SMC_ACT(reg_bank_addr, 0, reg_row_addr, 0), 0, tras_cycles - 1);
                    remaining_cycs = add_op_with_delay(prog, n_fake_hammer ? SMC_NOP() : SMC_PRE(reg_bank_addr, 0, 0), 0, trp_cycles - 5);
                }
            }

            prog.add_inst(SMC_ADDI(reg_cur_hammers, 1, reg_cur_hammers));
            prog.add_branch(Program::BR_TYPE::BL, reg_cur_hammers, reg_num_hammers, lbl_rh);


            // this subtracts min_elem from every non-zero element
            std::for_each(hammers_per_ref.begin(), hammers_per_ref.end(), [&](uint& a) {if (a > 0) a -= min_elem;});
        }

        
    } else { // cascaded_hammer == true
        for (int ind_row = 0; ind_row < rows_to_hammer.size(); ind_row++) {
            int row_id = rows_to_hammer[ind_row];

            if(num_hammers[ind_row] == 0) // do not hammer rows with 0 hammer count
                continue;

            // take into account fake_dummy_hammer for the dummy rows. Dummy rows follow the aggressor rows in the rows_to_hammer vector
            bool n_fake_hammer = ind_row >= num_aggressors ? fake_hammer | fake_dummy_hammer : fake_hammer;

            prog.add_inst(SMC_LI(row_id, reg_row_addr));
            prog.add_inst(SMC_LI(num_hammers[ind_row], reg_num_hammers));
            prog.add_inst(SMC_LI(0, reg_cur_hammers));

            string lbl_rh = createSMCLabel("ROWHAMMERING");
            prog.add_label(lbl_rh);

            prog.add_inst(SMC_ADDI(reg_cur_hammers, 0, reg_cur_hammers)); // Hasan: this is doing functionally nothing. 
                                                                          // I put it here to make the interleaved and cascaded hammering loops have the same latency. 
                                                                          // The reason is that these additional 4 cycles spent while the bank is precharged causes more bitflips (at least for Micron 'mic03').

            if(ind_row >= num_aggressors) { // hammering a dummy aggressor row
                // hammer the same dummy row in all dummy_banks
                const uint RRD_CYCLES = 5; // the spec says 6.4ns when activating rows in the same bank group. 5*1.5ns = 7.5ns
                for(auto dummy_bank_id : dummy_banks){
                    prog.add_inst(SMC_LI(dummy_bank_id, reg_bank_addr));
                    remaining_cycs = add_op_with_delay(prog, n_fake_hammer ? SMC_NOP() : SMC_ACT(reg_bank_addr, 0, reg_row_addr, 0), remaining_cycs, RRD_CYCLES - 1 - 4);
                }

                prog.add_inst(SMC_LI(aggr_bank_id, reg_bank_addr)); // reload the aggressors' bank ID
                remaining_cycs = add_op_with_delay(prog, n_fake_hammer ? SMC_NOP() : SMC_PRE(reg_bank_addr, 0, 1), tras_cycles - 5, 0); // precharge all banks

            } else { // hammering an aggressor row
                remaining_cycs = add_op_with_delay(prog, n_fake_hammer ? SMC_NOP() : SMC_ACT(reg_bank_addr, 0, reg_row_addr, 0), 0, tras_cycles - 1);
                remaining_cycs = add_op_with_delay(prog, n_fake_hammer ? SMC_NOP() : SMC_PRE(reg_bank_addr, 0, 0), 0, 0);
            }
            remaining_cycs = 0;
            prog.add_inst(SMC_ADDI(reg_cur_hammers, 1, reg_cur_hammers));
            prog.add_branch(Program::BR_TYPE::BL, reg_cur_hammers, reg_num_hammers, lbl_rh);
        }
    }

    reg_alloc.free_SMC_REG(reg_cur_hammers);
    reg_alloc.free_SMC_REG(reg_num_hammers);
    reg_alloc.free_SMC_REG(reg_row_addr);
}

uint scheduleHammersPerStep(std::vector<uint>& scheduled_rows_to_hammer, std::vector<uint>& scheduled_num_hammers, uint& scheduled_num_aggr,
        std::vector<uint>& rows_to_hammer, std::vector<uint>& aggr_hammers, uint& num_aggr,
        const uint hammer_bugdet, const uint cascaded_hammer) {

    scheduled_num_aggr = 0;
    uint scheduled_hammer_cnt = 0;
    uint remaining_hammer_cnt = hammer_bugdet;

    scheduled_rows_to_hammer.clear();
    scheduled_num_hammers.clear();

    if(cascaded_hammer) {
        while(remaining_hammer_cnt > 0) {
            if(rows_to_hammer.size() == 0)
                break; // break if no more rows to hammer although there is some hammer bugdet left

            scheduled_rows_to_hammer.push_back(rows_to_hammer[0]);

            if(num_aggr > 0)
                scheduled_num_aggr++;
            
            uint hammers_to_cur_row = aggr_hammers[0] > remaining_hammer_cnt ? remaining_hammer_cnt : aggr_hammers[0];
            scheduled_num_hammers.push_back(hammers_to_cur_row);

            remaining_hammer_cnt -= hammers_to_cur_row;
            aggr_hammers[0] -= hammers_to_cur_row;
            scheduled_hammer_cnt += hammers_to_cur_row;

            if (aggr_hammers[0] == 0){
                rows_to_hammer.erase(rows_to_hammer.begin());
                aggr_hammers.erase(aggr_hammers.begin());

                if(num_aggr > 0)
                    num_aggr--;
            }
        }
    } else {
        scheduled_rows_to_hammer = rows_to_hammer;
        scheduled_num_hammers = std::vector<uint>(scheduled_rows_to_hammer.size(), 0);
        scheduled_num_aggr = num_aggr;

        while(remaining_hammer_cnt > 0) {

            if (std::accumulate(aggr_hammers.begin(), aggr_hammers.end(), 0) == 0)
                break; // break if reached the target hammer count for each row

            // distribute remaining_hammer_cnt equally to all rows
            uint hammers_per_row = max(int(remaining_hammer_cnt/scheduled_rows_to_hammer.size()), 1);

            for (uint i = 0; i < aggr_hammers.size(); i++) {    
                uint cur_row_hammer_cnt = aggr_hammers[i] > hammers_per_row ? hammers_per_row : aggr_hammers[i];

                scheduled_num_hammers[i] += cur_row_hammer_cnt;
                aggr_hammers[i] -= cur_row_hammer_cnt;
                remaining_hammer_cnt -= cur_row_hammer_cnt;
                scheduled_hammer_cnt += cur_row_hammer_cnt;

                if (remaining_hammer_cnt == 0)
                    break;
            }
        }
    }

    return scheduled_hammer_cnt;
}

void hammer_aggressors(Program& prog, SoftMCRegAllocator& reg_alloc, const SMC_REG reg_bank_addr, const vector<uint>& aggressors,
                        const std::vector<uint>& hammers_per_aggressor, const uint aggr_bank_id, const std::vector<LogicalRowID> dummy_rows, const uint hammers_per_dummy, 
                        const bool hammer_dummies_independently, const bool hammer_dummies_before, const bool hammer_dummies_after, const std::vector<uint> dummy_banks,
                        const uint num_refs, const uint refs_per_loop, const bool trrref_sync, const bool cascaded_hammer_aggr, const bool cascaded_hammer_dummy, 
                        const bool fake_hammer = false, const bool fake_dummy_hammer = false, const bool fake_ref = false) {


    
    vector<uint> rows_to_hammer(aggressors);
    std::vector<uint> num_hammers(hammers_per_aggressor);

    if(!hammer_dummies_independently){
        rows_to_hammer.insert(rows_to_hammer.end(), dummy_rows.begin(), dummy_rows.end()); // add dummy row ids at the back

        for (uint i = 0; i < dummy_rows.size(); i++) // add dummy hammer count at the back, one hammer count per dummy
            num_hammers.push_back(hammers_per_dummy);
    }
    

    if(rows_to_hammer.size() < 1 && (!hammer_dummies_independently || dummy_rows.size() < 1))
        return; // nothing to hammer

    uint initial_free_regs = reg_alloc.num_free_regs();
    uint remaining_cycs = 0;

    SMC_REG reg_ref_it = reg_alloc.allocate_SMC_REG();
    SMC_REG reg_num_refs = reg_alloc.allocate_SMC_REG();

    prog.add_inst(SMC_LI(0, reg_ref_it));
    prog.add_inst(SMC_LI(num_refs, reg_num_refs));

    std::string lbl_hammer_loop = createSMCLabel("HAMMER_LOOP");
    prog.add_label(lbl_hammer_loop);

    if (trrref_sync) { 
        // when --trrref_sync is specified, instead of issuing all of the --refs_per_loop as a batch,
        // the REF commands are issued with tREFI intervals such that they are interleaved with the hammers after hammering for tREFI.

        uint total_aggr_hammers = std::accumulate(num_hammers.begin(), num_hammers.end(), 0);
        uint total_dummy_only_hammers = 0;
        if(hammer_dummies_independently)
            total_dummy_only_hammers = hammers_per_dummy*dummy_rows.size();

        uint dummy_hammer_phases = (hammer_dummies_after & hammer_dummies_before) ? 2 : 
                                    (hammer_dummies_after | hammer_dummies_before) ? 1 : 0;
        uint total_hammer_bugdet = total_aggr_hammers + total_dummy_only_hammers*dummy_hammer_phases;

        // std::cout << YELLOW_TXT << "[DEBUG] Total hammer bugdet:  " << total_hammer_bugdet << NORMAL_TXT << std::endl;

        uint step_hammer_bugdet = total_hammer_bugdet/refs_per_loop;

        uint remaining_total_aggr_hammers = total_aggr_hammers;

        auto vec_remaining_rows_to_hammer = rows_to_hammer;
        auto vec_remaining_aggr_hammers = num_hammers;
        uint remaining_num_aggr = aggressors.size();

        auto vec_remaining_dummies_to_hammer = dummy_rows;
        std::vector<uint> vec_remaining_dummy_hammers_before, vec_remaining_dummy_hammers_after;

        if(hammer_dummies_independently)
            for (uint i = 0; i < dummy_rows.size(); i++) { // add dummy hammer counts, one hammer count per dummy
                vec_remaining_dummy_hammers_before.push_back(hammers_per_dummy);
                vec_remaining_dummy_hammers_after.push_back(hammers_per_dummy);
            }

        for (uint i = 0; i < refs_per_loop; i++) {
            uint remaining_step_hammers = step_hammer_bugdet;
            if(i == (refs_per_loop - 1)) {
                // Performing step_hammer_bugdet*refs_per_loop times may result in fewer total hammers due to integer rounding when dividing total_hammer_bugdet/refs_per_loop above
                remaining_step_hammers = total_hammer_bugdet - (step_hammer_bugdet*(refs_per_loop - 1));
            }

            // determine what to hammer in this step
            std::vector<uint> step_row_to_hammer;
            std::vector<uint> step_num_hammers;
            uint step_num_aggr = 0;

            // std::cout << YELLOW_TXT << "[DEBUG] Ref step: " << i << NORMAL_TXT << std::endl;
            // std::cout << YELLOW_TXT << "[DEBUG] Performing " << remaining_step_hammers << " hammers in this step" << NORMAL_TXT << std::endl;

            if(remaining_step_hammers > 0 && hammer_dummies_independently && hammer_dummies_before){
                // std::cout << YELLOW_TXT << "[DEBUG] Scheduling only dummy hammers via --hammer_dummies_independently..." << NORMAL_TXT << std::endl;
                remaining_num_aggr = 0;
                remaining_step_hammers -= scheduleHammersPerStep(step_row_to_hammer, step_num_hammers, step_num_aggr, 
                    vec_remaining_dummies_to_hammer, vec_remaining_dummy_hammers_before, remaining_num_aggr, remaining_step_hammers, cascaded_hammer_dummy);

                // std::cout << YELLOW_TXT << "[DEBUG] Hammering only dummies via --hammer_dummies_independently" << NORMAL_TXT << std::endl;
                // std::cout << YELLOW_TXT << "[DEBUG] Scheduled row hammers:" << NORMAL_TXT << std::endl;

                // for(uint row_ind = 0; row_ind < step_row_to_hammer.size(); row_ind++){
                //     std::cout << YELLOW_TXT << "[DEBUG] RowID (hammers):" << step_row_to_hammer[row_ind] << " (" << step_num_hammers[row_ind] << ")" << NORMAL_TXT << std::endl;
                // }

                perform_hammers(prog, reg_alloc, reg_bank_addr, step_row_to_hammer, step_num_hammers, step_num_aggr, aggr_bank_id, cascaded_hammer_dummy, dummy_banks, fake_hammer, fake_dummy_hammer);
            }

            if (remaining_total_aggr_hammers > 0) {
                // std::cout << YELLOW_TXT << "[DEBUG] Scheduling aggressor and dummy hammers..." << NORMAL_TXT << std::endl;

                uint scheduled_hammers = scheduleHammersPerStep(step_row_to_hammer, step_num_hammers, step_num_aggr, 
                    vec_remaining_rows_to_hammer, vec_remaining_aggr_hammers, remaining_num_aggr, remaining_step_hammers, cascaded_hammer_aggr);

                remaining_step_hammers -= scheduled_hammers;
                remaining_total_aggr_hammers -= scheduled_hammers;

                // std::cout << YELLOW_TXT << "[DEBUG] Step num aggr: " << step_num_aggr << NORMAL_TXT << std::endl;
                // std::cout << YELLOW_TXT << "[DEBUG] Remaining num aggr: " << remaining_num_aggr << NORMAL_TXT << std::endl;
                // std::cout << YELLOW_TXT << "[DEBUG] Remaining step hammers: " << remaining_step_hammers << NORMAL_TXT << std::endl;
                // std::cout << YELLOW_TXT << "[DEBUG] Step hammer budget: " << step_hammer_bugdet << NORMAL_TXT << std::endl;
                // std::cout << YELLOW_TXT << "[DEBUG] Hammering aggressors and dummies" << NORMAL_TXT << std::endl;
                // std::cout << YELLOW_TXT << "[DEBUG] Scheduled row hammers:" << NORMAL_TXT << std::endl;

                // for(uint row_ind = 0; row_ind < step_row_to_hammer.size(); row_ind++){
                //     std::cout << YELLOW_TXT << "[DEBUG] RowID (hammers):" << step_row_to_hammer[row_ind] << " (" << step_num_hammers[row_ind] << ")" << NORMAL_TXT << std::endl;
                // }
                
                perform_hammers(prog, reg_alloc, reg_bank_addr, step_row_to_hammer, step_num_hammers, step_num_aggr, aggr_bank_id, cascaded_hammer_aggr, dummy_banks, fake_hammer, fake_dummy_hammer);
            }

            if(remaining_step_hammers > 0 && hammer_dummies_independently && hammer_dummies_after){
                // std::cout << YELLOW_TXT << "[DEBUG] Scheduling only dummy hammers via --hammer_dummies_independently..." << NORMAL_TXT << std::endl;
                remaining_num_aggr = 0;
                remaining_step_hammers -= scheduleHammersPerStep(step_row_to_hammer, step_num_hammers, step_num_aggr, 
                    vec_remaining_dummies_to_hammer, vec_remaining_dummy_hammers_after, remaining_num_aggr, remaining_step_hammers, cascaded_hammer_dummy);

                // std::cout << YELLOW_TXT << "[DEBUG] Hammering only dummies via --hammer_dummies_independently" << NORMAL_TXT << std::endl;
                // std::cout << YELLOW_TXT << "[DEBUG] Scheduled row hammers:" << NORMAL_TXT << std::endl;

                // for(uint row_ind = 0; row_ind < step_row_to_hammer.size(); row_ind++){
                //     std::cout << YELLOW_TXT << "[DEBUG] RowID (hammers):" << step_row_to_hammer[row_ind] << " (" << step_num_hammers[row_ind] << ")" << NORMAL_TXT << std::endl;
                // }

                perform_hammers(prog, reg_alloc, reg_bank_addr, step_row_to_hammer, step_num_hammers, step_num_aggr, aggr_bank_id, cascaded_hammer_dummy, dummy_banks, fake_hammer, fake_dummy_hammer);
            }

            assert(remaining_step_hammers == 0);
            
            perform_refresh(prog, reg_alloc, 1, 0, fake_ref);
        }
    } else {
        if(hammer_dummies_independently && hammer_dummies_before){
            std::vector<uint> dummy_hammers;
            for (uint i = 0; i < dummy_rows.size(); i++) // add dummy hammer counts, one hammer count per dummy
                dummy_hammers.push_back(hammers_per_dummy);

            perform_hammers(prog, reg_alloc, reg_bank_addr, dummy_rows, dummy_hammers, 0, aggr_bank_id, cascaded_hammer_dummy, dummy_banks, fake_hammer, fake_dummy_hammer);
        }

        perform_hammers(prog, reg_alloc, reg_bank_addr, rows_to_hammer, num_hammers, aggressors.size(), aggr_bank_id, cascaded_hammer_aggr, dummy_banks, fake_hammer, fake_dummy_hammer);

        if(hammer_dummies_independently && hammer_dummies_after){
            std::vector<uint> dummy_hammers;
            for (uint i = 0; i < dummy_rows.size(); i++) // add dummy hammer counts, one hammer count per dummy
                dummy_hammers.push_back(hammers_per_dummy);

            perform_hammers(prog, reg_alloc, reg_bank_addr, dummy_rows, dummy_hammers, 0, aggr_bank_id, cascaded_hammer_dummy, dummy_banks, fake_hammer, fake_dummy_hammer);
        }

        perform_refresh(prog, reg_alloc, refs_per_loop, 0, fake_ref);
    }

    prog.add_inst(SMC_ADDI(reg_ref_it, 1, reg_ref_it));
    prog.add_branch(prog.BR_TYPE::BL, reg_ref_it, reg_num_refs, lbl_hammer_loop);

    reg_alloc.free_SMC_REG(reg_num_refs);
    reg_alloc.free_SMC_REG(reg_ref_it);

    assert(reg_alloc.num_free_regs() == initial_free_regs);
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


// performing REF at nominal rate, i.e., a REF cmd is issued once every 7.8us
// dummy rows are hammered between the REF cmds
void hammer_dummies(Program& prog, SoftMCRegAllocator& reg_alloc, const uint bank_id, 
    const vector<LogicalRowID>& dummy_aggrs, const uint num_refs) {

    
    SMC_REG reg_num_refs = reg_alloc.allocate_SMC_REG();
    SMC_REG reg_issued_refs = reg_alloc.allocate_SMC_REG();

    SMC_REG reg_row_id = reg_alloc.allocate_SMC_REG();
    SMC_REG reg_bank_id = reg_alloc.allocate_SMC_REG();
    prog.add_inst(SMC_LI(bank_id, reg_bank_id));
  
    prog.add_inst(SMC_LI(num_refs, reg_num_refs));
    prog.add_inst(SMC_LI(0, reg_issued_refs));

    int remaining_cycs = 0;

    uint num_dummies = dummy_aggrs.size();
    uint cycs_hammer_dummies_once = (tras_cycles + trp_cycles)*num_dummies + 28;
    uint hammers_per_ref = floor((trefi_cycles - trfc_cycles)/cycs_hammer_dummies_once);

    if(remaining_cycs > 0) { 
        remaining_cycs = add_op_with_delay(prog, SMC_NOP(), remaining_cycs, 0);
    }

    std::string lbl_issue_refs = createSMCLabel("ISSUE_REFS");
    prog.add_label(lbl_issue_refs);

    add_op_with_delay(prog, SMC_REF(), 0, trfc_cycles - 1);

        // hammer the dummy rows
        SMC_REG reg_hammer_it = reg_alloc.allocate_SMC_REG();
        SMC_REG reg_hammers_per_ref = reg_alloc.allocate_SMC_REG();
        prog.add_inst(SMC_LI(0, reg_hammer_it));
        prog.add_inst(SMC_LI(hammers_per_ref, reg_hammers_per_ref));

        std::string lbl_hammer = createSMCLabel("AFTER_INIT_DUMMY_HAMMERING");
        prog.add_label(lbl_hammer);
            remaining_cycs = 0;
            for (uint dummy_row_id : dummy_aggrs) {
                prog.add_inst(SMC_LI(dummy_row_id, reg_row_id));

                remaining_cycs = add_op_with_delay(prog, SMC_ACT(reg_bank_id, 0, reg_row_id, 0), remaining_cycs, tras_cycles - 1);
                remaining_cycs = add_op_with_delay(prog, SMC_PRE(reg_bank_id, 0, 0), remaining_cycs, trp_cycles - 1);
            }

            add_op_with_delay(prog, SMC_ADDI(reg_hammer_it, 1, reg_hammer_it), 0, 0);
        prog.add_branch(prog.BR_TYPE::BL, reg_hammer_it, reg_hammers_per_ref, lbl_hammer);


        prog.add_inst(SMC_ADDI(reg_issued_refs, 1, reg_issued_refs));
    prog.add_branch(prog.BR_TYPE::BL, reg_issued_refs, reg_num_refs, lbl_issue_refs);

    reg_alloc.free_SMC_REG(reg_num_refs);
    reg_alloc.free_SMC_REG(reg_issued_refs);
    reg_alloc.free_SMC_REG(reg_bank_id);
    reg_alloc.free_SMC_REG(reg_row_id);
    reg_alloc.free_SMC_REG(reg_hammer_it);
    reg_alloc.free_SMC_REG(reg_hammers_per_ref);    
}


void performHammer(SoftMCPlatform& platform, const uint target_bank, const PhysicalRowID anchor_row, const std::string row_layout, 
                    const std::vector<uint> num_hammers, const std::vector<LogicalRowID> dummy_rows, const uint hammers_per_dummy, 
                    const bool hammer_dummies_independently, const bool hammer_dummies_before, const bool hammer_dummies_after, const std::vector<uint> dummy_banks,
                    const uint num_refs, const uint refs_per_loop, const bool trrref_sync, const bool cascaded_hammer_aggr, const bool cascaded_hammer_dummy, 
                    const bitset<512>& victim_data, const bitset<512>& aggr_data, 
                    const bool fake_hammer, const bool fake_dummy_hammer, const bool fake_ref) {


    //get victims and aggressor IDs
    auto aggr_ids = getRowIDsOfType(row_layout, anchor_row, 'A');
    auto victim_ids = getRowIDsOfType(row_layout, anchor_row, 'V');

    // in-place conversion from physical row ID to logical row ID 
    toLogicalRowIDs(aggr_ids);
    toLogicalRowIDs(victim_ids);

    // std::cout << "Logical Victim IDs: ";
    // for (auto v : victim_ids)
    //     std::cout << v << " ";
    // std::cout << std::endl;

    // std::cout << "Logical Aggr IDs: ";
    // for (auto a : aggr_ids)
    //     std::cout << a << " ";
    // std::cout << std::endl;

    Program prog;
    SoftMCRegAllocator reg_alloc(NUM_SOFTMC_REGS, reserved_regs);

    SMC_REG reg_bank_addr = reg_alloc.allocate_SMC_REG();
    SMC_REG reg_num_cols = reg_alloc.allocate_SMC_REG();
    prog.add_inst(SMC_LI(target_bank, reg_bank_addr));
    prog.add_inst(SMC_LI(NUM_COLS_PER_ROW*8, reg_num_cols));

    add_op_with_delay(prog, SMC_PRE(reg_bank_addr, 0, 1), 0, 0); // precharge all banks

    std::vector<LogicalRowID> rows_to_init;
    rows_to_init.reserve(victim_ids.size() + aggr_ids.size());
    rows_to_init.insert(rows_to_init.end(), victim_ids.begin(), victim_ids.end());
    rows_to_init.insert(rows_to_init.end(), aggr_ids.begin(), aggr_ids.end());

    std::vector<bitset<512>> data_patts;
    data_patts.reserve(victim_ids.size() + aggr_ids.size());
    data_patts.insert(data_patts.end(), victim_ids.size(), victim_data);
    data_patts.insert(data_patts.end(), aggr_ids.size(), aggr_data);

    init_row_data(prog, reg_alloc, reg_bank_addr, reg_num_cols, rows_to_init, data_patts);

    // hammer the dummy rows while refreshing to kick out victims/aggressors from the counter table
    // const uint AFTER_INIT_DUMMY_HAMMER_REFS = 3758;
    // hammer_dummies(prog, reg_alloc, target_bank, dummy_rows, AFTER_INIT_DUMMY_HAMMER_REFS);

    // std::vector<LogicalRowID> tmp_dummy_rows;
    // for(uint i = 0; i < 64; i++)
    //     tmp_dummy_rows.push_back(NUM_ROWS/2 + i*2);
    // hammer_dummies(prog, reg_alloc, target_bank, tmp_dummy_rows, AFTER_INIT_DUMMY_HAMMER_REFS);


    // perform the actual hammers
    hammer_aggressors(prog, reg_alloc, reg_bank_addr, aggr_ids, num_hammers, target_bank, dummy_rows, hammers_per_dummy, 
        hammer_dummies_independently, hammer_dummies_before, hammer_dummies_after, dummy_banks, num_refs, refs_per_loop, 
        trrref_sync, cascaded_hammer_aggr, cascaded_hammer_dummy, fake_hammer, fake_dummy_hammer, fake_ref);

    // issue DRAM reads to read back the victim data
    read_row_data(prog, reg_alloc, reg_bank_addr, reg_num_cols, victim_ids);

    reg_alloc.free_SMC_REG(reg_bank_addr);
    reg_alloc.free_SMC_REG(reg_num_cols);

    prog.add_inst(SMC_END());

    platform.execute(prog);
    #ifdef PRINT_SOFTMC_PROGS
    if(print_times > 0) {
        prog.pretty_print();
        print_times--;
    }
    #endif
}

//last_row_id is inclusive
void readBankRegion(Program& prog, SoftMCRegAllocator& reg_alloc, const SMC_REG reg_bank_addr, const SMC_REG reg_num_cols,
                    const uint first_row_id, const uint last_row_id) {

    int remaining_cycs = 0;

    SMC_REG reg_row_addr = reg_alloc.allocate_SMC_REG();
    SMC_REG reg_col_addr = reg_alloc.allocate_SMC_REG();
    SMC_REG reg_last_row_addr = reg_alloc.allocate_SMC_REG();

    prog.add_inst(SMC_LI(first_row_id, reg_row_addr));
    prog.add_inst(SMC_LI(last_row_id + 1, reg_last_row_addr));

    prog.add_inst(SMC_LI(1, RASR)); // Load 1 into RASR
    prog.add_inst(SMC_LI(8, CASR)); // Load 8 into CASR since each READ reads 8 columns
    
    std::string lbl_read_row = createSMCLabel("READ_ROW");
    prog.add_label(lbl_read_row);
        remaining_cycs = add_op_with_delay(prog, SMC_ACT(reg_bank_addr, 0, reg_row_addr, 1), remaining_cycs, trcd_cycles - 1);

        prog.add_inst(SMC_LI(0, reg_col_addr));
        std::string lbl_read_cols = createSMCLabel("READ_COLS_OF_A_ROW");
        prog.add_label(lbl_read_cols);
            uint unroll_factor = 16;
            for(uint i_unroll = 0; i_unroll < unroll_factor; i_unroll++){
                remaining_cycs = add_op_with_delay(prog, SMC_READ(reg_bank_addr, 0, reg_col_addr, 1, 0, 0), remaining_cycs, 3);
            }
            remaining_cycs = 0;
        prog.add_branch(prog.BR_TYPE::BL, reg_col_addr, reg_num_cols, lbl_read_cols);

        remaining_cycs = add_op_with_delay(prog, SMC_PRE(reg_bank_addr, 0, 0), remaining_cycs, trp_cycles - 1);
    prog.add_branch(prog.BR_TYPE::BL, reg_row_addr, reg_last_row_addr, lbl_read_row);


    reg_alloc.free_SMC_REG(reg_row_addr);
    reg_alloc.free_SMC_REG(reg_col_addr);
    reg_alloc.free_SMC_REG(reg_last_row_addr);
}

void pick_dummy_aggressors(vector<uint>& dummy_aggrs, const uint num_dummies, const LogicalRowID region_start) {

    if(region_start + num_dummies*DUMMY_ROW_DIST >= NUM_ROWS)
        std::cerr << "ERROR: cannot allocate " << num_dummies << " dummy rows starting from row ID " 
            << region_start << ". The minimum distance between two dummy rows is " << DUMMY_ROW_DIST << std::endl;


    dummy_aggrs.clear();

    for(uint i = 0; i < num_dummies; i++) {
        dummy_aggrs.push_back(region_start + i*DUMMY_ROW_DIST);
    }
}

// returns a vector of bitflips per the specified data chunk size
std::vector<uint32_t> count_bitflips_in_chunks(const std::vector<uint32_t>& bitflips, const uint bitflip_counting_granularity){

    if(bitflip_counting_granularity == 0)
        return std::vector<uint32_t>{(uint32_t)bitflips.size()};

    std::vector<uint32_t> cl_bitflips(ROW_SIZE/bitflip_counting_granularity, 0);

    for (auto bitflip_loc : bitflips){
        uint cl_ind = bitflip_loc/(bitflip_counting_granularity << 3); // bitflip_loc shows the bit location of the bit flip. This is why we multiply bitflip_counting_granularity by 8
        cl_bitflips[cl_ind]++;
    }

    assert(bitflips.size() == std::accumulate(cl_bitflips.begin(), cl_bitflips.end(), 0));

    return cl_bitflips;
}

// returns true if 'hammer_count' causes bitflips
bool testHammerCount(SoftMCPlatform& platform, const uint bank_id, const uint row_id, const uint hammer_count,
                        const bitset<512>& victim_data, const bitset<512>& aggr_data) {

    Program prog;
    SoftMCRegAllocator reg_alloc(NUM_SOFTMC_REGS, reserved_regs);

    SMC_REG reg_bank_addr = reg_alloc.allocate_SMC_REG();
    SMC_REG reg_num_cols = reg_alloc.allocate_SMC_REG();
    prog.add_inst(SMC_LI(bank_id, reg_bank_addr));
    prog.add_inst(SMC_LI(NUM_COLS_PER_ROW*8, reg_num_cols));

    add_op_with_delay(prog, SMC_PRE(reg_bank_addr, 0, 1), 0, 0); // precharge all banks

    // std::vector<uint> aggr_ids = {row_id - 1, row_id + 1};
    std::vector<uint> aggr_ids = {row_id + 1};
    std::vector<uint> aggr_hammers = std::vector<uint>(aggr_ids.size(), hammer_count);
    std::vector<bitset<512>> aggr_data_pattern = std::vector<bitset<512>>(aggr_ids.size(), aggr_data);

    init_row_data(prog, reg_alloc, reg_bank_addr, reg_num_cols, {row_id}, {victim_data});
    init_row_data(prog, reg_alloc, reg_bank_addr, reg_num_cols, aggr_ids, aggr_data_pattern);
    perform_hammers(prog, reg_alloc, reg_bank_addr, aggr_ids, aggr_hammers, 1, bank_id, false, {1}, false, false);
    read_row_data(prog, reg_alloc, reg_bank_addr, reg_num_cols, {row_id});

    prog.add_inst(SMC_END());

    platform.execute(prog);
    // prog.pretty_print();

    reg_alloc.free_SMC_REG(reg_bank_addr);
    reg_alloc.free_SMC_REG(reg_num_cols);

    char buf[ROW_SIZE*2];
    platform.receiveData(buf, ROW_SIZE);

    std::vector<uint> bitflips;
    collect_bitflips(bitflips, buf, victim_data);

    return bitflips.size() > 0;
}

uint findHCFirst(SoftMCPlatform& platform, const uint bank_id, const uint row_id, const bitset<512>& victim_data, const bitset<512>& aggr_data){

    const uint INIT_HAMMER_COUNT = 5000;
    const uint MAX_HAMMER_COUNT = 1500000;
    const uint HAMMER_COUNT_PRECISION = 500; // keep refining until the distance between the latest two tested hammer counts is less than or equal to this value

    uint hc_low = INIT_HAMMER_COUNT;
    uint hc_high = MAX_HAMMER_COUNT;
    uint t_hc = INIT_HAMMER_COUNT;
    while (true) {

        // std::cout << YELLOW_TXT << "[DEBUG] Testing hammer count: " << t_hc << NORMAL_TXT << std::endl;

        bool causes_bitflips = testHammerCount(platform, bank_id, row_id, t_hc, victim_data, aggr_data);

        if (causes_bitflips) {
            
            if (t_hc == INIT_HAMMER_COUNT)
                return 0;

            hc_high = t_hc;
            t_hc = (hc_high + hc_low)/2;

            if (hc_high - t_hc <= HAMMER_COUNT_PRECISION)
                return hc_high;
        } else {

            if ((MAX_HAMMER_COUNT - t_hc) <= HAMMER_COUNT_PRECISION*2)
                return 0;

            hc_low = t_hc;
            t_hc = (hc_high + hc_low)/2;

            if ((t_hc - hc_low) <= HAMMER_COUNT_PRECISION)
                return hc_high;
        }
    }

    return 0; // not reachable
}

std::vector<uint> gatherRunsWithTRR(SoftMCPlatform& platform, const uint num_runs, 
                const uint bank_id, const uint row_id, const uint hammer_count, 
                const bitset<512>& victim_data, const bitset<512>& aggr_data) {

    std::vector<uint> num_bitflips_per_run;
    num_bitflips_per_run.reserve(num_runs);
    int target_cell = -1;

    // std::vector<uint> aggr_ids = {row_id - 1, row_id + 1}; // double-sided does not work well with Micron's TRR
    std::vector<uint> aggr_ids = {row_id + 1};
    std::vector<uint> aggr_hammers = std::vector<uint>(aggr_ids.size(), hammer_count/2);
    std::vector<bitset<512>> aggr_data_pattern = std::vector<bitset<512>>(aggr_ids.size(), aggr_data);


    progresscpp::ProgressBar progress_bar(num_runs, 70, '#', '-');

    char buf[ROW_SIZE*2];
    for (uint run_id = 0; run_id < num_runs; run_id++) {

        // Generating SoftMC program
        Program prog;
        SoftMCRegAllocator reg_alloc(NUM_SOFTMC_REGS, reserved_regs);

        SMC_REG reg_bank_addr = reg_alloc.allocate_SMC_REG();
        SMC_REG reg_num_cols = reg_alloc.allocate_SMC_REG();
        prog.add_inst(SMC_LI(bank_id, reg_bank_addr));
        prog.add_inst(SMC_LI(NUM_COLS_PER_ROW*8, reg_num_cols));

        add_op_with_delay(prog, SMC_PRE(reg_bank_addr, 0, 1), 0, 0); // precharge all banks

        uint aggr_id = row_id + 1;

        init_row_data(prog, reg_alloc, reg_bank_addr, reg_num_cols, {row_id}, {victim_data});
        init_row_data(prog, reg_alloc, reg_bank_addr, reg_num_cols, aggr_ids, aggr_data_pattern);

        perform_hammers(prog, reg_alloc, reg_bank_addr, aggr_ids, aggr_hammers, 1, bank_id, false, {1}, false, false);
        perform_refresh(prog, reg_alloc, 1, 0, false);
        perform_hammers(prog, reg_alloc, reg_bank_addr, aggr_ids, aggr_hammers, 1, bank_id, false, {1}, false, false);
        
        read_row_data(prog, reg_alloc, reg_bank_addr, reg_num_cols, {row_id});

        prog.add_inst(SMC_END());
        // END - generating SoftMC program

        platform.execute(prog);

        platform.receiveData(buf, ROW_SIZE);

        std::vector<uint> bitflips;
        collect_bitflips(bitflips, buf, victim_data);

        if(target_cell == -1 && bitflips.size() > 0){
            target_cell = bitflips[0];
        }

        if(find(bitflips.begin(), bitflips.end(), target_cell) != bitflips.end())
            num_bitflips_per_run.push_back(1);
        else
            num_bitflips_per_run.push_back(0);

        ++progress_bar;
        progress_bar.display();
    }
    progress_bar.done();

    return num_bitflips_per_run;
}


// this program finds which REF commands perform TRR and synchronizes with TRR REFs such that the last REF issued
// prior to the program finishes is a TRR REF
// It is expected to use this with Samsung and Hynix modules where every 9th and 4th REF command, respectively, is a TRR REF
// Method:
//  1. Pick an arbitrary row and find its HC_first value. Assume the HC_first is N.
//  2. Initialize the row and hammer it (N*1.2)/2 times
//  3. Issue a REF command
//  4. Hammer the row (N*1.2)/2 more times
//  5. Collect the bitflips. Iterations w/o bitflips would indicate the TRR REFs
//  6. Repeat steps 2-5 to gather a good sample size
uint syncTRRREF(SoftMCPlatform& platform, const uint input_data_victims, const uint input_data_aggressors) {

    const uint TARGET_BANK = 1;
    const uint TARGET_ROW = 1300;

    bitset<512> victims_data = setup_data_pattern(input_data_victims);
    bitset<512> aggrs_data = setup_data_pattern(input_data_aggressors);

    std::cout << BLUE_TXT << "Finding HC_first for bank: " << TARGET_BANK << ", row: " << TARGET_ROW << "..." << NORMAL_TXT << std::endl;
    const uint hc_first = findHCFirst(platform, TARGET_BANK, TARGET_ROW, victims_data, aggrs_data);

    if (hc_first == 0) {
        std::cout << RED_TXT << "ERROR: HC_first is either smaller or bigger than the hammer count search range. Consider updating the range or try another potentially weaker/stronger row" << hc_first << NORMAL_TXT << std::endl;
        return 0;
    }
    
    std::cout << YELLOW_TXT << "HC_first is " << hc_first << NORMAL_TXT << std::endl;

    const float HC_MULT = 1.2f;
    uint test_hammers = hc_first * HC_MULT;

    const uint NUM_RUNS = 1024;
    std::cout << BLUE_TXT << "Gathering data to find TRR REFs..." << NORMAL_TXT << std::endl;
    std:vector<uint> num_bitflips_per_run = gatherRunsWithTRR(platform, NUM_RUNS, TARGET_BANK, TARGET_ROW, test_hammers, victims_data, aggrs_data);

    // for (uint i = 0; i < NUM_RUNS; i++) {
    //     std::cout << "Run #" << i << ": " << num_bitflips_per_run[i] << std::endl;
    // }

    // build a histogram of distances between runs with 0 bitflips
    int last_no_bitflip_run = -1;
    std::vector<uint> hist_zero_bitflips = std::vector<uint>(NUM_RUNS, 0);
    for (uint i = 0; i < NUM_RUNS; i++) {
        if (last_no_bitflip_run == -1) {
            if (num_bitflips_per_run[i] == 0)
                last_no_bitflip_run = i;
            continue;
        }

        if (num_bitflips_per_run[i] == 0) {
            hist_zero_bitflips[i - last_no_bitflip_run]++;
            last_no_bitflip_run = i;
        }
    }

    if (last_no_bitflip_run == -1) {
        std::cout << RED_TXT << "ERROR: Could not find a run without bitflips." << NORMAL_TXT << std::endl;
        return 0;
    }

    int most_frequent_distance = 0; // = std::max_element(hist_zero_bitflips.begin(), hist_zero_bitflips.end()) - hist_zero_bitflips.begin();

    // instead of using the most frequent distance, implementing the following logic because Micron modules can take the other aggressor 
    // that does not cause refresh on the target victim row (because Micron's TRR performs targeted refresh on one of the neighbors of an aggressor row but not both)
    const float PASS_THRESHOLD = 0.15f;
    for (uint i = 1; i < hist_zero_bitflips.size(); i++) {
        if (hist_zero_bitflips[i] >= ((NUM_RUNS/i)*PASS_THRESHOLD)){
            most_frequent_distance = i;
            break;
        }
    }

    std::cout << YELLOW_TXT << "The distance between TRR REFs: " << most_frequent_distance << NORMAL_TXT << std::endl;

    uint num_last_regular_refs = NUM_RUNS - last_no_bitflip_run - 1;
    uint refs_to_sync = (most_frequent_distance - num_last_regular_refs) % most_frequent_distance;
    std::cout << BLUE_TXT << "Issuing " << refs_to_sync << " REFs so that the last REF is a TRR REF..." << NORMAL_TXT << std::endl;

    issue_REFs(platform, refs_to_sync);

    std::cout << BLUE_TXT << "Completing synchronizing with TRR REFs" << NORMAL_TXT << std::endl;

    return most_frequent_distance;
}

// last_row_id is inclusive
void hammerBank(SoftMCPlatform& platform, const uint target_bank, const PhysicalRowID first_row_id, const PhysicalRowID last_row_id, 
        const std::vector<uint>& num_hammers, const uint num_ref_loops, const uint refs_per_loop, const bool trrref_sync, 
        const uint num_dummies, const bool hammer_dummies_independently, const bool hammer_dummies_before, const bool hammer_dummies_after, const std::vector<uint>& dummy_banks,
        const bool cascaded_hammer_aggr, const bool cascaded_hammer_dummy, const bool fake_hammer, const bool fake_dummy_hammer, const bool fake_ref,
        const std::string& row_layout, const uint input_data_victims, const uint input_data_aggressors, 
        const uint bitflip_counting_granularity,
        boost::filesystem::ofstream& out_file){

    std::vector<uint> victims_pos;
    for (uint i = 0; i < row_layout.size(); i++) {
        if (row_layout[i] == 'V')
            victims_pos.push_back(i);
    }

    char buf[ROW_SIZE*victims_pos.size()*2];
    bitset<512> victims_data = setup_data_pattern(input_data_victims);
    bitset<512> aggrs_data = setup_data_pattern(input_data_aggressors);

    std::map<uint, std::vector<uint>> num_bitflips_data;

    std::cout << YELLOW_TXT << "There are " << victims_pos.size() << " victim rows" << NORMAL_TXT << std::endl;

    // Setting up a progress bar
    uint total_iterations = last_row_id - first_row_id + 1;
    progresscpp::ProgressBar progress_bar(total_iterations, 70, '#', '-');

    std::vector<uint32_t> bitflips;
    std::vector<uint> num_bitflips_per_victim;
    num_bitflips_per_victim.reserve(victims_pos.size());

    std::vector<LogicalRowID> dummy_rows;
    LogicalRowID dummy_row_region_start = 3*NUM_ROWS/4;
    pick_dummy_aggressors(dummy_rows, num_dummies, dummy_row_region_start);

    // calculate the maximum ammount of hammers_per_dummy based on remaining ACTs after hammering the aggressor rows
    // For example, assuming we can perform 150 ACTs between two REFs with nominal refresh rate,
        // hammering two aggressors 20 times each would leave floor(110/16) ACTs for each of the 16 dummy rows
    uint hammers_per_dummy = 0;

    const uint TREFI = 7800; // assume tREFI = 7800ns
    const uint TRC = 50; // assume activate-precharge complete in 50ns
    const uint MAX_ACTS_PER_REF = (TREFI/TRC)*refs_per_loop; 

    uint dummy_hammer_phases = (hammer_dummies_after & hammer_dummies_before) ? 2 : 
                                    (hammer_dummies_after | hammer_dummies_before) ? 1 : 0;
    
    uint total_aggressor_ACTs = std::accumulate(num_hammers.begin(), num_hammers.end(), 0);
    if(total_aggressor_ACTs >= MAX_ACTS_PER_REF){
        std::cerr << YELLOW_TXT << "Warning: Activating more than the time between two REF commands permit under default refresh" << NORMAL_TXT << std::endl;
    } else {
        if (dummy_hammer_phases == 0)
            hammers_per_dummy = 0;
        else
            hammers_per_dummy = std::floor((MAX_ACTS_PER_REF-total_aggressor_ACTs)/(float)(num_dummies*dummy_hammer_phases));
    }

    std::cout << BLUE_TXT << "Hammering the dummy rows as many times as possible given the ACT count for the aggressors." << NORMAL_TXT << std::endl;
    std::cout << BLUE_TXT << "Number of dummies: " << num_dummies << NORMAL_TXT << std::endl;
    std::cout << BLUE_TXT << "Hammers per dummy (per bank): " << hammers_per_dummy << NORMAL_TXT << std::endl;
    std::cout << BLUE_TXT << "Total dummy hammers: " << hammers_per_dummy*num_dummies*dummy_banks.size() << NORMAL_TXT << std::endl;

    for(PhysicalRowID anchor_row = first_row_id; anchor_row <= last_row_id; anchor_row++){

        // std::cout << "Physical Row ID offset: " << anchor_row << std::endl; // DEBUG

        if(to_logical_row_id(anchor_row) >= dummy_row_region_start/2) {
            dummy_row_region_start = 0;
            pick_dummy_aggressors(dummy_rows, num_dummies, dummy_row_region_start);
        }

        performHammer(platform, target_bank, anchor_row, row_layout, num_hammers, dummy_rows, hammers_per_dummy, hammer_dummies_independently, hammer_dummies_before, hammer_dummies_after, dummy_banks,
            num_ref_loops, refs_per_loop, trrref_sync, cascaded_hammer_aggr, cascaded_hammer_dummy, victims_data, aggrs_data, fake_hammer, fake_dummy_hammer, fake_ref);

        platform.receiveData(buf, ROW_SIZE*victims_pos.size());
        // std::cout << "Succesffuly received the row data!" << std::endl; // DEBUG

        ulong total_bitflips = 0;
        num_bitflips_per_victim.clear();

        for (uint i_victim = 0; i_victim < victims_pos.size(); i_victim++) {
            // check for bitflips
            collect_bitflips(bitflips, buf + i_victim*ROW_SIZE, victims_data);

            total_bitflips += bitflips.size();

            auto bitflips_per_cl = count_bitflips_in_chunks(bitflips, bitflip_counting_granularity);
            num_bitflips_per_victim.insert(num_bitflips_per_victim.end(), bitflips_per_cl.begin(), bitflips_per_cl.end());

        }

        if(total_bitflips > 0) {
            out_file << std::setw(5) << anchor_row << ": ";

            for(auto num_bf : num_bitflips_per_victim)
                out_file << num_bf << " ";

            out_file << std::endl;
        }

        ++progress_bar;
        progress_bar.display();
        
    }

    progress_bar.done();
}


int main(int argc, char** argv)
{
    /* Program options */
    string out_filename = "./out.txt";
    bool append_output = false;

    uint target_bank = 1;
    vector<int> row_range{-1, -1};
    std::vector<uint> num_hammers = {70, 70};
    uint num_ref_loops = 8192;
    uint refs_per_loop = 1;
    bool trrref_sync = false;
    uint num_dummy_rows = 0;
    bool hammer_dummies_independently = false;
    bool hammer_dummies_before = false;
    bool hammer_dummies_after = false;
    std::vector<uint> dummy_banks{1};
    bool cascaded_hammer_aggr = false;
    bool cascaded_hammer_dummy = false;
    bool fake_hammer = false;
    bool fake_dummy_hammer = false;
    bool fake_ref = false;
    bool shift_refs = false;
    std::string row_layout = "VAVAV";
    uint input_data_victims = 2;
    uint input_data_aggressors = 1;

    uint bitflip_counting_granularity = 0; // When 0, counts and outputs the bitflips for each row. Otherwise based on the byte granularity provided with this parameter. E.g., 8 would report bitflips in every 8-byte chunk in the tested memory region.

    uint arg_log_phys_conv_scheme = 0;

    // try{
    options_description desc("RowHammerAttacker Options");
    desc.add_options()
        ("help,h", "Prints this usage statement.")
        ("out,o", value(&out_filename)->default_value(out_filename), "Specifies a path for the output file.")
        ("bank,b", value(&target_bank)->default_value(target_bank), "Specifies the address of the bank to perform the RowHammer attack on.")
        ("range", value<vector<int>>(&row_range)->multitoken(), "Specifies a range of row addresses (start and end values are both inclusive) to perform the RowHammer attack on. The range is set to cover the entire bank when --range is not provided.")
        ("row_layout", value(&row_layout)->default_value(row_layout), "The layout of victim (V) and aggressor (A) rows. E.g., VAV, VVVAVVV.")
        
        ("hammers_per_ref_loop", value<vector<uint>>(&num_hammers)->multitoken(), "Number of ACTs (hammers) to perform per aggressor row before each REF command, i.e., per refresh loop.")
        ("cascaded_hammer_aggr", bool_switch(&cascaded_hammer_aggr)->default_value(cascaded_hammer_aggr), "When specified, the aggressor rows are hammered in non-interleaved manner, e.g., one row is hammered N times and then the next row is hammered. Otherwise, the aggressor rows get activated only once, one after another N times.")
        ("fake_hammer", bool_switch(&fake_hammer)->default_value(fake_hammer), "When specified, the aggressor rows are not actually activated but the program pretends to activate them by spending the activation time.")

        ("num_ref_loops", value(&num_ref_loops)->default_value(num_ref_loops), "Specifies the number refresh loops to perform.")
        ("refs_per_loop", value(&refs_per_loop)->default_value(refs_per_loop), "Specifies the number of REF commands to issue at the end of each refresh loop")
        ("trrref_sync", bool_switch(&trrref_sync)->default_value(trrref_sync), "When specified, the program finds the period of TRR-induced refresh operations and synchronizes with them. Overwrites --refs_per_loop with the found period.")
        ("shift_refs", bool_switch(&shift_refs)->default_value(shift_refs), "When specified, a single REF is issued at the beginning of the experiment.")
        ("fake_ref", bool_switch(&fake_ref)->default_value(fake_ref), "When specified, the REF commands are not actually issued but the program pretends to issue REFs by spending the same amount of time as when issuing REFs.")

        ("num_dummy_rows", value(&num_dummy_rows)->default_value(num_dummy_rows), "Specifies the number of dummy rows.")
        ("hammer_dummies_independently", bool_switch(&hammer_dummies_independently)->default_value(hammer_dummies_independently), "When specified, the the dummy rows are always hammered after being done hammering the actual aggressor rows, w/ or w/o --cascaded.")
        ("hammer_dummies_before", bool_switch(&hammer_dummies_before)->default_value(hammer_dummies_before), "When specified, the the dummy rows are hammered before the actual aggressor rows. Can be specified together with --hammer_dummies_after to hammer the dummies both before and after.")
        ("hammer_dummies_after", bool_switch(&hammer_dummies_after)->default_value(hammer_dummies_after), "When specified, the the dummy rows are hammered after the actual aggressor rows. Can be specified together with --hammer_dummies_before to hammer the dummies both before and after.")
        ("dummy_banks", value<vector<uint>>(&dummy_banks)->multitoken(), "Specifies the bank ID to simultaneously hammer dummy rows from. By default, dummies are selected from bank 1.")
        ("cascaded_hammer_dummy", bool_switch(&cascaded_hammer_dummy)->default_value(cascaded_hammer_dummy), "When specified, the dummy rows are hammered in non-interleaved manner, e.g., one row is hammered N times and then the next row is hammered. Otherwise, the dummy rows get activated only once, one after another N times.")
        ("fake_dummy_hammer", bool_switch(&fake_dummy_hammer)->default_value(fake_dummy_hammer), "When specified, the dummy aggressor rows are not actually activated but the program pretends to activate them by spending the activation time.")
        
        ("input_victims", value(&input_data_victims)->default_value(input_data_victims), "Input data pattern for victim rows. 0: random, 1: all ones, 2: all zeros, 3: colstripe (0101), 4: inverse colstripe (1010), 5: checkered (0101, 1010), 6: inverse checkered (1010, 0101)")
        ("input_aggressors", value(&input_data_aggressors)->default_value(input_data_aggressors), "Input data pattern for aggressor rows. 0: random, 1: all ones, 2: all zeros, 3: colstripe (0101), 4: inverse colstripe (1010), 5: checkered (0101, 1010), 6: inverse checkered (1010, 0101)")
        ("log_phys_scheme", value(&arg_log_phys_conv_scheme)->default_value(arg_log_phys_conv_scheme), "Specifies how to convert logical row IDs to physical row ids and the other way around. Pass 0 (default) for sequential mapping, 1 for the mapping scheme typically used in Samsung chips.")
        ("bitflip_counting_granularity", value(&bitflip_counting_granularity)->default_value(bitflip_counting_granularity), "When 0, counts and outputs the bitflips for each row. Otherwise based on the byte granularity provided with this parameter. E.g., 8 would report bitflips in every 8-byte chunk in the tested memory region.")
        
        ("append", bool_switch(&append_output)->default_value(append_output), "When specified, the output is appended to the --out file (if it exists). Otherwise the --out file is cleared.")
        ;

    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
        cout << desc << endl;
        return 0;
    }

    notify(vm);

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

    // make sure row_pattern contains only V(v) or A(a)
    if(!std::regex_match(row_layout, std::regex("^[VvAa]+$"))) {
        std::cerr << RED_TXT << "ERROR: --row_layout should contain only 'V' and 'A' characters. Provided: " << row_layout << NORMAL_TXT << std::endl;
        exit(-3);
    }

    // make sure row_layout is uppercase
    for (auto& c : row_layout) c = toupper(c);

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
    srand(0);
  
    auto t_prog_started = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed;
    bool check_time;

    if (trrref_sync) {
        uint trrref_dist = syncTRRREF(platform, input_data_victims, input_data_aggressors);

        if (trrref_dist == 0) {
            std::cout << RED_TXT << "ERROR: Could not find the period of the TRR REF operations!" << NORMAL_TXT << std::endl;
        }

        refs_per_loop = trrref_dist;
    }

    if(shift_refs)
        issue_REFs(platform, 1);

    hammerBank(platform, target_bank, row_range[0], row_range[1], num_hammers, num_ref_loops, refs_per_loop, trrref_sync, 
        num_dummy_rows, hammer_dummies_independently, hammer_dummies_before, hammer_dummies_after, dummy_banks,
        cascaded_hammer_aggr, cascaded_hammer_dummy, fake_hammer, fake_dummy_hammer, fake_ref, row_layout, input_data_victims, input_data_aggressors, 
        bitflip_counting_granularity, out_file);


    std::cout << "The test has finished!" << endl;

    out_file.close();

    return 0;
}
