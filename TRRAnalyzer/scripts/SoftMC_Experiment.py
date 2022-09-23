#!/usr/bin/env python3

# SingleTest and Experiment classes to use for storing experiment data

from math import floor
import copy
import pandas as pd
import numpy as np

class SingleTest:
    def __init__(self, path):
        self.path = path
        self.configs = self.__extractConfigParams(path)

    def __repr__(self):
        return f"confs: {self.configs}"

    def __extractConfigParams(self, path):
        confs = dict()

        dfile = open(self.path, 'r')

        while True:
            line = dfile.readline().strip()

            if line == "--- END OF HEADER ---":
                break

            if line == "": # end of the file
                break

            key, val = line.split("=")
            confs[key] = val

        dfile.close()

        return confs
    
    def getConfig(self, config_name):
        if config_name in self.configs:
            return self.configs[config_name]

        return ''

    def __parseTRRAnalyzerDataLine(self, line):
        tokens = line.split()
        row_id = int(tokens[2].replace(':', ''))
        num_bitflips = int(tokens[3].replace(':', ''))
        loc_bitflips = [int(x.replace(',', '')) for x in tokens[4:]]
        return row_id, num_bitflips, loc_bitflips

    def __parseTRRAnalyzerData(self, file_handle):
    
        bitflip_samples = dict()
        bitflip_locations = dict()
        dfs = dict()

        expr_count = 0
        colName = 'NumBitflips'

        lines = file_handle.readlines()

        for i, line in enumerate(lines):

            if line.startswith('Victim row'):
                row_id, num_bitflips, loc_bitflips = self.__parseTRRAnalyzerDataLine(line)
                bitflip_samples.setdefault(row_id,[]).append(num_bitflips)
                bitflip_locations.setdefault(row_id,[]).append(loc_bitflips)

            if line.startswith('Total bitflips') or i == (len(lines) - 1):
                
                for row_id in bitflip_samples.keys():
                    if row_id in dfs:
                        pd.concat([dfs[row_id], pd.DataFrame(bitflip_samples[row_id], columns=[colName])], axis=1)
                    else:
                        dfs[row_id] = pd.DataFrame(bitflip_samples[row_id], columns=[colName])

                bitflip_samples = dict()
                expr_count += 1 # this is for parsing multiple test results in a file. Currently, we parse only the most recent (i.e., latest) test results in the file
                colName = 'NumBitflips-' + str(expr_count)

        self.data = dfs
        self.bitflip_locations = bitflip_locations

    def parseTestData(self):
        dfile = open(self.path, 'r')
        self.__parseTRRAnalyzerData(dfile)
        dfile.close()

    def configsDigest(self, excluded_configs = []):
        filtered_configs = copy.deepcopy(self.configs)

        for c in excluded_configs:
            if c in filtered_configs.keys():
                del filtered_configs[c]

        return str(filtered_configs)
    
    def compareConfig(self, other_expr, configs_to_comp):

        for c in configs_to_comp:
            assert c in self.configs.keys(), f"ERROR: The provided configuration parameters '{c}' does not exist in {self.path}"
            assert c in other_expr.configs.keys(), f"ERROR: The provided configuration parameters '{c}' does not exist in {other_expr.path}"

            if self.configs[c] != other_expr.configs[c]:
                return False

        return True
        
    
    def convertToTRR(self, colname):

        for victim in self.data.keys():
            self.data[victim][colname] = self.data[victim][colname] == 0

    def listToVector(self, mdata, CHIP_NUM):
        label_data = max(mdata, key=len)
        vector_data = []
        for data in mdata:
            single_vector = []
            for label in label_data:
                if label in data:
                    single_vector.append(floor((label%512)/(512/CHIP_NUM)))
                else:
                    single_vector.append(99999)
            vector_data.append(single_vector)
        return label_data, vector_data

    def getChipIts(self, CHIP_NUM):
        bitflip_locations = []
        victims = list(self.bitflip_locations.keys())
        for it in range(len(self.bitflip_locations[victims[0]])):
            bitflip_list = []
            for victim in victims:
                bitflip_list += self.bitflip_locations[victim][it]
            bitflip_locations.append(bitflip_list)
        bitflips, bitflip_its = self.listToVector(bitflip_locations, CHIP_NUM)
        return bitflips, bitflip_its

    def printPotentialTRRTargets(self, colname, row_layout, trr_range = 1):
        victim_locs = []

        for i, c in enumerate(row_layout):
            if c == 'R' or c == 'r' or c == 'U' or c == 'u':
                victim_locs.append(i)

        # combine the colname column of all victims into a single dataframe
        combined_df = pd.DataFrame()
        for victim in self.data.keys():
            combined_df[victim] = self.data[victim][colname]

        trr_targets = combined_df.apply(lambda x: self.findTRRTarget(x, victim_locs, trr_range), axis=1)

        return trr_targets.to_csv(header=False)


    def findTRRTarget(self, victims_trr_data, victim_locs, trr_range):

        # return None if no victims were TRR'd
        if not any(victims_trr_data):
            return None

        # consider each aggressor and victim row as potential TRR targets and check if the trr data matches
        for i in range(0, victim_locs[-1] + 1):

            must_TRR_rows = []

            # check rows on the left
            cur_row = max(i - trr_range, 0)

            for ii in range(cur_row, i):
                if ii in victim_locs:
                    must_TRR_rows.append(ii)

            # check rows on the right
            cur_row = i + 1
            last_row = min(i + trr_range, victim_locs[-1])

            for ii in range(cur_row, last_row + 1):
                if ii in victim_locs:
                    must_TRR_rows.append(ii)

            is_match = True

            for ind in must_TRR_rows:
                if victims_trr_data[victims_trr_data.keys()[victim_locs.index(ind)]] == False: # probably it is not a good idea to index keys() since the code depends on the order of the keys
                    is_match = False
                    break

            if len(must_TRR_rows) != victims_trr_data.value_counts()[True]:
                is_match = False

            if is_match:
                return i

        return "No TRR Candidate, normal refresh?"

    # Returns true if the test configuration matches all parameters in 'conf'
    def matches(self, conf):
        if all((c[1] == 'none' if c[0] not in self.configs.keys() else c in self.configs.items()) or (c[1].lower() == 'false' and c[0] not in self.configs.keys()) for c in conf.items()):
            return True

        return False

    # returns a list of iteration IDs that satisfy the condition specified by function 'df_apply', which goes to the DataFrame.apply()
    # An example df_apply() -> lambda x: (x == 1).all(), True if all elements in a row are 0
    def get_iterations_satisfying(self, df_apply):
        
        res = self.data.apply(df_apply, axis='columns')

        return res.index[res == True].tolist()

    # Should be called after convertToTRR() is called
    def get_iterations_with_common_REF(self):

        row_ids = list(self.data)

        if len(row_ids) == 0:
            return ''

        res = self.data[row_ids[0]].apply(lambda x: (x == 1).all(), axis='columns')

        for row_id in row_ids:
            res = res & self.data[row_id].apply(lambda x: (x == 1).all(), axis='columns')

        return res.index[res == True].tolist()


    # Should be called after convertToTRR() is called
    def get_iterations_with_single_REF(self):

        row_ids = list(self.data)

        if len(row_ids) == 0:
            return ''
        
        all_TRRs = self.data[row_ids[0]].apply(lambda x: (x == 1).all(), axis='columns')

        for row_id in row_ids:
            all_TRRs = all_TRRs & self.data[row_id].apply(lambda x: (x == 1).all(), axis='columns')

        single_REF_indices = []

        for row_id in row_ids:
            tmp_df = ~all_TRRs & self.data[row_id].apply(lambda x: (x == 1).all(), axis='columns')

            single_REF_indices.append(tmp_df.index[tmp_df == True].tolist())
            single_REF_indices.append(-1)

        return single_REF_indices