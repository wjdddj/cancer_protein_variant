import re
import warnings
import itertools
from class_variant import variant
from class_protein import protein

class peptide:
    
    def __init__(self, peptide_profile):
        self.seq, self.start_pos, self.length, self.count, self.n_miss, self.transcript, self.pepID, self.mutation = peptide_profile
        #self.transcript, self.pepID, self.mutation = peptide_profile[5:8]
        return
    
    # find variants associated with the peptide
    # shift variants refer to peptide start position in the reference protein
    def find_variants(self, variants, refSeq):
        protVariants = variants[self.transcript]
        self.peptVariants = []
        for eachVar in protVariants:
            var = variant(self.transcript, eachVar)
            var.get_seq(refSeq)
            if self.start_pos <= var.pos <= self.start_pos + self.length - 1:
                if not var.is_RKP():
                    #self.peptVariants.append(var.var3)
                    self.peptVariants.append(var)
        return
    
    # organize by position
    def organize_variants(self):
        positions = set([item.pos for item in self.peptVariants])
        self.peptVarDict = {}
        for key in positions:
            self.peptVarDict[key] = []
        
        for eachVar in self.peptVariants:
            self.peptVarDict[eachVar.pos].append(eachVar)
        return
    
    # var_set a set of variant objects
    def multi_mutate(self, var_set):
        newSeq = self.seq
        for var in var_set:
            if var.newAA1 == '*': # if nonsense mutation, no need to continue mutate after the stop
                break
            else:
                newSeq = protein.mutate(newSeq, var)
        return newSeq   
    
    # generate all combinations of mutations
    def generate_multi_mut_peptide(self):
        max_mut = len(self.peptVarDict)
        var_sets = []
        for i in range(max_mut):
            ids = list(itertools.combinations(self.peptVarDict.keys(), i+1))
            for idx in ids:
                candidates = [self.peptVarDict[key] for key in idx]
                var_sets.extend(list(itertools.product(*candidates)))
        
        self.mut_pept = {}
        for var_set in var_sets:
            var_set_list = list(var_set)
            var_set_list.sort(key = lambda x:x.pos)
            var1_set = [var.var1 for var in var_set_list]
            shifted_var_set = [variant.shift_var(var.var1, self.start_pos) for var in var_set_list]
            shifted_var_set = [variant(self.pepID, var) for var in shifted_var_set]
            # shifted_var_set.sort(key = lambda x:x.pos) # sort variants based on ascending positions
            self.mut_pept['_'.join(var1_set)] = self.multi_mutate(shifted_var_set)
        return
    
    def generate_single_mut_peptide(self):
        self.mut_pept = {}
        for var in self.peptVariants:
            shifted_var = variant.shift_var(var.var1, self.start_pos)
            shifted_var = variant(self.pepID, shifted_var)
            self.mut_pept[var.var1] = protein.mutate(self.seq, shifted_var)
        return
    
    def annotate(self):
        self.annotated_profile = []
        for key in self.mut_pept.keys():
            item = self.mut_pept[key]
            mut_profile = [item, self.start_pos, len(item), self.count, self.n_miss, self.transcript, self.pepID+'_'+key, self.mutation+':'+':'.join(key.split('_'))]
            self.annotated_profile.append(mut_profile)
        return
    
    def exclude_len(self, min_len = 6, max_len = 300):
        self.annotated_profile = [item for item in self.annotated_profile if min_len <= item[2] <= max_len]
        return
    