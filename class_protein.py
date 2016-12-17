import re
import warnings
from class_variant import variant

class protein:
    
    def __init__(self, header, sequence):
        
        self.head = header
        self.seq = sequence
        self.length = self.get_len(self.seq)
        self.rm_X()
        
        self.transcript = self.get_transcript(self.head)
        self.protein = self.get_protein(self.head)
        self.gene = self.get_gene(self.head)
        self.name = self.get_name(self.head)
        self.uniprotAC = self.get_uniprotAC(self.head)
        self.id = self.transcript + '_REF'
        self.mutation = 'REF'
    
    @staticmethod
    def get_transcript(head):
        transcript = re.search('ENST[0-9]+', head).group(0)
        return transcript
    
    @staticmethod
    def get_protein(head):
        protein = re.search('ENSP[0-9]+', head).group(0)
        return protein
    
    @staticmethod
    def get_gene(head):
        gene = re.search('ENSG[0-9]+', head).group(0)
        return gene
    
    @staticmethod
    def get_name(head):
        name = re.search('(?<=Gene_Name:).*', head).group(0)
        return name
    
    @staticmethod
    def get_uniprotAC(head):
        uniprotAC = re.search('(?<=UniProtKBAC:).*', head).group(0)
        return uniprotAC
    
    @staticmethod
    def get_len(seq):
        length = len(seq)
        return length
    
    #@staticmethod
    #def get_id(ID, idType, mapTab):
        
    #    return
    
    ## trypsin digestion generator, cut after K or R when not before P.
    ## returns substring of digested peptide, peptide start position, peptide length, peptide No (order from the start)
    #  ['peptide seq', 
    #  'start position in reference', 
    #  'peptide length', 
    #  'order from the start']
    @staticmethod
    def trypsin(bases):
        #subs = []
        pos = [1]
        count = 0
        while bases:
            k, r = re.search('K[^P]', bases), re.search('R[^P]', bases)
            if k is None and r is None:
                cut = len(bases)
            elif k is None:
                cut = r.start() + 1
            elif r is None:
                cut = k.start() + 1   
            else:
                cut = min(k.start(), r.start()) + 1
            pos.append(cut + pos[-1])
            sub = bases[:cut]
            sub_len = cut
            bases = bases[cut:]
            count = count + 1
            # subs.append(sub)
            pep_profile = [sub, pos[-2], sub_len, count]
            yield pep_profile
    
    
    # output structure for each peptide is 
    # [str(peptide seq), 
    #  int(start position in reference), 
    #  int(peptide length), 
    #  int(order from the start), 
    #  int(number of miss cleavage)]
    @staticmethod
    def get_trypsin_profile(seq, max_miss = 0):
        trypsin_out = list(protein.trypsin(seq))
        trypsin_profile = []
        for n_miss in range(max_miss + 1):
            for i in range(len(trypsin_out) - n_miss):
                sub, pos, sub_len, count = trypsin_out[i]
                for j in range(n_miss):
                    sub = sub + trypsin_out[i+j+1][0]
                    sub_len = sub_len + trypsin_out[i+j+1][2]
                trypsin_profile.append([sub, pos, sub_len, count, n_miss])
        return trypsin_profile
    
    
    # input: var, a variant object
    # input: seq, a string object
    @staticmethod
    def mutate(seq, var):
        var_len = len(var.refAA1)
        if var.pos > len(seq) or var.pos < 0:
            warnings.warn('variant position is out of range! No change at this time!\n')
            print var.transcript, var.var1
            newSeq = seq
        else:
            posAA = seq[var.pos-1 : var.pos-1+var_len]
            if posAA != var.refAA1:
                warnings.warn('Reference Sequence has different AA from the variant call! No change at this time! Please double check!\n')
                print var.transcript, var.var1
                newSeq = seq
            else:
                if re.search('\\*', var.newAA1):
                    stop_site = re.search('\\*', var.newAA1).span()[0]
                    newSeq = seq[:var.pos-1] + var.newAA1[:stop_site] # stop_gained
                else:
                    newSeq = ''.join(seq[:var.pos-1] + var.newAA1 + seq[var.pos-1+var_len:])
        return newSeq      
    
    # remove leading "X" from reference sequence.
    def rm_X(self):
        if self.seq[0] == 'X':
            self.seq = self.seq[1:]
        return
    
    def add_trypsin_profile(self, max_miss = 0):
        self.trypsin_profile = self.get_trypsin_profile(self.seq, max_miss = max_miss)
        return
    
    # For reference sequence
    # append last three items to trypsin_profile to make it
    # [str(peptide seq), 
    #  int(start position in reference), 
    #  int(peptide length), 
    #  int(order from the start), 
    #  int(number of miss cleavage),
    #  str(ENST), 
    #  str(pepID), 
    #  str(variant type),
    #  str(Gene Name),
    #  str(UniProtKBAC)]
    def annotate_pep(self):
        for item in self.trypsin_profile:
            item.append(self.transcript) # Ensembl transcript id
            item.append(self.transcript + 
                        '_TRP' + 
                        str(item[1]).zfill(6) + 
                        str(item[2]).zfill(4) + 
                        'REF') # pepID
            item.append('REF') # variant information reference to full length protein
            item.append(self.name)
            item.append(self.uniprotAC)
            # self.trypsin_profile[i].append('Ref') # variant information reference to digested peptide
        return
    
    # For mutant sequence, mutations that affecting trypsinization. 
    # mainly used to append to trypsin_profile only novel peptides resulted from trypsinization
    # [str(peptide seq), 
    #  int(start position in reference), 
    #  int(peptide length), 
    #  int(order from the start), 
    #  int(number of miss cleavage),
    #  str(ENST), 
    #  str(pepID), 
    #  str(variant type),
    #  str(Gene Name),
    #  str(UniProtKBAC)]
    def annotate_pepvar(self, var, max_miss = 0):
        pep_profile_seq = set([item[0]+str(item[1]) for item in self.trypsin_profile]) # use both seq and position to id novel peptide
        newSeq = self.mutate(self.seq, var)
        if newSeq != self.seq:
            pep_profile_var = protein.get_trypsin_profile(newSeq, max_miss = max_miss)
            for item in pep_profile_var:
                if item[0]+str(item[1]) not in pep_profile_seq:
                    item.append(self.transcript) # Ensembl transcript id
                    item.append(self.transcript + 
                                '_TRP' + 
                                str(item[1]).zfill(6) + 
                                str(item[2]).zfill(4) + 
                                var.var1) # pepID
                    item.append(var.var1) # variant information reference to full length protein defined by ENST ID
                    item.append(self.name)
                    item.append(self.uniprotAC)
                    self.trypsin_profile.append(item)
        return
    
    def exclude_len(self, min_len = 6, max_len = 300):
        self.trypsin_profile = [item for item in self.trypsin_profile if min_len <= item[2] <= max_len]
        return
    
    def add_mutation(self, var):
        self.seq = protein.mutate(self.seq, var)
        self.id = self.id + "_" + var.var1
        self.mutation = var.var1
        self.length = len(self.seq)
        return
    
    
    
    
    
    