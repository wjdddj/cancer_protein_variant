import re
from Bio import SeqUtils

class variant:
    
    def __init__(self, transcript, var):
        self.transcript = transcript
        vtype = variant.var_type(var)
        if vtype == 'var3':
            self.refAA3 = variant.get_refAA(var)
            self.newAA3 = variant.get_newAA(var)
            self.pos = variant.get_pos(var)
            self.var3 = ''.join([self.refAA3, str(self.pos), self.newAA3])
            self.update_var1()
        elif vtype == 'var1':
            self.var1 = var
            self.refAA1 = variant.get_refAA(var)
            self.newAA1 = variant.get_newAA(var)
            self.pos = variant.get_pos(var)
            self.var1 = ''.join([self.refAA1, str(self.pos), self.newAA1])
            self.update_var3()            

    @staticmethod
    def get_refAA(var):
        #aa_var = self.var
        refAA = re.search('[A-Za-z\\*]*(?=[0-9])', var).group(0)
        return refAA
    
    @staticmethod
    def get_newAA(var):
        #aa_var = self.var
        newAA = re.search('(?<=[0-9])[A-Za-z\\*]*$', var).group(0)
        return newAA

    @staticmethod
    def get_pos(var):
        #aa_var = self.var
        pos = int(re.search('[0-9]+', var).group(0))
        return pos

    @staticmethod
    def var_type(var):
        newAA = variant.get_newAA(var)
        refAA = variant.get_refAA(var)
        if (len(newAA)%3 == 0 or len(refAA)%3 == 0) and (re.search('[a-z]', newAA) or re.search('[a-z]', refAA)):
            return 'var3'
        elif (re.search('[^a-z]', newAA) and re.search('[^a-z]', refAA)):
            return 'var1'
        else:
            return None
    
    ## convert aa variant annotated by three letters to a single letter
    @staticmethod
    def var3_to_var1(var3):
        refAA3 = variant.get_refAA(var3)
        newAA3 = variant.get_newAA(var3)
        pos = variant.get_pos(var3)
        refAA1 = "*" if refAA3 == "*" else SeqUtils.seq1(refAA3)
        # refAA1 = SeqUtils.seq1(refAA3)
        newAA1 = "*" if newAA3 == "*" else SeqUtils.seq1(newAA3)
        # newAA1 = SeqUtils.seq1(newAA3)
        var1 = ''.join([refAA1, str(pos), newAA1])
        return var1

    ## convert aa variant annotated by a single letter to three letters
    @staticmethod
    def var1_to_var3(var1):
        refAA1 = variant.get_refAA(var1)
        newAA1 = variant.get_newAA(var1)
        pos = variant.get_pos(var1)
        refAA3 = "*" if refAA1 == "*" else SeqUtils.seq3(refAA1)
        # refAA1 = SeqUtils.seq1(refAA3)
        newAA3 = "*" if newAA1 == "*" else SeqUtils.seq3(newAA1)
        # newAA1 = SeqUtils.seq1(newAA3)
        var3 = ''.join([refAA3, str(pos), newAA3])
        return var3
    
    # input: var, a string object. e.g. "A15T"
    # output: shifted_var, a string object
    @staticmethod
    def shift_var(var, start_pos):
        pos = variant.get_pos(var)
        new_pos = str(pos - start_pos + 1)
        shifted_var = re.sub('[0-9]+', new_pos, var)
        return shifted_var
    
    def update_var1(self):
        self.var1 = variant.var3_to_var1(self.var3) 
        self.refAA1 = variant.get_refAA(self.var1)
        self.newAA1 = variant.get_newAA(self.var1) 
        return

    def update_var3(self):
        self.var3 = variant.var1_to_var3(self.var1) 
        self.refAA3 = variant.get_refAA(self.var3)
        self.newAA3 = variant.get_newAA(self.var3) 
        return
    
    ##
    def get_seq(self, refSeq):
        self.seq = refSeq[self.transcript]['seq']
        return
        
    ## examine whether the variant will affect trypsin digestion
    def is_RKP(self):
        if hasattr(self, 'seq'):
            seq = self.seq
            pos = self.pos
            length = len(self.seq)
            if re.search('[RK]', self.var1):
                if pos == length:
                    return False
                elif seq[pos] == 'P':
                    return False
                else:
                    return True
            elif re.search('P', self.var1):
                if re.search('[RK]', seq[pos-2]):
                    return True
                else:
                    return False
            else:
                return False
        else:
            print "please run variant_obj.get_seq() before this step!\n"
            return
    
    def update_shift(self, start_pos):
        self.var1 = variant.shift_var(self.var1, start_pos)
        self.var3 = variant.shift_var(self.var3, start_pos)
        self.pos = variant.get_pos(self.var3)
        self.seq = None
        return

    
    
    
    
    
    
    
    
    
    
    