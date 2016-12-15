import pandas as pd
import getopt, sys, vcf, re
from Bio import SeqIO
from time import localtime, strftime

def help_print():
    print "Obtain reference protein sequence from uniprot database using parsed snpEff txt file."
    print "-i, --input      str     output .txt file from parse_vcf_to_variants.py."
    print "-d, --database   str     uniprot reference fasta database"
    print "-o, --output     str     a .fa file containing reference sequences ready to be called by createVarPepDB.py."
    return

def my_fa_filter(records, va_enst):
    for record in records:
        fields = record.description.split(' ')
        transcript = re.search('ENST[0-9]+', fields[4].split(':')[1]).group(0)
        if transcript in va_enst:
            yield record

def get_refSeq(vcf_file, uniprotDB_file, refSeq_file):
    variant = pd.read_table(vcf_file, sep='\t')
    va_enst = set(variant['Feature_ID'])
    
    # uniprotDB_file = '/Users/jwang/Documents/Project6_Proteomics/databases/20160930_Ensembl_GRCh37/pep/Homo_sapiens.GRCh37.pep.all.fa'
    input_handle = open(uniprotDB_file, "rU")
    entr = SeqIO.parse(input_handle, 'fasta')
    SeqIO.write(my_fa_filter(entr, va_enst), refSeq_file, 'fasta')
    return
    
def main():
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "i:o:d:h", ["input=", "output=", "database=", "help"])
    except getopt.GetoptError as err:
        print str(err)
        help_print()
        sys.exit()
    
    for o,a in opts:
        if o in ("-i", "--input"):
            vcf_file = a
        elif o in ("-d", "--database"):
            uniprotDB_file = a        
        elif o in ("-h", "--help"):
            help_print()
            sys.exit()
    
    if len(sys.argv[1:]) > 0:
        #now = strftime("%Y%m%d-%H%M%S", localtime())
        #refSeq_file = "uniprot_refSeq_" + now + ".fa"
        get_refSeq(vcf_file, uniprotDB_file, refSeq_file)
    else:
        help_print()
    
if __name__=="__main__":
    main()
