import getopt, sys, vcf
from time import localtime, strftime

def help_print():
    print "Parse SnpEff output .vcf file into a .txt file."
    print "-i, --input    str     output .vcf file from SnpEff."
    print "-o, --output   str     a .txt file containing variants ready to be called by createVarPepDB.py."
    
def parse_SnpEFF(vcf_file, variant_file):
    vcf_reader = vcf.Reader(open(vcf_file))
    outfile = open(variant_file, 'w')
    header = ['CHROM', 
              'POS', 
              'ALT', 
              'Annotation', 
              'Putative_impact', 
              'Gene_Name', 
              'Gene_ID', 
              'Feature_type', 
              'Feature_ID', 
              'Transcript_biotype', 
              'Rank', 
              'HGVS.c', 
              'HGVS.p', 
              'cDNA_position', 
              'CDS_position', 
              'Protein_Position', 
              'Distance_to_Feature', 
              'ErrorsWarnings']
    first_row = '\t'.join(header) + '\n'
    outfile.write(first_row)
    for record in vcf_reader:
        INFO = record.INFO['ANN']
        CHROM = record.CHROM
        POS = record.POS
        for info in INFO:
            info = info.split('|')
            info.insert(0, str(CHROM))
            info.insert(1, str(POS))
            newline = '\t'.join(info) + '\n'
            outfile.write(newline)

def main():
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "i:o:h", ["input=", "output=", "help"])
    except getopt.GetoptError as err:
        print str(err)
        help_print()
        sys.exit()
    
    for o,a in opts:
        if o in ("-i", "--input"):
            vcf_file = a
        elif o in ("-o", "--output"):
            variant_file = a        
        elif o in ("-h", "--help"):
            help_print()
            sys.exit()
    
    if len(sys.argv[1:]) > 0:
        #now = strftime("%Y%m%d-%H%M%S", localtime())
        #variant_file = "variant_" + now + ".txt"
        parse_SnpEFF(vcf_file, variant_file)
    else:
        help_print()
    
if __name__=="__main__":
    main()

