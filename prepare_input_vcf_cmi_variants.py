import getopt, sys, re
from time import localtime, strftime

def help_print():
    print "Parse cmi variant output given by get_cmi_variant.R. Format for SnpEFF."
    print "-i, --input    str     output .txt file from get_cmi_variant.R."
    print "-o, --output   str     a .vcf file ready to be called by SnpEFF."
    
def prep_SnpEFF(db_variant_file, SnpEFF_vcf_file):
    with open(SnpEFF_vcf_file, 'w') as outfile, open(db_variant_file, 'r') as infile:
        next(infile)
        outfile.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for line in infile:
            line = line[:-1]
            fields = line.split('\t')
            CHROM = re.search('([0-9XY]+)$', fields[5]).group(0)
            POS, REF, ALT = fields[6:9]
            newline = "\t".join([CHROM, POS, "", REF, ALT, "", "", ""]) + '\n'
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
            db_variant_file = a
        elif o in ("-o", "--output"):
            SnpEFF_vcf_file = a        
        elif o in ("-h", "--help"):
            help_print()
            sys.exit()
    
    if len(sys.argv[1:]) > 0:
        #now = strftime("%Y%m%d-%H%M%S", localtime())
        #SnpEFF_vcf_file = "ngs_result_" + now + ".vcf"
        prep_SnpEFF(db_variant_file, SnpEFF_vcf_file)
    else:
        help_print()
    
if __name__=="__main__":
    main()
    
