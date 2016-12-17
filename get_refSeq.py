import pandas as pd
import getopt, sys, vcf, re
from Bio import SeqIO
from time import localtime, strftime

def help_print():
    print "Obtain reference protein sequence from Ensembl database using parsed snpEff txt file."
    print "-i, --input      str     output .txt file from parse_vcf_to_variants.py."
    print "-d, --database   str     Ensembl protein reference fasta database."
    print "-m, --map_tab    str     uniprot ID mapping table to convert ensembl transcript ID to uniprot accession id."
    print "-o, --output     str     a .fa file containing reference sequences ready to be called by createVarPepDB.py."
    return

def read_mapTab(map_tab_file, va_enst):
    map_tab = {}
    tmp_tab = {}
    with open(map_tab_file, 'r') as infile:
        for line in infile:
            line = line[:-1]
            fields = line.split('\t')
            if fields[1] == 'Ensembl_TRS':
                if fields[2] in va_enst:
                    map_tab[fields[2]] = {}
                    map_tab[fields[2]]['uniprotID'] = fields[0]
                    map_tab[fields[2]]['Gene_Name'] = ""
                    tmp_tab[fields[0]] = fields[2]
                    
    map_tab_values = set([value['uniprotID'] for value in map_tab.values()])
    with open(map_tab_file, 'r') as infile:
        for line in infile:
            line = line[:-1]
            fields = line.split('\t')
            if fields[0] in map_tab_values:
                if fields[1] == 'Gene_Name':
                    key = tmp_tab[fields[0]]
                    map_tab[key]['Gene_Name'] = fields[2]
        
    return map_tab

def my_fa_filter(records, va_enst, map_tab):
    for record in records:
        fields = record.description.split(' ')
        transcript = re.search('ENST[0-9]+', fields[4].split(':')[1]).group(0)
        if transcript in va_enst:
            if transcript in map_tab.keys():
                uniprotID = map_tab[transcript]['uniprotID']
                record.description = record.description + " UniProtKBAC:" + uniprotID
                gene_name = map_tab[transcript]['Gene_Name']
                record.description = record.description + " Gene_Name:" + gene_name
                yield record
            
def get_refSeq(vcf_file, enspDB_file, map_tab_file, refSeq_file):
    
    variant = pd.read_table(vcf_file, sep='\t', low_memory = False)
    va_enst = []
    for enst in variant['Feature_ID']:
        if re.search('ENST[0-9]+', enst):
            va_enst.append(re.search('ENST[0-9]+', enst).group(0))
    va_enst = set(va_enst)
    map_tab = read_mapTab(map_tab_file, va_enst)
    
    # enspDB_file = '/Users/jwang/Documents/Project6_Proteomics/databases/20160930_Ensembl_GRCh37/pep/Homo_sapiens.GRCh37.pep.all.fa'
    input_handle = open(enspDB_file, "rU")
    entr = SeqIO.parse(input_handle, 'fasta')
    SeqIO.write(my_fa_filter(entr, va_enst, map_tab), refSeq_file, 'fasta')
    return
    
def main():
    try: 
        opts, args = getopt.getopt(
            sys.argv[1:], 
            "i:o:d:m:h", 
            ["input=", 
             "output=", 
             "database=", 
             "map_tab=",
             "help"]
        )
    except getopt.GetoptError as err:
        print str(err)
        help_print()
        sys.exit()
    
    for o,a in opts:
        if o in ("-i", "--input"):
            vcf_file = a
        elif o in ("-o", "--output"):
            refSeq_file = a
        elif o in ("-d", "--database"):
            enspDB_file = a
        elif o in ("-m", "--map_tab"):
            map_tab_file = a
        elif o in ("-h", "--help"):
            help_print()
            sys.exit()
    
    if len(sys.argv[1:]) > 0:
        #now = strftime("%Y%m%d-%H%M%S", localtime())
        #refSeq_file = "ensp_refSeq_" + now + ".fa"
        get_refSeq(vcf_file, enspDB_file, map_tab_file, refSeq_file)
    else:
        help_print()
    
if __name__=="__main__":
    main()
