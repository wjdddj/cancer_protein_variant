from class_variant import variant
from class_protein import protein
from class_peptide import peptide
import cmi_var_functions as cfunc
import time, sys, getopt
from Bio import SeqIO, SeqUtils
from time import localtime, strftime

##
def print_help():
    print "Create fasta database for mutant proteins or mutant trypsinized peptides."
    print "USAGE: python createVarPepDB.py -s protein -o xxx [uniprot_refSeqxxx.txt] [variantxxx.txt]"
    print "       python createVarPepDB.py -s peptide -o xxx -n2 -m6 -x144 [uniprot_refSeqxxx.txt] [variantxxx.txt]"
    print "\n"
    print "OPTIONS: "
    print "-s, --setting     str     type of output fasta database. [protein|peptide]"
    print "-o, --output_dir  str     output directory."
    print "-n, --max_miss    int     maximum number of miss cleavage site allowed for trypsinization."
    print "-m, --min_len     int     mininum length of reported peptides."
    print "-x, --max_len     int     maximum length of reported peptides"
    print "-w, --mut_only    boo     when 1, the output database only contains mutants. [0|1]"
    print "-h, --help"  

##    
def input_argv(argv):
    try: 
        opts, args = getopt.getopt(
            argv, 
            "s:o:n:m:x:w:h",
            ["setting=", 
             "output==",
             "max_miss=", 
             "min_len=", 
             "max_len=", 
             "mut_only=", 
             "help"]
        )
    except getopt.GetoptError as err:
        print str(err)
        print_help()
        sys.exit()
        
    for opt, arg in opts:
        if opt in ("-s", "--setting"):
            setting = arg
        elif opt in ("-o", "--output_dir"):
            output_dir = arg
        elif opt in ("-n", "--max_miss"):
            max_miss = int(arg)
        elif opt in ("-m", "--min_len"):
            min_len = int(arg)
        elif opt in ("-x", "--max_len"):
            max_len = int(arg)
        elif opt in ("-w", "--mut_only"):
            mut_only = int(arg)
        elif opt in ("-h", "--help"):
            print_help()
            sys.exit()
        else:
            assert False
    
    if len(args) != 2:
        assert False
        print "please input refSeq and variants file."
        print_help()
        sys.exit()
    else:
        refSeq_file = args[0]
        variants_file = args[1]
    
    if setting == 'peptide':
        return [setting, refSeq_file, variants_file, output_dir, mut_only, max_miss, min_len, max_len]
    elif setting == 'protein':
        return [setting, refSeq_file, variants_file, output_dir, mut_only]
    
##
def build_pepDB(refSeq_file, variants_file, output_dir, mut_only, max_miss, min_len, max_len):
    
    #max_miss = 2
    #min_len = 6
    #max_len = 144
    
    variants_read = cfunc.read_cmi_var(variants_file)
    refSeq_read = cfunc.read_uniprotFa(refSeq_file)
    
    # filter variants to focus on the ones only associated with the longest transcript
    refSeq = refSeq_read
    # refSeq = cfunc.get_longest_records(refSeq)
    variants = {}
    for key, var in variants_read.iteritems():
        if key in refSeq.keys():
            variants[key] = var    
    
    # obtain trypsinized peptide profiles considering mutations in R, K, P that affect trypsin digestion
    print 'Annotating reference peptides...'
    pep_profile = []
    len_protein = 0
    for key in variants.keys():
        eachProtein = refSeq[key]
        protein1 = protein(eachProtein['header'], eachProtein['seq'])
        protein1.add_trypsin_profile(max_miss = max_miss)
        protein1.annotate_pep()
        matchVariants = variants[protein1.transcript]
        for var in matchVariants:
            variant1 = variant(protein1.transcript, var)
            variant1.get_seq(refSeq)
            if variant1.is_RKP(): # check whether mutation will affect trypsin digestion
                protein1.annotate_pepvar(variant1, max_miss = max_miss)
        protein1.exclude_len(min_len = min_len, max_len = max_len)
        pep_profile.extend(protein1.trypsin_profile)
        
        len_protein+= 1
        len_profile = len(pep_profile)
        if len_protein % 100 == 0:
            print 'Finished %d reference proteins.' % (len_protein)
            print 'Annotated %d reference peptides.' % (len_profile)
    print '\nFinished %d reference proteins.' % (len_protein)
    print 'Annotated %d reference peptides.\n' % (len_profile)
    
    # generate all peptides with single mutations at trypsinized peptide level
    print 'Annotating mutant peptides...'
    mut_pept = []
    len_peptide = 0
    for pep in pep_profile:
        peptide1 = peptide(pep)
        if peptide1.mutation == 'REF': # exclude mutants that affect trypsin digestion
            peptide1.find_variants(variants, refSeq)
            peptide1.generate_single_mut_peptide()
            peptide1.annotate()
            peptide1.exclude_len(min_len = min_len, max_len = max_len)
            mut_pept.extend(peptide1.annotated_profile)
        
        len_peptide+= 1
        len_profile = len(mut_pept)
        if len_peptide % 10000 == 0:
            print 'Finished %d reference peptides.' % (len_peptide)
            print 'Annotated %d mutations.' % (len_profile)
    print '\nFinished %d reference peptides.' % (len_peptide)
    print 'Annotated %d mutations.\n' % (len_profile)
    
    # remove duplicates due to the same nonsense mutation on different miss cleavage yielded peptides
    print 'Removing replicates...'
    uniqIDs = set()
    uniq_mut_pept = []
    for pep in mut_pept:
        uniqID = '_'.join([pep[0], str(pep[1]), str(pep[2])]) #unique peptide defined by seq+position+length
        if uniqID not in uniqIDs:
            uniq_mut_pept.append(pep)
        uniqIDs.update([uniqID])
        
    len_profile = len(uniq_mut_pept)
    print 'Removed %d mutant peptides due to nonsense mutation on different miss cleavage yielded peptides.' % (len(mut_pept) - len_profile)
    print 'Finally, annotated %d peptides.' % (len_profile)
    
    pep_profile.extend(uniq_mut_pept)
    
    if mut_only == 1:
        print "Export only mutant peptides..."
        pep_profile = [pep for pep in pep_profile if pep[7] != 'REF']
                
    records = []
    for pep in pep_profile:
        records.append(cfunc.pep_to_SeqRecord(pep))
    now = strftime("%Y%m%d-%H%M%S", localtime())
    out_db_file = output_dir + "mutPeptideDB_" + now + ".fa"
    SeqIO.write(records, out_db_file, 'fasta')
    
    return


def build_protDB(refSeq_file, variants_file, output_dir, mut_only):
    
    variants_read = cfunc.read_cmi_var(variants_file)
    refSeq = cfunc.read_uniprotFa(refSeq_file)
    
    # filter variants according to refSeq.
    variants = {}
    for key, var in variants_read.iteritems():
        if key in refSeq.keys():
            variants[key] = var      
    
    # Annotate proteins
    print 'Annotating proteins and mutants...'
    prot_profile = []
    len_protein = 0
    len_mutant = 0
    for key in variants.keys():
        eachProtein = refSeq[key]
        protein1 = protein(eachProtein['header'], eachProtein['seq'])
        prot_profile.append(protein1)
        len_protein+= 1
        matchVariants = variants[protein1.transcript]
        for var in matchVariants:
            mprotein = protein(eachProtein['header'], eachProtein['seq'])
            variant1 = variant(protein1.transcript, var)
            mprotein.add_mutation(variant1)
            prot_profile.append(mprotein)
            len_mutant+= 1
    len_profile = len(prot_profile)
    print '\nFinished %d reference proteins.' % (len_protein)
    print '\nFinished %d mutant proteins.' % (len_mutant)
    print 'Annotated %d peptides.\n' % (len_profile)    
    
    if mut_only == 1:
        print "Export only mutant proteins..."
        prot_profile = [prot for prot in prot_profile if prot.mutation != 'REF']
    
    records = []
    for prot in prot_profile:
        records.append(cfunc.protein_to_SeqRecord(prot))
    now = strftime("%Y%m%d-%H%M%S", localtime())
    out_db_file = output_dir + "mutProteinDB_" + now + ".fa"
    SeqIO.write(records, out_db_file, 'fasta')
    return

def main():
    config = input_argv(sys.argv[1:])
    if config[0] == 'protein':
        refSeq_file, variants_file, output_dir, mut_only = config[1:]
        build_protDB(refSeq_file, variants_file, output_dir, mut_only)
    elif config[0] == 'peptide':
        refSeq_file, variants_file, output_dir, mut_only, max_miss, min_len, max_len = config[1:]
        #print refSeq_file, variants_file, mut_only, max_miss, min_len, max_len
        #print type(max_miss)
        build_pepDB(refSeq_file, variants_file, output_dir, mut_only, max_miss, min_len, max_len)
    return

if __name__=="__main__":
    main()