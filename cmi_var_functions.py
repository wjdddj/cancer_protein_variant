import re
import warnings
import mygene
from class_protein import protein
from Bio import SeqIO, SeqUtils
import pandas as pd
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from Bio.SeqRecord import SeqRecord

##
def read_cmi_var(cmi_var):
    var_oi = pd.read_table(cmi_var, sep = '\t', low_memory=False)
    print "Total: %d annotated variants." % (var_oi.shape[0])
    
    var_oi_filter = var_oi.loc[
        pd.notnull(var_oi.Protein_Position) &
        pd.isnull(var_oi.ErrorsWarnings) &
        (var_oi.Transcript_biotype == 'protein_coding') &
        (var_oi.Putative_impact.isin(('HIGH', 'MODERATE'))) & # not LOW, or MODIFIER
        (var_oi.Annotation.isin(('missense_variant', 'stop_gained')))
    ]
    print "After filtering: %d annotated variants." % (var_oi_filter.shape[0])
    counts = var_oi_filter.Annotation.value_counts()
    print "missense: %d; nonsense: %d" % (counts['missense_variant'], counts['stop_gained'])
    print "\n"
    
    var_oi_byFeature = var_oi_filter.groupby('Feature_ID')
    variants = {}
    for name, group in var_oi_byFeature:
        variants[name] = list(set(group['HGVS.p']))
    return variants

## 
def read_uniprotFa(uniprotFa):
    handle = open(uniprotFa, 'r')
    refSeq = {}
    for record in SeqIO.parse(handle, "fasta"):
        enst_id = re.search('ENST[0-9]+', record.description).group(0)
        # ensp_id = re.search('ENSP[0-9]+', record.description).group(0)
        # ensg_id = re.search('ENSG[0-9]+', record.description).group(0)
        header = record.description
        seq = str(record.seq)
        # length = len(seq)
        # refSeq[enst_id] = {'transcript': enst_id, 'protein': ensp_id, 'gene': ensg_id, 'seq': seq, 'length': length}
        refSeq[enst_id] = {'header': header, 'seq': seq}
    handle.close()
    return refSeq

##
def get_longest_records(refSeq):
    #genes = []
    long_proteins = {}
    refSeq_long = {}
    for key in refSeq.keys():
        record = refSeq[key]
        protein1 = protein(record['header'], record['seq'])
        if protein1.gene not in long_proteins.keys():
            #genes.append(protein1.gene)
            long_proteins[protein1.gene] = protein1
            refSeq_long[protein1.gene] = record
        else:
            if protein1.length > long_proteins[protein1.gene].length:
                long_proteins[protein1.gene] = protein1
                refSeq_long[protein1.gene] = record
    
    refSeq_long_t = {}
    for key in refSeq_long.keys():
        record = refSeq_long[key]
        protein1 = protein(record['header'], record['seq'])
        refSeq_long_t[protein1.transcript] = record
        
    return refSeq_long_t

##
def convert_id(enst, field):
    mg = mygene.MyGeneInfo()
    ids = mg.query(
        enst, 
        scopes='ensembl.transcript',
        # fields = ['uniprot', 'symbol', 'name', 'summary']
        fields = field
    )
    return ids['hits'][0][field]

##
def pep_to_SeqRecord(pep):
    seq = Seq(pep[0], generic_protein)
    record = SeqRecord(seq, id = pep[6])
    record.annotations['Mutation'] = pep[7].replace('REF:', '')
    record.annotations['nMissCleavage'] = pep[4]
    record.annotations['ENST'] = str(pep[5])
    record.annotations['Position'] = int(pep[1])
    record.description = '|'.join([str(k)+':'+str(v) for k,v in record.annotations.items()])
    return record

##
def protein_to_SeqRecord(protein):
    seq = Seq(protein.seq, generic_protein)
    record = SeqRecord(seq, id = protein.id)
    record.annotations['ENST'] = protein.transcript
    record.annotations['ENSG'] = protein.gene
    record.annotations['ENSP'] = protein.protein    
    record.annotations['Mutation'] = protein.mutation
    record.description = '|'.join([str(k)+':'+str(v) for k,v in record.annotations.items()])
    return record

##
#def export_lookup():
#    return



















