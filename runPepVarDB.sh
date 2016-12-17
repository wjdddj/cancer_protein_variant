ts(){
    date +"%Y%m%d-%H%M%S"
}

current_dir=`pwd`
output_dir="fastaDB_"`ts`"/"
tmp_dir=$output_dir"temp/"
snpEff=$HOME"/snpEff/snpEff.jar"
get_cmi_var=$current_dir"/get_cmi_variant.R"
prep_input_vcf=$current_dir"/prepare_input_vcf_cmi_variants.py"
parse_snpEff_vcf=$current_dir"/parse_vcf_to_variants.py"
get_refSeq=$current_dir"/get_refSeq.py"
ensemblDB=$HOME"/Documents/Project6_Proteomics/databases/20160930_Ensembl_GRCh37/pep/Homo_sapiens.GRCh37.pep.all.fa"
map_tab_file=$HOME"/Users/jwang/Documents/Project6_Proteomics/databases/20160915_Uniprot_ID_mapping/HUMAN_9606_idmapping.dat"
createDB=$current_dir"/createVarPepDB.py"
Genome="GRCh37.75"

mkdir $output_dir
mkdir $tmp_dir
# obtain cmi variants
echo "obtaining cmi variants from current database..."
db_variant_file=$tmp_dir"ngs_result_"`ts`".txt"
Rscript $get_cmi_var $db_variant_file

# prepare and run snpEff
echo "preparing snpEff input vcf file..."
in_vcf_file=$tmp_dir"ngs_result_"`ts`".vcf"
python $prep_input_vcf -i $db_variant_file -o $in_vcf_file
out_vcf_file=$tmp_dir"cmi_"$Genome"_"`ts`".ann.vcf"
echo "running snpEff..."
java -Xmx4g -jar $snpEff -v $Genome $in_vcf_file > $out_vcf_file

# parse snpEff output .vcf file
echo "parsing snpEff output vcf file..."
variant_file=$tmp_dir"variant_"`ts`".txt"
python $parse_snpEff_vcf -i $out_vcf_file -o $variant_file

# obtain reference uniprot protein sequences for the annotated variants
echo "obtaining reference uniprot sequences..."
refSeq_file=$tmp_dir"ensembl_refSeq_"`ts`".fa"
python $get_refSeq -i $variant_file -d $ensemblDB -o $refSeq_file -m $map_tab_file

# create mutant peptide fasta database
echo "creating mutant database..."
python $createDB -s peptide -o $output_dir -n2 -m6 -x144 -w1 $refSeq_file $variant_file










