rm(list = ls())
source('~/R_modules/caris_basic/caris_basic.R')
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))

args = commandArgs(trailingOnly = TRUE)
if(length(args) == 0)
    stop('Please provide output file name and directory!\n')
    
output_file = args[1]

# now <- format(Sys.time(), '%Y%m%d-%H%M%S')

#### obtain CMI variants
NGS_proteinChange <- get_sting_query(
  'bioinfo',
  "select distinct test, ORG_Result, Result, calculatedconclusion, proteinChange, chrom, position, Ref, Alt 
  from bioinfo.nextgen_all_test_results 
  where technology = 'NGS Q3';
  "
)
# head(NGS_proteinChange)
# table(NGS_proteinChange$ORG_Result)
# table(NGS_proteinChange$Result)

# NGS_proteinChange[is.na(NGS_proteinChange$chrom)[1:10], ]
# sapply(NGS_proteinChange, function(x) sum(is.na(x)))

NGS_proteinChange <- NGS_proteinChange[-which(
  is.na(NGS_proteinChange$proteinChange)
  | is.na(NGS_proteinChange$chrom)
  | is.na(NGS_proteinChange$position)
  | is.na(NGS_proteinChange$Ref)
  | is.na(NGS_proteinChange$Alt)
), ]

NGS_proteinChange$calculatedconclusion <- ifelse(
  is.na(NGS_proteinChange$calculatedconclusion), 
  '', NGS_proteinChange$calculatedconclusion
)
NGS_proteinChange <- unique(NGS_proteinChange)

NGS_proteinChange  <- NGS_proteinChange %>%
  filter(grepl('Pathogenic', ORG_Result))

# write.table(NGS_proteinChange, file = paste0('ngs_result_', now, '.txt'), sep = '\t', row.names = F, quote = F)

write.table(NGS_proteinChange, file = output_file, sep = '\t', row.names = F, quote = F)
