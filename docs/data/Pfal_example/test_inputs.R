
library(argparse)
library(stringr)
library(rmarkdown)

wd = '~/Documents/Github/MHap-Analysis/docs/data/Pfal_example/'
wd = gsub('/$', '', wd)
# Tools or functions directory
fd = '~/Documents/Github/MHap-Analysis/docs/functions_and_libraries'
fd = gsub('/$', '', fd)
# Reference files directory
rd = '~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/'
rd = gsub('/$', '', rd)


output =  "MHap_test"

# Starting vcf file
cigar_paths = "null"
cigar_paths = if(cigar_paths == 'null'){NULL}else{file.path(wd, cigar_paths)}
# Starting vcf file

cigar_files = "cigar_tables"
cigar_files = if(cigar_files == 'null'){NULL}else{file.path(wd, cigar_files)}

# Starting vcf file
ampseq_jsonfile = "null"
ampseq_jsonfile = if(ampseq_jsonfile == 'null'){NULL}else{file.path(wd, ampseq_jsonfile)}

# Starting vcf file
ampseq_excelfile = "null"
ampseq_excelfile = if(ampseq_excelfile == 'null'){NULL}else{file.path(wd, args$ampseq_excelfile)}

# Pattern to diferentiate between samples of interest and controls
sample_id_pattern = "^ID"

# csv table with markers information
markers = 'markers.csv'
markers = if(markers == 'null'){NULL}else{file.path(rd, markers)}

# Minimum abundance to call an allele
min_abd = 10

# Minimum ratio to call minor alleles
min_ratio = 0.1


off_target_formula = "dVSITES_ij>=0.3"
off_target_formula = gsub('"',"",off_target_formula)

off_target_formula = gsub('&'," & ",off_target_formula, ignore.case = TRUE)
off_target_formula = gsub('\\|'," \\| ",off_target_formula, ignore.case = TRUE)

if(grepl("\\w>\\d",off_target_formula)){
  patterns = str_extract_all(off_target_formula, "\\w>\\d")[[1]]
  
  for(pattern in patterns){
    
    replacement = gsub('>',' > ',pattern)
    off_target_formula = gsub(pattern,
                              replacement,
                              off_target_formula, ignore.case = TRUE)
  }
  
}

if(grepl("\\w<\\d",off_target_formula)){
  patterns = str_extract_all(off_target_formula, "\\w<\\d")[[1]]
  
  for(pattern in patterns){
    
    replacement = gsub('<',' < ',pattern)
    off_target_formula = gsub(pattern,
                              replacement,
                              off_target_formula, ignore.case = TRUE)
  }
  
}

off_target_formula = gsub('>=', " >= ", off_target_formula, ignore.case = TRUE)
off_target_formula = gsub('<=', " <= ", off_target_formula, ignore.case = TRUE)
off_target_formula = gsub('==', " == ", off_target_formula, ignore.case = TRUE)
off_target_formula = gsub('!=', " != ", off_target_formula, ignore.case = TRUE)

off_target_formula = gsub('\\+', " \\+ ", off_target_formula, ignore.case = TRUE)
off_target_formula = gsub('-', " - ", off_target_formula, ignore.case = TRUE)
off_target_formula = gsub('\\*', " \\* ", off_target_formula, ignore.case = TRUE)
off_target_formula = gsub('/', " / ", off_target_formula, ignore.case = TRUE)


print(paste0('off_target_formula: ', off_target_formula))


flanking_INDEL_formula = "flanking_INDEL==TRUE&h_ij>=0.66"
flanking_INDEL_formula = gsub('"',"",flanking_INDEL_formula)

flanking_INDEL_formula = gsub('&'," & ",flanking_INDEL_formula, ignore.case = TRUE)
flanking_INDEL_formula = gsub('\\|'," \\| ",flanking_INDEL_formula, ignore.case = TRUE)

if(grepl("\\w>\\d",flanking_INDEL_formula)){
  patterns = str_extract_all(flanking_INDEL_formula, "\\w>\\d")[[1]]
  
  for(pattern in patterns){
    
    replacement = gsub('>',' > ',pattern)
    flanking_INDEL_formula = gsub(pattern,
                                  replacement,
                                  flanking_INDEL_formula, ignore.case = TRUE)
  }
  
}

if(grepl("\\w<\\d",flanking_INDEL_formula)){
  patterns = str_extract_all(flanking_INDEL_formula, "\\w<\\d")[[1]]
  
  for(pattern in patterns){
    
    replacement = gsub('<',' < ',pattern)
    flanking_INDEL_formula = gsub(pattern,
                                  replacement,
                                  flanking_INDEL_formula, ignore.case = TRUE)
  }
  
}

flanking_INDEL_formula = gsub('>='," >= ",flanking_INDEL_formula, ignore.case = TRUE)
flanking_INDEL_formula = gsub('<='," <= ",flanking_INDEL_formula, ignore.case = TRUE)
flanking_INDEL_formula = gsub('=='," == ",flanking_INDEL_formula, ignore.case = TRUE)
flanking_INDEL_formula = gsub('!='," != ",flanking_INDEL_formula, ignore.case = TRUE)

flanking_INDEL_formula = gsub('\\+'," \\+ ",flanking_INDEL_formula, ignore.case = TRUE)
flanking_INDEL_formula = gsub('-'," - ",flanking_INDEL_formula, ignore.case = TRUE)
flanking_INDEL_formula = gsub('\\*'," \\* ",flanking_INDEL_formula, ignore.case = TRUE)
flanking_INDEL_formula = gsub('/'," / ",flanking_INDEL_formula, ignore.case = TRUE)


print(paste0('flanking_INDEL_formula: ', flanking_INDEL_formula))


PCR_errors_formula = "h_ij>=0.66&h_ijminor>=0.66"
PCR_errors_formula = gsub('"',"",PCR_errors_formula)

PCR_errors_formula = gsub('&'," & ",PCR_errors_formula, ignore.case = TRUE)
PCR_errors_formula = gsub('\\|'," \\| ",PCR_errors_formula, ignore.case = TRUE)

if(grepl("\\w>\\d",PCR_errors_formula)){
  patterns = str_extract_all(PCR_errors_formula, "\\w>\\d")[[1]]
  
  for(pattern in patterns){
    
    replacement = gsub('>',' > ',pattern)
    PCR_errors_formula = gsub(pattern,
                              replacement,
                              PCR_errors_formula, ignore.case = TRUE)
  }
  
}

if(grepl("\\w<\\d",PCR_errors_formula)){
  patterns = str_extract_all(PCR_errors_formula, "\\w<\\d")[[1]]
  
  for(pattern in patterns){
    
    replacement = gsub('<',' < ',pattern)
    PCR_errors_formula = gsub(pattern,
                              replacement,
                              PCR_errors_formula, ignore.case = TRUE)
  }
  
}

PCR_errors_formula = gsub('>='," >= ", PCR_errors_formula, ignore.case = TRUE)
PCR_errors_formula = gsub('<='," <= ", PCR_errors_formula, ignore.case = TRUE)
PCR_errors_formula = gsub('=='," == ", PCR_errors_formula, ignore.case = TRUE)
PCR_errors_formula = gsub('!='," != ", PCR_errors_formula, ignore.case = TRUE)

PCR_errors_formula = gsub('\\+'," \\+ ", PCR_errors_formula, ignore.case = TRUE)
PCR_errors_formula = gsub('-'," - ", PCR_errors_formula, ignore.case = TRUE)
PCR_errors_formula = gsub('\\*'," \\* ", PCR_errors_formula, ignore.case = TRUE)
PCR_errors_formula = gsub('/'," / ", PCR_errors_formula, ignore.case = TRUE)


print(paste0('PCR_errors_formula: ', PCR_errors_formula))

PerformanceReport = TRUE

# sample_ampl_rate
sample_ampl_rate = 0.5

# locus_ampl_rate
locus_ampl_rate  = 0.5

Drug_Surveillance_Report = TRUE
Variants_of_Interest_Report = FALSE

# Reference gff and fasta files, require to translate DNA cigar formats to aminoacid sequences

ref_gff = "PlasmoDB-59_Pfalciparum3D7.gff"
ref_gff = if(ref_gff == 'null'){NULL}else{file.path(rd, ref_gff)}

ref_fasta = "PlasmoDB-59_Pfalciparum3D7_Genome.fasta"
ref_fasta = if(ref_fasta == 'null'){NULL}else{file.path(rd, ref_fasta)}


reference_alleles = "drugR_alleles.csv"
reference_alleles = if(reference_alleles == 'null'){NULL}else{file.path(rd, reference_alleles)}


# gene_names
gene_names = "PfDHFR,PfMDR1,PfDHPS,PfKelch13,PF3D7_1447900"
if(gene_names == 'null'){
  gene_names = NULL
}else{
  gene_names = strsplit(gene_names, ',')[[1]]
}

# gene_ids
gene_ids = "PF3D7_0417200,PF3D7_0523000,PF3D7_0810800,PF3D7_1343700,PF3D7_1447900"
if(gene_ids == 'null'){
  gene_ids = NULL
}else{
  gene_ids = strsplit(gene_ids, ',')[[1]]
}


# ibd_thres
ibd_thres = 0.99
ibd_thres = if(ibd_thres == 'null'){NULL}else{as.numeric(ibd_thres)}

pairwise_relatedness_table = 'pairwise_relatedness.csv'
pairwise_relatedness_table = if(pairwise_relatedness_table == 'null'){NULL}else{file.path(wd, pairwise_relatedness_table)}

nChunks = 500

# metadata
metadata_file = "Pfal_metadata.csv"
metadata_file = if(metadata_file == 'null'){NULL}else{file.path(wd, metadata_file)}

# join metadata by
join_by = "Sample_id"
join_by = if(join_by == 'null'){NULL}else{join_by}

# pop
Variable1 = "Subnational_level2"
Variable1 = if(Variable1 == 'null'){NULL}else{Variable1}

# temporal_population
Variable2 = "Quarter_of_Collection"
Variable2 = if(Variable2 == 'null'){NULL}else{Variable2}


Longitude = "Longitude"
Longitude = if(Longitude == 'null'){NULL}else{Longitude}

Latitude = "Latitude"
Latitude = if(Latitude == 'null'){NULL}else{Latitude}


na_var_rm = FALSE
na_hap_rm = FALSE

drugs = "Artemisinin,Chloroquine,Pyrimethamine,Sulfadoxine,Lumefantrine,Mefloquine"
if(drugs == 'null'){
  drugs = NULL
}else{
  drugs = strsplit(drugs, ',')[[1]]
}

include_all_drug_markers = TRUE

var_filter = 'Subnational_level2;Keep;Municipality.1,Municipality.2,Municipality.3/Quarter_of_Collection;Keep;2021-Q1,2021-Q2,2021-Q3,2021-Q4,2022-Q1,2022-Q2'
#var_filter = 'null'
if(var_filter == 'null'){
  var_filter = NULL
}else{
  var_filter = gsub('\\.', ' ', var_filter)
  var_filter = strsplit(var_filter, '/')[[1]]
}
var_filter


# parallel
parallel = TRUE

# nTasks

nTasks = 1

# Task_id
Task_id = 0



poly_quantile = 0.75

poly_formula = "NHetLoci >= 1 & Fws < 1"

hap_color_palette = 'random'

# ampseq_object name
# if(nTasks > 1){
#   ampseq_object_name = paste0(output, '_ampseq_Chunk', Task_id)
#
# }else{
#   ampseq_object_name = paste0(output, '_ampseq')
#
# }

# R image name

if(nTasks > 1){
  imagename = paste0('Chunks/',output, '_ampseq_Chunk', Task_id, '.RData')

}else{
  imagename = paste0(output, '_ampseq.RData')

}





#ibd_step = 'merge'






# {
#   "vmem": 32,
#   "cores": 8,
#   "h_rt1": 05:00:00,
#   "h_rt2": 20:00:00,
#   "h_rt3": 05:00:00,
#   "nTasks": 50,
#   
#   "wd": "/gsap/garage-protistvector/",
#   "fd": "/gsap/garage-protistvector/",
#   "rd": "/gsap/garage-protistvector/",
#   
#   "cigar_paths": "/gsap/garage-protistvector/",
#   "cigar_files": "/gsap/garage-protistvector/",
#   "ampseq_jsonfile": "/gsap/garage-protistvector/",
#   "ampseq_excelfile": "/gsap/garage-protistvector/",
#   
#   "output": "MHap_test",
#   
#   "sample_id_pattern": ".",
#   "markers": NaN,
#   "min_abd": NaN,
#   "min_ratio": NaN,
#   "PerformanceReport": false,
#   
#   "sample_ampl_rate": NaN,
#   "locus_ampl_rate": NaN,
#   
#   "ref_gff": "genes.gff",
#   "ref_fasta": "ref_fasta.fasta",
#   
#   "ibd_thres": 0.99,
#   "parallel": true,
#   "ibd_ncol": 4,
#   "pop_levels":null,
#   
#   "nchunks": 500,
#   
#   "metadata": "Pviv_ColPerVen_metadata.csv", 
#   "join_by": "Sample_id",
#   "geographic_population": "Country",
#   "temporal_population": "Country"
#   
# }
# 
# 
# 
# metadata %<>% mutate(
#   Latitude = case_when(
#     Subnational_level2 == 'Municipality 1' ~ -76.65806,
#     Subnational_level2 == 'Municipality 2' ~ -77.03116,
#     Subnational_level2 == 'Municipality 3' ~ -77.883,
#     Subnational_level2 == 'Municipality 4' ~ -78.764725,
#     Subnational_level2 == 'Municipality 5' ~ -67.92389
#   ),
#   Longitude = case_when(
#     Subnational_level2 == 'Municipality 1' ~ 5.69222,
#     Subnational_level2 == 'Municipality 2' ~ 3.8801,
#     Subnational_level2 == 'Municipality 3' ~ 2.567,
#     Subnational_level2 == 'Municipality 4' ~ 1.806667,
#     Subnational_level2 == 'Municipality 5' ~ 3.86528
#   )
# )
# 
# write.csv(metadata, file.path(wd, 'Pfal_metadata.csv'), quote = F, row.names = F)
#  
# library(tmaptools)
# library(sf)
# 
# data(NLD_prov)
# 
# 
# origin_data <- NLD_prov %>% 
#   st_set_geometry(NULL) %>% 
#   select(name, origin_native, origin_west, origin_non_west) %>% 
#   gather(key=origin, value=perc, origin_native, origin_west, origin_non_west, factor_key=TRUE)
# 
# plot_cols <- get_brewer_pal("Dark2", 3)
# 
# grobs2 <- lapply(split(origin_data, origin_data$name), function(x) {
#   ggplotGrob(ggplot(x, aes(x="", y=-perc, fill=origin)) +
#                geom_bar(width=1, stat="identity") +
#                coord_polar("y", start=0) +
#                scale_y_continuous(expand=c(0,0)) +
#                scale_fill_manual(values=plot_cols) +
#                theme_ps(plot.axes = FALSE))
# })
# 
# 
# tm_shape(NLD_prov) +
#   tm_polygons(group = "Provinces") +
#   tm_symbols(shape="name", 
#              shapes= grobs2, 
#              border.lwd = 0,
#              legend.shape.show = FALSE, 
#              legend.size.is.portrait = TRUE, 
#              shapes.legend = 22)+
#   tm_add_legend(type="fill", 
#                 col=plot_cols, 
#                 labels=c("Native", "Western", "Non-western"), 
#                 title="Origin")
