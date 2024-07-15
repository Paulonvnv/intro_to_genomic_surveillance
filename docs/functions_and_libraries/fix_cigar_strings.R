
cigar_files = file.path("~/Documents/Github/intro_to_genomic_surveillance/docs/data/Pfal_example/cigar_tables",
                        list.files("~/Documents/Github/intro_to_genomic_surveillance/docs/data/Pfal_example/cigar_tables"))

for(file in cigar_files){
  cigar_file = read.table(file,
                          header = T,
                          check.names = F)
  
  
  cigar_table_cols = names(cigar_file)
  
  for(amplicon in markers$amplicon){
    
    cigar_table_cols[grepl(paste0(amplicon, '\\.'), cigar_table_cols)] = gsub(paste0(amplicon, '\\.'), paste0(amplicon, ','), cigar_table_cols[grepl(paste0(amplicon, '\\.'), cigar_table_cols)])
    
  }
  
  cigar_table_cols = gsub('I\\.', 'I=', cigar_table_cols)
  cigar_table_cols = gsub('D\\.', 'D=', cigar_table_cols)
  
  names(cigar_file) = cigar_table_cols
  
  
  write.table(cigar_file, 
              file,
              sep = '\t',
              quote = F,
              row.names = F
  )
  
}


