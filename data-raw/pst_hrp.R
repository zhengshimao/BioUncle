## code to prepare `hrp_gene.fasta` dataset goes here
# pst_hrp <- refdb::refdb_import_NCBI(query = "Pseudomonas syringae pv. tomato hrp gene ",
#                                     full = TRUE,
#                                     max_seq_length = 1000,
#                                     seq_bin = 200,
#                                     verbose = TRUE
# )
# pst_hrp <- pst_hrp %>% dplyr::select(id,sequence,species,length) %>%
#   tidyr::unite(col = "description",species,length,sep = " ",remove = TRUE, na.rm = FALSE)
pst_hrp <- read_bio_fasta("data-raw/pst_hrp_gene.fasta")
usethis::use_data_row(pst_hrp, overwrite = TRUE)
