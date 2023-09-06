# setwd("/Users/analisamilkey/Documents/projects/SMC-memory-efficient/SMC/test2")

filenames <- Sys.glob("*log")
gene_number <- 0

# read in data
for (file_name in filenames) {
  if (startsWith(file_name, "g")) {
    gene_number <- gene_number + 1
    name <- paste("gene", gene_number, sep="")
    name <- read.table(file = file_name, sep = '\t', header = TRUE)
    if (gene_number == 1) {
      concat <- name
      n_gene_iter <- nrow(name)
    }
    else {
      concat <- merge(name, concat, by = "iter")
    }
  }
  else {
    species <- read.table(file = file_name, sep = '\t', header = TRUE)
    # not sure how to deal with extra particles in species trees - for now, remove extra iterations
    if (n_gene_iter == 0) {
      exit()
    }
    species<-species[1:(n_gene_iter),]
    concat <- merge(species, concat, by = "iter")
  }
}

write.table(concat, file='output.log', quote=FALSE, sep='\t', col.names = NA)
