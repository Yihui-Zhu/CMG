#usage: Genotype(sample_info_file_name, VCF_dir); e.g., Genotype("sample_information.csv", "1120_frq.count_files/")

Genotype <-
  function(sample_info_file_name, VCF_dir) {

library(dplyr)
sample_information = read.csv(sample_info_file_name, header=T)

options(warn=-1)

cn = colnames(sample_information)
# n, number of samples
n = length(sample_information[,1])

listOfDataFrames <- vector(mode = "list", length = n)

# combine all people's data into combined_snp_file
for (k in 1:n){
  listOfDataFrames[[k]] <- read.delim(paste(VCF_dir, sample_information[k,2], sep = ""), header=FALSE)[-1,]
  listOfDataFrames[[k]]$V7 = NA
  if (length(which(listOfDataFrames[[k]]$V3 == 3)) > 0){
    listOfDataFrames[[k]]$V7[which(listOfDataFrames[[k]]$V3 == 3)] = toString(listOfDataFrames[[k]]$V1[which(listOfDataFrames[[k]]$V3 == 3) + 1])
    listOfDataFrames[[k]] = listOfDataFrames[[k]][-(which(listOfDataFrames[[k]]$V3 == 3) + 1),]
  }
}
combined_snp_file = bind_rows(listOfDataFrames)

# rename colnames
colnames(combined_snp_file) = c("CHROM","LOC","N_alleles","N_chrom","REF","ALT","ALT2")


n1 = dim(unique(combined_snp_file[,1:2]))[1]
n2 = n + dim(unique(combined_snp_file[,1:2]))[2]

M = matrix(0, nrow=n1, ncol=n2)
M[,1] = unique(combined_snp_file[,1:2])[,1]
M[,2] = unique(combined_snp_file[,1:2])[,2]

for (k in 1:length(combined_snp_file$REF)){
  ref_temp = substr(toString(combined_snp_file$REF[k]), start = 1, stop = 1)
  combined_snp_file$REF[k] = paste(ref_temp,ref_temp, sep = "")
}

temp = vector(mode = 'character',length = n1)
temp = unique(combined_snp_file[,c(1,2,5)])[,3]

# function: genotype
# arguments: ref, reference allele; alt, alternative allele; alt2, second alternative allele. e.g. ref = C:0, alt = A:1, alt2 = G:1.
# output: String of genotype, e.g. "AG"
genotype = function(ref, alt, alt2){
  ref = toString(ref)
  alt = toString(alt)
  alt2 = toString(alt2)
  n_ref = substr(ref, start = 3, stop = 3)
  n_alt = substr(alt, start = 3, stop = 3)
  n_alt2 = substr(alt2, start = 3, stop = 3)
  g_ref = substr(ref, start = 1, stop = 1)
  g_alt = substr(alt, start = 1, stop = 1)
  g_alt2 = substr(alt2, start = 1, stop = 1)
  if (alt2 != "NA"){
    result = paste(g_alt, g_alt2, sep = "")
  } else {
    if (n_ref == 1){
      result = paste(g_ref, g_alt, sep = "")
    } else { 
      result = paste(g_alt, g_alt, sep = "")
    }
  }
  return(result)
}

# Assemble matrix M
for (k in 1:n){
  listOfDataFrames[[k]]
  for (p in 1:n1){
    idx = which(as.character(listOfDataFrames[[k]]$V1) == M[p,1] & as.character(listOfDataFrames[[k]]$V2) == M[p,2])
    if (length(idx) > 0){
      M[p, k+2] = genotype(listOfDataFrames[[k]][idx,5], listOfDataFrames[[k]][idx,6], listOfDataFrames[[k]][idx,7])
    } else {
      M[p, k+2] = temp[p]
    }
  }
}

Genotype_finalize = as.data.frame(M)

colnames(Genotype_finalize)[1] = "CHROM"
colnames(Genotype_finalize)[2] = "LOC"
for (k in 1:n){
  colnames(Genotype_finalize)[2+k] = toString(sample_information[k,1])
}

# output csv file, i.e., Genotype_finalize.csv
write.csv(Genotype_finalize, "Genotype_finalize.csv")
print("Done!")

}






