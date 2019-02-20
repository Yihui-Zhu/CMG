SNP_WGS <- function(x,y){
  
  Sample_Chrom = as.character(x$CHROM)
  Sample_Pos = as.character(x$LOC)
  Sample_REF_Nucleotide = substr(x$REF, start=1, stop=1)
  Sample_REF_Count = substr(x$REF, start=3, stop=3)
  Sample_ALT_Nucleotide = substr(x$ALT, start=1, stop=1)
  Sample_ALT_Count = substr(x$ALT, start=3, stop=3)
  Sample_genotype = as.character(x$N_alleles)
  A = cbind(Sample_Chrom, Sample_Pos, Sample_genotype, Sample_REF_Nucleotide, Sample_REF_Count, Sample_ALT_Nucleotide, Sample_ALT_Count)
  head(A)
  dim(combined_snp_file)
  head(combined_snp_file)
  Combine_Chrom = as.character(combined_snp_file$CHROM)
  Combine_Pos = as.character(combined_snp_file$LOC)
  Combine_REF_Nucleotide = substr(combined_snp_file$REF, start=1, stop=1)
  Combine_REF_Count = substr(combined_snp_file$REF, start=3, stop=3)
  Combine_ALT_Nucleotide = substr(combined_snp_file$ALT, start=1, stop=1)
  Combine_ALT_Count = substr(combined_snp_file$ALT, start=3, stop=3)
  Combine_Nalleles = as.character(combined_snp_file$N_alleles)
  B = cbind(Combine_Chrom, Combine_Pos, Combine_Nalleles, Combine_REF_Nucleotide, Combine_REF_Count, Combine_ALT_Nucleotide, Combine_ALT_Count)
  head(B)
  
  dim(combined_snp_file)
  C = unique(combined_snp_file)
  dim(C)
  
  head(combined_snp_file)
  Combine_Chrom = as.character(combined_snp_file$CHROM)
  Combine_Pos = as.character(combined_snp_file$LOC)
  Combine_REF_Nucleotide = substr(combined_snp_file$REF, start=1, stop=1)
  Combine_REF_Count = substr(combined_snp_file$REF, start=3, stop=3)
  Combine_ALT_Nucleotide = substr(combined_snp_file$ALT, start=1, stop=1)
  Combine_ALT_Count = substr(combined_snp_file$ALT, start=3, stop=3)
  Combine_Nalleles = as.character(combined_snp_file$N_alleles)
  B = cbind(Combine_Chrom, Combine_Pos, Combine_Nalleles, Combine_REF_Nucleotide, Combine_REF_Count, Combine_ALT_Nucleotide, Combine_ALT_Count)
  head(B)
  D = unique(B)
  dim(D)
  head(D)
  
  z = rep(0,dim(D)[1])
  E= cbind(D, z)
  head(E)
  dim(E)
  
  m = dim(A)[1]
  
  for (i in 1:m){
    E[which(A[i,1]==E[,1] & A[i,2]==E[,2] &A[i,3]==E[,3] &A[i,4]==E[,4] &A[i,5]==E[,5] &A[i,6]==E[,6] &A[i,7]==E[,7]) ,8] =1
  }
  summary(E)
  head(E)
  head(A)
  
  for (i in 1:m){
    E[which(A[i,1]==E[,1] & A[i,2]==E[,2] & E[,8]!=1), 8] = 2
  }
  head(E)
  summary(E)
  dim(E)
  
  n = dim(E)[1]
  genotype = c()
  
  for (j in 1:n){
    if (E[j,8] == 0){
      E[j,8] = paste(E[j,4],E[j,4],sep="") 
    }
    if (E[j,8] == 1){
      if (E[j,5] ==1){
        E[j,8] = paste(E[j,4], E[j,6], sep="")
      } 
      if (E[j,7] ==2){
        E[j,8] = paste(E[j,6], E[j,6], sep="")
      }
    }
    if (E[j,8] == 2){
      E[j,8] = NA
    }
  }
  E=na.omit(E)
  
  I = as.character(y)
  
  colnames(E)[8] <- I
  K = paste (I,"_snp_datastory",sep = "")
  K = as.data.frame(K)
  K = E
  E <- unique( E[ ,c(1,2,8)] )
  
  name_snp =paste(I,"_snp.csv",sep="")
  
  print(dim(E))
  print(summary(E))
  print(K)
  print(name_snp)
  write.csv(E,file=name_snp)
}
