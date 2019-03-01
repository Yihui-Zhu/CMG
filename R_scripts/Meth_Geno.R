#usage: Meth_Geno(sample_info_file_name, VCF_dir, Methylation_file); 
#e.g., Meth_Geno(Genotype_finalize, "Meth_Smooth_Adjust.csv")

Meth_Geno <-
  function(Genotype_finalize, Methylation_file) {
    
library(ggplot2)
library(ggcorrplot)
library(reshape2)

Methylation = read.csv(Methylation_file)

col_idx = vector(length = ncol(Methylation))
col_idx[1:3] = TRUE

for (i in 4:ncol(Methylation)){
  if (colnames(Methylation)[i] %in% colnames(Genotype_finalize)){
    col_idx[i] = TRUE
  }
}
Methylation = Methylation[,col_idx]


Meth_N = nrow(Methylation)
Geno_N = nrow(Genotype_finalize)

mymat <- matrix(0, nrow=Meth_N, ncol=Geno_N)
rownames(mymat) <- c(rownames(Methylation))
colnames(mymat) <- c(rownames(Genotype_finalize))
for (i in 1:Meth_N){
  rownames(mymat)[i] = paste("M", toString(Methylation[i,1]), toString(Methylation[i,2]), toString(Methylation[i,3]), sep = "_")
  for (j in 1:Geno_N){
    if (i == 1){
      colnames(mymat)[j] = paste("G", toString(Genotype_finalize[j,1]), toString(Genotype_finalize[j,2]), sep = "_")
    }
    if (toString(Genotype_finalize[j,1]) == toString(Methylation[i,1])){
      if ((as.integer(toString(Genotype_finalize[j,2])) >= as.integer(Methylation[i,2])) 
          & (as.integer(toString(Genotype_finalize[j,2])) <= as.integer(Methylation[i,3]))){
        mymat[i,j] = 1
      }
    }
  }
}

rowSums(mymat)
write.csv(mymat,"Meth_Geno_count.csv")

## get the p value
mymat_p <- matrix(1, nrow=Meth_N, ncol=Geno_N)
rownames(mymat_p) <- c(rownames(Methylation))
colnames(mymat_p) <- c(rownames(Genotype_finalize))

for (i in 1:Meth_N){
  rownames(mymat_p)[i] = paste("M", toString(Methylation[i,1]), toString(Methylation[i,2]), toString(Methylation[i,3]), sep = "_")
  for (j in 1:Geno_N){
    if (i == 1){
      colnames(mymat_p)[j] = paste("G", toString(Genotype_finalize[j,1]), toString(Genotype_finalize[j,2]), sep = "_")
    }
    if (toString(Genotype_finalize[j,1]) == toString(Methylation[i,1])){
      if ((as.integer(toString(Genotype_finalize[j,2])) >= as.integer(Methylation[i,2])) 
          & (as.integer(toString(Genotype_finalize[j,2])) <= as.integer(Methylation[i,3]))){
        if (length(unique(unlist(Genotype_finalize[j,(3:ncol(Genotype_finalize))]) )) > 1){
          mymat_p[i, j] = summary(aov(as.numeric(Methylation[i,(4:ncol(Methylation))]) ~ unlist(Genotype_finalize[j,(3:ncol(Genotype_finalize))])           ))[[1]][["Pr(>F)"]][1]
        }
      }
    }
  }
}

write.csv(mymat_p,"Meth_Geno_p_value.csv")

## Ploting
melted_cormat <- melt(mymat_p)
head(melted_cormat)

melted_cormat_sig <- melted_cormat[ which(melted_cormat$value < 0.05), ]
dim(melted_cormat_sig)

pdf("Heapmap_cis_sig.pdf", width = 20, height = 20, family = "Helvetica")
ggplot(data = melted_cormat_sig, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "red", high = "blue", mid = "white", 
                       midpoint = 0.05,  space = "Lab", 
                       name="p value") +
  geom_tile(colour="gray80") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 12, hjust = 1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(size = 20),axis.title.x = element_text(size = 16),axis.title.y = element_text(size = 16))+
  ylab("DMRs from WGS")+ggtitle("WGS and WGBS")+xlab("SNPs from WGBS")+
  coord_fixed()+
  theme( plot.margin = margin(1, 1, 1, 1, "cm")  )
dev.off()
print("Done Meth_Geno !")
}
