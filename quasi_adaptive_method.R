
library(data.table)
args <- commandArgs(trailingOnly = TRUE)

infn1 = args[1] #Intersnps results-based on N significant loci from primary GWAS.
infn2 = args[2] #Conditional analysis results. Output from GCTA software was included in Rscript with necessary columns including rsid of candidate SNP with corresponding conditional P-value
outfn = args[3] #output file name, Signal name: secondary_signal, tersiary_signal, signal_of_4th, ...

#The results if Intersnps remain constant when determining secondary signal, tersiaty signal and multiple independent signals from higher order.
#The results of conditional analysis are changed depending on which signal is determined.



boundery <- 5.00E-08
r2=0.3 #optimal value for r2
wd=1 #wd=1 weights based on distance-step-wise-strong
k=5
alpha=0.05 #alpha is a type one error rate



print("Reading data...")

interSNPs <- fread(infn1)
GCTA <- fread(infn2)

colnames(GCTA) <- c("Chr","candidate_SNP","bp","refA","freq","b","se","p","n","freq_geno","bC","bC_se","pC")

interSNPs <- interSNPs[,c(1:8,18)]
colnames(interSNPs) <- c("chr","index","index_allele1","index_allele2","candidate_SNP","candidate_SNP_allele1","candidate_SNP_allele2","D","R2")


data1 <- merge(interSNPs,GCTA,by.x ="candidate_SNP",by.y= "candidate_SNP")
data2 <- data1[,c(2:5,1,6:9,11:21)]
data <- data2[-which(is.na(data2$pC)),]


#Assign pre-weight based on r2 (w_r)

data$w_r <- (1-abs(data$R2-0.3)-0.3)/(1-0.3)

#Assign pre-weight based on distance (w_d)

data$w_d[abs(data$D) > 0 & abs(data$D) <= 1000] <- 1
data$w_d[abs(data$D) > 1000 & abs(data$D) <= 10000] <- 0.5
data$w_d[abs(data$D) > 10000 & abs(data$D) <= 50000] <- 0.25
data$w_d[abs(data$D) > 50000 & abs(data$D) <= 100000] <- 0.125
data$w_d[abs(data$D) > 100000 & abs(data$D)  <= 500000] <- 0.0625
data$w_d[abs(data$D) > 500000 & abs(data$D) <= 2500000] <- 0.03125
data$w_d[abs(data$D) > 2500000 & abs(data$D) <= 5000000] <- 0.015625
#data$w_d[abs(data$D) > 2500000 & abs(data$D) <= max(abs(data$D))] <- 0.015625

#Combine w_d and w_r to get the final weight (w)
data$mwd <- ((data$w_d^k)*data$w_r)^(1/(k+1))

#m is the number of multiple tests (number of the candidate SNPs across N independent loci).
m <- dim(data)[1]
print(m)
#W is the optimal weight
S <- sum(data$mwd)
data$W <- ((data$mwd)*m)/S


data$SNP_specifc_alpha <- 1-(1-alpha)^{data$W/m}

data <- data[order(data$chr,data$bp),]
write.table(data,"all_SNP_specifc_alpha.txt", row.names=F, sep="\t", quote=F)

last.result <- data[data$pC <= data$SNP_specifc_alpha,]
################################################################################################################
cat("Determine one",outfn,"per locus.\n")


signal <- data.frame(matrix(rep(NA,length(unique(last.result$index))*dim(last.result)[2]),length(unique(last.result$index)),dim(last.result)[2]))
colnames(signal) <- colnames(last.result)
rownames(signal) <- unique(last.result$index)


for (i in unique(last.result$index)){
A <- last.result[which(last.result$index==i),]
B <- A[which(A[,"pC"]==min(A[,"pC"])),]
if (dim(B)[1]==1){signal[i,] <- A[which(A[,"pC"]==min(A[,"pC"])),]}
{signal[i,] <-  B[which(B[,"D"]==min(B[,"D"])),]
}
}


final_signal <- signal[order(signal$chr,signal$bp),]

write.table(final_signal,paste0(outfn,".txt"), row.names=F, sep="\t", quote=F)

###########################################################################################################

helper2<-subset(data,pC<= boundery)
SNPs_establishedMethod2<-unique(helper2$index)


new2<-subset(last.result,!is.element(last.result$index,SNPs_establishedMethod2))
SNPs_newMethod_2<-unique(new2$index) # these are the index SNPs for which new secondary SNPs were found
hits2<-length(SNPs_newMethod_2)



cat(length(unique(last.result$index)),"regions",outfn,"with conditional p smaller than SNP_specifc_alpha threshold.",hits2,
"regions identified exclusively by quasi-adaptive method.\n")































