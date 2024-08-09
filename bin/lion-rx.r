# Rx_identifier
# Based on the ratio of X chromosome-derived shotgun sequencing data to the autosomal coverage to establish the probability of an XX or XY karyotype for ancient samples.
# Author: Chuan-Chao Wang
# Contact: wang@shh.mpg.de, chuan-chao_wang@hms.harvard.edu; Department of Genetics, Harvard Medical School; Department of Archaeogenetics, Max Planck Institute for the Science of Human History.
# Date: 30 Jan, 2016
# example:samtools view -q 30 -b Ajv52.hs37d5.fa.merged.bam > Ajv52_q30.bam
        # samtools index Ajv52_q30.bam
        # samtools idxstats Ajv52_q30.bam > Ajv52.idxstats
        # Rscript Rx_identifier.r Ajv52 > Ajv52.Rx
		
# Original script published under a Creative Commons Attribution License (CC BY 4.0)
# DOI: 10.1371/journal.pone.0163019.s003

# Modified for lions by M.G. Campana, 31 Jan 2024
# Using samtools coverage output instead of idxstats
# Need to add in the Y apparently...

args=(commandArgs(TRUE))
FILENAME=as.character(args[1])

idxstats<-read.table(FILENAME,header=F,row.names=1)
c1 <- c(as.numeric(idxstats[,2][1:19]), as.numeric(idxstats[,2][97]))
c2 <- c(as.numeric(idxstats[,3][1:19]), as.numeric(idxstats[,3][97]))
total_ref <- sum(c1)
total_map <- sum(c2)

if (total_map < 1000) {print("Fewer than 1000 reads. Quitting.")
} else {
  
LM <- lm(c1~c2)
summary(LM)  
  
Rt1 <- (idxstats[1,3]/total_map)/(idxstats[1,2]/total_ref)
Rt2 <- (idxstats[2,3]/total_map)/(idxstats[2,2]/total_ref)
Rt3 <- (idxstats[3,3]/total_map)/(idxstats[3,2]/total_ref)
Rt4 <- (idxstats[4,3]/total_map)/(idxstats[4,2]/total_ref)
Rt5 <- (idxstats[5,3]/total_map)/(idxstats[5,2]/total_ref)
Rt6 <- (idxstats[6,3]/total_map)/(idxstats[6,2]/total_ref)
Rt7 <- (idxstats[7,3]/total_map)/(idxstats[7,2]/total_ref)
Rt8 <- (idxstats[8,3]/total_map)/(idxstats[8,2]/total_ref)
Rt9 <- (idxstats[9,3]/total_map)/(idxstats[9,2]/total_ref)
Rt10 <- (idxstats[10,3]/total_map)/(idxstats[10,2]/total_ref)
Rt11 <- (idxstats[11,3]/total_map)/(idxstats[11,2]/total_ref)
Rt12 <- (idxstats[12,3]/total_map)/(idxstats[12,2]/total_ref)
Rt13 <- (idxstats[13,3]/total_map)/(idxstats[13,2]/total_ref)
Rt14 <- (idxstats[14,3]/total_map)/(idxstats[14,2]/total_ref)
Rt15 <- (idxstats[15,3]/total_map)/(idxstats[15,2]/total_ref)
Rt16 <- (idxstats[16,3]/total_map)/(idxstats[16,2]/total_ref)
Rt17 <- (idxstats[17,3]/total_map)/(idxstats[17,2]/total_ref)
Rt18 <- (idxstats[18,3]/total_map)/(idxstats[18,2]/total_ref)
Rt19 <- (idxstats[19,3]/total_map)/(idxstats[19,2]/total_ref)

tot <- c(Rt19/Rt1,Rt19/Rt2,Rt19/Rt3,Rt19/Rt4,Rt19/Rt5,Rt19/Rt6,Rt19/Rt7,Rt19/Rt8,Rt19/Rt9,Rt19/Rt10,Rt19/Rt11,Rt19/Rt12,Rt19/Rt13,Rt19/Rt14,Rt19/Rt15,Rt19/Rt16,Rt19/Rt17,Rt19/Rt18)
Rx <- mean(tot)
cat("Rx :",Rx,"\n")
confinterval <- 1.96*(sd(tot)/sqrt(18))
CI1 <- Rx-confinterval
CI2 <- Rx+confinterval
cat("95% CI :",CI1, CI2,"\n")

if (CI1 > 0.8) {print ("Sex assignment:The sample should be assigned as Female")
} else if (CI2 < 0.6) {print ("Sex assignment:The sample should be assigned as Male")
} else if (CI1 > 0.6 & CI2 > 0.8) {print ("Sex assignment:The sample is consistent with XX but not XY")
} else if (CI1 < 0.6 & CI2 < 0.8) {print ("Sex assignment:The sample is consistent with XY but not XX")
} else print ("Sex assignment:The sample could not be assigned")

print ("***It is important to realize that the assignment is invalid, if there is no correlation between the number of reference reads and that of the mapped reads***")

}
