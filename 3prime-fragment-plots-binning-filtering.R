'''
3 prime fragments - normalisation and filtering

The script will import 3 prime fragments tables and added aditional columns with binned and normalised values. Each one of them will be ploted in 1 nucleotide, 10 nt and 30 nt resolution.
Set all paths to "*_same.tab" tables from 3 prime fragments results.
'''

# function will return binned vector
binning <- function(vector, bin_size) {
  counter = 0
  sum = 0.0
  binned_vector <- 0
  for (i in 1:length(vector)) {
    if (i >= length(vector) - 1)
      bin_size = counter
    if (counter < bin_size) {
      sum = sum + vector[i]
      counter = counter + 1
    } else {
      if (sum == 0) {
        norm = 0
      } else {
        norm = (sum / bin_size)
      }
      for (j in 0:bin_size) {
        binned_vector[i + j - bin_size] <- norm
      }
      sum = vector[i]
      counter = 1
    }
  }
  return (binned_vector)
}


### library Ule_Nov3 ###
tpi_48 <- read.table("../Ule_Nov3-Library1-4/TPI-PTC48/Ule_Nov3_NoIndex_L007_R1_001._CCGG_same.tab", header = TRUE)
tpi_100 <- read.table("../Ule_Nov3-Library1-4/TPI-PTC100/Ule_Nov3_NoIndex_L007_R1_001._TTAG_same.tab", header = TRUE)
tpi_120 <- read.table("../Ule_Nov3-Library1-4/TPI-PTC120/Ule_Nov3_NoIndex_L007_R1_001._AATG_same.tab", header = TRUE)
tpi_160 <- read.table("../Ule_Nov3-Library1-4/TPI-PTC160/Ule_Nov3_NoIndex_L007_R1_001._GGTC_same.tab", header = TRUE)
tpi_189 <- read.table("../results/Ule_Nov3-Library1-4/TPI-PTC189/Ule_Nov3_NoIndex_L007_R1_001._AACC_same.tab", header = TRUE)
tpi_smg5 <- read.table("../results/Ule_Nov3-Library1-4/TPI-SMG5/Ule_Nov3_NoIndex_L007_R1_001._GGCA_same.tab",header=TRUE)
tpi_wt <- read.table("../results/Ule_Nov3-Library1-4/TPI-WT/Ule_Nov3_NoIndex_L007_R1_001._ACCT_same.tab",header=TRUE)

#normalisation
tpi_48$normalised <- tpi_48$cDNA / sum(tpi_48$cDNA)
tpi_100$normalised <- tpi_100$cDNA / sum(tpi_100$cDNA)
tpi_120$normalised <- tpi_120$cDNA / sum(tpi_120$cDNA)
tpi_160$normalised <- tpi_160$cDNA / sum(tpi_160$cDNA)
tpi_189$normalised <- tpi_189$cDNA / sum(tpi_189$cDNA)
tpi_smg5$normalised <- tpi_smg5$cDNA / sum(tpi_smg5$cDNA)
tpi_wt$normalised <- tpi_wt$cDNA / sum(tpi_wt$cDNA)

# 10nt binning
tpi_48$bin10 <- binning(tpi_48$normalised, 10)
tpi_100$bin10 <- binning(tpi_100$normalised, 10)
tpi_120$bin10 <- binning(tpi_120$normalised, 10)
tpi_160$bin10 <- binning(tpi_160$normalised, 10)
tpi_189$bin10 <- binning(tpi_189$normalised, 10)
tpi_smg5$bin10 <- binning(tpi_smg5$normalised, 10)
tpi_wt$bin10 <- binning(tpi_wt$normalised, 10)

# 30 nt binning
tpi_48$bin30 <- binning(tpi_48$normalised, 30)
tpi_100$bin30 <- binning(tpi_100$normalised, 30)
tpi_120$bin30 <- binning(tpi_120$normalised, 30)
tpi_160$bin30 <- binning(tpi_160$normalised, 30)
tpi_189$bin30 <- binning(tpi_189$normalised, 30)
tpi_smg5$bin30 <- binning(tpi_smg5$normalised, 30)
tpi_wt$bin30 <- binning(tpi_wt$normalised, 30)

#filter TPI-WT
tpi_48$bin30wt <- tpi_48$bin30 - tpi_wt$bin30
tpi_48$bin30wt[tpi_48$bin30wt < 0] <- 0
tpi_100$bin30wt <- tpi_100$bin30 - tpi_wt$bin30
tpi_100$bin30wt[tpi_100$bin30wt < 0] <- 0
tpi_120$bin30wt <- tpi_120$bin30 - tpi_wt$bin30
tpi_120$bin30wt[tpi_120$bin30wt < 0] <- 0
tpi_160$bin30wt <- tpi_160$bin30 - tpi_wt$bin30
tpi_160$bin30wt[tpi_160$bin30wt < 0] <- 0
tpi_189$bin30wt <- tpi_189$bin30 - tpi_wt$bin30
tpi_189$bin30wt[tpi_189$bin30wt < 0] <- 0
tpi_smg5$bin30wt <- tpi_smg5$bin30 - tpi_wt$bin30
tpi_smg5$bin30wt[tpi_smg5$bin30wt < 0] <- 0

# ploting TPI - 48
par(mfrow=c(4,1))
plot(tpi_48$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-48 - 1nt")
plot(tpi_48$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-48 - 10nt")
plot(tpi_48$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-48 - 30nt")
plot(tpi_48$bin30wt, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-48 - 30nt - filtered by WT")

# ploting TPI - 100
par(mfrow=c(4,1))
plot(tpi_100$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-100 - 1nt")
plot(tpi_100$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-100 - 10nt")
plot(tpi_100$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-100 - 30nt")
plot(tpi_100$bin30wt, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-100 - 30nt - filtered by WT")

# ploting TPI - 120
par(mfrow=c(4,1))
plot(tpi_120$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-120 - 1nt")
plot(tpi_120$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-120 - 10nt")
plot(tpi_120$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-120 - 30nt")
plot(tpi_120$bin30wt, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-120 - 30nt - filtered by WT")

# ploting TPI - 160
par(mfrow=c(4,1))
plot(tpi_160$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160 - 1nt")
plot(tpi_160$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160 - 10nt")
plot(tpi_160$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160 - 30nt")
plot(tpi_160$bin30wt, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160 - 30nt - filtered by WT")

# ploting TPI - 189
par(mfrow=c(4,1))
plot(tpi_189$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-189 - 1nt")
plot(tpi_189$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-189 - 10nt")
plot(tpi_189$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-189 - 30nt")
plot(tpi_189$bin30wt, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-189 - 30nt - filtered by WT")

# ploting TPI - SMG5
par(mfrow=c(4,1))
plot(tpi_smg5$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-SMG5 - 1nt")
plot(tpi_smg5$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-SMG5 - 10nt")
plot(tpi_smg5$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-SMG5 - 30nt")
plot(tpi_smg5$bin30wt, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-SMG5 - 30nt - filtered by WT")

# ploting TPI - WT
par(mfrow=c(3,1))
plot(tpi_wt$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-WT - 1nt")
plot(tpi_wt$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-WT - 10nt")
plot(tpi_wt$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-WT - 30nt")


# save tables
write.table(tpi_48, file= "../Ule_Nov3-Library1-4/Ule_Nov3-Library1-4-TPI-48.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_100, file= "../Ule_Nov3-Library1-4/Ule_Nov3-Library1-4-TPI-100.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_120, file= "../Ule_Nov3-Library1-4/Ule_Nov3-Library1-4-TPI-120.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_160, file= "../Ule_Nov3-Library1-4/Ule_Nov3-Library1-4-TPI-160.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_189, file= "../Ule_Nov3-Library1-4/Ule_Nov3-Library1-4-TPI-189.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_smg5, file= "../Ule_Nov3-Library1-4/Ule_Nov3-Library1-4-TPI-SMG5.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_wt, file= "../Ule_Nov3-Library1-4/Ule_Nov3-Library1-4-TPI-WT.tab", quote=FALSE, row.names= FALSE, sep="\t")


#####################################################

### library Ule_Nov4 ###
tpi_48 <- read.table("../Ule_Nov4-Library1-4/TPI-PTC48/Ule_Nov4_NoIndex_L008_R1_001._CCGG_same.tab", header = TRUE)
tpi_100 <- read.table("../Ule_Nov4-Library1-4/TPI-PTC100/Ule_Nov4_NoIndex_L008_R1_001._TTAG_same.tab", header = TRUE)
tpi_120 <- read.table("../Ule_Nov4-Library1-4/TPI-PTC120/Ule_Nov4_NoIndex_L008_R1_001._AATG_same.tab", header = TRUE)
tpi_160 <- read.table("../Ule_Nov4-Library1-4/TPI-PTC160/Ule_Nov4_NoIndex_L008_R1_001._GGTC_same.tab", header = TRUE)
tpi_189 <- read.table("../Ule_Nov4-Library1-4/TPI-PTC189/Ule_Nov4_NoIndex_L008_R1_001._AACC_same.tab", header = TRUE)
tpi_smg5 <- read.table("../Ule_Nov4-Library1-4/TPI-SMG5/Ule_Nov4_NoIndex_L008_R1_001._GGCA_same.tab",header=TRUE)
tpi_wt <- read.table("../Ule_Nov4-Library1-4/TPI-WT/Ule_Nov4_NoIndex_L008_R1_001._ACCT_same.tab",header=TRUE)

#normalisation
tpi_48$normalised <- tpi_48$cDNA / sum(tpi_48$cDNA)
tpi_100$normalised <- tpi_100$cDNA / sum(tpi_100$cDNA)
tpi_120$normalised <- tpi_120$cDNA / sum(tpi_120$cDNA)
tpi_160$normalised <- tpi_160$cDNA / sum(tpi_160$cDNA)
tpi_189$normalised <- tpi_189$cDNA / sum(tpi_189$cDNA)
tpi_smg5$normalised <- tpi_smg5$cDNA / sum(tpi_smg5$cDNA)
tpi_wt$normalised <- tpi_wt$cDNA / sum(tpi_wt$cDNA)

# 10nt binning
tpi_48$bin10 <- binning(tpi_48$normalised, 10)
tpi_100$bin10 <- binning(tpi_100$normalised, 10)
tpi_120$bin10 <- binning(tpi_120$normalised, 10)
tpi_160$bin10 <- binning(tpi_160$normalised, 10)
tpi_189$bin10 <- binning(tpi_189$normalised, 10)
tpi_smg5$bin10 <- binning(tpi_smg5$normalised, 10)
tpi_wt$bin10 <- binning(tpi_wt$normalised, 10)

# 30 nt binning
tpi_48$bin30 <- binning(tpi_48$normalised, 30)
tpi_100$bin30 <- binning(tpi_100$normalised, 30)
tpi_120$bin30 <- binning(tpi_120$normalised, 30)
tpi_160$bin30 <- binning(tpi_160$normalised, 30)
tpi_189$bin30 <- binning(tpi_189$normalised, 30)
tpi_smg5$bin30 <- binning(tpi_smg5$normalised, 30)
tpi_wt$bin30 <- binning(tpi_wt$normalised, 30)

#filter TPI-WT
tpi_48$bin30wt <- tpi_48$bin30 - tpi_wt$bin30
tpi_48$bin30wt[tpi_48$bin30wt < 0] <- 0
tpi_100$bin30wt <- tpi_100$bin30 - tpi_wt$bin30
tpi_100$bin30wt[tpi_100$bin30wt < 0] <- 0
tpi_120$bin30wt <- tpi_120$bin30 - tpi_wt$bin30
tpi_120$bin30wt[tpi_120$bin30wt < 0] <- 0
tpi_160$bin30wt <- tpi_160$bin30 - tpi_wt$bin30
tpi_160$bin30wt[tpi_160$bin30wt < 0] <- 0
tpi_189$bin30wt <- tpi_189$bin30 - tpi_wt$bin30
tpi_189$bin30wt[tpi_189$bin30wt < 0] <- 0
tpi_smg5$bin30wt <- tpi_smg5$bin30 - tpi_wt$bin30
tpi_smg5$bin30wt[tpi_smg5$bin30wt < 0] <- 0

# ploting TPI - 48
par(mfrow=c(4,1))
plot(tpi_48$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-48 - 1nt")
plot(tpi_48$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-48 - 10nt")
plot(tpi_48$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-48 - 30nt")
plot(tpi_48$bin30wt, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-48 - 30nt - filtered by WT")

# ploting TPI - 100
par(mfrow=c(4,1))
plot(tpi_100$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-100 - 1nt")
plot(tpi_100$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-100 - 10nt")
plot(tpi_100$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-100 - 30nt")
plot(tpi_100$bin30wt, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-100 - 30nt - filtered by WT")

# ploting TPI - 120
par(mfrow=c(4,1))
plot(tpi_120$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-120 - 1nt")
plot(tpi_120$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-120 - 10nt")
plot(tpi_120$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-120 - 30nt")
plot(tpi_120$bin30wt, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-120 - 30nt - filtered by WT")

# ploting TPI - 160
par(mfrow=c(4,1))
plot(tpi_160$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160 - 1nt")
plot(tpi_160$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160 - 10nt")
plot(tpi_160$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160 - 30nt")
plot(tpi_160$bin30wt, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160 - 30nt - filtered by WT")

# ploting TPI - 189
par(mfrow=c(4,1))
plot(tpi_189$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-189 - 1nt")
plot(tpi_189$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-189 - 10nt")
plot(tpi_189$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-189 - 30nt")
plot(tpi_189$bin30wt, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-189 - 30nt - filtered by WT")

# ploting TPI - SMG5
par(mfrow=c(4,1))
plot(tpi_smg5$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-SMG5 - 1nt")
plot(tpi_smg5$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-SMG5 - 10nt")
plot(tpi_smg5$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-SMG5 - 30nt")
plot(tpi_smg5$bin30wt, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-SMG5 - 30nt - filtered by WT")

# ploting TPI - WT
par(mfrow=c(3,1))
plot(tpi_wt$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-WT - 1nt")
plot(tpi_wt$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-WT - 10nt")
plot(tpi_wt$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-WT - 30nt")

# save tables
write.table(tpi_48, file= "../Ule_Nov4-Library1-4/Ule_Nov4-Library1-4-TPI-48.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_100, file= "../Ule_Nov4-Library1-4/Ule_Nov4-Library1-4-TPI-100.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_120, file= "../Ule_Nov4-Library1-4/Ule_Nov4-Library1-4-TPI-120.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_160, file= "../Ule_Nov4-Library1-4/Ule_Nov4-Library1-4-TPI-160.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_189, file= "../Ule_Nov4-Library1-4/Ule_Nov4-Library1-4-TPI-189.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_smg5, file= "../Ule_Nov4-Library1-4/Ule_Nov4-Library1-4-TPI-SMG5.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_wt, file= "../Ule_Nov4-Library1-4/Ule_Nov4-Library1-4-TPI-WT.tab", quote=FALSE, row.names= FALSE, sep="\t")


########################
### library Ule_Nov5 ###
########################

tpi_160_ACCT <- read.table("../Ule_Nov1-Library5/Ule_Nov1_NoIndex_L005_R1_001._ACCT_same.tab", header = TRUE)
tpi_160_CCGG <- read.table("../Ule_Nov1-Library5/Ule_Nov1_NoIndex_L005_R1_001._CCGG_same.tab", header = TRUE)
tpi_160_TTAG <- read.table("../Ule_Nov1-Library5/Ule_Nov1_NoIndex_L005_R1_001._TTAG_same.tab", header = TRUE)
tpi_160_AATG <- read.table("../Ule_Nov1-Library5/Ule_Nov1_NoIndex_L005_R1_001._AATG_same.tab", header = TRUE)
tpi_160_GGTC <- read.table("../Ule_Nov1-Library5/Ule_Nov1_NoIndex_L005_R1_001._GGTC_same.tab", header = TRUE)
tpi_160_AACC <- read.table("../Ule_Nov1-Library5/Ule_Nov1_NoIndex_L005_R1_001._AACC_same.tab", header = TRUE)

#normalisation
tpi_160_ACCT$normalised <- tpi_160_ACCT$cDNA / sum(tpi_160_ACCT$cDNA)
tpi_160_CCGG$normalised <- tpi_160_CCGG$cDNA / sum(tpi_160_CCGG$cDNA)
tpi_160_TTAG$normalised <- tpi_160_TTAG$cDNA / sum(tpi_160_TTAG$cDNA)
tpi_160_AATG$normalised <- tpi_160_AATG$cDNA / sum(tpi_160_AATG$cDNA)
tpi_160_GGTC$normalised <- tpi_160_GGTC$cDNA / sum(tpi_160_GGTC$cDNA)
tpi_160_AACC$normalised <- tpi_160_AACC$cDNA / sum(tpi_160_AACC$cDNA)
#tpi_wt$normalised <- tpi_wt$cDNA / sum(tpi_wt$cDNA) ??? #which wild type?

# 10nt binning
tpi_160_ACCT$bin10 <- binning(tpi_160_ACCT$normalised, 10)
tpi_160_CCGG$bin10 <- binning(tpi_160_CCGG$normalised, 10)
tpi_160_TTAG$bin10 <- binning(tpi_160_TTAG$normalised, 10)
tpi_160_AATG$bin10 <- binning(tpi_160_AATG$normalised, 10)
tpi_160_GGTC$bin10 <- binning(tpi_160_GGTC$normalised, 10)
tpi_160_AACC$bin10 <- binning(tpi_160_AACC$normalised, 10)
#tpi_wt$bin10 <- binning(tpi_wt$normalised, 10)

# 30 nt binning
tpi_160_ACCT$bin30 <- binning(tpi_160_ACCT$normalised, 30)
tpi_160_CCGG$bin30 <- binning(tpi_160_CCGG$normalised, 30)
tpi_160_TTAG$bin30 <- binning(tpi_160_TTAG$normalised, 30)
tpi_160_AATG$bin30 <- binning(tpi_160_AATG$normalised, 30)
tpi_160_GGTC$bin30 <- binning(tpi_160_GGTC$normalised, 30)
tpi_160_AACC$bin30 <- binning(tpi_160_AACC$normalised, 30)
#tpi_wt$bin30 <- binning(tpi_wt$normalised, 30)

#filter TPI-WT
#????

# ploting TPI - 160_ACCT
par(mfrow=c(3,1))
plot(tpi_160_ACCT$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160_ACCT - 1nt")
plot(tpi_160_ACCT$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160_ACCT - 10nt")
plot(tpi_160_ACCT$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160_ACCT - 30nt")

# ploting TPI - 160_CCGG
par(mfrow=c(3,1))
plot(tpi_160_CCGG$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160_CCGG - 1nt")
plot(tpi_160_CCGG$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160_CCGG - 10nt")
plot(tpi_160_CCGG$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160_CCGG - 30nt")

# ploting TPI - 160_TTAG
par(mfrow=c(3,1))
plot(tpi_160_TTAG$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160_TTAG - 1nt")
plot(tpi_160_TTAG$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160_TTAG - 10nt")
plot(tpi_160_TTAG$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160_TTAG - 30nt")

# ploting TPI - 160_AATG
par(mfrow=c(3,1))
plot(tpi_160_AATG$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160_AATG - 1nt")
plot(tpi_160_AATG$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160_AATG - 10nt")
plot(tpi_160_AATG$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160_AATG - 30nt")

# ploting TPI - 160_GGTC
par(mfrow=c(3,1))
plot(tpi_160_GGTC$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160_GGTC - 1nt")
plot(tpi_160_GGTC$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160_GGTC - 10nt")
plot(tpi_160_GGTC$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160_GGTC - 30nt")

# ploting TPI - 160_AACC
par(mfrow=c(3,1))
plot(tpi_160_AACC$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160_AACC - 1nt")
plot(tpi_160_AACC$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160_AACC - 10nt")
plot(tpi_160_AACC$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-160_AACC - 30nt")

# save tables
write.table(tpi_160_AACC, file="../Ule_Nov1-Library5/Ule_Nov1-Library5-TPI-160_AACC.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_160_AATG, file="../Ule_Nov1-Library5/Ule_Nov1-Library5-TPI-160_AATG.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_160_ACCT, file="../Ule_Nov1-Library5/Ule_Nov1-Library5-TPI-160_ACCT.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_160_CCGG, file="../Ule_Nov1-Library5/Ule_Nov1-Library5-TPI-160_CCGG.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_160_GGTC, file="../Ule_Nov1-Library5/Ule_Nov1-Library5-TPI-160_GGTC.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_160_TTAG, file="../Ule_Nov1-Library5/Ule_Nov1-Library5-TPI-160_TTAG.tab", quote=FALSE, row.names= FALSE, sep="\t")


########################
### library Ule_Nov6 ###
########################

tpi_SMG5_AATG <- read.table("../Ule_Nov2-Library6/Library6-SMG5_NM_015327/Ule_Nov2_NoIndex_L006_R1_001._AATG_same.tab", header = TRUE)
tpi_SMG5_GGTC <- read.table("../Ule_Nov2-Library6/Library6-SMG5_NM_015327/Ule_Nov2_NoIndex_L006_R1_001._GGTC_same.tab", header = TRUE)
tpi_ACTB <- read.table("../Ule_Nov2-Library6/Library6-ACTB_NM_001101/Ule_Nov2_NoIndex_L006_R1_001._AACC_same.tab", header = TRUE)

#normalisation
tpi_SMG5_AATG$normalised <- tpi_SMG5_AATG$cDNA / sum(tpi_SMG5_AATG$cDNA)
tpi_SMG5_GGTC$normalised <- tpi_SMG5_GGTC$cDNA / sum(tpi_SMG5_GGTC$cDNA)
tpi_ACTB$normalised <- tpi_ACTB$cDNA / sum(tpi_ACTB$cDNA)

# 10nt binning
tpi_SMG5_AATG$bin10 <- binning(tpi_SMG5_AATG$normalised, 10)
tpi_SMG5_GGTC$bin10 <- binning(tpi_SMG5_GGTC$normalised, 10)
tpi_ACTB$bin10 <- binning(tpi_ACTB$normalised, 10)

# 10nt binning
tpi_SMG5_AATG$bin30 <- binning(tpi_SMG5_AATG$normalised, 30)
tpi_SMG5_GGTC$bin30 <- binning(tpi_SMG5_GGTC$normalised, 30)
tpi_ACTB$bin30 <- binning(tpi_ACTB$normalised, 30)


# ploting SMG5_AATG
par(mfrow=c(3,1))
plot(tpi_SMG5_AATG$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "SMG5(NM_015327) AATG - 1nt")
plot(tpi_SMG5_AATG$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "SMG5(NM_015327) AATG - 10nt")
plot(tpi_SMG5_AATG$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "SMG5(NM_015327) AATG - 30nt")

# ploting SMG5_GGTC
par(mfrow=c(3,1))
plot(tpi_SMG5_GGTC$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "SMG5(NM_015327) GGTC - 1nt")
plot(tpi_SMG5_GGTC$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "SMG5(NM_015327) GGTC - 10nt")
plot(tpi_SMG5_GGTC$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "SMG5(NM_015327) GGTC - 30nt")

# ploting ACTB (NM_001101)
par(mfrow=c(3,1))
plot(tpi_ACTB$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "ACTB(NM_001101) - 1nt")
plot(tpi_ACTB$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "ACTB(NM_001101) - 10nt")
plot(tpi_ACTB$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "ACTB(NM_001101) - 30nt")


# save tables
write.table(tpi_SMG5_AATG, file="../Ule_Nov2-Library6/Ule_Nov2-Library6-SMG5(NM_015327)-AATG.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_SMG5_GGTC, file="../Ule_Nov2-Library6/Ule_Nov2-Library6-SMG5(NM_015327)-GGTC.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_ACTB, file="../Ule_Nov2-Library6/Ule_Nov2-Library6-ACTB(NM_001101).tab", quote=FALSE, row.names= FALSE, sep="\t")


