#3 prime fragments - normalisation and filtering
#The script will import 3 prime fragments tables and added aditional columns with binned and normalised values. Each one of them will be ploted in 1 nucleotide, 10 nt and 30 nt resolution.
#Set all paths to "*_same.tab" tables from 3 prime fragments results.

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

#################
### library 1 ###
#################
tpi_48 <- read.table("/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/laneU1406B_L002_R1_001._TTGT_same.tab", header = TRUE)
tpi_100 <- read.table("/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/laneU1406B_L002_R1_001._CAAT_same.tab", header = TRUE)
tpi_120 <- read.table("/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/laneU1406B_L002_R1_001._ACCT_same.tab", header = TRUE)
tpi_160 <- read.table("/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/laneU1406B_L002_R1_001._GGCG_same.tab", header = TRUE)
tpi_189 <- read.table("/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/laneU1406B_L002_R1_001._CCGG_same.tab", header = TRUE)
tpi_smg5 <- read.table("/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/laneU1406B_L002_R1_001._TTAG_same.tab",header=TRUE)
tpi_wt <- read.table("/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/laneU1406B_L002_R1_001._GGTT_same.tab",header=TRUE)

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
#tpi_smg5$bin30wt <- tpi_smg5$bin30 - tpi_wt$bin30
#tpi_smg5$bin30wt[tpi_smg5$bin30wt < 0] <- 0

pdf("Library1-figures.pdf")

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
#plot(tpi_smg5$bin30wt, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-SMG5 - 30nt - filtered by WT")

# ploting TPI - WT
par(mfrow=c(3,1))
plot(tpi_wt$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-WT - 1nt")
plot(tpi_wt$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-WT - 10nt")
plot(tpi_wt$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-WT - 30nt")


# save tables
write.table(tpi_48, file= "/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/Library1-TPI-48.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_100, file= "/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/Library1-TPI-100.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_120, file= "/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/Library1-TPI-120.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_160, file= "/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/Library1-TPI-160.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_189, file= "/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/Library1-TPI-189.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_smg5, file= "/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/Library1-TPI-SMG5.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_wt, file= "/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/Library1-TPI-WT.tab", quote=FALSE, row.names= FALSE, sep="\t")

#####################################################
#################
### library 2 ###
#################
tpi_48 <- read.table("/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/laneU1406B_L002_R1_001._TGGC_same.tab", header = TRUE)
tpi_100 <- read.table("/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/laneU1406B_L002_R1_001._GGTC_same.tab", header = TRUE)
tpi_120 <- read.table("/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/laneU1406B_L002_R1_001._AACC_same.tab", header = TRUE)
tpi_160 <- read.table("/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/laneU1406B_L002_R1_001._CCAC_same.tab", header = TRUE)
tpi_189 <- read.table("/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/laneU1406B_L002_R1_001._CGGA_same.tab", header = TRUE)
tpi_smg5 <- read.table("/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/laneU1406B_L002_R1_001._GGCA_same.tab",header=TRUE)
tpi_wt <- read.table("/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/laneU1406B_L002_R1_001._AATG_same.tab",header=TRUE)

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
#tpi_smg5$bin30wt <- tpi_smg5$bin30 - tpi_wt$bin30
#tpi_smg5$bin30wt[tpi_smg5$bin30wt < 0] <- 0

dev.off()
pdf("Library2-figures.pdf")

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
#plot(tpi_smg5$bin30wt, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-SMG5 - 30nt - filtered by WT")

# ploting TPI - WT
par(mfrow=c(3,1))
plot(tpi_wt$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-WT - 1nt")
plot(tpi_wt$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-WT - 10nt")
plot(tpi_wt$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-WT - 30nt")

# save tables
write.table(tpi_48, file= "/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/Library2-TPI-48.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_100, file= "/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/Library2-TPI-100.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_120, file= "/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/Library2-TPI-120.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_160, file= "/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/Library2-TPI-160.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_189, file= "/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/Library2-TPI-189.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_smg5, file= "/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/Library2-TPI-SMG5.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_wt, file= "/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/Library2-TPI-WT.tab", quote=FALSE, row.names= FALSE, sep="\t")


#################
### library 3 ###
#################

tpi_48 <- read.table("/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/laneU1406B_L002_R1_001._TTAA_same.tab", header = TRUE)
tpi_100 <- read.table("/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/laneU1406B_L002_R1_001._ATTT_same.tab", header = TRUE)
tpi_120 <- read.table("/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/laneU1406B_L002_R1_001._CCTT_same.tab", header = TRUE)
tpi_160 <- read.table("/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/laneU1406B_L002_R1_001._TATT_same.tab", header = TRUE)
tpi_189 <- read.table("/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/laneU1406B_L002_R1_001._GCGT_same.tab", header = TRUE)
tpi_smg5 <- read.table("/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/laneU1406B_L002_R1_001._AAGT_same.tab",header=TRUE)
tpi_wt <- read.table("/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/laneU1406B_L002_R1_001._AATA_same.tab",header=TRUE)

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
#tpi_smg5$bin30wt <- tpi_smg5$bin30 - tpi_wt$bin30
#tpi_smg5$bin30wt[tpi_smg5$bin30wt < 0] <- 0

dev.off()
pdf("Library3-figures.pdf")

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
#plot(tpi_smg5$bin30wt, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-SMG5 - 30nt - filtered by WT")

# ploting TPI - WT
par(mfrow=c(3,1))
plot(tpi_wt$normalised, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-WT - 1nt")
plot(tpi_wt$bin10, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-WT - 10nt")
plot(tpi_wt$bin30, type = "h", ylab="normalised counts", xlab = "position", main = "TPI-WT - 30nt")

# save tables
write.table(tpi_48, file= "/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/Library3-TPI-48.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_100, file= "/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/Library3-TPI-100.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_120, file= "/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/Library3-TPI-120.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_160, file= "/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/Library3-TPI-160.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_189, file= "/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/Library3-TPI-189.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_smg5, file= "/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/Library3-TPI-SMG5.tab", quote=FALSE, row.names= FALSE, sep="\t")
write.table(tpi_wt, file= "/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/2014.07.16-data/results/Library3-TPI-WT.tab", quote=FALSE, row.names= FALSE, sep="\t")



