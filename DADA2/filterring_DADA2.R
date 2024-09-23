
# 1. fastq fileをzshにてPrimerのと整合性を確認 (sample metadataとprimer F/Rをよく確認) → adapterはcutaduptにて除去すること
# 2. dada2にて各fastq fileのplotQualityProfile
# 3. primerを含む末端をTrimming


ls(name = "package:dada2")

# ----- May 10, 2024 at 11:45:35 -----
dir <- getwd()
path <- paste(dir,"PRJNA643596",sep = "/" ) # .fastq.gz fileを格納するdirectory
list.files(path)


fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# 不要なITS,Culturomeを削除
fnFs <- fnFs[-81]
fnFs_v1 <- fnFs[-c(1:6)]
fnFs_v2 <- fnFs_v1[-c(11:44)]
fnFs <- fnFs_v2

fnRs <- fnRs[-81]
fnRs_v1 <- fnRs[-c(1:6)]
fnRs_v2 <- fnRs_v1[-c(11:44)]
fnRs <- fnRs_v2
rm(fnFs_v1,fnRs_v1,fnFs_v2,fnRs_v2)

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


# 各fastq fileのQuality check (末端でQualityが落ちていないかを確認)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz")) 
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz")) 
names(filtFs) <- sample.names
names(filtRs) <- sample.names


## parameterによって結果が異なるのは以下codeから
# matrix型のデータ
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(225,225), # Primer(20or21bp)のおよその位置でtruncLen
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

# EpiとEndで区別
Epi_out <- out[-c(1:9,11:21),]
End_out <- out[c(1:21),]
End_out <- End_out[-10,]

# 40168_2022_1234_MOESM3_ESM_TableS1_aと一致
sum(Epi_out[,1])
sum(End_out[,1])

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Quality scoreの上昇と共にerror率は低下しているかどうか
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


# sampleの予測, dereplicaredをinput?? dada2 algorithmによってreadをdenoising
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

# F/Rのmerge
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)


getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

sum(track[,6])


#####
R_data <- read.csv("Supplementary Table 3- ASV table and Wilcoxon test of the epiphytic bacteria_relative abundance_v1.csv")







