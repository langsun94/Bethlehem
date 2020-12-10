## --------------------------------------------------------------------
R.Version() 
#R version 3.4.3
library(dada2)
packageVersion("dada2")
#1.9.0
path <- "/Volumes/efs/CANR/Ansci/D\'Amico\ Lab/Langsun_las17015/2.\ Bethlehem\ project/Bethlehem2020/Bacteria/fastq"


## --------------------------------------------------------------------
#fastq files in matched order
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


## --------------------------------------------------------------------
#Examine quality profiles of forward and reverse reads
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])


## --------------------------------------------------------------------
#Perform filtering and trimming
filt_path <- file.path(path, "filtered") 
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)


## --------------------------------------------------------------------
#Learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


## --------------------------------------------------------------------
#Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names


## --------------------------------------------------------------------
#Apply the core sequence-variant inference algorithm
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool="pseudo")
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool="pseudo")
dadaFs[[1]]


## --------------------------------------------------------------------
#Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


## --------------------------------------------------------------------
#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))


## --------------------------------------------------------------------
#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
save(seqtab.nochim, file = "seqtab.nochim.RData") #save with data of the output data
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)


## --------------------------------------------------------------------
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
save(track, file = "track.RData") 


## --------------------------------------------------------------------
taxa_silva <- assignTaxonomy(seqtab.nochim,"/Volumes/efs/CANR/Ansci/D\'Amico\ Lab/Langsun_las17015/Bacteria\ reference/silva_nr_v132_train_set.fa.gz", minBoot=80)
#genus.species <- addSpecies(taxa_silva, "/Volumes/efs/CANR/Ansci/D\'Amico\ Lab/Langsun_las17015/Bacteria\ reference/silva_species_assignment_v132.fa.gz", allowMultiple=TRUE)
save(taxa_silva, file = "genus.RData")


## --------------------------------------------------------------------
#only maintain bacteria, and remove mitochondria and chloroplast
taxa<-as.data.frame(taxa_silva)
taxa1<-taxa[taxa$Kingdom %in% "Bacteria",]
taxa2<-taxa1[!taxa1$Family %in% "Mitochondria",]
taxa3<-taxa2[!taxa2$Order %in% "Chloroplast",]

#save taxa4 with exact sequences
taxa4<-as.matrix(taxa3)
save(taxa4, file = "taxa.sequence.RData")

## --------------------------------------------------------------------
#import metadata
samdf=read.table("/Volumes/efs/CANR/Ansci/D\'Amico\ Lab/Langsun_las17015/2.\ Bethlehem\ project/Bethlehem2020/Bacteria/metadata.txt",header = T)
rownames(samdf) <- samdf$Barcode

## --------------------------------------------------------------------
#create phyloseq
library(phyloseq)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa4))

#save exact sequences and short names to ASV
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

save(ps, file = "phyloseq.RData")
