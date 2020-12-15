## --------------------------------------------------------------------
#load pyloseq and metadata

load("/Volumes/efs/CANR/Ansci/D\'Amico\ Lab/Langsun_las17015/2.\ Bethlehem\ project/Bethlehem2020/Bacteria/dada2/phyloseq.RData")


## --------------------------------------------------------------------
#libraries
library(ape)
library(phyloseq)
library(ggplot2)
library(vegan)
library(dplyr)
library(tibble)
library(ggdendro)
library(plyr)
library(gridExtra)
library(pals)
library(reshape)
library(ggsignif)
library(pairwiseAdonis)
library(indicspecies)

## --------------------------------------------------------------------
#dairy alpha

ps.dairy<-subset_samples(ps,source=="dairy")

p = plot_richness(ps.dairy, x="sampletype", measures=c("Shannon", "Simpson"))
p$data$sampletype <- as.character(p$data$sampltype)
p$data$sampletype <- factor(p$data$sampletype, 
                            levels=c("rawmilk","filtered_milk","milk_before_ripening",
                                     "milk_stored_overnight","ripened_milk","mixed_milk",
                                     "whey","cheese_curd","cheese"))

p1 = p + geom_boxplot() + 
  theme(axis.title=element_blank(),axis.text.x = element_text(angle = 270),legend.position = "none",panel.background = element_rect(fill = 'white', colour = 'black'))


## --------------------------------------------------------------------
#dairy beta 

prop.dairy <- transform_sample_counts(ps.dairy, function(otu) otu/sum(otu))

ord.nmds.bray <- ordinate(prop.dairy, method="NMDS", distance="bray")

p = plot_ordination(prop.dairy, ord.nmds.bray, color="sampletype")
p$data$sampletype <- factor(p$data$sampletype, 
                            levels=c("rawmilk","filtered_milk","milk_before_ripening",
                                     "milk_stored_overnight","ripened_milk","mixed_milk",
                                     "whey","cheese_curd","cheese"))

p2 = p + theme(legend.text=element_text(size=8),legend.title = element_blank(),panel.background = element_rect(fill = 'white', colour = 'black'), panel.grid=element_blank())


## --------------------------------------------------------------------
#dairy bar plot

merge.dairy <- merge_samples(ps.dairy, "sampletype")
merge2.dairy <- tax_glom(merge.dairy, "Genus")
prop.merge.dairy <- transform_sample_counts(merge2.dairy, function(otu) otu/sum(otu))
melt.dairy <- psmelt(prop.merge.dairy)

melt.dairy1<-melt.dairy
melt.dairy1$Genus <- as.character(melt.dairy1$Genus)
melt.dairy1$Genus[melt.dairy1$Abundance < 0.01] <- "< 1% abund."
melt.dairy1$Genus[is.na(melt.dairy1$Genus)] <- "unclassified_bacteria"
melt.dairy1$Genus <- reorder(melt.dairy1$Genus, -melt.dairy1$Abundance)

melt.dairy1$Sample <- factor(melt.dairy1$Sample, 
                    c("rawmilk","filtered_milk","milk_before_ripening",
                      "milk_stored_overnight","ripened_milk","mixed_milk",
                      "whey","cheese_curd","cheese"))

library(randomcoloR)
d = length(unique(melt.dairy1$Genus))
palette <- distinctColorPalette(d)

p.s3 <- ggplot(melt.dairy1, aes(x=Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette) +
  theme(axis.text.x = element_text(vjust=0.5, size=12, angle = 270),axis.title = element_text(size=12),legend.text=element_text(size=12),legend.key.size = unit(0.4, "cm"), legend.title = element_blank(),title = element_text(size=12,face = "bold"),panel.background = element_blank(), panel.grid=element_blank(), axis.line=element_line('black')) + guides(fill = guide_legend(ncol = 2)) +  ylab("Relative Abundance (Genera > 1%) \n") + xlab("Dairy samples") 


## --------------------------------------------------------------------
#dairy top20 ASV


top20 <- names(sort(taxa_sums(ps.dairy), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps.dairy, function(otu) otu/sum(otu))
ps.top20 <- prune_taxa(top20, ps.top20)

top20.melt <- psmelt(ps.top20) 

top20.new <- top20.melt
top20.new$Abundance_per <- (top20.melt$Abundance)*100

top20.new$sampletype <- factor(top20.new$sampletype, 
                             c("rawmilk","filtered_milk","milk_before_ripening",
                               "milk_stored_overnight","ripened_milk","mixed_milk",
                               "whey","cheese_curd","cheese"))

heatmap_top20 <- ggplot(data = top20.new, mapping = aes(x = sample,y = OTU,fill = Abundance_per)) + 
  geom_tile() + 
  scale_fill_viridis_c(option = 'D', trans = "log10",  breaks = c(0.01, 0.1, 1, 10)) + 
  facet_grid(Family~sampletype,switch = "both", scales = "free", space = "free") + 
  labs(fill = "% relative abundance") +
  theme(axis.text.x = element_blank(),axis.text.y = element_text(size=8),axis.title = element_blank(),legend.title = element_text(size = 8),legend.text = element_text(size = 8),strip.placement = "outside", strip.background = element_blank(),strip.text.x= element_text(angle = 270),strip.text.y.left = element_text(angle = 0))



cheese.heat<-subset(top20.new,top20.new$sampletype=="cheese")
cheese.heat$section <- factor(cheese.heat$section,c("core","middlesection","rindsection","rind"))
cheese.heat$day <- factor(cheese.heat$day,c("0","4","7","14","21","44","60"))


cheese.heatmap <- ggplot(data = cheese.heat, mapping = aes(x = day,y = OTU,fill = Abundance_per)) + 
  geom_tile() + 
  scale_fill_viridis_c(option = 'D', trans = "log10",  breaks = c(0.01, 0.1, 1, 10)) + 
  facet_grid(Family~section,switch = "both", scales = "free", space = "free") + 
  labs(fill = "% relative abundance") +
  theme(axis.text.y = element_text(size=8),axis.title = element_blank(),legend.title = element_text(size = 8),legend.text = element_text(size = 8),strip.placement = "outside", strip.background = element_blank(),strip.text.x= element_text(angle = 270),strip.text.y.left = element_text(angle = 0))

## --------------------------------------------------------------------
#dairy cheese bar

ps.cheese <- subset_samples(ps,sampletype=="cheese")
merge.cheese <- tax_glom(ps.cheese, "Genus")
prop.merge.cheese <- transform_sample_counts(merge.cheese, function(otu) otu/sum(otu))
melt.cheese <- psmelt(prop.merge.cheese)

melt.cheese1<-melt.cheese
melt.cheese1$Genus <- as.character(melt.cheese1$Genus)
melt.cheese1$Genus[melt.cheese1$Abundance < 0.01] <- "< 1% abund."
melt.cheese1$Genus[is.na(melt.cheese1$Genus)] <- "unclassified_bacteria"
melt.cheese1$Genus <- reorder(melt.cheese1$Genus, -melt.cheese1$Abundance)


palette <- c(#"all others"="grey",
             "< 1% abund."="antiquewhite",
             #"Acinetobacter"="#ED2E49",
             #"Ruminococcaceae_UCG-005"="#45A6D1",
             #"Corynebacterium_1"="#468423",
             "Lactococcus"="#885EA5",
             #"Kocuria"="#EA7C52",
             "Lactobacillus"="#F3D01F",
             "Brachybacterium" = "#B5664C",
             #"Staphylococcus"="#EB65AF",
             "Brevibacterium"= "#6EAB97",
             "Streptomyces"="#8FCC14",
             #"Pseudonocardia"="#6A8ACA",
             "Leuconostoc"="#EC9E8E",
             "Streptococcus"="#D8A52F",
             "Chryseobacterium"="#C1D82F",
             "Enterococcus"="#14C20C",
             "Actinomyces"="#0CC2BF",
             "Stenotrophomonas"="#837AF0"
             )

melt.cheese1$day <- factor(melt.cheese1$day,c("0","4","7","14","21","44","60"))
melt.cheese1$section <- factor(melt.cheese1$section,c("core","middlesection","rindsection","rind"))

p.s4 <- ggplot(melt.cheese1, aes(x=day, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette) +
  facet_grid(section~.,switch = "both", scales = "free", space = "free") + 
  theme(axis.text.x = element_text(vjust=0.5, size=12),axis.title = element_text(size=12),legend.text=element_text(size=12),legend.key.size = unit(0.4, "cm"), legend.title = element_blank(),title = element_text(size=12,face = "bold"),panel.background = element_blank(), panel.grid=element_blank(), axis.line=element_line('black'))  +  ylab("Relative Abundance (Genera > 1%) \n") + xlab("Cheese samples") 



## --------------------------------------------------------------------
#common and unique asv between cheese sections

merge.cheese <- merge_samples(ps.cheese, "section")
f.cheese = filter_taxa(merge.cheese, function(x) mean(x) > 0, TRUE)
df.cheese<-as.data.frame(otu_table(f.cheese))
df<-df.cheese

a <- apply(df, 1, function(x) paste(colnames(df)[which(x > 0)], collapse = ", "))

listInput <- list(core = c("ASV1", "ASV2", "ASV3", "ASV4", "ASV5", "ASV6", "ASV7", "ASV8", "ASV9", "ASV10", "ASV11", "ASV14", "ASV16", "ASV18", "ASV19", "ASV20", "ASV21", "ASV22", "ASV24", "ASV25", "ASV30", "ASV31", "ASV34", "ASV35", "ASV36", "ASV38", "ASV39", "ASV40", "ASV42", "ASV48", "ASV52", "ASV54", "ASV58", "ASV59", "ASV60", "ASV61", "ASV62", "ASV63", "ASV68", "ASV70", "ASV71", "ASV75", "ASV76", "ASV78", "ASV80", "ASV82", "ASV90", "ASV93", "ASV94", "ASV96", "ASV99", "ASV101", "ASV112", "ASV137", "ASV154", "ASV159", "ASV160", "ASV164", "ASV190", "ASV203", "ASV204", "ASV211", "ASV215", "ASV217", "ASV286", "ASV322", "ASV328", "ASV378", "ASV406", "ASV408", "ASV438", "ASV515", "ASV566", "ASV573", "ASV609", "ASV624", "ASV641", "ASV703", "ASV743", "ASV778", "ASV814", "ASV918", "ASV1105", "ASV1165", "ASV1244", "ASV1280", "ASV1297", "ASV2122", "ASV2635"), 
                  middlesection = c("ASV1", "ASV2", "ASV3", "ASV4", "ASV5", "ASV6", "ASV7", "ASV8", "ASV10", "ASV11", "ASV14", "ASV16", "ASV18", "ASV19", "ASV20", "ASV22", "ASV25", "ASV30", "ASV31", "ASV35", "ASV38", "ASV39", "ASV40", "ASV42", "ASV48", "ASV49", "ASV54", "ASV58", "ASV59", "ASV60", "ASV61", "ASV62", "ASV63", "ASV68", "ASV70", "ASV71", "ASV76", "ASV78", "ASV80", "ASV82", "ASV90", "ASV93", "ASV94", "ASV104", "ASV137", "ASV154", "ASV159", "ASV160", "ASV164", "ASV165", "ASV190", "ASV203", "ASV204", "ASV211", "ASV215", "ASV295", "ASV344", "ASV362", "ASV364", "ASV366", "ASV378", "ASV406", "ASV475", "ASV492", "ASV515", "ASV609", "ASV611", "ASV703", "ASV782", "ASV785", "ASV1165", "ASV1791", "ASV2205"), 
                  rindsection = c("ASV1", "ASV2", "ASV3", "ASV4", "ASV5", "ASV6", "ASV7", "ASV10", "ASV11", "ASV14", "ASV15", "ASV16", "ASV19", "ASV20", "ASV21", "ASV22", "ASV24", "ASV25", "ASV26", "ASV28", "ASV30", "ASV31", "ASV32", "ASV35", "ASV36", "ASV38", "ASV39", "ASV40", "ASV42", "ASV46", "ASV48", "ASV54", "ASV55", "ASV57", "ASV58", "ASV59", "ASV60", "ASV61", "ASV62", "ASV63", "ASV64", "ASV66", "ASV67", "ASV68", "ASV69", "ASV70", "ASV71", "ASV76", "ASV78", "ASV79", "ASV80", "ASV82", "ASV86", "ASV90", "ASV93", "ASV94", "ASV98", "ASV103", "ASV107", "ASV115", "ASV120", "ASV122", "ASV137", "ASV154", "ASV155", "ASV159", "ASV160", "ASV164", "ASV183", "ASV186", "ASV190", "ASV191", "ASV203", "ASV204", "ASV211", "ASV217", "ASV244", "ASV366", "ASV378", "ASV406", "ASV506", "ASV513", "ASV515", "ASV537", "ASV564", "ASV609", "ASV617", "ASV703", "ASV743", "ASV766", "ASV774", "ASV782", "ASV785", "ASV803", "ASV883", "ASV946", "ASV1065", "ASV1165", "ASV1458", "ASV1829"),
                  rind=c("ASV1", "ASV2", "ASV3", "ASV4", "ASV5", "ASV6", "ASV7", "ASV8", "ASV10", "ASV11", "ASV14", "ASV15", "ASV16", "ASV17", "ASV18", "ASV19", "ASV20", "ASV21", "ASV22", "ASV24", "ASV25", "ASV26", "ASV28", "ASV30", "ASV31", "ASV32", "ASV35", "ASV36", "ASV38", "ASV39", "ASV40", "ASV42", "ASV45", "ASV46", "ASV47", "ASV48", "ASV54", "ASV58", "ASV59", "ASV60", "ASV62", "ASV64", "ASV65", "ASV66", "ASV68", "ASV70", "ASV71", "ASV78", "ASV79", "ASV80", "ASV81", "ASV82", "ASV86", "ASV89", "ASV91", "ASV93", "ASV94", "ASV108", "ASV111", "ASV137", "ASV146", "ASV154", "ASV155", "ASV159", "ASV160", "ASV164", "ASV166", "ASV178", "ASV190", "ASV195", "ASV199", "ASV200", "ASV201", "ASV203", "ASV204", "ASV211", "ASV217", "ASV221", "ASV227", "ASV239", "ASV246", "ASV262", "ASV269", "ASV280", "ASV297", "ASV344", "ASV366", "ASV378", "ASV490", "ASV515", "ASV537", "ASV573", "ASV609", "ASV624", "ASV645", "ASV653", "ASV703", "ASV743", "ASV756", "ASV782", "ASV785", "ASV819", "ASV918", "ASV1065", "ASV1090", "ASV1165", "ASV1244", "ASV1432", "ASV1686", "ASV2063", "ASV2205"))

library(UpSetR)
upset(fromList(listInput),order.by = "freq")


## --------------------------------------------------------------------
#cluster dendrogram

prop.cheese<-transform_sample_counts(ps.cheese, function(otu) otu/sum(otu))
sample_names(prop.cheese)<-(sample_data(prop.cheese))$sample
bray <- phyloseq::distance(prop.cheese, method = "bray")

cheese.hclust<- hclust(bray, method="average")
plot(cheese.hclust)






