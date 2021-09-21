## --------------------------------------------------------------------
#load pyloseq and metadata
#bacteria
load("/Volumes/efs/CANR/Ansci/D\'Amico\ Lab/Langsun_las17015/2.\ Bethlehem\ project/Bethlehem2020/Bacteria/dada2/phyloseq.RData")

#fungi
load("/Volumes/efs/CANR/Ansci/D\'Amico\ Lab/Langsun_las17015/2.\ Bethlehem\ project/Bethlehem2020/fungi/dada2/phyloseq.RData")
#remove sample with 0 read
ps<-prune_samples(sample_sums(ps)>0, ps)

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
library(dendextend)
library(cowplot)
library(scales)
library(RColorBrewer)
library(randomcoloR)

## --------------------------------------------------------------------
#dairy alpha

ps.dairy<-subset_samples(ps,source=="dairy")
ps.dairy1<-subset_samples(ps.dairy, sampletype!="cheese")
ps.dairy2<-subset_samples(ps.dairy, day=="60")
ps.dairy<-merge_phyloseq(ps.dairy1,ps.dairy2)


p = plot_richness(ps.dairy, x="sampletype", measures=c("Shannon"))
p$data$sampletype <- as.character(p$data$sampletype)
p$data$sampletype[p$data$sampletype=="cheese"] <- "cheese_d60"
p$data$sampletype <- factor(p$data$sampletype, 
                            levels=c("rawmilk","filtered_milk","milk_before_ripening",
                                     "milk_stored_overnight","ripened_milk","mixed_milk",
                                     "whey","cheese_curd","cheese_d60"))

p4 = p + geom_boxplot() + theme_bw() + 
  theme(axis.title=element_blank(),axis.text.x = element_text(angle = -45),legend.position = "none",panel.background = element_rect(fill = 'white', colour = 'black'))


## --------------------------------------------------------------------
#dairy beta 

prop.dairy <- transform_sample_counts(ps.dairy, function(otu) otu/sum(otu))

ord.pcoa.bray <- ordinate(prop.dairy, method="PCoA", distance="bray")

p = plot_ordination(prop.dairy, ord.pcoa.bray, color="sampletype")
p$data$sampletype <- as.character(p$data$sampletype)
p$data$sampletype[p$data$sampletype=="cheese"] <- "cheese_d60"
p$data$sampletype <- factor(p$data$sampletype, 
                            levels=c("rawmilk","filtered_milk","milk_before_ripening",
                                     "milk_stored_overnight","ripened_milk","mixed_milk",
                                     "whey","cheese_curd","cheese_d60"))

p5 = p + theme_bw() + 
  theme(legend.title=element_blank(),legend.text = element_text(size=10)) + 
  ylab("PCoA 1 (13.6%)") + 
  xlab("PCoA 2 (33.7%)")


#save(p1, file = "p1.RData") 
#save(p2, file = "p2.RData") 
#save(p3, file = "p3.RData") 
#save(p4, file = "p4.RData") 
#save(p5, file = "p5.RData") 
#save(p6, file = "p6.RData") 

p4 = p4 + theme(text = element_text(size = 7))
p5 = p5 +theme(text = element_text(size = 7),
               legend.text = element_text(size=5),
               legend.key.size = unit(0.3, "cm"),
               legend.margin=margin(0,0,0,0))

p6[["layers"]][[3]][["mapping"]][["size"]] = 1.5
p6[["layers"]][[2]][["mapping"]][["size"]] = 1.5
p6[["layers"]][[1]][["mapping"]][["size"]] = 0.3

p3[["layers"]][[3]][["mapping"]][["size"]] = 1.5
p3[["layers"]][[2]][["mapping"]][["size"]] = 1.5
p3[["layers"]][[1]][["mapping"]][["size"]] = 0.3

fig4 <- plot_grid(p1, p2, p3, p4, p5, p6,
                  labels = c("A", "B","C","D","E","F"),
                  ncol = 3, nrow = 2, rel_widths = c(1,2,1.1),
                  label_size = 9)

ggsave(fig4, file="~/Desktop/fig4.eps",height = 5, width = 6.86,units = "in",dpi = 300)



## --------------------------------------------------------------------
#dairy bar plot
merge.dairy <- merge_samples(ps.dairy, "sampletype")
merge2.dairy <- tax_glom(merge.dairy, "Genus")
prop.merge.dairy <- transform_sample_counts(merge2.dairy, function(otu) otu/sum(otu))
melt.dairy <- psmelt(prop.merge.dairy)

melt.dairy1<-melt.dairy
melt.dairy1$Genus <- as.character(melt.dairy1$Genus)
melt.dairy1$Genus[melt.dairy1$Abundance < 0.02] <- "All others (1-2%)"
melt.dairy1$Genus[melt.dairy1$Abundance < 0.01] <- "<1% abun."
melt.dairy1$Genus[is.na(melt.dairy1$Genus)] <- "unclassified_bacteria"
#for fungi
#melt.dairy1$Genus <- gsub('g__', '', melt.dairy1$Genus)
melt.dairy1$Genus <- reorder(melt.dairy1$Genus, -melt.dairy1$Abundance)

melt.dairy1$Sample[melt.dairy1$Sample == "cheese"] <- "cheese_d60"
melt.dairy1$Sample <- factor(melt.dairy1$Sample, 
                    c("rawmilk","filtered_milk","milk_before_ripening",
                      "milk_stored_overnight","ripened_milk","mixed_milk",
                      "whey","cheese_curd","cheese_d60"))


d = length(unique(melt.dairy1$Genus))
palette <- distinctColorPalette(d)

p2 <- ggplot(melt.dairy1, aes(x=Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = palette) + theme_bw() +
  theme(axis.text.x = element_text(vjust=0.5, size=11, angle = 270),
        axis.title.x = element_blank(),
        legend.text=element_text(size=10),
        legend.title = element_blank()) + 
  guides(fill = guide_legend(ncol = 2)) +
  ylab("Relative Abundance")

#save(p1, file = "p1.RData") 
#save(p2, file = "p2.RData") 
p1 = p1 + theme(text = element_text(size=6),
                axis.text.x = element_text(vjust=0.5,size = 6, angle = -90),
                axis.title.x = element_blank(),
                legend.text=element_text(size=5),
                legend.key.size = unit(0.4, "cm"), 
                legend.title = element_blank(),
                #title = element_text(size=12,face = "bold")
                )
p2 = p2 + theme(text = element_text(size=6),
                axis.text.x = element_text(vjust=0.5,size = 6, angle = -90),
                axis.title.x = element_blank(),
                legend.text=element_text(size=5),
                legend.key.size = unit(0.4, "cm"), 
                legend.title = element_blank(),
                #title = element_text(size=12,face = "bold")
                )
legend1 <- get_legend(p1)
legend2 <- get_legend(p2)

fig5<-plot_grid(p1+theme(legend.position = "none"),
                legend1,
                p2+theme(legend.position = "none"),
                labels = c("A", "","B"),
                legend2,nrow = 2,ncol = 2,rel_widths = c(1,0.7))

ggsave(fig5, file="~/Desktop/fig5.eps",height = 6, width = 6.86,units = "in",dpi = 300)


## --------------------------------------------------------------------
#cheese top20 ASV

ps.cheese<-subset_samples(ps,sampletype=="cheese")

top <- names(sort(taxa_sums(ps.cheese), decreasing=TRUE))[1:20]
ps.top <- transform_sample_counts(ps.cheese, function(otu) otu/sum(otu))
ps.top <- prune_taxa(top, ps.top)


dend <- as.dendrogram(hclust(dist(t(otu_table(ps.top)))))
dend_data <- dendro_data(dend)
segment_data <- with(
  segment(dend_data), 
  data.frame(x = y, y = x, xend = yend, yend = xend))

#height settings
asv_pos_table <- with(
  dend_data$labels, 
  data.frame(y_center = x, OTU = as.character(label), height = 1))

asv_axis_limits <- with(
  asv_pos_table, 
  c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))
) + 0.1 * c(-1, 1)

#dendrogram
p.dendro2=ggplot(segment_data) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  scale_x_reverse(expand = c(0, 0)) + 
  scale_y_continuous(breaks = asv_pos_table$y_center, 
                     labels = asv_pos_table$OTU, 
                     limits = asv_axis_limits, 
                     expand = c(0, 0)) +
  labs(x = "", y = "", colour = "", size = "") +
  theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.line = element_blank())

#heatmap dataframe
top.melt <- psmelt(ps.top) 

top.new <- top.melt
top.new$Abundance_per <- (top.melt$Abundance)*100

top.new$section <- factor(top.new$section,c("core","middlesection","rindsection","rind"))
top.new$day <- factor(top.new$day,c("0","4","7","14","21","44","60"))

#dataframe for heatmap
top.new1<-left_join(top.new,asv_pos_table,by="OTU")

#bacteria
#heatmap dairy databse 
species<-read.table(file="/Volumes/efs/CANR/Ansci/D\'Amico\ Lab/Langsun_las17015/2.\ Bethlehem\ project/Bethlehem2020/Bacteria/db_dairy/asv.txt",header = T)
asv_pos_table<-left_join(asv_pos_table,species,by="OTU")

#fungi
species<-as.data.frame(tax_table(ps.top))
species$Species <- paste(species$Genus,species$Species)
species$Species <- gsub('g__', '', species$Species)
species$Species <- gsub('s__', '', species$Species)
species$OTU <- rownames(species)
asv_pos_table<-left_join(asv_pos_table,species,by="OTU")


myPalette <- RColorBrewer::brewer.pal(5, "Set1")

cheese.heatmap2 <- ggplot(data = top.new1, mapping = aes(x = day,y = y_center,fill = Abundance_per,height=height)) + 
  geom_tile() + 
  scale_fill_viridis_c(option = 'D', trans = "log10",  breaks = c(0.01, 0.1, 1, 10)) + 
  facet_grid(.~section,switch = "both", scales = "free", space = "free") + 
  labs(fill = "% relative abundance") +
  scale_y_continuous(breaks = asv_pos_table[, "y_center"], 
                     labels = asv_pos_table$Species,
                     position = "right",
                     limits = asv_axis_limits, 
                     expand = c(0, 0)) +
  theme(axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=270),
        axis.text.y = element_text(colour = myPalette[asv_pos_table$Phylum]),
        axis.title = element_blank(),
        strip.placement = "outside", 
        strip.background = element_blank())

cheese.heatmap1 = cheese.heatmap1 + 
  theme(axis.text.y = element_text(size = 5.5),
        legend.title = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 5),
        axis.text = element_text(size = 5.5),
        strip.text.x = element_text(size = 6))

cheese.heatmap2 = cheese.heatmap2 + 
  theme(axis.text.y = element_text(size = 5.5),
        legend.title = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 5),
        axis.text = element_text(size = 5.5),
        strip.text.x = element_text(size = 6))

p1 = plot_grid(p.dendro1,p.dendro2,labels = c("A","B"),nrow = 2,align = "hv")
p2 = plot_grid(cheese.heatmap1,cheese.heatmap2,labels = c("",""),nrow = 2,align = "hv")


#save(p.dendro1, file = "p.dendro1.RData") 
#save(p.dendro2, file = "p.dendro2.RData") 
#save(cheese.heatmap1, file = "cheese.heatmap1.RData") 
#save(cheese.heatmap2, file = "cheese.heatmap2.RData") 
fig6<-plot_grid(p1,p2,
                align="h",axis = "tb",
                ncol = 2, rel_widths = c(0.3,1.2))
ggsave(fig6, file="~/Desktop/fig6.tiff",height = 6, width = 6.86,units = "in",dpi = 300)

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
sample_names(prop.cheese)<-paste((sample_data(prop.cheese))$section,(sample_data(prop.cheese))$day)

#unifrac
set.seed(1)
unifrac <- phyloseq::distance(prop.cheese, method = "unifrac",weighted=T)
cheese.hclust2<- hclust(unifrac, method="average")
hcd <- as.dendrogram(cheese.hclust2)

labels(hcd) = sub("middlesection", " ms", labels(hcd))
labels(hcd) = sub("rindsection", " rs", labels(hcd))
labels(hcd) = sub("rind", " rind", labels(hcd))
labels(hcd) = sub("core", " core", labels(hcd))

dend1 <- hcd %>%
  set("branches_k_color", k=4) %>% set("labels_cex", 0.6) %>% set("branches_lwd", 0.6)
p3=ggplot(dend1, horiz = TRUE,theme = theme_minimal())+ 
  theme(axis.text.y = element_blank(),axis.title = element_blank())

dend2 <- hcd %>%
  set("branches_k_color", k=2) %>% set("labels_cex", 0.6) %>% set("branches_lwd", 0.6)
p6=ggplot(dend2, horiz = TRUE,theme = theme_minimal())+ 
  theme(axis.text.y = element_blank(),axis.title = element_blank())




