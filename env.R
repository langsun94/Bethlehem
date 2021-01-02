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
library(reshape2)
library(ggsignif)
library(pairwiseAdonis)
library(indicspecies)
library(cowplot)
library(DESeq2)

## --------------------------------------------------------------------
#environment alpha

ps.env<-subset_samples(ps,source=="direct"|source=="indirect")

p = plot_richness(ps.env, x="room", measures=c("Shannon"),color = "room")
p$data$room <- as.character(p$data$room)
p$data$room <- factor(p$data$room, levels=c("milking","cheesemaking","ripening"))

p3 = p + geom_violin() + theme_bw() +
  theme(axis.title=element_blank(),axis.text.x = element_text(angle = -45),legend.position = "none") +  
  geom_signif(comparisons = list(c("cheesemaking", "ripening"),c("milking", "cheesemaking"),c("milking", "ripening")), map_signif_level=TRUE,step_increase = 0.1,color = "black") 


## --------------------------------------------------------------------
#environment beta 

prop.env <- transform_sample_counts(ps.env, function(otu) otu/sum(otu))

ord.pcoa.bray <- ordinate(prop.env, method="PCoA", distance="bray")

p = plot_ordination(prop.env, ord.pcoa.bray, color="room")
p$data$room <- as.character(p$data$room)
p$data$room <- factor(p$data$room, levels=c("milking","cheesemaking","ripening"))

p4 = p + theme_bw() + 
  theme(legend.title=element_text(size=12),legend.text = element_text(size=11)) + 
  ylab("PCoA 1 (13.4%)") + 
  xlab("PCoA 2 (17.1%)")

#save(p1, file = "p1.RData") 
#save(p2, file = "p2.RData") 
#save(p3, file = "p3.RData") 
#save(p4, file = "p4.RData") 
#fig1 <- plot_grid(p1, p2,p3,p4,
#                  labels = c("A", "B","C","D"),
#                  ncol = 2, nrow = 2, rel_widths = c(1, 1.5))
#ggsave(fig1, file="~/Desktop/fig1.jpeg",height = 8.5, width = 10)


## --------------------------------------------------------------------
#environment beta permanova

set.seed(1)
# Calculate bray curtis distance matrix
bray <- phyloseq::distance(prop.env, method = "jaccard")
# make a data frame from the sample_data
sampledf <- data.frame(sample_data(prop.env))
# Adonis test
adonis(bray ~ room, data = sampledf)
#pairwise
pair.adonis = pairwise.adonis(otu_table(prop.env),sampledf$room)


## --------------------------------------------------------------------
#environment bar plot

merge.env <- merge_samples(ps.env, "sampletype")
merge2.env <- tax_glom(merge.env, "Genus")
prop.merge.env <- transform_sample_counts(merge2.env, function(otu) otu/sum(otu))
melt.env <- psmelt(prop.merge.env)

melt.env1<-melt.env
melt.env1$Genus <- as.character(melt.env1$Genus)
melt.env1$Genus[melt.env1$Abundance < 0.05] <- "All others (1-5%)"
melt.env1$Genus[melt.env1$Abundance < 0.01] <- "<1% abun."
melt.env1$Genus[is.na(melt.env1$Genus)] <- "unclassified"
#for fungi
#melt.env1$Genus <- gsub('g__', '', melt.env1$Genus)
melt.env1$Genus <- reorder(melt.env1$Genus, -melt.env1$Abundance)

melt.env1$room <- as.character(melt.env1$room)
melt.env1$room[melt.env1$room == 3] <- "ripening"
melt.env1$room[melt.env1$room == 1] <- "cheesemaking"
melt.env1$room[melt.env1$room == 2] <- "milking"
melt.env1$room<-factor(melt.env1$room, c("milking","cheesemaking","ripening"))


library(randomcoloR)
d = length(unique(melt.env1$Genus))
palette <- distinctColorPalette(d)

p1 <- ggplot(melt.env1, aes(x=Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette) + theme_bw() +
  theme(axis.text.x = element_text(vjust=0.5, size=11, angle = -90),
        axis.title.x = element_blank(),
        legend.text=element_text(size=10),
        #legend.key.size = unit(0.4, "cm"), 
        legend.title = element_blank(),
        #title = element_text(size=12,face = "bold")
        ) + 
  guides(fill = guide_legend(ncol = 2)) + 
  facet_grid(.~room,switch = "both", scales = "free", space = "free") +
  ylab("Relative Abundance") 

#save(p1, file = "p1.RData") 
#save(p2, file = "p2.RData") 
#legend2 <- get_legend(p2)
#fig2<-plot_grid(p1+theme(legend.position = "none"),
#                legend1,
#                p2+theme(legend.position = "none"),
#                labels = c("A", "","B"),
#                legend2,nrow = 2,ncol = 2,rel_widths = c(1,0.4))
#
#ggsave(fig2, file="~/Desktop/fig2.jpeg",height = 12, width = 17)
## --------------------------------------------------------------------
#environment prevalent genera

preva.env<-filter_taxa(prop.merge.env, function(otu) sum(otu > 0.01) >= (0.3*length(otu)), TRUE)

ps.milking<-subset_samples(prop.merge.env, room == "2") # 2==milking
preva.milk<-filter_taxa(ps.milking, function(otu) sum(otu > 0.01) >= (0.5*length(otu)), TRUE)

ps.cheesemak<-subset_samples(prop.merge.env, room == "1") # 1==cheesemaking
preva.cheesemak<-filter_taxa(ps.cheesemak, function(otu) sum(otu > 0.01) >= (0.5*length(otu)), TRUE)

ps.aging<-subset_samples(prop.merge.env, room == "3") # 3==aging
preva.aging<-filter_taxa(ps.aging, function(otu) sum(otu > 0.01) >= (0.5*length(otu)), TRUE)


## --------------------------------------------------------------------
#environment indicator species

#genera
merge2.env.f<-prune_taxa(taxa_sums(merge2.env)>0, merge2.env)
room = (as.data.frame(sample_data(merge2.env.f)))$room
otu<-as.data.frame(otu_table(merge2.env.f))

indic = multipatt(otu, room, func = "r.g", control = how(nperm=999))
summary(indic)
sign2=as.data.frame(indic[["sign"]])
sign2$OTU = rownames(sign2)

prop.merge.env <- transform_sample_counts(merge2.env.f, function(otu) otu/sum(otu))
melt.env <- psmelt(prop.merge.env)

#milking barn
sign.melt <- left_join(sign2, melt.env, by = "OTU")
milk2 = subset(sign.melt,sign.melt$p.value < 0.05)
milk2 = subset(milk2,milk2$room == "2")
milk2 = subset(milk2,milk2$s.2 == "1")
milk2 = subset(milk2,milk2$Abundance > 0.01)

#ripening cellar
sign.melt <- left_join(sign2, melt.env, by = "OTU")
ripen2 = subset(sign.melt,sign.melt$p.value < 0.05)
ripen2 = subset(ripen2,ripen2$room == "3")
ripen2 = subset(ripen2,ripen2$s.3 == "1")
ripen2 = subset(ripen2,ripen2$Abundance > 0.01)

## --------------------------------------------------------------------
#environment differential abundance analysis

#genera

otu_table(merge2.env.f) <- otu_table(merge2.env.f) + 1

#milking vs. cheesemaking
group1 = subset_samples(merge2.env.f,room != "3")
deseq = phyloseq_to_deseq2(group1, ~ room)
deseq = DESeq(deseq, test="Wald", fitType="parametric")
res = results(deseq, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab1 = cbind(as(sigtab, "data.frame"), as(tax_table(merge2.env.f)[rownames(sigtab), ], "matrix"))

#milking vs. ripening
group2 = subset_samples(merge2.env.f,room != "1")
deseq = phyloseq_to_deseq2(group2, ~ room)
deseq = DESeq(deseq, test="Wald", fitType="parametric")
res = results(deseq, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab2 = cbind(as(sigtab, "data.frame"), as(tax_table(merge2.env.f)[rownames(sigtab), ], "matrix"))

#cheesemaking vs. ripening
group3 = subset_samples(merge2.env.f,room != "2")
deseq = phyloseq_to_deseq2(group3, ~ room)
deseq = DESeq(deseq, test="Wald", fitType="parametric")
res = results(deseq, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab3 = cbind(as(sigtab, "data.frame"), as(tax_table(merge2.env.f)[rownames(sigtab), ], "matrix"))

##### milking
sigtab1$OTU = rownames(sigtab1)
sigtab2$OTU = rownames(sigtab2)
milking = inner_join(subset(sigtab2,sigtab2$log2FoldChange < 0),subset(sigtab1,sigtab1$log2FoldChange > 0), by="OTU")
milking = subset(melt.env1,OTU %in% milking$OTU)

#### cheesemaking
sigtab1$OTU = rownames(sigtab1)
sigtab3$OTU = rownames(sigtab3)
cheesemaking = inner_join(subset(sigtab3,sigtab3$log2FoldChange < 0),subset(sigtab1,sigtab1$log2FoldChange < 0), by="OTU")
cheesemaking = subset(melt.env1,OTU %in% cheesemaking$OTU)

#### ripening
sigtab2$OTU = rownames(sigtab2)
sigtab3$OTU = rownames(sigtab3)
ripening = inner_join(subset(sigtab2,sigtab2$log2FoldChange > 0),subset(sigtab3,sigtab3$log2FoldChange > 0), by="OTU")
ripening = subset(melt.env1,OTU %in% ripening$OTU)

## --------------------------------------------------------------------
#source  beta 

prop.ps <- transform_sample_counts(ps, function(otu) otu/sum(otu))

ord.nmds.bray <- ordinate(prop.ps, method="NMDS", distance="bray")

p = plot_ordination(prop.ps, ord.nmds.bray, color="room",shape="source")
p$data$room <- as.character(p$data$room)
p$data$room <- factor(p$data$room, levels=c("milking","cheesemaking","ripening"))

p2 = p + theme(legend.text=element_text(size=8),panel.background = element_rect(fill = 'white', colour = 'black'), panel.grid=element_blank())


## --------------------------------------------------------------------
#direct vs.indirect beta box plot

# calc distances
bray <- phyloseq::distance(prop.ps, method = "bray")
wu.m = melt(as.matrix(bray))

# remove self-comparisons
wu.m = wu.m %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor, as.character)

# get sample data (S4 error OK and expected)
sd = sample_data(prop.ps) %>%
  select(Barcode, source) %>%
  mutate_if(is.factor,as.character)

# combined distances with sample data
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(wu.m, sd, by = "Var1")

colnames(sd) = c("Var2", "Type2")
wu.sd = left_join(wu.sd, sd, by = "Var2")

#only keep comparisons to dairy
wu.f = subset(wu.sd,wu.sd$Type2=="dairy")
wu.f = subset(wu.f,wu.f$Type1 !="dairy")

#plot
p = ggplot(wu.f, aes(x = Type1, y = value, color=Type1)) +
  theme_bw() +
  geom_point() +
  theme(legend.position = "none") +
  ylab("Bray-curtis distance") +
  xlab("Disimilarity to dairy samples")


























