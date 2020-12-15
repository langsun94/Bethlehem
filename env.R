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
library(reshape2)
library(ggsignif)
library(pairwiseAdonis)
library(indicspecies)

## --------------------------------------------------------------------
#environment alpha

ps.env<-subset_samples(ps,source=="direct"|source=="indirect")

p = plot_richness(ps.env, x="room", measures=c("Shannon", "Simpson"),color = "room")
p$data$room <- as.character(p$data$room)
p$data$room <- factor(p$data$room, levels=c("milking","cheesemaking","ripening"))

p1 = p + geom_violin() + 
  theme(axis.title=element_blank(),axis.text.x = element_text(angle = 270),legend.position = "none",panel.background = element_rect(fill = 'white', colour = 'black')) +  
  geom_signif(comparisons = list(c("milking", "cheesemaking"),c("milking", "ripening")), map_signif_level=TRUE,step_increase = 0.1,color = "black") 

save(p1, file = "fig1a.RData") 

## --------------------------------------------------------------------
#environment beta 

prop.env <- transform_sample_counts(ps.env, function(otu) otu/sum(otu))

ord.nmds.bray <- ordinate(prop.env, method="NMDS", distance="bray")

p = plot_ordination(prop.env, ord.nmds.bray, color="room")
p$data$room <- as.character(p$data$room)
p$data$room <- factor(p$data$room, levels=c("milking","cheesemaking","ripening"))

p2 = p + theme(legend.text=element_text(size=8),legend.title = element_blank(),legend.position = c(0.8, 0.9),panel.background = element_rect(fill = 'white', colour = 'black'), panel.grid=element_blank(),legend.background = element_rect(fill="white",colour ="black",size=0.5, linetype="solid"))

save(p2, file = "fig1b.RData") 


## --------------------------------------------------------------------
#environment beta permanova

set.seed(1)
# Calculate bray curtis distance matrix
bray <- phyloseq::distance(prop.env, method = "bray")
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
melt.env1$Genus[melt.env1$Abundance < 0.01] <- "< 1% abund."
melt.env1$Genus[is.na(melt.env1$Genus)] <- "unclassified_bacteria"
melt.env1$Genus <- reorder(melt.env1$Genus, -melt.env1$Abundance)

melt.env1$room <- as.character(melt.env1$room)
melt.env1$room[melt.env1$room == 3] <- "ripening"
melt.env1$room[melt.env1$room == 1] <- "cheesemaking"
melt.env1$room[melt.env1$room == 2] <- "milking"
melt.env1$room<-factor(melt.env1$room, c("milking","cheesemaking","ripening"))

library(randomcoloR)
d = length(unique(melt.env1$Genus))
palette <- distinctColorPalette(d)

p.s1 <- ggplot(melt.env1, aes(x=Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette) +
  theme(axis.text.x = element_text(vjust=0.5, size=12, angle = 270),axis.title = element_text(size=12),legend.text=element_text(size=12),legend.key.size = unit(0.4, "cm"), legend.title = element_blank(),title = element_text(size=12,face = "bold"),panel.background = element_blank(), panel.grid=element_blank(), axis.line=element_line('black')) + guides(fill = guide_legend(ncol = 2)) + facet_grid(.~room,switch = "both", scales = "free", space = "free")+  ylab("Relative Abundance (Genera > 1%) \n") + xlab("Environmental samples") 

save(p.s1, file = "figs1.RData") 


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
#environment prevalent genera

merge.env <- merge_samples(ps.env, "sampletype")
merge2.env <- tax_glom(merge.env, "Genus")
prop.merge.env <- transform_sample_counts(merge2.env, function(otu) otu/sum(otu))
melt.env <- psmelt(prop.merge.env)

melt.env1<-melt.env
melt.env1$Genus <- as.character(melt.env1$Genus)
melt.env1$Genus[melt.env1$Abundance < 0.01] <- "< 1% abund."
melt.env1$Genus[is.na(melt.env1$Genus)] <- "unclassified_bacteria"
melt.env1$Genus[(melt.env1$Genus!="Acinetobacter")&
                  (melt.env1$Genus!="Ruminococcaceae_UCG-005")&
                  (melt.env1$Genus!="Corynebacterium_1")&
                  (melt.env1$Genus!="Lactococcus")&
                  (melt.env1$Genus!="Kocuria")&
                  (melt.env1$Genus!="Lactobacillus")&
                  (melt.env1$Genus!="Brachybacterium")&
                  (melt.env1$Genus!="Staphylococcus")&
                  (melt.env1$Genus!="Brevibacterium")&
                  (melt.env1$Genus!="Streptomyces")&
                  (melt.env1$Genus!="Pseudonocardia")&
                  (melt.env1$Genus!="< 1% abund.")&
                  (melt.env1$Genus!="unclassified_bacteria")] <- "all others"
melt.env1$Genus <- reorder(melt.env1$Genus, melt.env1$Abundance)


melt.env1$room <- as.character(melt.env1$room)
melt.env1$room[melt.env1$room == 3] <- "ripening"
melt.env1$room[melt.env1$room == 1] <- "cheesemaking"
melt.env1$room[melt.env1$room == 2] <- "milking"
melt.env1$room<-factor(melt.env1$room, c("milking","cheesemaking","ripening"))

palette <- c("all others"="grey",
             "< 1% abund."="antiquewhite",
             "Acinetobacter"="#ED2E49",
             "Ruminococcaceae_UCG-005"="#45A6D1",
             "Corynebacterium_1"="#468423",
             "Lactococcus"="#885EA5",
             "Kocuria"="#EA7C52",
             "Lactobacillus"="#F3D01F",
             "Brachybacterium" = "#B5664C",
             "Staphylococcus"="#EB65AF",
             "Brevibacterium"= "#6EAB97",
             "Streptomyces"="#8FCC14",
             "Pseudonocardia"="#6A8ACA")

p3 <- ggplot(melt.env1, aes(x=Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette) +
  theme(axis.text.x = element_text(vjust=0.5, size=12, angle = 270),axis.title = element_text(size=12),legend.text=element_text(size=12),legend.key.size = unit(0.4, "cm"), legend.title = element_blank(),title = element_text(size=12,face = "bold"),panel.background = element_blank(), panel.grid=element_blank(), axis.line=element_line('black')) + facet_grid(.~room,switch = "both", scales = "free", space = "free")+  ylab("Relative Abundance (Genera > 1%) \n") + xlab("Environmental samples") 

save(p3, file = "fig2a.RData") 


## --------------------------------------------------------------------
#environment indicator species

room = (as.data.frame(sample_data(ps.env)))$room
otu<-as.data.frame(otu_table(ps.env))
indic = multipatt(otu, room, func = "r.g", control = how(nperm=9999))
summary(indic)


## --------------------------------------------------------------------
#source  beta 

prop.ps <- transform_sample_counts(ps, function(otu) otu/sum(otu))

ord.nmds.bray <- ordinate(prop.ps, method="NMDS", distance="bray")

p = plot_ordination(prop.ps, ord.nmds.bray, color="room",shape="source")
p$data$room <- as.character(p$data$room)
p$data$room <- factor(p$data$room, levels=c("milking","cheesemaking","ripening"))

p2 = p + theme(legend.text=element_text(size=8),panel.background = element_rect(fill = 'white', colour = 'black'), panel.grid=element_blank())


## --------------------------------------------------------------------
#source  beta box plot

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


























