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

p1$layers
p1$layers[[3]] <- NULL
p1 = p1 + 
  theme(text = element_text(size=10),axis.text.x = element_text(angle = -15)) +
  geom_signif(comparisons = list(c("cheesemaking", "ripening"),
                                 c("milking", "cheesemaking"),
                                 c("milking", "ripening")), 
              map_signif_level=TRUE,step_increase = 0.17,color = "black",
              size = 0.3,textsize = 1.7)

p3$layers
p3$layers[[3]] <- NULL
p3 = p3 + 
  theme(text = element_text(size=10),axis.text.x = element_text(angle = -15)) + 
  geom_signif(comparisons = list(c("cheesemaking", "ripening"),
                                 c("milking", "cheesemaking"),
                                 c("milking", "ripening")), 
              map_signif_level=TRUE,step_increase = 0.17,color = "black",
              size = 0.3,textsize = 1.7)

p2 = p2 + theme(legend.title=element_text(size=11),
                legend.text = element_text(size=10))
p4 = p4 + theme(legend.title=element_text(size=11),
                legend.text = element_text(size=10))


fig1 <- plot_grid(p1, p2,p3,p4,
                  labels = c("A", "B","C","D"),
                  ncol = 2, nrow = 2, rel_widths = c(1, 1.5),label_size = 12)



ggsave(fig1, file="~/Desktop/fig1.eps",height = 5, width = 6.86,units = "in",dpi = 300)


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

p1 = p1  + theme(text = element_text(size=6),
                axis.text.x = element_text(vjust=0.5,size = 6, angle = -90),
           axis.title.x = element_blank(),
           legend.text=element_text(size=4),
           legend.key.size = unit(0.4, "cm"), 
           legend.title = element_blank(),
           #title = element_text(size=12,face = "bold")
           )

p2 = p2 + theme(text = element_text(size=6),
                axis.text.x = element_text(vjust=0.5,size = 6, angle = -90),
                axis.title.x = element_blank(),
                legend.text=element_text(size=4),
                legend.key.size = unit(0.4, "cm"), 
                legend.title = element_blank(),
                #title = element_text(size=12,face = "bold")
                )
  
legend1 <- get_legend(p1)
fig2<-plot_grid(p1+theme(legend.position = "none"),
                legend1,
                p2+theme(legend.position = "none"),
                labels = c("A", "","B"),
                legend2,nrow = 2,ncol = 2,rel_widths = c(1,0.6),label_size = 12)

ggsave(fig2, file="~/Desktop/fig2.eps",height = 7, width = 6.86,units = "in",dpi = 300)
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
#differential abundance analysis: deseq2

#transform
ps.2rm = subset_samples(ps.env,room != "milking")
ps.2rm.f<-prune_taxa(taxa_sums(ps.2rm)>0, ps.2rm)
otu_table(ps.2rm.f) <- otu_table(ps.2rm.f) + 1

#run deseq2
deseq = phyloseq_to_deseq2(ps.2rm.f, ~room)
deseq = DESeq(deseq, test="Wald", fitType="parametric")

#result
res = results(deseq, cooksCutoff = FALSE)
res.df = as.data.frame(res@listData)
# get a vector of significant and non-significant rows
res$effect.sig[abs(res$log2FoldChange) <= 2] <- "non.sig"
res$effect.sig[res$log2FoldChange > 2] <- "sig.pos"
res$effect.sig[res$log2FoldChange < -2] <- "sig.neg"


sigtab = cbind(as(res, "data.frame"), as(tax_table(ps.2rm.f)[rownames(res), ], "matrix"))  
sigtab = sigtab[which(sigtab$padj < 0.01), ]

sig.gen = as.character(unique(sigtab$Genus[sigtab$effect.sig!="non.sig"]))
x.sig.gen = sigtab[sigtab$Genus %in% sig.gen, ]
x.sig.gen = subset(x.sig.gen,x.sig.gen$Genus !="NA")

x.sig.gen$Class <- gsub('c__', '', x.sig.gen$Class)
x.sig.gen$Genus <- gsub('g__', '', x.sig.gen$Genus)
#plot
p4 = ggplot(d, aes(x = -log2FoldChange,y = Class)) +
  geom_jitter(aes(color = effect.sig),size = 0.2,position = position_jitter(width=0.1, height=0)) +
  scale_color_manual(values = c("non.sig" = "grey", "sig.neg" = "red","sig.pos" = "blue")) + 
  xlim(-11.5, 11.5) +
  theme_bw() + 
  theme(text = element_text(size=6),
        axis.title = element_blank(),legend.position = "none",
        plot.title = element_text(hjust = 0.5,size = 6)) +
  ggtitle("Milking barn : Cheesemaking room")
  #facet_grid(Class~.,switch = "both", scales = "free", space = "free")

#save(p1,p2,p3,p4,p5,p6, file ="fig3.RData")

fig3<-plot_grid(p1,p2,p3,p4,p5,p6,
                 labels = c("A","","","B","",""),
                 nrow = 2,ncol = 3,
                 label_size = 9)


ggsave(fig3,file="~/Desktop/fig3.eps",height = 5, width = 6.875,units = "in",dpi = 300)

## --------------------------------------------------------------------
#source  beta 

prop.ps <- transform_sample_counts(ps, function(otu) otu/sum(otu))

ord.nmds.bray <- ordinate(prop.ps, method="NMDS", distance="bray")

p = plot_ordination(prop.ps, ord.nmds.bray, color="room",shape="source")
p$data$room <- as.character(p$data$room)
p$data$room <- factor(p$data$room, levels=c("milking","cheesemaking","ripening"))

p2 = p + theme(legend.text=element_text(size=8),panel.background = element_rect(fill = 'white', colour = 'black'), panel.grid=element_blank())

## --------------------------------------------------------------------
#differential abundance analysis: compositional  
#centered log ratio (clr) transform
#compare using aldex.2
library(ALDEx2)
ps.2rm = subset_samples(ps.env,room != "ripening")

df = as.data.frame(otu_table(ps.2rm))
map = as.data.frame(sample_data(ps.2rm))
taxa = as.data.frame(tax_table(ps.2rm))
rownames(df) = map$sample
conds <- as.character(map$room)
df = t(df)

# apply aldex2
x.b.m_c <- aldex(df, conds, mc.samples=128, test="t", effect=TRUE, 
                 include.sample.summary=FALSE, denom="all", verbose=FALSE)
x.b.c_r <- aldex(df, conds, mc.samples=128, test="t", effect=TRUE,
                 include.sample.summary=FALSE, denom="all", verbose=FALSE)
x.b.m_r <- aldex(df, conds, mc.samples=128, test="t", effect=TRUE,
                 include.sample.summary=FALSE, denom="all", verbose=FALSE)
x.f.m_c <- aldex(df, conds, mc.samples=128, test="t", effect=TRUE,
                 include.sample.summary=FALSE, denom="all", verbose=FALSE)
x.f.c_r <- aldex(df, conds, mc.samples=128, test="t", effect=TRUE,
                 include.sample.summary=FALSE, denom="all", verbose=FALSE)
x.f.m_r <- aldex(df, conds, mc.samples=128, test="t", effect=TRUE,
                 include.sample.summary=FALSE, denom="all", verbose=FALSE)

x.all = x.f.m_c
# get a vector of significant and non-significant rows
x.all$effect.sig[abs(x.all$effect) <= 1] <- "non.sig"
x.all$effect.sig[x.all$effect > 1] <- "sig.pos"
x.all$effect.sig[x.all$effect < -1] <- "sig.neg"

# select only dignificant genera in table for plotting
x.all$otu = rownames(x.all)
taxa$otu = rownames(taxa)
x.all.taxa = left_join(x.all,taxa,by="otu")

sig.gen = as.character(unique(x.all.taxa$Genus[x.all.taxa$effect.sig!="non.sig"]))
x.sig.gen = x.all.taxa[x.all.taxa$Genus %in% sig.gen, ]
x.sig.gen = subset(x.sig.gen,x.sig.gen$Genus !="NA")
#create plot
p4 = ggplot(x.sig.gen, aes(x = effect,y = Genus)) +
  geom_jitter(aes(color = effect.sig),position = position_jitter(width=0.1, height=0)) +
  scale_color_manual(values = c("non.sig" = "grey", "sig.pos" = "red","sig.neg" = "blue")) + 
  xlim(-2.3, 2.3) +
  theme_bw() + 
  theme(axis.title = element_blank(),legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Milking barn : Cheesemaking room")


p6 = p6 + theme(text = element_text(size=6),
                plot.title = element_text(size = 6))


figs1<-plot_grid(p1,p2,p3,p4,p5,p6,
                labels = c("A","","","B","",""),
                nrow = 2,ncol = 3,
                label_size = 9)


ggsave(figs1,file="~/Desktop/figs1.eps",height = 4, width = 6,units = "in",dpi = 300)
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


























