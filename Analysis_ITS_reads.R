library("phyloseq")
library("ggplot2")
library("readxl")
library("microbiome")
library("vegan")
library("dplyr")
library("tibble")
library("fantaxtic")
library("decontam")
library("tidyr")
library("ggpubr")
library("vegan")
library("phyloseq")
library("tidyverse")
library("patchwork")
library("agricolae")
library("FSA")
library("ggplot2")
library("kableExtra")
library("tibble")
library("plotly")
#########################################Creating phyloseq object
taxtab <- readRDS("path to taxonomic table .rds file")
ASVtab <- readRDS("path to ASV table .rds file")
samples_df <- read_excel("path to metadata file")
head(samples_df)
samples_df <- samples_df %>% 
  tibble::column_to_rownames("Sample_ID")
samples = sample_data(samples_df)
taxa.print <- taxtab 
rownames(taxa.print) <- NULL
head(taxa.print)
pso <- phyloseq(otu_table(ASVtab, taxa_are_rows=FALSE), 
               tax_table(taxtab), sample_data(samples))
pso
dna <- Biostrings::DNAStringSet(taxa_names(pso))
names(dna) <- taxa_names(pso)
pso <- merge_phyloseq(pso, dna)
taxa_names(pso) <- paste0("ASV", seq(ntaxa(pso)))
pso
tax_table(pso)   
set.seed(1782) # set seed for analysis reproducibility

#########################################Remove single and double tones
pso_no_double <- physeq_filtered <- filter_taxa(pso, function(x) sum(x) > 2, TRUE)

#########################################Removal of potential contaminants
df <- as.data.frame(sample_data(pso_no_double))
df$LibrarySize <- sample_sums(pso_no_double)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

sample_data(pso_no_double)$is.neg <- sample_data(pso)$Sample_or_Control == "Control"

contamdf.freq2 <- isContaminant(pso_no_double, method="frequency", conc="quant_reading", neg="is.neg")
table(contamdf.freq2$contaminant)
pso.noncontam <- prune_taxa(!contamdf.freq2$contaminant, pso_no_double)

#########################################Resample on even depth
pso.rare <-rarefy_even_depth(pso.onlybac, rngseed=1 , replace=F, verbose=TRUE)

#########################################Safe new taxonomy and ASV table as csv
pso.tax <- tax_table(pso.rare)
pso.otu <- otu_table(pso.rare)

write.csv(pso.tax, file = "C:/file")
write.csv(pso.otu, file = "C:/file")

##################################Calculate alpha diversity
richness <- estimate_richness(pso.rare)
richness

##################################Calculate beta diversity for height comparison
veganotu = function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}
#transform to relative abundances
pso.rare.rel = transform_sample_counts(pso.rare, function(x){x / sum(x)})

#remove samples with dust influence
pso.rare.rel.backgr <- remove_samples (samples = "B13", pso.rare.rel)
pso.rare.rel.backgr <- remove_samples (samples = "B17", pso.rare.rel.backgr)
pso.rare.rel.backgr <- remove_samples (samples = "B18", pso.rare.rel.backgr)
pso.rare.rel.backgr <- remove_samples (samples = "A13", pso.rare.rel.backgr)
pso.rare.rel.backgr <- remove_samples (samples = "A14", pso.rare.rel.backgr)
pso.rare.rel.backgr <- remove_samples (samples = "A16", pso.rare.rel.backgr)
pso.rare.rel.backgr <- remove_samples (samples = "A17", pso.rare.rel.backgr)
pso.rare.rel.backgr <- remove_samples (samples = "A18", pso.rare.rel.backgr)

# export data from phyloseq to vegan-compatible object
pso.rare.rel.backgr.vegan <- veganotu(pso.rare.rel.backgr)
pso.rare.rel.backgr.df <- data.frame(sample_data(pso.rare.rel.backgr))

#Bray-Curtis dissimilarities
backgr_BC <- vegdist(pso.rare.rel.backgr.vegan, method = "bray") 
backgr_BC

#NMDS plot
backgr_nmds.bray <- ordinate(pso.rare.rel.backgr, method = "NMDS", distance = "bray")
backgr_nmds.bray
allGroupsColors<- c('#ADD8E6','#007FFF')

p1 <- plot_ordination(pso.rare.rel.backgr, backgr_nmds.bray, color = "Height")
p1 +  scale_color_manual(values = allGroupsColors) + geom_point(size = 1)

backgr_betadisp <- betadisper(backgr_BC, pso.rare.rel.backgr.df$Height, type = "centroid")
backgr_betadisp$distances

boxplot(backgr_betadisp,
        col = c('#ADD8E6','#007FFF'), 
        main = "", 
        xlab = "", 
        ylab = "", 
        cex.axis = 1.5, 
        ylim = c(0.3, 0.75), 
        horizontal = TRUE, 
        at = c(2,1))

anova(backgr_betadisp)
permutest(backgr_betadisp)
TukeyHSD(backgr_betadisp)
adonis2(backgr_BC ~ Height, data = pso.rare.rel.backgr.df)


##################################Calculate beta diversity for dust vs. background comparison

#remove samples
pso.rare.rel.dust <- remove_samples (samples = "A04", pso.rare.rel)
pso.rare.rel.dust <- remove_samples (samples = "A05", pso.rare.rel.dust)
pso.rare.rel.dust <- remove_samples (samples = "A19", pso.rare.rel.dust)
pso.rare.rel.dust <- remove_samples (samples = "A22", pso.rare.rel.dust)
pso.rare.rel.dust <- remove_samples (samples = "A24", pso.rare.rel.dust)
pso.rare.rel.dust <- remove_samples (samples = "A25", pso.rare.rel.dust)
pso.rare.rel.dust <- remove_samples (samples = "A26", pso.rare.rel.dust)
pso.rare.rel.dust <- remove_samples (samples = "A27", pso.rare.rel.dust)
pso.rare.rel.dust <- remove_samples (samples = "A28", pso.rare.rel.dust)
pso.rare.rel.dust <- remove_samples (samples = "A31", pso.rare.rel.dust)
pso.rare.rel.dust <- remove_samples (samples = "B04", pso.rare.rel.dust)
pso.rare.rel.dust <- remove_samples (samples = "B05", pso.rare.rel.dust)
pso.rare.rel.dust <- remove_samples (samples = "B19", pso.rare.rel.dust)
pso.rare.rel.dust <- remove_samples (samples = "B22", pso.rare.rel.dust)
pso.rare.rel.dust <- remove_samples (samples = "B23", pso.rare.rel.dust)
pso.rare.rel.dust <- remove_samples (samples = "B24", pso.rare.rel.dust)
pso.rare.rel.dust <- remove_samples (samples = "B25", pso.rare.rel.dust)
pso.rare.rel.dust <- remove_samples (samples = "B26", pso.rare.rel.dust)
pso.rare.rel.dust <- remove_samples (samples = "B27", pso.rare.rel.dust)
pso.rare.rel.dust <- remove_samples (samples = "B28", pso.rare.rel.dust)

# export data from phyloseq to vegan-compatible object
pso.rare.rel.dust.vegan <- veganotu(pso.rare.rel.dust)
pso.rare.rel.dust.df <- data.frame(sample_data(pso.rare.rel.dust))

#Bray-Curtis dissimilarities
dust_BC <- vegdist(pso.rare.rel.dust.vegan, method = "bray") 
dust_BC

#NMDS plot
dust_nmds.bray <- ordinate(pso.rare.rel.dust, method = "NMDS", distance = "bray")
dust_nmds.bray
allGroupsColors<- c('#ADD8E6', '#FCB2BF','#007FFF','#FF0000')

p2 <- plot_ordination(pso.rare.rel.dust, dust_nmds.bray, color = "Dust")
p2 +  scale_color_manual(values = allGroupsColors) + geom_point(size = 1)

dust_betadisp <- betadisper(dust_BC, pso.rare.rel.dust.df$Dust, type = "centroid")
dust_betadisp$distances

boxplot(dust_betadisp,
        col = c('#ADD8E6', '#FCB2BF','#007FFF','#FF0000'), 
        main = "", 
        xlab = "", 
        ylab = "", 
        cex.axis = 1.5, 
        ylim = c(0.3, 0.75), 
        horizontal = TRUE, 
        at = c(2,1))

anova(dust_betadisp)
permutest(dust_betadisp)
TukeyHSD(dust_betadisp)
adonis2(dust_BC ~ Height, data = pso.rare.rel.dust.df)


#########################################plotting barcharts
ps.rare.phylum <- phyloseq::tax_glom(ps.rare, "Order", NArm=FALSE)
ps.rare.phylum_top11 <- get_top_taxa(ps.rare.phylum, 11, relative = TRUE, discard_other = FALSE,
                                other_label = "Other")

taxic <- as.data.frame(ps.rare.phylum_top11@tax_table)  
colourCount = length(unique(taxic$Order))  
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
bars_phylum <- ps.rare.phylum_top11 %>%                     
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  psmelt() %>%                         
  arrange(Order)

ggplot(bars_phylum, 
       aes(x = Date, y = Abundance, fill = Order)) + 
  facet_grid(Height~.) +
  theme(legend.title = element_text(size=12),
        legend.text = element_text(size=10))+
  geom_bar(stat = "identity") + scale_fill_manual("Order", values = getPalette(colourCount)) +
  theme(axis.title.x = element_blank()) + 
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance \n")  +
  ggtitle("Order")+
  theme_classic()

#########################################statistical testing
#Shapiro-Wilk for normal distribution
shapiro.test(df)
#Levene´s test  homogeneity of variances
leveneTest(df)
#two-sided t-test
t.test(x, y = NULL, alternative = "two.sided", paired = FALSE, conf.level = 0.95)
#Wilcox rank sum test
wilcox.test(x, y = NULL, alternative = "two.sided", paired = FALSE, conf.level = 0.95)
