# R scripts used for data analyses of ARGs from a tropical metagenome
# importing important r packages
#for data loading and sorting
library(readxl); library(dplyr); library(tidyverse);library(reshape2)
#for data visualization
library(ggplot2); library(ggpubr); library(ggalluvial); library(ggVennDiagram);library(patchwork);
library(gplots)
#for diversity analyses
library(vegan)
#for BART analysis
library(vegan)
#############################################################################################################
#loading datasets
setwd("C:/Users/aisabwe/Desktop/Ben")
args_data<- read_excel("Data.xlsx", sheet = "ARGs", na = "")
site_data<- read_excel("Data.xlsx", sheet = "metadata2", na = "")
sample_data<- read_excel("Data.xlsx", sheet = "metadata", na = "")
metagenome_data_quality<- read_excel("Data.xlsx", sheet = "data_quality", na = "")
MGEs_level1<- read_excel("Data.xlsx", sheet = "MGEs_1", na = "")
MGEs_level2<- read_excel("Data.xlsx", sheet = "MGEs_2", na = "")
MGEs_level3<- read_excel("Data.xlsx", sheet = "MGEs_3", na = "")
bacteria<- read_excel("Data.xlsx", sheet = "OTUs", na = "")
#############################################################################################################
# Create the box plots of other water characteristics and comparisons
source("water_chemisty_analysis.R")#loading a script that does comparisons

Fig1<-alt_temper+temp_boxplot
  
plot2<-ggarrange(pH_boxplot, TDS_boxplot,EC_boxplot, ORP_boxplot,
                 TP_boxplot,TN_boxplot,NO3N_boxplot,NH3N_boxplot,
                 ncol = 2, nrow = 4)
#############################################################################################################
# Comparison of ARGs and bacteria richness between low and high water temperature
# Richness estimation
bact_richness <- specnumber(t(bacteria[, -c(1, 13,14,15,16,17,18,19)]))
arg_richness <- specnumber(t(args_data[, -c(1,2,3,4, 16)]))
mge_richness <- specnumber(t(MGEs_level3[, -c(1,13)]))
combined_bact_arg_richness <- cbind(bact_richness, arg_richness,mge_richness)
rich_by_sample<-cbind(combined_bact_arg_richness,sample_data)
mat_richness<-rich_by_sample[rich_by_sample$ Substrate == "Mat", ]
sed_richness<-rich_by_sample[rich_by_sample$ Substrate == "Sediment", ]

#Box plots of bacterial and ARGs richness
bac_boxplot <- ggbarplot(rich_by_sample, x = "Location", y = "bact_richness",
                         add = c("mean_se", "point", size =4),
                         color = "Substrate", fill = "Substrate", alpha = 0.5,
                         palette = c("#009E73", "#E7B800"))

Plot3 <- bac_boxplot + labs(y = expression("Bacterial richness"))
#############################################################################################################
ARGboxplot <- ggviolin(rich_by_sample, "Location", "arg_richness",
                       fill = "Location", palette = c("#009E73", "#E7B800")) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 3) +
  xlab("") +theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "right") +
  ylim(0, 2000) +facet_grid(. ~ Substrate) +labs(fill = "Location")+theme_bw()
Plot4<-ARGboxplot+ labs(y = expression("ARGs richness"))
#############################################################################################################
#Sorting data for taxonomic composition
source("sorting_data_by_kingdom.R")#loading script that sort data by kingdom
AllKingPlot<-ggplot(merged_kingdom_data, aes(x = Sample, y=relative_abundance, 
                                             stratum = Kingdom, alluvium = Kingdom, fill = Kingdom)) +
  geom_flow(stat = "alluvium") +geom_stratum() +theme_minimal() 
Plot5<-AllKingPlot
#############################################################################################################
source("sorting data by bacterial phyla.R")#loading script that sort bacterial data by phylum
BactPhylPlot<-ggplot(merged_bact_phylum_data, aes(x = Sample, y=phylum_relative_abundance, 
                                                  stratum = Phylum, alluvium = Phylum, fill = Phylum)) +
  geom_flow(stat = "alluvium") +geom_stratum() +theme_minimal() 
Plot6<-BactPhylPlot
#############################################################################################################
#Phylum-level taxonomic composition of bacterial community
BactPhylPlot2<-ggplot(merged_bact_phylum_data2, aes(x = Sample, y=phylum_relative_abundance, 
                                                    stratum = Phylum_Grouped, alluvium = Phylum_Grouped, 
                                                    fill = Phylum_Grouped)) +geom_flow(stat = "alluvium") +
  geom_stratum() +theme_minimal()+theme(legend.position = "right")
Plot7<-BactPhylPlot2+ 
  labs(y = expression("Relative abundance"))
#############################################################################################################
#Bacterial community structure in PCOA
bact.hel<-decostand(t(bacteria[, -c(1, 13,14,15,16,17,18,19)]), "hellinger")
bactBC<-vegdist(bact.hel)
pcoaBact = cmdscale(bactBC, k=3, eig=T) 
points = as.data.frame(pcoaBact$points) 
PCOaPoints <- data.frame(Sample = row.names(points), points)
colnames(points) = c("x", "y", "z") 
eigBact = pcoaBact$eig
PCOaPoints_sample_data <- merge(PCOaPoints, sample_data, by.x = 0, by.y = "Sample", all = TRUE)
# Make the figure with ggplot2
BactPCoA = ggplot(PCOaPoints_sample_data, aes(x=V1, y=V2,shape=Location, color=Substrate)) +
  geom_point(size=5) + 
  labs(x=paste("PCoA 1 (", format(100 * eigBact[1] / sum(eigBact), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigBact[2] / sum(eigBact), digits=4), "%)", sep=""))+
  theme_bw()+scale_shape_manual(values=c(15,16))+
  scale_color_manual(values = c("#009E73", "#E7B800"))+
  theme(legend.position = "right")
BactANOSIM_ByLoc<- anosim(bactBC, sample_data$Location)
BactANOSIM_BySub<- anosim(bactBC, sample_data$Substrate)
AddANOSIM<-BactPCoA+ annotate("text", x = 0.2, y = -0.3, label = "p-value (Loc): 0.363", size = 4)
AddANOSIM2<-AddANOSIM+ annotate("text", x = 0.2, y = -0.35, label = "p-value (Sub): 0.002", size = 4)
Plot8<-AddANOSIM2
#############################################################################################################
#ARGs community structure in PCOA
ARG.hel<-decostand(t(args_data[, -c(1,2,3,4, 16)]), "hellinger")
ARGBC<-vegdist(ARG.hel)
pcoaARG = cmdscale(ARGBC, k=3, eig=T) 
points2 = as.data.frame(pcoaARG$points) 
PCOaPoints2 <- data.frame(Sample = row.names(points2), points2)
colnames(points2) = c("x", "y", "z") 
eigARGs = pcoaARG$eig
PCOaPoints_sample_data2 <- merge(PCOaPoints2, sample_data, by.x = 0, by.y = "Sample", all = TRUE)
# Make the figure with ggplot2
ARGPCoA = ggplot(PCOaPoints_sample_data2, aes(x=V1, y=V2,shape=Location, color=Substrate)) +
  geom_point(size=5) + 
  labs(x=paste("PCoA 1 (", format(100 * eigARGs[1] / sum(eigARGs), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigARGs[2] / sum(eigARGs), digits=4), "%)", sep=""))+
  theme_bw()+scale_shape_manual(values=c(15,16))+
  scale_color_manual(values = c("#009E73", "#E7B800"))+
  theme(legend.position = "right")
ARGANOSIM_ByLoc<- anosim(ARGBC, sample_data$Location)#the result was written manually on the figure
ARGANOSIM_BySub<- anosim(ARGBC, sample_data$Substrate)#the result was written manually on the figure
AddANOSIM<-ARGPCoA+ annotate("text", x = 0.2, y = 0.10, label = "p-value (Loc): 0.309 ", size = 4)
AddANOSIM2<-AddANOSIM+ annotate("text", x = 0.2, y = 0.08, label = "p-value (Sub): 0.023", size = 4)
Plot9<-AddANOSIM2
#############################################################################################################
source("sorting data by antibiotic class and MGEs by taxonomies.R")#loading script that sort bacterial data by phylum
#Drug Class ARGs
color_codes_args <- c("Aminoglycoside" = "#D02028",  
                      "Diaminopyrimidine" = "#F5A720",  
                      "Fluoroquinolone" = "#C49A6C",
                      "Glycopeptide" = "#7553A2",
                      "MDR" = "#85C441",
                      "MLSB" = "#427638",
                      "Peptide" = "#984B9D",
                      "Rifamycin" = "#8B582B",
                      "Sulfonamide" = "#558DCA",
                      "Tetracycline" = "#6DC7B5",
                      "etc" = "#CCCBCA",
                      "Beta−lactam" = "#F3766E",
                      "Amphenicol" = "#603913") 
argClassPlot<-ggplot(merged_args_class_data, aes(x = Sample, y=class_relative_abundance, 
                                                 stratum = Drug_class, alluvium = Drug_class, 
                                                 fill = Drug_class)) +geom_flow(stat = "alluvium") +
  geom_stratum() +theme_minimal()+theme(legend.position = "right")
Plot10<-argClassPlot+ 
  labs(y = expression("Relative abundance of drug class"))+
  scale_fill_manual(values = color_codes_args)
#############################################################################################################
#MGEs by taxonomy
mgeTaxPlot<-ggplot(merged_mge_data, aes(x = Sample, y=Value, 
                                        stratum = Taxonomy, alluvium = Taxonomy, 
                                        fill = Taxonomy)) +geom_flow(stat = "alluvium") +
  geom_stratum() +theme_minimal()+theme(legend.position = "right")
Plot11<-mgeTaxPlot+ 
  labs(y = expression("Relative abundance of MGEs"))

#MGEs by taxonomy
source("comparing_MGEs.R")
mge_boxplot <- ggarrange(transp_boxplot,plasmid_boxplot, insertion_boxplot,
                         integrase_boxplot, ncol = 1, nrow = 4)
Plot11 <- mge_boxplot

#############################################################################################################
#Shared and unique ARGs across the two habitat types
bugarama_group <- subset(args_data, select = c(Short_name,Drug_class, Resistance_mechanism, 
                                               BM4a, BM4b, BM1c, BS1a,BS2c))
bugarama_group1<- bugarama_group %>%
  mutate(total = rowSums(select(., BM4a, BM4b, BM1c, BS1a,BS2c)))
all_bug<-bugarama_group1[, c("Short_name", "Drug_class","Resistance_mechanism", "total")]

gisenyi_group <- subset(args_data, select = c(Short_name,Drug_class, Resistance_mechanism, 
                                              GM4a, GM8b, GM1b, GS1b,GS6b,GS3a))
gisenyi_group1<- gisenyi_group %>%
  mutate(total = rowSums(select(.,GM4a, GM8b, GM1b, GS1b,GS6b,GS3a)))
all_gis<-gisenyi_group1[, c("Short_name", "Drug_class","Resistance_mechanism", "total")]
############

combined_args<-merge(all_bug, all_gis, by= "Short_name")
unique_bug <- combined_args[combined_args$total.y == 0,]
unique_gis <- combined_args[combined_args$total.x == 0,]

#search all args present in bugarama and gisenyi separately
gis_args<-combined_args[combined_args$total.y> 0,]
bug_args<-combined_args[combined_args$total.x> 0,]

shared_args<-combined_args[combined_args$total.x> 0 & combined_args$total.y>0,]

total_args<-nrow(unique_bug)+nrow(unique_gis)+nrow(shared_args)
#creating a venn diagram
data_for_venn <- list(Bugarama = row.names(bug_args), 
                      Gisenyi = row.names(gis_args))
venn_args<-ggVennDiagram(data_for_venn , color = "black", lwd = 1, lty = 1)+
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
Plot12<-venn_args
#############################################################################################################
#taxonomic composition of unique ARGs 
#Subset unique ARGs from original ARGs dataset
unique_bug_args<- all_bug[all_bug$Short_name %in% unique_bug$Short_name, ]
unique_gis_args<- all_gis[all_gis$Short_name %in% unique_gis$Short_name, ]

unique_bug_args$Location <- "Bugarama"
unique_gis_args$Location <- "Gisenyi"

combined_args_tax <- dplyr::bind_rows(unique_bug_args, unique_gis_args)
drug_class_by_location <- combined_args_tax %>%
  group_by(Location, Drug_class) %>%
  summarise(total= sum(total),.groups = "drop")
total_by_location <- aggregate(total ~ Location+  Drug_class, drug_class_by_location, sum)
SharedPlot<-ggplot(total_by_location, aes(x = Location, y=total, 
                                          stratum = Drug_class, alluvium = Drug_class, 
                                          fill = Drug_class)) +geom_flow(stat = "alluvium") +
  geom_stratum() +theme_minimal()+theme(legend.position = "right")
Plot13<-SharedPlot+ 
  labs(y = expression("Absolute abundance"))+
  scale_fill_manual(values = color_codes_args)+coord_flip()
#############################################################################################################
#ARGs number per mode
eff_pump<-nrow(filter(args_data, Resistance_mechanism == "efflux pump"))
ant_targ_alter<-nrow(filter(args_data, Resistance_mechanism == "antibiotic target alteration"))
ant_inact<-nrow(filter(args_data, Resistance_mechanism == "antibiotic inactivation"))
ant_targ_repl<-nrow(filter(args_data, Resistance_mechanism == "antibiotic target replacement"))
ant_targ_prot<-nrow(filter(args_data, Resistance_mechanism == "antibiotic target protection"))
tot_mode<-nrow(args_data)
percentages_mode <- data.frame(
  Category = c("efflux pump", "antibiotic target alteration", "antibiotic inactivation", 
               "antibiotic target replacement","antibiotic target protection"),
  Percentage = c(eff_pump, ant_targ_alter, ant_inact, ant_targ_repl, ant_targ_prot) / tot_mode * 100)

labs3 <- paste0(percentages_mode$Category, " (", round(percentages_mode$Percentage, 1), "%)")

pieModeAll<-ggdonutchart(percentages_mode, "Percentage", label = labs3,
                         lab.pos = "in", lab.font = "white",
                         fill = "Category", color = "white")

Plot14<-pieModeAll
#############################################################################################################
#Search for efflux pump elements
keyword1 <- "ABC"; keyword2 <- "RND"; keyword3 <- "MFS"; keyword4 <- "SMR";keyword5 <- "MATE" 
arg_eff_pump<-filter(args_data, Resistance_mechanism == "efflux pump")
ABC_data <- filter(arg_eff_pump, grepl(keyword1, Description))
RND_data <- filter(arg_eff_pump, grepl(keyword2, Description))
MFS_data <- filter(arg_eff_pump, grepl(keyword3, Description))
SMR_data <- filter(arg_eff_pump, grepl(keyword4, Description))
MATE_data <- filter(arg_eff_pump, grepl(keyword5, Description))
#Percentages
tot_ef_pump<-nrow(arg_eff_pump)
ABC<-nrow(ABC_data); RND<-nrow(RND_data); MFS<-nrow(MFS_data); 
SMR<-nrow(SMR_data); MATE<-nrow(MATE_data)
Others<-tot_ef_pump-(ABC+RND+MFS+SMR+MATE)

percentages_eff_pumps <- data.frame(
  Category = c("ABC", "MATE", "MFS", "Others", "RND","SMR"),
  Percentage = c(ABC, MATE, MFS, Others, RND, SMR) / tot_ef_pump * 100)

# Create the donut pie chart
labs <- paste0(percentages_eff_pumps$Category, " (", round(percentages_eff_pumps$Percentage, 1), "%)")

pieEffPlot<-ggdonutchart(percentages_eff_pumps, "Percentage", label = labs,
                         lab.pos = "in", lab.font = "white",
                         fill = "Category", color = "white",
                         palette = c("#00AFBB", "#E7B800", "#FC4E07","#4981BF","#009E73", "#D5FF7F") )

Plot15<-pieEffPlot
#############################################################################################################
#Predictors of ARGs patterns
# Search potential human pathogens
bacteria_genus<-filter(bacteria, Genus == "g__")
kw1 <- "Escherichia"; kw2 <- "Salmonella"; kw3 <- "Staphylococcus"; kw4 <- "Streptococcus";
kw5 <- "Neisseria";kw6 <- "Clostridium";kw7 <- "Helicobacter";kw8 <- "Mycobacterium";
kw9 <- "Bordetella"; kw10 <- "Vibrio"; kw11 <- "Listeria"; kw12 <- "Campylobacter";
kw13 <- "Shigella"; kw14 <- "Yersinia"; kw15 <- "Legionella"; kw16 <- "Pseudomonas";
kw17 <- "Enterococcus"

# filtering potential pathogens
Escherichia <- filter(bacteria, grepl(kw1, Genus));Salmonella <- filter(bacteria, grepl(kw2, Genus));
Staphylococcus <- filter(bacteria, grepl(kw3, Genus)); Streptococcus <- filter(bacteria, grepl(kw4, Genus));
Neisseria <- filter(bacteria, grepl(kw5, Genus)); Clostridium <- filter(bacteria, grepl(kw6, Genus)); 
Helicobacter <- filter(bacteria, grepl(kw7, Genus)); Mycobacterium <- filter(bacteria, grepl(kw8, Genus));
Bordetella <- filter(bacteria, grepl(kw9, Genus)); Vibrio <- filter(bacteria, grepl(kw10, Genus));
Listeria <- filter(bacteria, grepl(kw11, Genus)); Campylobacter <- filter(bacteria, grepl(kw12, Genus));
Shigella <- filter(bacteria, grepl(kw13, Genus)); Yersinia <- filter(bacteria, grepl(kw14, Genus));
Legionella <- filter(bacteria, grepl(kw15, Genus)); Pseudomonas <- filter(bacteria, grepl(kw16, Genus));
Enterococcus <- filter(bacteria, grepl(kw17, Genus))
###########################################################################################################
pathogens<-rbind(Escherichia,Salmonella,Staphylococcus,Streptococcus,Neisseria,Clostridium,
                 Helicobacter, Mycobacterium, Bordetella, Vibrio, Listeria, Campylobacter,
                 Shigella, Yersinia, Legionella, Pseudomonas, Enterococcus)
all_OTUs<-nrow(bacteria)
pathogens_OTUs<-nrow(pathogens)
non_pathogens_OTUs<-all_OTUs-pathogens_OTUs

percentages_path <- data.frame(
  Cat= c("Pathogens", "Non Pathogens"),
  Perc = c(pathogens_OTUs,non_pathogens_OTUs) / all_OTUs * 100)

# Create the donut pie chart
lab2 <- paste0(percentages_path$Cat, " (", round(percentages_path$Perc, 1), "%)")

piePathog<-ggdonutchart(percentages_path, "Perc", label = lab2,
                        lab.pos = "in", lab.font = "white",
                        fill = "Cat", color = "white",
                        palette = c("#4981BF","#FC4E07"))

Plot16<-piePathog
################################################################################################################
# Create 
path_sum<-colSums(pathogens[,-c(1, 13,14,15,16,17,18,19)])
bact_sum<-colSums(bacteria[,-c(1, 13,14,15,16,17,18,19)])
path_by_bact<-cbind(path_sum,bact_sum,sample_data)
tot_path_reads <- sum(path_by_bact$path_sum)
tot_bact_reads <- sum(path_by_bact$bact_sum)
percentages_path_reads <- data.frame(
  Cat= c("Pathogens", "Non Pathogens"),
  Perc = c(tot_path_reads,tot_bact_reads-tot_path_reads) / tot_bact_reads * 100)

# Create the donut pie chart
lab3 <- paste0(percentages_path_reads$Cat, " (", round(percentages_path_reads$Perc, 1), "%)")

piePathogReads<-ggdonutchart(percentages_path_reads, "Perc", label = lab3,
                        lab.pos = "in", lab.font = "white",
                        fill = "Cat", color = "white",
                        palette = c("#4981BF","#FC4E07"))

Plot17<-piePathogReads
################################################################################################################

#Create a plots of pathogens per sample type
path_richness <- specnumber(t(pathogens[, -c(1, 13,14,15,16,17,18,19)]))
path_by_sample<-cbind(path_richness,sample_data)

#Box plots of pathogen richness
path_boxplot <- ggviolin(path_by_sample, "Location", "path_richness",
                        fill = "Location", palette = c("#009E73", "#E7B800")) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 3) +
  facet_grid(. ~ Substrate) +
  theme_bw()
Plot17 <- path_boxplot + labs(y = expression("pathogen richness"))
################################################################################################################
#Creating a heatmap at genus level of pathogens
aggregated_data <- aggregate(. ~ Genus, data = pathogens[,-c(1, 13,14,15,16,17,19)], sum)
aggregated_data2 <- aggregated_data[-1]  
row.names(aggregated_data2) <- aggregated_data$Genus
htmp<-heatmap.2(as.matrix(aggregated_data2), scale = "column", col = bluered(100), 
          dendrogram = c("column"),
          trace = "none", density.info = "none")
Plot18<-htmp
################################################################################################################
#Bubble plot of the 20 most dominant pathogens
#search for 20 dominant pathogens
path_row_sum <- pathogens %>% mutate(Sum = BM4a+BM4b+BM1c+BS1a+BS2c+GM1b+GM4a+GM8b+GS1b+GS6b+GS3a)
path_row_sum_order <- path_row_sum[order(-path_row_sum$Sum), ]
top30_dominant <- path_row_sum_order [1:20,]
melted_top20_dominant <- pivot_longer(top30_dominant,cols = c(BM4a,BM4b,BM1c,BS1a,BS2c,GM1b,GM4a,GM8b,GS1b,GS6b,GS3a),
                                      names_to = "Variable", values_to = "Value")

bbPlot<-ggplot(melted_top20_dominant, aes(x = Variable, y = Species, size = log1p(Value), color=Order)) +
  geom_point()  +
  scale_size(range = c(2, 10)) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

Plot19<-bbPlot
################################################################################################################
#barplot with averages and sd
head(path_by_bact)
PathReads<-ggbarplot(path_by_bact, x = "Location", y = "path_sum",
          add = c("mean_se", "point", size = 2),
          color = "Location", fill = "Location", alpha = 0.5,
          palette = c("#009E73", "#E7B800"))
  
Plot20<-PathReads

################################################################################################################
# Bayesian hierarchical model (BHM)
# Load required libraries
library(brms)
library(ggplot2)
library(tidybayes)
library(dplyr)
library(tidyr)

library(dplyr)
library(tibble)

# First, let's convert bact_args_mge to a tibble and add the row names as a column
bact_args_mge <- bact_args_mge %>%
  rownames_to_column(var = "Sample")

# Now, let's combine the datasets
combined_data <- bact_args_mge %>%
  left_join(path_by_bact, by = "Sample") %>%
  mutate(path_ratio = path_sum / bact_sum)

# View the first few rows of the combined dataset
head( combined_for_bhm)

model_new <- brm(
  formula = arg_richness ~ bact_richness + 
    (1 | Substrate) + (0 + bact_richness | Substrate) + 
    (1 | Location) + (0 + bact_richness | Location),
  data = combined_for_bhm,
  family = gaussian(),
  prior = c(
    prior(normal(0, 10), class = "Intercept"),
    prior(normal(0, 5), class = "b"),
    prior(normal(0, 5), class = "sd", lb = 0),
    prior(normal(0, 5), class = "sigma")
      ),
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 25
  ),
  chains = 4,
  iter = 8000,
  warmup = 4000
)

library(tidyverse)
library(ggplot2)
library(brms)

# Extract posterior samples
post_samples <- as_draws(model_new)

# Combine all chains into a single data frame
post_samples_df <- post_samples %>%
  map_df(~ as.data.frame(.), .id = "chain")

# Extract population-level effects
pop_effects <- post_samples_df %>%
  select(b_Intercept, b_bact_richness) %>%
  pivot_longer(everything(), names_to = "Parameter", values_to = "Effect") %>%
  group_by(Parameter) %>%
  summarise(
    Effect_mean = mean(Effect),
    Effect_lower = quantile(Effect, 0.025),
    Effect_upper = quantile(Effect, 0.975)
  ) %>%
  mutate(Parameter = case_when(
    Parameter == "b_Intercept" ~ "Intercept",
    Parameter == "b_bact_richness" ~ "Bacterial richness",
    TRUE ~ Parameter
  ))

pop_effects
# Extract group-level effects (standard deviations)
group_effects <- post_samples_df %>%
  select(starts_with("sd_")) %>%
  pivot_longer(everything(), names_to = "Group", values_to = "Effect") %>%
  group_by(Group) %>%
  summarise(
    Effect_mean = mean(Effect),
    Effect_lower = quantile(Effect, 0.025),
    Effect_upper = quantile(Effect, 0.975)
  ) %>%
  mutate(Group = case_when(
    str_detect(Group, "Substrate") ~ "Substrate",
    str_detect(Group, "Location") ~ "Location",
    TRUE ~ Group
  ))

group_effects


#####

# View the first few rows of the combined dataset
head( combined_for_bhm)

model_2<- brm(
  formula = arg_richness ~ mge_richness + 
    (1 | Substrate) + (0 + mge_richness | Substrate) + 
    (1 | Location) + (0 + mge_richness | Location),
  data = combined_for_bhm,
  family = gaussian(),
  prior = c(
    prior(normal(0, 10), class = "Intercept"),
    prior(normal(0, 5), class = "b"),
    prior(normal(0, 5), class = "sd", lb = 0),
    prior(normal(0, 5), class = "sigma")
  ),
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 25
  ),
  chains = 4,
  iter = 8000,
  warmup = 4000
)

library(tidyverse)
library(ggplot2)
library(brms)

# Extract posterior samples
post_samples2 <- as_draws(model_2)

# Combine all chains into a single data frame
post_samples_df2 <- post_samples2 %>%
  map_df(~ as.data.frame(.), .id = "chain")

# Extract population-level effects
pop_effects2 <- post_samples_df2 %>%
  select(b_Intercept, b_mge_richness) %>%
  pivot_longer(everything(), names_to = "Parameter", values_to = "Effect") %>%
  group_by(Parameter) %>%
  summarise(
    Effect_mean = mean(Effect),
    Effect_lower = quantile(Effect, 0.025),
    Effect_upper = quantile(Effect, 0.975)
  ) %>%
  mutate(Parameter = case_when(
    Parameter == "b_Intercept" ~ "Intercept",
    Parameter == "b_mge_richness" ~ "MGE richness",
    TRUE ~ Parameter
  ))

pop_effects2
# Extract group-level effects (standard deviations)
group_effects2 <- post_samples_df2 %>%
  select(starts_with("sd_")) %>%
  pivot_longer(everything(), names_to = "Group", values_to = "Effect") %>%
  group_by(Group) %>%
  summarise(
    Effect_mean = mean(Effect),
    Effect_lower = quantile(Effect, 0.025),
    Effect_upper = quantile(Effect, 0.975)
  ) %>%
  mutate(Group = case_when(
    str_detect(Group, "Substrate") ~ "Substrate",
    str_detect(Group, "Location") ~ "Location",
    TRUE ~ Group
  ))

group_effects2
conditional_effects(model_2)

library(brms)
library(ggplot2)
library(dplyr)
# Generate predictions
predictions <- predict(model_new)
# Combine predictions with observed data
plot_data1 <-  combined_for_bhm %>%
  mutate(predicted = predictions[, "Estimate"],
         lower = predictions[, "Q2.5"],
         upper = predictions[, "Q97.5"])
# Create the plot
ggplot(plot_data1, aes(x = predicted, y = arg_richness)) +
  geom_point(alpha = 0.6) +
  geom_errorbar(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "Predicted AMR gene richness",
       y = "Observed AMR gene richness",
       title = "Pre") +
  theme_minimal() +
  coord_equal()+ylim(500,1800)+ xlim(500,1800)

# Add correlation coefficient
cor_coef <- cor(plot_data1$predicted, plot_data$arg_richness)
Predplot1<-ggplot(plot_data1, aes(x = predicted, y = arg_richness)) +
  geom_point(alpha = 0.6) +
  geom_errorbar(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "Predicted AMR gene richness",
       y = "Observed AMR gene richness",
       title = "PRED",
       subtitle = paste("Correlation coefficient:", round(cor_coef, 3))) +
  theme_minimal() +
  coord_equal()

predictions <- predict(model_2)
# Combine predictions with observed data
plot_data <-  combined_for_bhm %>%
  mutate(predicted = predictions[, "Estimate"],
         lower = predictions[, "Q2.5"],
         upper = predictions[, "Q97.5"])
# Create the plot
ggplot(plot_data, aes(x = predicted, y = arg_richness)) +
  geom_point(alpha = 0.6) +
  geom_errorbar(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "Predicted AMR gene richness",
       y = "Observed AMR gene richness",
       title = "args - pred") +
  theme_minimal() +
  coord_equal()+ylim(500,1800)+ xlim(500,1800)

# Add correlation coefficient
cor_coef <- cor(plot_data$predicted, plot_data$arg_richness)
Predplot2<-ggplot(plot_data, aes(x = predicted, y = arg_richness)) +
  geom_point(alpha = 0.6) +
  geom_errorbar(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "Predicted AMR gene richness",
       y = "Observed AMR gene richness",
       title = "arg pred",
       subtitle = paste("Correlation coefficient:", round(cor_coef, 3))) +
  theme_minimal() +
  coord_equal()


######
# Bact rich model
ce <- conditional_effects(model_new, re_formula = NULL)
library(ggplot2)
library(dplyr)

# Extract the data from the conditional effects object
ce_data <- ce[[1]]  # Assuming you're interested in the first (and probably only) effect

# Create the plot
Modl1Plot<-ggplot() +
  # Add the conditional effects
  geom_ribbon(data = ce_data, aes(x = bact_richness, ymin = lower__, ymax = upper__, fill = "BHM CI"), 
              alpha = 0.3) +
  geom_line(data = ce_data, aes(x = bact_richness, y = estimate__, color = "BHM"), 
            size = 1) +
  
  # Add the original data points
  geom_point(data =  combined_for_bhm, aes(x = bact_richness, y = arg_richness), 
             alpha = 0.5) +
  
  # Add OLS regression line
  geom_smooth(data =  combined_for_bhm, aes(x = bact_richness, y = arg_richness, color = "OLS"), 
              method = "lm", se = FALSE) +
  
  # Customize the plot
  labs(x = "Bacterial richness", y = "AMR gene richness", 
       title = "Cond and OLS",
       color = "Model", fill = "Interval") +
  scale_color_manual(values = c("BHM" = "blue", "OLS" = "red"),
                     labels = c("BHM" = "Bayesian Hierarchical Model", "OLS" = "Ordinary Least Squares")) +
  scale_fill_manual(values = c("BHM CI" = "blue"),
                    labels = c("BHM CI" = "BHM 95% Credible Interval")) +
  theme_minimal() +
  theme(legend.position = "bottom")+ylim(0,1700)


#####
# MGE rich model
ce <- conditional_effects(model_2, re_formula = NULL)
library(ggplot2)
library(dplyr)

# Extract the data from the conditional effects object
ce_data <- ce[[1]]  # Assuming you're interested in the first (and probably only) effect

# Create the plot
Modl2Plot<-ggplot() +
  # Add the conditional effects
  geom_ribbon(data = ce_data, aes(x = mge_richness, ymin = lower__, ymax = upper__, fill = "BHM CI"), 
              alpha = 0.3) +
  geom_line(data = ce_data, aes(x = mge_richness, y = estimate__, color = "BHM"), 
            size = 1) +
  
  # Add the original data points
  geom_point(data =  combined_for_bhm, aes(x = mge_richness, y = arg_richness), 
             alpha = 0.5) +
  
  # Add OLS regression line
  geom_smooth(data =  combined_for_bhm, aes(x = mge_richness, y = arg_richness, color = "OLS"), 
              method = "lm", se = FALSE) +
  
  # Customize the plot
  labs(x = "MGE richness", y = "AMR gene richness", 
       title = "cond and OLS",
       color = "Model", fill = "Interval") +
  scale_color_manual(values = c("BHM" = "blue", "OLS" = "red"),
                     labels = c("BHM" = "Bayesian Hierarchical Model", "OLS" = "Ordinary Least Squares")) +
  scale_fill_manual(values = c("BHM CI" = "blue"),
                    labels = c("BHM CI" = "BHM 95% Credible Interval")) +
  theme_minimal() +
  theme(legend.position = "bottom")+ylim(0,1700)


PredPlots<-ggarrange(Modl1Plot,Modl2Plot, 
                     ncol = 2, nrow = 1,common.legend = TRUE,
                     legend = "right")


random_effects <- ranef(model_new)
substrate_sd <- sd(random_effects$Substrate[,,"bact_richness"])
location_sd <- sd(random_effects$Location[,,"bact_richness"])

print(paste("Substrate SD:", substrate_sd))
print(paste("Location SD:", location_sd))
total_variance <- substrate_sd^2 + location_sd^2
substrate_proportion <- substrate_sd^2 / total_variance
location_proportion <- location_sd^2 / total_variance

print(paste("Proportion of variance explained by Substrate:", substrate_proportion))
print(paste("Proportion of variance explained by Location:", location_proportion))
library(ggplot2)

ggplot(data.frame(effect = c(ranef(model_new)$Substrate[,,"bact_richness"], 
                             ranef(model_new)$Location[,,"bact_richness"]),
                  group = c(rep("Substrate", nrow(ranef(model_new)$Substrate)),
                            rep("Location", nrow(ranef(model_new)$Location)))),
       aes(x = group, y = effect)) +
  geom_boxplot() +
  labs(title = "Random Effects by Group",
       y = "Effect on bacterial richness slope")+theme_minimal()

#####
random_effects <- ranef(model_2)
substrate_sd <- sd(random_effects$Substrate[,,"mge_richness"])
location_sd <- sd(random_effects$Location[,,"mge_richness"])

print(paste("Substrate SD:", substrate_sd))
print(paste("Location SD:", location_sd))
total_variance <- substrate_sd^2 + location_sd^2
substrate_proportion <- substrate_sd^2 / total_variance
location_proportion <- location_sd^2 / total_variance

print(paste("Proportion of variance explained by Substrate:", substrate_proportion))
print(paste("Proportion of variance explained by Location:", location_proportion))

ggplot(data.frame(effect = c(ranef(model_2)$Substrate[,,"mge_richness"], 
                             ranef(model_2)$Location[,,"mge_richness"]),
                  group = c(rep("Substrate", nrow(ranef(model_2)$Substrate)),
                            rep("Location", nrow(ranef(model_2)$Location)))),
       aes(x = group, y = effect)) +
  geom_boxplot() +
  labs(title = "Random Effects by Group",
       y = "Effect on bacterial richness slope")+theme_minimal()


