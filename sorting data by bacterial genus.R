#sorting data by bacterial phyla 
filtered_bact_data <- bacteria %>%
  filter(Kingdom == "k__Bacteria")
#defining sample names
samples <- c("BM4a", "BM4b", "BM1c", "BS1a","BS2c", "GM1b", "GM4a", "GM8b", "GS1b", "GS6b", "GS3a")
#sorting out absolute abundance of phyla in bacteria kingdom
sorted_by_Genus <- filtered_bact_data  %>%
  group_by(Genus) %>%
  summarize(across(all_of(samples), sum, na.rm = TRUE))

sorted_by_Genus_trans<- t(sorted_by_Genus)
head(sorted_by_Genus_trans)
sorted_by_GenusTot <- sorted_by_Genus%>%
  mutate(Overall_Sum = rowSums(.[,-1]))

top_10_genus <- sorted_by_GenusTot %>%
  arrange(desc(Overall_Sum)) %>%
  head(11)

row_label <- "2"
column_label <- "Overall_Sum"

top_10_genus <- top_10_genus[!(rownames(top_10_genus) %in% row_label), !(colnames(top_10_genus) %in% column_label)]



















