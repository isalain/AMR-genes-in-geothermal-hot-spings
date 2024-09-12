#sorting data by bacterial phyla 
filtered_bact_data <- bacteria %>%
  filter(Kingdom == "k__Bacteria")
#defining sample names
samples <- c("BM4a", "BM4b", "BM1c", "BS1a","BS2c", "GM1b", "GM4a", "GM8b", "GS1b", "GS6b", "GS3a")
#sorting out absolute abundance of phyla in bacteria kingdom
sorted_by_phylum <- filtered_bact_data  %>%
  group_by(Phylum) %>%
  summarize(across(all_of(samples), sum, na.rm = TRUE))

sample_bact_phyl_values <- sorted_by_phylum%>%
  pivot_longer(cols = -Phylum, names_to = "Sample", values_to = "Value")

# Calculate the total abundance for each sample
sample_totals_bact_phyl <- sample_bact_phyl_values%>%
  group_by(Sample) %>%
  summarise(total_abundance_phylum = sum(Value))

# Join the total abundances with the original data
phylum_data_relative <- sample_bact_phyl_values %>%
  left_join(sample_totals_bact_phyl, by = "Sample") %>%
  mutate(phylum_relative_abundance = (Value / total_abundance_phylum) * 100)

merged_bact_phylum_data <- merge(phylum_data_relative, sample_data, by = "Sample", all.x = TRUE)

phylum_data_relative2 <- phylum_data_relative %>%
  mutate(Phylum_Grouped = ifelse(phylum_relative_abundance < 0.99, "Others", Phylum))%>%
  group_by(Sample, Phylum_Grouped) %>%
  summarise(Value = sum(Value),
            total_abundance_phylum = sum(total_abundance_phylum),
            phylum_relative_abundance = sum(phylum_relative_abundance),
            .groups = "drop")

merged_bact_phylum_data2 <- merge(phylum_data_relative2, sample_data, by = "Sample", all.x = TRUE)


































