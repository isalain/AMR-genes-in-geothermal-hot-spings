#sorting data by Kingdom 
sort_by_kingdom <- bacteria %>%
  group_by(Kingdom) 
samples <- c("BM4a", "BM4b", "BM1c", "BS1a","BS2c", "GM1b", "GM4a", "GM8b", "GS1b", "GS6b", "GS3a") 
sorted_by_kingdom <- bacteria %>%
  group_by(Kingdom) %>%
  summarize(across(all_of(samples), sum, na.rm = TRUE))

sample_kingdom_values <- sorted_by_kingdom%>%
  pivot_longer(cols = -Kingdom, names_to = "Sample", values_to = "Value")

# Calculate the total abundance for each sample
sample_totals <- sample_kingdom_values%>%
  group_by(Sample) %>%
  summarise(total_abundance = sum(Value))

# Join the total abundances with the original data
kingdom_data_relative <- sample_kingdom_values %>%
  left_join(sample_totals, by = "Sample") %>%
  mutate(relative_abundance = (Value / total_abundance) * 100)
merged_kingdom_data <- merge(kingdom_data_relative, sample_data, by = "Sample", all.x = TRUE)

