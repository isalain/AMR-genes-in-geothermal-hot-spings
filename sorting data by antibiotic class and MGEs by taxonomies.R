#sorting data by antibiotic class and MGEs by taxonomies
samples <- c("BM4a", "BM4b", "BM1c", "BS1a","BS2c", "GM1b", "GM4a", "GM8b", "GS1b", "GS6b", "GS3a")
#sorting out absolute abundance ARG drug class
sorted_by_class <- args_data  %>%
  group_by(Drug_class) %>%
  summarize(across(all_of(samples), sum, na.rm = TRUE))

sample_arg_class_values <- sorted_by_class%>%
  pivot_longer(cols = -Drug_class, names_to = "Sample", values_to = "Value")

# Calculate the total abundance for each sample
sample_totals_args_class <- sample_arg_class_values%>%
  group_by(Sample) %>%
  summarise(total_abundance_class = sum(Value))

# Join the total abundances with the original data
args_class_data_relative <- sample_arg_class_values %>%
  left_join(sample_totals_args_class, by = "Sample") %>%
  mutate(class_relative_abundance = (Value / total_abundance_class) * 100)

merged_args_class_data <- merge(args_class_data_relative, sample_data, by = "Sample", all.x = TRUE)


####MGEs####

#sorting MGEs data
#sorting out absolute abundance ARG drug class
samples <- c("BM4a", "BM4b", "BM1c", "BS1a","BS2c", "GM1b", "GM4a", "GM8b", "GS1b", "GS6b", "GS3a")
sorted_by_tax <- MGEs_level1%>%
  group_by(Taxonomy) %>%
  summarize(across(all_of(samples), sum, na.rm = TRUE))
sample_mge_values <- sorted_by_tax%>%
  pivot_longer(cols = -Taxonomy, names_to = "Sample", values_to = "Value")
merged_mge_data <- merge(sample_mge_values, sample_data, by = "Sample", all.x = TRUE)















