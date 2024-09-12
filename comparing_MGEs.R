# Comparing MGEs
# Box plots of transposase
plasmid_data <- MGEs_level1%>% 
  filter(Taxonomy == "plasmid")
integrase_data <- MGEs_level1%>% 
  filter(Taxonomy == "integrase")
transposase_data <- MGEs_level1%>% 
  filter(Taxonomy == "transposase")
insertion_data <- MGEs_level1%>% 
  filter(Taxonomy == "insertion_element_IS91")
#making a new dataframe of MGEs
combined_MGEs <- rbind(plasmid_data, integrase_data, transposase_data,insertion_data)
combined_MGEs <- combined_MGEs[, -which(names(combined_MGEs) == "Tax_detail")]
combined_MGEs_sample <- combined_MGEs%>%
  pivot_longer(cols = -Taxonomy, names_to = "Sample", values_to = "Value")
combined_MGEs_sample_info<-merge(combined_MGEs_sample, sample_data, by = 'Sample')

#filtering combined data
filtered_plasmid_data <- combined_MGEs_sample_info %>%
  filter(Taxonomy == "plasmid")
filtered_trans_data <- combined_MGEs_sample_info %>%
  filter(Taxonomy == "transposase")
filtered_insertion_data <- combined_MGEs_sample_info %>%
  filter(Taxonomy == "insertion_element_IS91")
filtered_integrase_data <- combined_MGEs_sample_info %>%
  filter(Taxonomy == "integrase")

#filtering combined data

transp_boxplot <-ggbarplot(filtered_trans_data, x = "Substrate", y = "Value",
          add = c("mean_se", "point", size = 2),
          color = "Location", fill = "Location", alpha = 0.5,
          palette = c("#009E73", "#E7B800"))+ 
  labs(y = expression("Transposase"))+coord_flip()

plasmid_boxplot <- ggbarplot(filtered_plasmid_data, x = "Substrate", y = "Value",
          add = c("mean_se", "point", size = 2),
          color = "Location", fill = "Location", alpha = 0.5,
          palette = c("#009E73", "#E7B800"))+ 
  labs(y = expression("Plasmid"))+coord_flip()
 
insertion_boxplot <- ggbarplot(filtered_insertion_data, x = "Substrate", y = "Value",
                             add = c("mean_se", "point", size = 2),
                             color = "Location", fill = "Location", alpha = 0.5,
                             palette = c("#009E73", "#E7B800"))+ 
  labs(y = expression("IS91"))+coord_flip()

integrase_boxplot <- ggbarplot(filtered_integrase_data, x = "Substrate", y = "Value",
                               add = c("mean_se", "point", size = 2),
                               color = "Location", fill = "Location", alpha = 0.5,
                               palette = c("#009E73", "#E7B800"))+ 
  labs(y = expression("Integrase"))+coord_flip()

  