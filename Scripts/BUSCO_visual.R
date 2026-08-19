## The following code was developed by Laurel C. Diaz on 8/11/2026
## This code is for visualizing the BUSCO output for the Orbicella genome annotation project

## Load libraries 
library(ggplot2)
library(dplyr)
library(tidyverse)
library(here)

busco_results <- data.frame(
  species = c("O. franksi", "O.annularis", "O.faveolata", "C. salae"),
  complete_single = c(95.8, 91.0, 89.1, 95.1),
  complete_duplicated = c(1.5, 4.9, 1.4, 1.8),
  fragmented = c(0.5, 0.6, 0.6, 0.8), 
  missing = c(2.2, 3.5, 8.9, 2.3)
)

busco_long <- busco_results %>%
  pivot_longer(cols = -species, 
               names_to = "Category", 
               values_to = "Percentage")

busco_long <- busco_long %>%
  mutate(vjust_val = case_when(
    Category == "complete_single"      ~ 0.5,    # middle
    Category == "complete_duplicated"  ~ 0.5,    # middle
    Category == "fragmented"           ~ 0.9,    # lower
    Category == "missing"              ~ -0.5     # middle
  ))

busco_plot <- ggplot(busco_long, aes(x = species, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%"),
                vjust = vjust_val),
            position = position_stack(vjust = 0.5),
            fontface = "bold", 
            size = 5,
            color = "black") +
  scale_fill_manual(values = c("#56B4E9", "#009E73", "#F0E442", "#CC79A7"),
                    labels = c("Complete (duplicated)", "Complete (single)",
                               "Fragmented", "Missing")) +
  coord_flip() +
  labs(title = "BUSCO Assessment Results",
       x = "Species",
       y = "Percentage (%)",
       fill = "Category") +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 12, face = "italic"),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12))
busco_plot

ggsave(plot = busco_plot, here("Outputs", "busco_plot.png"), width = 12, height = 8)
