library(devtools)
library(GSVA)
library(dplyr)
install_github("lydiaMyr/ImmuCellAI-mouse@main")
library(ImmuCellAImouse)
library(pracma)
library(tidyr)
library(ggplot2)
library(cowplot)

#the input counts file must be a tab delimited .txt file
#load and manipulate raw counts file
data <- read.table("counts_file.txt")
colnames(data) <- data[1,]
data <- data[-1,]

# Convert all count columns to numeric
data[-1] <- lapply(data[-1], as.numeric)

# Group by gene names and calculate the sum for each group(useful to consolidate counts data 
#from tools that counts isoforms of same genes inidividually )
data <- data %>%
  group_by_at(1) %>%
  summarise(across(everything(), sum, na.rm = TRUE))

#(Optional, if need to run in group mode)create a vector with the group details to incorporate in counts data
groups <- rep(c("Groups", "Control", "Treatment1", "Treatment2"), c(1, 10, 4, 4))
#(Optional, if need to run in group mode)add a row in data with group info
data <- rbind(groups, data)
#ensure the counts data is a data frame object
data <- as.data.frame(data)
#rownames should be the gene names/gene IDs and the data sheet should not contain gene details
rownames(data) <- data$geneID
data <- data[, -1]

#run ImmuCellAI_mouse in group mode
result <- ImmuCellAI_mouse(data, "rnaseq", 1, 0)
#run ImmucellAI_mouse in individual mode (uncomment)
#result <- ImmucellAI_mouse(data, "rnaseq", 0, 0)

#save results for downstream analysis or future reference
saveRDS(result, "result.rds")
write.csv(result$abundance, "Abundance.csv")

#optional- Load the results
result<- readRDS("result.rds")

#Downstream analysis and visualization of individual and group results
#applicable only if run in group mode
abundance <- as.data.frame(result$abundance)
group <- as.data.frame(result$group_result)
abundance$Group <- rep(c("Sham", "D7", "D70"), c(10, 4, 4))
abundance$Sample <- rownames(abundance)
group$Group <- rownames(group)
group_means <- abundance %>%
  group_by(Group) %>%
  summarise(across(B_cell:Infiltration_score, mean, na.rm = TRUE))

abundance_long <- gather(abundance, key = "Cell_Type", value = "Abundance", -Sample, -Group)
group_long <- group %>%
  pivot_longer(cols = -Group, names_to = "Cell_Type", values_to = "Mean_Abundance")


p_value <- group_long %>% filter(Group == "p value")
p_value$p_value <- p_value$Mean_Abundance
p_values <- p_value %>%
  mutate(p_value_label = paste("p =", p_value))
group_long <- group_long %>% filter(Group != "p value")
cell_types <- unique(abundance_long$Cell_Type)

group_long$Group <- factor(group_long$Group, levels = c("Sham", "D7", "D70"))
abundance_long$Group <- factor(abundance_long$Group, levels = c("Sham", "D7", "D70"))

#create box plots for all celltypes in a single graph
All_plots <- ggplot() +
  geom_boxplot(data = abundance_long, aes(x = Group, y = Abundance, fill = Group)) +
  geom_errorbar(data = group_long, 
                aes(x = Group, ymin = Mean_Abundance, ymax = Mean_Abundance), 
                width = 0.3, color = "red", size = 0.8) + # Custom mean lines
  facet_wrap(~ Cell_Type, scales = "free_y") +
  geom_text(data = p_values, aes(x = Inf, y = Inf, label = p_value_label), 
            hjust = 2, vjust = 1.5, size = 2.5, color = "red") + 
  scale_fill_manual(values = c("Sham" = "yellow", "D7" = "darkgreen", "D70" = "orange")) + # Custom colors
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),   # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    panel.border = element_blank(),       # Remove panel border
    plot.background = element_rect(fill = "white"),  # Set plot background to white
    panel.background = element_rect(fill = "white"), # Set panel background to white
    axis.line.x.bottom = element_line(color = "black", size = 0.5), # Add x-axis line at the bottom
    axis.line.y.left  = element_line(color = "black", size = 0.5), # Add y-axis line on the left
    axis.ticks = element_line(color = "black", size = 0.5), # Add major ticks
    axis.ticks.length = unit(0.3, "cm"), # Adjust length of major ticks
  ) +
  labs(title = "Boxplots and Means of Cell Types by Treatment- n(Sham)=10", y = "ImmuCellAI-Mouse.Score", x = "Treatment")

#save the combined plot as PNG and SVG
ggsave(filename = "plots/ImmuCellAI_result_sham-n(10).png", plot = All_plots, width = 12, height = 8, dpi = 300)
ggsave(filename = "plots/ImmuCellAI_result_sham-n(10).svg", plot = All_plots, device = "svg", width = 12, height = 8, dpi = 300)


#create an empty list to store all indivdual plot data
plot_list <- list()

#for loop to make plots for individual cell types and store in list
for (i in cell_types) {
  plot_list[[i]] <-  ggplot() +
    geom_boxplot(data = subset(abundance_long, Cell_Type == i), 
                 aes(x = Group, y = Abundance, fill = Group)) +
    geom_errorbar(data = subset(group_long, Cell_Type == i), 
                  aes(x = Group, ymin = Mean_Abundance, ymax = Mean_Abundance), 
                  width = 0.3, color = "red", size = 0.8) + # Custom mean lines
    geom_text(data = subset(p_values, Cell_Type == i), 
              aes(x = Inf, y = Inf, label = p_value_label), 
              hjust = 2.25, vjust = 1.5, size = 4, color = "red") +
    scale_fill_manual(values = c("Sham" = "yellow", "D7" = "darkgreen", "D70" = "orange")) + # Custom colors
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),   # Remove major grid lines
      panel.grid.minor = element_blank(),   # Remove minor grid lines
      panel.border = element_blank(),       # Remove panel border
      plot.background = element_rect(fill = "white"),  # Set plot background to white
      panel.background = element_rect(fill = "white"), # Set panel background to white
      axis.line.x.bottom = element_line(color = "black", size = 0.5), # Add x-axis line at the bottom
      axis.line.y.left  = element_line(color = "black", size = 0.5), # Add y-axis line on the left
      axis.ticks = element_line(color = "black", size = 0.5), # Add major ticks
      axis.ticks.length = unit(0.3, "cm"), # Adjust length of major ticks
    ) +
    labs(title = i, y = "ImmuCellAI-Mouse.Score", x = "Treatment")
}


# Loop through each plot in plot_list and save as SVG
for (i in names(plot_list)) {
  # Define the filename using the cell type name
  file_name <- paste0("plots/plot_", i, ".svg")
  
  # Save each plot as a SVG file
  ggsave(filename = file_name, plot = plot_list[[i]], device = "svg", width = 8, height = 6, dpi = 300)
}


# Loop through each plot in plot_list and save as PNG
for (i in names(plot_list)) {
  # Define the filename using the cell type name
  file_name <- paste0("plots/plot_", i, ".png")
  
  # Save each plot as a PNG file
  ggsave(filename = file_name, plot = plot_list[[i]], device = "png", width = 8, height = 6, dpi = 300)
}

