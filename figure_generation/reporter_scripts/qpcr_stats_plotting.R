library(dplyr) #(v1.2.0)
library(ggplot2) #(v4.0.2)
library(ggpubr) #(v0.6.3)
library(stats) #(v4.5.2)
library(RColorBrewer) #(v1.1-3)
library(tidyr) #(v1.3.2)

#Load in data and generating summary/descriptive statistics
data <- read.csv("qpcr_data.csv")
data$group <- factor(data$group, levels = c("control", "1", "2")) #switch group names as needed

summary_stats <- data %>%
  group_by(group) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value, na.rm = TRUE),
    n    = n(),
    se   = sd / sqrt(n)
  )
print(summary_stats)

#Perform one-way ANOVA with Tukey post-hoc (if ANOVA is significant)
anova_model <- aov(value ~ group, data = data)
summary(anova_model)

tukey <- TukeyHSD(anova_model)
print(tukey)

#Convert TukeyHSD results to a data frame
tukey_df <- as.data.frame(tukey$group)

#Add rownames as a column
tukey_df$comparison <- rownames(tukey_df)

#Split comparison into group1 and group2
groups <- do.call(rbind, strsplit(tukey_df$comparison, "-"))
tukey_df$group1 <- groups[,1]
tukey_df$group2 <- groups[,2]

#Rename p-value column to p.adj
names(tukey_df)[names(tukey_df) == "p adj"] <- "p.adj"

#Convert p-values to stars
tukey_df$p.adj <- symnum(tukey_df$p.adj, 
                         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                         symbols = c("****","***","**","*","ns"))

#Scale mean values relative to control
ref_mean <- summary_stats$mean[summary_stats$group == "control"]

summary_stats_2 <- summary_stats %>%
  mutate(
    ScaledMean = mean / ref_mean,
    ScaledSE = se / ref_mean
  )

#Select comparisons and their significance statistics to be displayed on plot
my_comparisons <- c("control-1", "control-2") #swap for whichever comparisons you want to be plotted with significance stars/ns

tukey_df_subset <- tukey_df %>%
  filter(comparison %in% my_comparisons)

tukey_df_subset$y.position <- c(10, 11) #height of significance stars/ns on plot

#Plotting the scaled mean from summary_stats_2 for all samples and selected comparisons with their significance statistics
ggplot(summary_stats_2, aes(x = group, y = ScaledMean)) +
  geom_bar(stat = "identity", color = "#000000", fill = "#000000", width = 0.4) +
  geom_errorbar(aes(ymin = ScaledMean - ScaledSE, ymax = ScaledMean + ScaledSE),
                width = 0.2, linewidth = 0.4) +
  stat_compare_means(
    method = "anova",
    label.y = NA
  ) +
  stat_pvalue_manual(
    tukey_df_subset,
    label = "p.adj",       
    tip.length = 0.01,     
    bracket.size = 0.4
  ) +
    labs(
    title = "Title",
    x = "Sample",
    y = "Relative Expression" 
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth = 0.5),
    axis.text = element_text(color = "black"),
    #axis.title = element_text(face = "bold")
  ) +
  theme(
    axis.title.x = element_text(size = 8, family = "Arial", margin = margin(t = 10)),  #axis titles
    axis.title.y = element_text(size = 8, family = "Arial", margin = margin(r = 10)),
    axis.text  = element_text(size = 8, family = "Arial", color = "black", margin = margin(t = 0)), #tick labels
    legend.text = element_text(size = 8, family = "Arial"),
    legend.title = element_text(size = 8, family = "Arial"),
    plot.title = element_text(size = 10, hjust = 0.5, margin = margin(b = 0))
  ) +
  scale_y_continuous(
    limits = c(0, 15),           #range of y-axis
    breaks = seq(0, 15, by = 3),  #specifying tick marks
    expand = c(0, 0)
  )