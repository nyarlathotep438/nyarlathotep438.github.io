# 2026年3月20日 研究日志
绘制可用的会议展示用图

↓根据之前生成的生存分析结果绘制KM生存曲线


```R
library(readr)
library(tidyr)
library(here)
library(ggplot2)
library(dplyr)
library(patchwork)
library(gridExtra)
library(tibble)
library(scales)
options(repr.plot.width=8.8, repr.plot.height=6, repr.plot.res=600)

# Basic information
PatientsInSGL <- 98897
PatientsInGLP <- 98902

# Input data(Additional annotations have been remove before)
data <- read.csv(here("Outcome_7_Result_b_KM_graph.csv"))
colnames(data) <- c("Time", 
                    "Surv1", "Lower1", "Upper1",
                    "Surv2", "Lower2", "Upper2")

# transfer to long format
plot_data <- data %>%
  pivot_longer(
    cols = -Time,
    names_to = c(".value", "Cohort"),
    names_pattern = "([A-Za-z]+)(\\d)"
  ) %>%
  mutate(
    Cohort = paste0("Cohort", Cohort)
  )

plot_data <- plot_data %>%
  mutate(Cohort = recode(Cohort,
                         "Cohort1" = "SGLT2i",
                         "Cohort2" = "GLP1RA"))

# Check data frame
if (nrow(plot_data) == 0 || all(is.na(plot_data$Lower))) {
  stop("Data conversion error: No confidence interval data was detected after conversion. Please check the original CSV structure.")
}

# Plot
# Color
color_SGLT2i <- "#8E00FA"  
fill_SGLT2i <- adjustcolor(color_SGLT2i, alpha.f = 0.2) 

color_GLP1RA <- "#117833" 
fill_GLP1RA <- adjustcolor(color_GLP1RA, alpha.f = 0.2) 

# label data
label_data <- data.frame(
  Cohort = c("SGLT2i", "GLP1RA"),
  Time = c(1000, 1000),  # fixed in x=1000
  Surv = c(
    approx(
      x = plot_data$Time[plot_data$Cohort == "SGLT2i"],
      y = plot_data$Surv[plot_data$Cohort == "SGLT2i"],
      xout = 1000
    )$y,
    approx(
      x = plot_data$Time[plot_data$Cohort == "GLP1RA"],
      y = plot_data$Surv[plot_data$Cohort == "GLP1RA"],
      xout = 1000
    )$y
  )
) %>%
  mutate(
    Surv_offset = ifelse(Cohort == "GLP1RA", Surv - 0.05, Surv + 0.05)
  )

# Create number table
number_table <- data %>%
  filter(Time %in% c(250, 500, 750, 1000)) %>%
  select(Surv1, Surv2) %>%
  `row.names<-`(as.character(c(250, 500, 750, 1000))) %>%
  rbind(data.frame(Surv1 = 1, Surv2 = 1, row.names = "0")) %>%
  slice(c(n(),1:n()-1)) %>%
  mutate(Surv1 = Surv1 * PatientsInSGL, Surv2 = Surv2 * PatientsInGLP) %>%
  mutate(across(everything(), round)) 

Ready_for_plot <- t(number_table) %>%
  as.data.frame(check.names = FALSE) %>%
  rownames_to_column("Drug") %>%
  mutate(Drug = recode(Drug, 
                       "Surv1" = "SGL T2i",
                       "Surv2" = "GLP-1RA")) %>%
  pivot_longer(
    cols = -Drug,
    names_to = "Time",
    values_to = "Value"
  ) %>%
  mutate(Value_num = as.numeric(Value)) %>%
  group_by(Drug) %>%
  mutate(
    Diff = first(Value_num) - Value_num,
    # table format
    Value_formatted = sprintf(
      "%s (%s)", 
      format(Value_num, big.mark = ",", trim = TRUE),  
      format(round(Diff), big.mark = ",", trim = TRUE) 
    )
  ) %>%
  ungroup() %>%
  select(Drug, Time, Value_formatted) %>%
  pivot_wider(
    names_from = Time,
    values_from = Value_formatted
  ) %>%
  column_to_rownames("Drug")

# Main plot----
main_plot <- ggplot(plot_data, aes(x = Time, y = Surv)) +
  geom_ribbon(
    aes(ymin = Lower, ymax = Upper, fill = Cohort), 
    alpha = 0.2, 
    linetype = "dashed",
    size = 0.3,
    color = "gray50"
  ) +
  geom_step(aes(color = Cohort), linewidth = 1.2) +
  
  # y axis scale
  scale_y_continuous(
    limits = c(0.969, 1.0),
    breaks = seq(0.97, 1.0, by = 0.01),  # identify y-axis range
    expand = c(0, 0),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  
  # x axis scale
  scale_x_continuous(
    expand = c(0, 0),
    breaks = function(limits) {
      # ensure 0 exist
      breaks <- unique(c(0, seq(0, max(limits), by = 250)))  # x-axis range: 250days
      breaks[breaks <= max(limits)]  
    }
  ) +
  
  labs(
    x = "Days", 
    y = "Survival Probability",
    title = "KM Curve(Mortality)"
  ) +
  
  scale_color_manual(
    name = NULL,
    values = c("SGLT2i" = color_SGLT2i, "GLP1RA" = color_GLP1RA),
    labels = c("SGLT2i" = "SGLT2i", "GLP1RA" = "GLP-1RA")
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("SGLT2i" = fill_SGLT2i, "GLP1RA" = fill_GLP1RA),
    labels = c("SGLT2i" = "SGLT2i", "GLP1RA" = "GLP-1RA")
  ) +
  
  annotate("text", 
           x = Inf, y = Inf, 
           label = "P = 0.0003", 
           hjust = 1.1, 
           vjust = 1.5,
           size = 5) +
  
  theme_bw(base_size = 14) +
  theme(
    legend.position = c(0.02, 0.02),
    legend.justification = c(0, 0),
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.key.size = unit(0.8, "cm"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.line = element_blank(),
    # use scale line
    axis.ticks = element_line(color = "black", size = 0.5),  # scale line size
    axis.ticks.length = unit(0.2, "cm"),  # scale line length
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(margin = margin(r = 10)),
    plot.margin = margin(10, 20, 10, 20)
  ) +
  coord_cartesian(clip = "off") +
  
  # show 0 label
  annotate("text", 
           x = 0, y = -Inf, 
           label = "", 
           vjust = 1.5, 
           hjust = 0.5,
           size = 4)

main_plot

# Table plot----
table_data <- as.data.frame(Ready_for_plot)

x_min <- -150
x_max <- 1100
x_breaks <- c(0, 250, 500, 750, 1000)

table_plot <- ggplot() +
  # title
  ggtitle("All-cause mortality") +
  
  scale_x_continuous(
    breaks = x_breaks,
    limits = c(x_min, x_max),
    expand = c(0, 0)
  ) +
  # Data
  annotate("text", 
           x = x_breaks,
           y = 1, 
           label = as.character(table_data["SGL T2i", ]),
           size = 4.5,
           fontface = "bold",
           hjust = 0.5) +
  
  annotate("text", 
           x = x_breaks,
           y = 0, 
           label = as.character(table_data["GLP-1RA", ]),
           size = 4.5,
           fontface = "bold",
           hjust = 0.5) +
  
  # Drug name
  annotate("text",
           x = -100,
           y = 1,
           label = "SGL T2i",
           size = 4.5,
           fontface = "bold",
           hjust = 1) +
  
  annotate("text",
           x = -100,
           y = 0,
           label = "GLP-1RA",
           size = 4.5,
           fontface = "bold",
           hjust = 1) +
  
  # Set y axis limit
  ylim(-0.5, 1.5) +
  
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    plot.margin = margin(20, 50, 10, 100),  
    plot.title = element_text(              
      size = 16, 
      face = "bold", 
      hjust = 0, 
      vjust = 1.5,
      margin = margin(b = 3)
    )
  )

table_plot

# Final plot(merge main_plot & table_plot)
final_plot <- (
  # adjust main plot
  main_plot +
    scale_x_continuous(
      breaks = c(0, 250, 500, 750, 1000),
      expand = c(0, 0)
    ) +
    theme(
      plot.margin = margin(0, 20, 0, 50)  # Control margin
    )
) / (
  # adjust table plot
  ggplot() +
    # Table data
    annotate("text", 
             x = c(0, 250, 500, 750, 1000),
             y = 1, 
             label = as.character(table_data["SGL T2i", ]),
             size = 4.5,
             fontface = "bold",
             hjust = 0.5) +
    annotate("text", 
             x = c(0, 250, 500, 750, 1000),
             y = 0, 
             label = as.character(table_data["GLP-1RA", ]),
             size = 4.5,
             fontface = "bold",
             hjust = 0.5) +
    # Drug name
    annotate("text",
             x = -73, 
             y = 1,
             label = "SGLT2i",
             size = 4.5,
             fontface = "bold",
             hjust = 1) +
    annotate("text",
             x = -74,
             y = 0,
             label = "GLP-1RA",
             size = 4.5,
             fontface = "bold",
             hjust = 1) +
    # Table title
    annotate("text",
             x = -160,  
             y = 1.5,    
             label = "All-cause mortality",
             size = 6,
             fontface = "bold",
             hjust = 0) + 
    # Set location
    ylim(-0.5, 1.5) +
    scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(clip = "off", xlim = c(-20, 1100)) +  # Limit shown region
    theme_void() +
    theme(plot.margin = margin(0, 20, 10, 20))  # balance margin
) + 
  plot_layout(heights = c(3, 1))

# show final result
final_plot

# Output
formats <- c("pdf", "png")
for(fmt in formats) {
  ggsave(paste0("KMFinal_plot_ver2.", fmt),
         plot = final_plot,
         dpi = ifelse(fmt=="png", 1200, 300),
         width=8.8, height=6)
}
```

*今天翻来覆去搞这个图真是耗费心神，今天就这样吧，只有再看看怎么给博客添加图片*
