############################################################
# Cox regression forest plot: progression endpoint
# GSE72970 source fields: pfs and pfs censored
############################################################

library(ggplot2)
library(dplyr)
library(svglite)

df <- data.frame(
  Group = c("Low", "High"),
  Patients = c(62, 62),
  Events = c(55, 59),
  HR = c(1.00, 1.232154807),
  lower = c(NA, 0.850574418),
  upper = c(NA, 1.784917859),
  P = c(NA, 0.269569937)
)

df$Group <- factor(df$Group, levels = c("High", "Low"))

df <- df %>%
  mutate(
    HR_CI = ifelse(
      Group == "Low",
      "1.00 (Reference)",
      sprintf("%.2f (%.2f–%.2f)", HR, lower, upper)
    ),
    P_label = ifelse(is.na(P), "—", sprintf("%.3f", P))
  )

p <- ggplot(df, aes(y = Group)) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.5,
    color = "#808080"
  ) +
  geom_errorbarh(
    data = df %>% filter(Group == "High"),
    aes(xmin = lower, xmax = upper),
    height = 0,
    linewidth = 0.9,
    color = "#A63D40"
  ) +
  geom_point(
    data = df %>% filter(Group == "High"),
    aes(x = HR),
    size = 3.2,
    shape = 16,
    color = "#A63D40"
  ) +
  geom_point(
    data = df %>% filter(Group == "Low"),
    aes(x = HR),
    size = 2.6,
    shape = 16,
    color = "#555555"
  ) +
  scale_x_log10(
    limits = c(0.7, 3.2),
    breaks = c(0.75, 1, 1.5, 2, 3),
    labels = c("0.75", "1", "1.5", "2", "3")
  ) +
  labs(
    x = "Hazard ratio (95% CI)",
    y = NULL,
    title = "Progression-free interval",
    subtitle = "High: HR 1.232 (0.851–1.785)   P = 0.270"
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 6.3, color = "#262626", hjust = 1),
    axis.title.x = element_text(size = 8.5),
    axis.text.x = element_text(size = 7.5, color = "#262626"),
    axis.text.y = element_text(size = 8.5, color = "#262626"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(
  "Cox_PFI_forest_plot.pdf", p,
  width = 3.2, height = 2.3, device = cairo_pdf
)

ggsave(
  "Cox_PFI_forest_plot.svg", p,
  width = 3.2, height = 2.3, device = svglite
)

ggsave(
  "Cox_PFI_forest_plot_600dpi.png", p,
  width = 3.2, height = 2.3, dpi = 600, bg = "white"
)

table_df <- df %>%
  mutate(`HR (95% CI)` = HR_CI, `Cox P` = P_label) %>%
  select(Group, Patients, Events, `HR (95% CI)`, `Cox P`)

write.csv(table_df, "Cox_PFI_results_table.csv", row.names = FALSE)

