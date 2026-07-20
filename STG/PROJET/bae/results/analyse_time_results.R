library(dplyr)
library(tidyr)
library(ggplot2)

df <- test3

colnames(df) <- c("filename", "Gamma", "tau", 
                  "Benders_master", "Benders_subproblem",
                  "Benders_Matrix_master", "Benders_Matrix_subproblem",
                  "Augmented_master", "Augmented_subproblem", 
                  "HOG_master", "HOG_subproblem")

df_long <- df %>%
  pivot_longer(
    cols = ends_with("master") | ends_with("subproblem"),
    names_to = c("Method", "Component"),
    names_pattern = "(.*)_(master|subproblem)",
    values_to = "Time"
  )

# Proportion of time spent in each
df_summary <- df %>%
  mutate(
    ratio_Benders = Benders_master / (Benders_master + Benders_subproblem),
    ratio_Benders_matrix = Benders_Matrix_master / (Benders_Matrix_master + Benders_Matrix_subproblem),
    ratio_Augmented = Augmented_master / (Augmented_master + Augmented_subproblem),
    ratio_HOG = HOG_master / (HOG_master + HOG_subproblem)
  )

cat("Average proportion of time spent in the Master :\n")
print(summary(df_summary %>% select(starts_with("ratio"))))

# Fig.01. Execution time
p1 <- ggplot(df_long, aes(x = Method, y = Time, fill = Component)) +
  geom_boxplot(alpha = 0.8) +
  theme_minimal() +
  labs(title = "Execution time: Master vs. subproblem",
       x = "Benders variant",
       y = "Time (second)",
       fill = "Component") +
  theme(legend.position = "bottom")

# Fig.02. Relative to Gamma
p2 <- ggplot(df_long, aes(x = Gamma, y = Time, color = Component)) +
  geom_smooth(method = "loess", se = FALSE, size = 1.2) +
  facet_wrap(~ Method) +
  theme_bw() +
  labs(title = "Time evolution relatve to Gamma",
       x = expression(Gamma),
       y = "Mean time (second)")

# Fig.03. Impact of Tau
p3 <- ggplot(df_long, aes(x = tau, y = Time, color = Component)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.2) +
  facet_wrap(~ Method) +
  theme_bw() +
  labs(title = expression("Impact of the critical path threshold (" *tau* ") on execution times"),
       x = expression(tau),
       y = "Time (seconds)",
       color = "Component")

# Fig.04. Interaction between Gamma and Tau
df_long$tau_group <- cut(df_long$tau, breaks = 3, labels = c("Tau Faible", "Tau Moyen", "Tau Élevé"))

p4 <- ggplot(df_long, aes(x = Gamma, y = Time, fill = Component)) +
  geom_area(position = "stack", stat = "summary", fun = mean, alpha = 0.8) +
  facet_grid(tau_group ~ Method) +
  theme_minimal() +
  labs(title = expression("Time allocation: Interaction between " * Gamma * " et " * tau),
       x = expression(Gamma),
       y = "Cumulative mean time (seconds)")

# Display
print(p1)
print(p2)
print(p3)
print(p4)

