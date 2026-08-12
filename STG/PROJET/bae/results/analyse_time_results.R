library(dplyr)
library(tidyr)
library(ggplot2)

df <- datas

df <- df %>%
  rename(
    limit_dfs = limit_number_paths,
    nb_path = nb_path_to_select,
    #Benders_BA_master         = time_master_BA,
    #Benders_BA_subproblem     = time_subproblem_BA,
    #Benders_KC_master         = time_master_KC,
    #Benders_KC_subproblem     = time_subproblem_KC,
    #Benders_KCRDK_master      = time_master_KCRDK,
    #Benders_KCRDK_subproblem  = time_subproblem_KCRDK,
    Benders_KCRDKL_master   = time_master_KCRDKL,
    Benders_KCRDKL_subproblem = time_subproblem_KCRDKL,
    #Benders_KCU_master        = time_master_KCU,
    #Benders_KCU_subproblem    = time_subproblem_KCU,
    #Benders_KCUD_master       = time_master_KCUD,
    #Benders_KCUD_subproblem   = time_subproblem_KCUD,
    #Benders_OG_master        = time_master_KCHOG,
    #Benders_OG_subproblem    = time_subproblem_KCHOG,
    Benders_OGL_master        = time_master_KCHOGL,
    Benders_OGL_subproblem    = time_subproblem_KCHOGL
  )

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
    #ratio_BA = Benders_BA_master / (Benders_BA_master + Benders_BA_subproblem),
    #ratio_KC = Benders_KC_master / (Benders_KC_master + Benders_KC_subproblem),
    #ratio_KCU = Benders_KCU_master / (Benders_KCU_master + Benders_KCU_subproblem),
    #ratio_KCUD = Benders_KCUD_master / (Benders_KCUD_master + Benders_KCUD_subproblem),
    #ratio_KCRDK = Benders_KCRDK_master / (Benders_KCRDK_master + Benders_KCRDK_subproblem),
    ratio_KCRDKL = Benders_KCRDKL_master / (Benders_KCRDKL_master + Benders_KCRDKL_subproblem),
    #ratio_OG = Benders_OG_master / (Benders_OG_master + Benders_OG_subproblem),
    ratio_OGL = Benders_OGL_master / (Benders_OGL_master + Benders_OGL_subproblem)
  )

cat("Average proportion of time spent in the Master:\n")
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
df_long$tau_group <- cut(df_long$tau, breaks = 3, labels = c("Low tau", "Medium tau", "High tau"))

df_ranking <- df_long %>%
  group_by(Gamma, tau_group, Method) %>%
  summarise(Total_Time = sum(Time), .groups = "drop") %>%
  group_by(tau_group, Method) %>%
  summarise(Mean_Time = mean(Total_Time), .groups = "drop") %>%
  group_by(tau_group) %>%
  mutate(
    Rank = rank(Mean_Time),
    Label = ifelse(Rank == 1,
                   sprintf("Rang 1\n(%.3f s)", Mean_Time),
                   sprintf("Rang %d\n(%.3f s)", Rank, Mean_Time))
  )

p4 <- ggplot(df_long, aes(x = Gamma, y = Time, fill = Component)) +
  geom_area(position = "stack", stat = "summary", fun = mean, alpha = 0.8) +
  geom_text(
    data = df_ranking,
    aes(x = 50, y = 2.5, label = Label),
    inherit.aes = FALSE,
    size = 3.5,
    fontface = "bold",
    color = "gray20",
    vjust = 1
  ) +
  facet_grid(tau_group ~ Method) +
  theme_minimal() +
  labs(
    title = expression("Time allocation relative to " * Gamma),
    x = expression(Gamma),
    y = "Cumulative mean time (seconds)"
  ) +
  theme(legend.position = "bottom")

# Fig.05. Temps cumulé par rapport à nb_path (Aires empilées)
p5 <- ggplot(df_long, aes(x = nb_path, y = Time, fill = Component)) +
  geom_area(position = "stack", stat = "summary", fun = mean, alpha = 0.8) +
  facet_wrap(~ Method) +
  theme_minimal() +
  labs(
    title = "Cumulative time by nb_path",
    x = "Number of paths (nb_path)",
    y = "Cumulative average time (s)"
  ) +
  theme(legend.position = "bottom")

# Fig.06. Temps cumulé par rapport à limit_dfs (Aires empilées)
p6 <- ggplot(df_long, aes(x = limit_dfs, y = Time, fill = Component)) +
  geom_area(position = "stack", stat = "summary", fun = mean, alpha = 0.8) +
  facet_wrap(~ Method) +
  theme_minimal() +
  labs(
    title = "Cumulative time by limit_dfs",
    x = "DFS limit (limit_dfs)",
    y = "Cumulative average time (s)"
  ) +
  theme(legend.position = "bottom")

# Fig.07. Évolution lissée selon nb_path (Courbes séparées pour voir les croisements)
p7 <- ggplot(df_long, aes(x = nb_path, y = Time, color = Component)) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.2) +
  facet_wrap(~ Method) +
  theme_bw() +
  labs(
    title = "Smoothed trend based on nb_path",
    x = "Number of paths (nb_path)",
    y = "Mean time (s)",
    color = "Component"
  ) +
  theme(legend.position = "bottom")

# Fig.08. Évolution lissée selon limit_dfs (Courbes séparées)
p8 <- ggplot(df_long, aes(x = limit_dfs, y = Time, color = Component)) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.2) +
  facet_wrap(~ Method) +
  theme_bw() +
  labs(
    title = "Smoothed trend based on limit_dfs",
    x = "DFS limit (limit_dfs)",
    y = "Mean time (s)",
    color = "Component"
  ) +
  theme(legend.position = "bottom")

# Display
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
print(p8)
