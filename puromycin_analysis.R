# Enzyme Kinetics Analysis: Puromycin

# 1. Installation & Setup
# install.packages(c("ggplot2", "dplyr"))
library(tidyverse)

# Load built-in dataset
data("Puromycin")

# 2. Michaelis-Menten Plot (Nonlinear Regression)
plot_mm <- ggplot(data = Puromycin, mapping = aes(x = conc, y = rate, color = state)) +
geom_point(size = 2) + geom_smooth(method = "nls", formula = y ~ (Vmax * x) / (Km + x), 
start = list(Vmax = 200, Km = 0.1), se = FALSE) + theme_light() + labs(title = "Enzyme Kinetics: Puromycin", 
x = "Substrate Concentration (ppm)", y = "Reaction Velocity (counts/min/min)", color = "Treatment State")

# Save Michaelis-Menten plot
ggsave("puromycin_michaelis_menten.png", plot = plot_mm, width = 7, height = 5, dpi = 300)

# 3. Lineweaver-Burk Plot (Linear Regression)
plot_lb <- ggplot(data = Puromycin, mapping = aes(x = 1/conc, y = 1/rate, color = state)) + 
geom_point(size = 3) + geom_smooth(method = "lm", se = FALSE, fullrange = TRUE) + theme_light() + 
labs(title = "Lineweaver-Burk Plot: Puromycin", x = "1 / Substrate Concentration (1/ppm)", 
y = "1 / Reaction Velocity (min/counts)", color = "Treatment State")

# Save Lineweaver-Burk plot
ggsave("puromycin_lineweaver_burk.png", plot = plot_lb, width = 7, height = 5, dpi = 300)

# 4. Statistical Summary & Quantile Analysis
# Calculate Means and Standard Deviations
untreated_mean <- Puromycin %>% filter(state == "untreated") %>% pull(conc) %>% mean(na.rm = TRUE)
untreated_sd <- Puromycin %>% filter(state == "untreated") %>% pull(conc) %>% sd(na.rm = TRUE)

treated_mean <- Puromycin %>% filter(state == "treated") %>% pull(conc) %>% mean(na.rm = TRUE)
treated_sd <- Puromycin %>% filter(state == "treated") %>% pull(conc) %>% sd(na.rm = TRUE)

untreated_mean_rate <- Puromycin %>% filter(state == "untreated") %>% pull(rate) %>% mean(na.rm = TRUE)
untreated_sd_rate <- Puromycin %>% filter(state == "untreated") %>% pull(rate) %>% sd(na.rm = TRUE)

treated_mean_rate <- Puromycin %>% filter(state == "treated") %>% pull(rate) %>% mean(na.rm = TRUE)
treated_sd_rate <- Puromycin %>% filter(state == "treated") %>% pull(rate) %>% sd(na.rm = TRUE)

# 5. Calculate 99.5% Confidence Interval: Rate & Concentration (Assuming Normal Distribution)
untreated_conc_michaelis <- qnorm(0.9975, untreated_mean, untreated_sd)
treated_conc_michaelis <- qnorm(0.9975, treated_mean, treated_sd)
michaelis_rate_treated <- qnorm(0.9975, treated_mean_rate, treated_sd_rate)
michaelis_rate_untreated <- qnorm(0.9975, untreated_mean_rate, untreated_sd_rate)

# 6. Boxplot for Rate & Concentration

plot_conc_box <- ggplot(Puromycin, aes(x = state, y = conc, fill = state)) + geom_boxplot(alpha = 0.7) + 
theme_light() + labs(title = "Puromycin Substrate Concentration by Treatment State", x = "Treatment State", y = "Substrate Concentration")

ggsave("puromycin_conc_boxplot.png", plot = plot_conc_box, width = 7, height = 5, dpi = 300)

plot_rate_box <- ggplot(Puromycin, aes(x = state, y = rate, fill = state)) + geom_boxplot(alpha = 0.7) + theme_light() + 
labs(title = "Puromycin Reaction Rate by Treatment State", x = "Treatment State",y = "Reaction Velocity")

ggsave("puromycin_rate_boxplot.png", plot = plot_rate_box, width = 7, height = 5, dpi = 300)

plot_rate_box
plot_conc_box

# 6. Parameters for Vmax & Km 
treated_model <- nls(rate ~ (Vmax * conc) / (Km + conc), data = filter(Puromycin, state == "treated"), start = list(Vmax = 200, Km = 0.1))
untreated_model <- nls(rate ~ (Vmax * conc) / (Km + conc), data = filter(Puromycin, state == "untreated"), start = list(Vmax = 200, Km = 0.1))
treated_params <- coef(treated_model)
untreated_params <- coef(untreated_model)
treated_Km <- treated_params["Km"]
treated_Vmax <- treated_params["Vmax"]
untreated_params <- coef(untreated_model)
untreated_Km <- untreated_params["Km"]
untreated_Vmax <- untreated_params["Vmax"]
cat("untreated Km:", untreated_Km,"\n")
cat("untreated Vmax:", untreated_Vmax,"\n")
cat("treated Km:", treated_Km,"\n")
cat("treated Vmax:", treated_Vmax,"\n")

parameter_df <- tibble(state = c("treated", "untreated", "treated", "untreated"), parameter = c("Vmax", "Vmax", "Km", "Km"), 
value = c(treated_Vmax, untreated_Vmax, treated_Km, untreated_Km))

plot_parameters <- ggplot(parameter_df, aes(x = state, y = value, fill = state)) + geom_col(width = 0.6) + facet_wrap(~ parameter, scales = "free_y") + theme_light() + 
labs( title = "Estimated Michaelis-Menten Parameters", x = "Treatment State", y = "Estimated Value")

plot_parameters
ggsave("puromycin_vmax_km_parameters.png", plot = plot_parameters, width = 7, height = 5, dpi = 300)
