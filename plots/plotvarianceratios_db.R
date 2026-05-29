
library(tidyverse)

set.seed(26) # seed for jitter - just the year

plot_df <- NULL

for (p in 5*(2:4)) {
  
for (density_level in 1:3) {
    
for (beta in 10^(0:2)){
  
if (file.exists(paste0("../sim_results/g", p, "d", density_level, "beta", beta, ".rData"))) {  
  
load(paste0("../sim_results/g", p, "d", density_level, "beta", beta, ".rData"))

# times are in minutes!
out_data %>% mutate(var = sd^2, exact = if_else(var==0, 1, 0), fail = if_else(is.na(var), 1, 0), time = 60*as.numeric(time), Vt = var*time) -> plot_data

plot_data %>% select(-n_samp, -time, -mean, -sd, -var) %>% 
  pivot_wider(names_from = method, values_from = c(Vt, exact, fail, na_count)) %>% mutate(rat = Vt_BD/Vt_GW) -> plot_df_p

plot_df <- rbind(plot_df, plot_df_p)
}
}
}
}

exact_val <- 1e9
fail_val <- 1e10
# clean up some fail values
plot_df %>% mutate(exact_BD = if_else(is.na(exact_BD), 0, exact_BD)) %>% 
  mutate(rat = if_else((exact_BD == 1) & (exact_GW == 1), 1, rat)) %>%
  mutate(rat = if_else((exact_BD == 0) & (exact_GW == 1), exact_val, rat)) %>% 
  mutate(rat = if_else((fail_BD == 1) & (fail_GW == 1), 1, rat)) %>%
  mutate(rat = if_else((fail_BD == 1) & (fail_GW == 0), fail_val, rat)) %>%
  mutate(fail = fail_BD + 2*fail_GW, exact = exact_GW + 2*exact_BD) -> plot_df_clean

p <- ggplot(plot_df_clean, aes(x = as.factor(p), y = rat)) +
  geom_violin(colour = "dodgerblue", fill = "dodgerblue", alpha = 0.2, linewidth = 0.2) + 
  geom_boxplot(colour = "dodgerblue", outliers = FALSE) +
  geom_jitter(aes(shape = interaction(as.factor(exact), as.factor(fail))), width = 0.25, colour = "steelblue") + 
  scale_y_log10() + 
  #coord_cartesian(ylim = c(1e-3, 1e8)) + 
  geom_hline(yintercept = 1) + xlab("p") + ylab("Variance ratio") +
  scale_shape_manual(values = c(`0.0` = 16, `0.1` = 4, `1.0` = 5)) + 
  theme_bw() + guides(shape = "none") +
  facet_grid(dl ~ beta, 
    labeller = labeller(
      beta = as_labeller(function(x) {parse(text = paste0("beta==", x))}, label_parsed),
      dl = c("1"="Lower density", "2"="Medium density", "3"="Higher density")
    )
  )

print(p)

ggsave("varianceratios_db.pdf", p, height = 9, width = 9)

