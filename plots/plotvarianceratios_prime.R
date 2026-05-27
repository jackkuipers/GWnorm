
library(tidyverse)
library(lubridate)
library(latex2exp)

set.seed(26) # seed for jitter - just the year

plot_df <- NULL

for (p in 5*(2:4)) {
  
load(paste0("../sim_results/g",p,"prime1.rData"))

# times are in minutes!
out_data %>% mutate(var = sd^2, time = 60*as.numeric(time), Vt = var*time) -> plot_data

plot_data %>% group_by(method, seed) %>% summarise(Vt = mean(Vt), tau = mean(tau), p = p, tp = tau/choose(p,2)) %>% 
  pivot_wider(names_from = method, values_from = Vt) %>% mutate(rat = if_else(is.na(BD), 1e6, BD/GW)) -> plot_df_p

plot_df <- rbind(plot_df, plot_df_p)
}

p <- ggplot(plot_df, aes(x = as.factor(p), y = rat)) +
  geom_violin(colour = "dodgerblue", fill = "dodgerblue", alpha = 0.2, linewidth = 0.2) + 
  geom_boxplot(colour = "dodgerblue", outliers = FALSE) +
  geom_jitter(aes(color = tp, shape = is.na(BD)), width = 0.25) + scale_y_log10() + 
  coord_cartesian(ylim = c(1e-2, 1e8)) + 
  scale_colour_gradient(
    limits = c(0, 0.4),
    name = expression(tilde(tau))) + guides(shape = "none") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4)) + 
  geom_hline(yintercept = 1) + xlab("p") + ylab("Variance ratio") +
  theme_bw()

print(p)

ggsave("varianceratios_prime.pdf", p, height = 5, width = 7)

plot_df %>% group_by(p) %>% summarise(BD_NA = mean(is.na(BD)), GW_NA = mean(is.na(GW)))

# for size 20
plot_data %>% group_by(method) %>% summarise(na_frac = mean(na_count)/40)

plot_df %>% filter(p == 20) %>% group_by(is.na(BD)) %>% summarise(mtp = mean(tp))

plot_df %>% filter(p == 20, is.na(BD)) %>% pull(tp) %>% mean()

plot_df %>% group_by(p) %>% summarise(better = mean(rat > 1, na.rm=TRUE))
