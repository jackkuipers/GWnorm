
library(tidyverse)
library(scales)   # for palette functions

load("../sim_results/g10prime6.rData")

out_data %>% mutate(var = sd^2, time = as.numeric(time)) -> plot_data

plot_data %>% group_by(method, n_samp) %>% 
  summarise(m_var = exp(median(log(var))), m_time = exp(median(log(time)))) -> plot_data_m

ggplot(plot_data, aes(x = time, y = var, color = method, group = interaction(method, seed))) +
  geom_line(alpha = 0.4) + scale_x_log10() + scale_y_log10() + theme_bw() +
  geom_line(
    data = plot_data_m,
    aes(x = m_time, y = m_var, color = method, group = method),
    linewidth = 1.5
  )

# Choose a palette for each method
method_palettes <- list(
  GW = colorRampPalette(brewer_pal(palette = "Blues")(9)),
  BD = colorRampPalette(brewer_pal(palette = "Reds")(9))
)

# colour per tau
# Create color vector
unique_methods <- unique(plot_data$method)
# Build colour vector
col_vec <- c()

for (m in unique_methods) {
  
  tau_m <- sort(unique(plot_data$tau[plot_data$method == m]))
  n <- length(tau_m)
  
  cols_m <- c(method_palettes[[m]](5)[4], method_palettes[[m]](2*n)[-(1:n)])
  
  names(cols_m) <- c(m, paste(m, tau_m, sep = "."))
  
  col_vec <- c(col_vec, cols_m)
}

av_tau <- mean(plot_data$tau)

method_labs <- c("BDgraph", "GWnorm")

p <- ggplot(plot_data, aes(x = time, y = var)) +
  geom_line(aes(color = interaction(method, tau), group = interaction(method, seed)), linewidth = 0.5, alpha = 0.5) +
  scale_color_manual(values = col_vec, 
                     breaks = paste(unique_methods, ceiling(av_tau), sep = "."),      # only show method colors in legend
                     labels = method_labs, name = "Method"
  ) + xlab("Time (s)") + ylab("Variance") +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim = c(0.04, 1.4)) +
  theme_bw() + 
  geom_line(
    data = plot_data_m,
    aes(x = m_time, y = m_var, color = method, group = method),
    linewidth = 1.5, linetype = 2
  )


print(p)

ggsave("converge10.pdf", p, height = 5, width = 7)

lm(log(m_var) ~ log(m_time), filter(plot_data_m, method == "BD"))
lm(log(m_var) ~ log(m_time), filter(plot_data_m, method == "GW"))
