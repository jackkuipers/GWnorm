
library(tidyverse)

res <- NULL

for (seednumber in 101:150) {#<- 101
  top_i <- 7
  
  fileName <- paste0("../iris_out/res", top_i, "seed", seednumber, ".Rdata")
  if (file.exists(fileName)) {
    load(file = fileName)
    res <- rbind(res, resall)
  }
}

res %>% mutate(log_t = log(time), log_d = log(rmse)) %>%
  mutate(method = case_when(
    method == "BD_rj" ~ "BD",
    method == "BDmpl_rj" ~ "BD.mpl",
    TRUE ~ method
  )) %>% filter(method != "BDmpl_bd") -> plot_df

#plot_df %>% group_by(method, is) %>% summarise(m_log_t = mean(log_t), m_log_d = mean(log_d)) %>%
#  mutate(time = exp(m_log_t), rmse = exp(m_log_d)) -> plot_summary_df

plot_df %>% group_by(method, is) %>% summarise(m_t = mean(time), m_d = sqrt(mean(rmse^2))) %>%
  mutate(time = m_t, rmse = m_d, log_t = log2(time), log_d = log2(rmse)) -> plot_summary_df

ggplot() + aes(x = time, y = rmse, colour = method) + geom_point(data = plot_df, alpha = 0.2) + scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2') + geom_point(data = plot_summary_df, size = 4) + geom_line(data = plot_summary_df)

ggplot(plot_df, aes(x = as.factor(is), y = rmse, colour = method)) + geom_boxplot() + 
  scale_y_continuous(trans = 'log2') + geom_point(position = position_jitterdodge(), alpha = 0.2)

ggplot() + aes(x = time, y = rmse, colour = method) + 
  geom_point(data = plot_df, alpha = 0.2) + scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2') +
  #scale_y_continuous(trans = 'log2', breaks = c(2^(-6:-1)),
  #             labels = c(expression(frac(1,64)), expression(frac(1,32)), expression(frac(1,16)), "1/8", "1/4", "1/2"))+ 
  geom_point(data = plot_summary_df, size = 4) + geom_line(data = plot_summary_df)

ggplot() + aes(x = time, y = rmse, colour = method) + 
  geom_point(data = plot_df, alpha = 0.2) + scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2', breaks = c(4^(-4:-1)),
               labels = c(expression(frac(1,256)),expression(frac(1,64)), expression(frac(1,16)), "1/4"))+ 
  geom_point(data = plot_summary_df, size = 4, shape = 3, stroke = 2) + 
  geom_line(data = plot_summary_df)


df_test <- subset(plot_summary_df, method == "WWA" & time > 10)
fit <- lm(log2(rmse) ~ log2(time), data = df_test)
a <- coef(fit)[1]
b <- coef(fit)[2]
pred_df <- data.frame(
  time = seq(
    min(plot_summary_df$time),
    max(plot_summary_df$time),
    length.out = 200
  )
)
pred_df$rmse <- 2^a * pred_df$time^b
pred_df$method <- "WWA"

p <- ggplot() + aes(x = time, y = rmse, colour = method) + 
  geom_point(data = plot_df, alpha = 0.33) + scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2', breaks = c(4^(-4:-1)),
                     labels = c(expression(frac(1,256)),expression(frac(1,64)), expression(frac(1,16)), "1/4"))+ 
  geom_point(data = plot_summary_df, size = 4, aes(shape = method), stroke = 2) + 
  geom_line(data = plot_summary_df) + theme_bw() + xlab("Time") + ylab("RMSE") +
  geom_line(data = pred_df, linetype = "dashed", linewidth = 1.5) +
  scale_color_manual(
    values = c(
      "BD" = "red",
      "BD.mpl" = "firebrick",
      "WWA" = "darkgreen"
    ),
    name = "Method"
  ) + 
  scale_shape_manual(
    values = c(
      "BD"     = 3,
      "BD.mpl" = 16,
      "WWA"    = 4
    ),
    name = "Method"
  ) +
  theme(legend.key.height = unit(1, "cm"))

print(p)

ggsave("iris_plot.pdf", p, height = 5, width = 7)
