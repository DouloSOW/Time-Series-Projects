# Question 3(a): ideal low-pass weights, cutoff at periods of 32 quarters
# w_0 = 1/16,  w_j = sin(j * pi/16) / (pi * j)  for j != 0

library(ggplot2)

omega_bar <- 2 * pi / 32
j <- -40:40
w <- ifelse(j == 0, omega_bar / pi, sin(j * omega_bar) / (pi * j))
weights <- data.frame(j = j, w = w)

p_w <- ggplot(weights, aes(j, w)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_segment(aes(xend = j, yend = 0), colour = "#2166AC", linewidth = 0.7) +
  geom_point(colour = "#2166AC", size = 1.6) +
  scale_x_continuous(breaks = seq(-40, 40, by = 8)) +
  labs(
    x = "j (quarters)",
    y = expression(w[j]),
    title = "Ideal low-pass filter weights (periods > 32 quarters)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
    panel.grid.minor = element_blank()
  )

print(p_w)

# ggsave("weights_3a.png", p_w, width = 7, height = 4, dpi = 300)