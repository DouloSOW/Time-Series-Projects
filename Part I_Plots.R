# Question 2(h): spectra of Delta eta_t, Delta tau_t and Delta y_t
# Model: y_t = tau_t + eta_t, tau_t = tau_{t-1} + v_t
# Spectral density> f(w) = (1/2pi) * sum_h gamma(h) e^{-ihw}, w in [-pi, pi]

library(ggplot2)
library(patchwork)

# ---------- spectra ----------
# Delta eta: |1 - e^{-iw}|^2 * s2eta / 2pi ; Delta tau: white noise
# v and eta independent, so the spectra add as showed
spectra_df <- function(w, s2eta, s2v) {
  deta <- s2eta * (2 - 2 * cos(w)) / (2 * pi)
  dtau <- rep(s2v / (2 * pi), length(w))
  data.frame(w = w, deta = deta, dtau = dtau, dy = deta + dtau)
}

w <- seq(-pi, pi, length.out = 500)

# ---------- shared styling ----------
pi_scale <- scale_x_continuous(breaks = c(-pi, -pi/2, 0, pi/2, pi), labels = expression(-pi, -pi/2, 0, pi/2, pi))


#cutoff <- geom_vline(xintercept = c(-pi/16, pi/16), linetype = "dotted", colour = "grey45")

common_theme <- theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(margin = margin(r = 6))
  )

axis_labels <- labs(x = expression(omega), y = "Spectral density")

# ---------- Panel 1: components at the 2(d) values ----------
s <- spectra_df(w, s2eta = 1, s2v = 0.5)
comp <- rbind(
  data.frame(w = w, f = s$dy, series = "dy"),
  data.frame(w = w, f = s$deta, series = "deta"),
  data.frame(w = w, f = s$dtau, series = "dtau")
)
comp$series <- factor(comp$series, levels = c("dy", "deta", "dtau"))

p1 <- ggplot(comp, aes(w, f, colour = series)) +
  #cutoff +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(
    name = NULL,
    breaks = c("dy", "deta", "dtau"),
    values = c(dy = "black", deta = "#B2182B", dtau = "#2166AC"),
    labels = expression(Delta * y[t], Delta * eta[t], Delta * tau[t])
  ) +
  pi_scale +
  axis_labels +
  labs(title = expression(paste("Components:  ", sigma[eta]^2, " = 1,  ", sigma[v]^2, " = 0.5"))) +
  common_theme

# ---------- Panel 2: Delta y, varying the noise variance ----------
s2eta_grid <- c(0.25, 1, 4)
vary_eta <- do.call(rbind, lapply(s2eta_grid, function(s2) {
  data.frame(w = w, f = spectra_df(w, s2, 0.5)$dy, par = factor(s2))
}))

p2 <- ggplot(vary_eta, aes(w, f, colour = par)) +
  #cutoff +
  geom_line(linewidth = 1.1) +
  scale_colour_viridis_d(name = expression(sigma[eta]^2), end = 0.85) +
  pi_scale +
  axis_labels +
  labs(title = expression(paste(Delta * y[t], ":  varying ", sigma[eta]^2, "  (", sigma[v]^2, " = 0.5)"))) +
  common_theme

# ---------- Panel 3: Delta y, varying the trend-shock variance ----------
s2v_grid <- c(0.1, 0.5, 2)
vary_v <- do.call(rbind, lapply(s2v_grid, function(s2) {
  data.frame(w = w, f = spectra_df(w, 1, s2)$dy, par = factor(s2))
}))

p3 <- ggplot(vary_v, aes(w, f, colour = par)) +
  #cutoff +
  geom_line(linewidth = 1.1) +
  scale_colour_viridis_d(name = expression(sigma[v]^2), end = 0.85) +
  pi_scale +
  axis_labels +
  labs(title = expression(paste(Delta * y[t], ":  varying ", sigma[v]^2, "  (", sigma[eta]^2, " = 1)"))) +
  common_theme

# ---------- combine and show ----------
print(p1 | p2 | p3)


# ggsave("spectra_2h.png", p1 | p2 | p3, width = 13, height = 4.5, dpi = 300)