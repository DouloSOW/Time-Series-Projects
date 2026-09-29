###############################################################################
# Core PCE inflation (excl. food & energy) — plotting and MA estimation
# Data   : FRED series DPCCRV1Q225SBEA (BEA), quarterly, % change at annual rate
# Refs open sources that I use:
##: https://nicolarighetti.github.io/Time-Series-Analysis-With-R/plot-time-series.html
##https://www.econometrics-with-r.org/14.8-niib.html
# Sample convention: first observation = y_0, last = y_T  (R index 1, ..., n)
###############################################################################


# =============================================================================
#  Packages
# =============================================================================
library(readxl)     # read .xlsx files
library(dplyr)      # data manipulation (provides dplyr::lag)
library(ggplot2)    # plotting
library(patchwork)  # combine several ggplots into one figure
library(lmtest)     # coeftest(): SEs, z-stats and p-values for arima fits
library(tidyr)     

# =============================================================================
#  Import the data
# -----------------------------------------------------------------------------
# Reads the "Quarterly" sheet of the FRED
# =============================================================================
sample <- read_excel("DPCCRV1Q225SBEA_PCE_excludingFoodAndEnergy.xlsx",
                     sheet = "Quarterly")
sample$observation_date <- as.Date(sample$observation_date)


# =============================================================================
# (1). Plot the series
# =============================================================================

# 1.1 Window to shade as the COVID-19 period
covid_start <- as.Date("2019-01-01")
covid_end   <- as.Date("2021-12-31")

# 1.2 Full-sample plot
plot1_series <- ggplot(sample, aes(x = observation_date, y = DPCCRV1Q225SBEA)) +
  annotate("rect", xmin = covid_start, xmax = covid_end,
           ymin = -Inf, ymax = Inf, fill = "grey70", alpha = 0.35) +
  annotate("text", x = covid_start, y = Inf, label = "COVID-19",
           hjust = -0.1, vjust = 1.5, size = 3.5, colour = "grey30") +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = 2, linetype = "dashed", colour = "firebrick") +
  geom_line(colour = "steelblue", linewidth = 0.6) +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  labs(title    = "Core PCE Price Index (excl. food & energy)",
       subtitle = "Quarterly, % change at annual rate; shaded area = COVID-19 period; dashed = 2% target",
       x = NULL, y = "Percent",
       caption  = "Source: BEA via FRED (DPCCRV1Q225SBEA)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

plot1_series

# 1.3 Zoomed plot (2010 onward)

plot1_series_covidZoom <- plot1_series +
  coord_cartesian(xlim = as.Date(c("2010-01-01", NA))) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")

plot1_series_covidZoom

# 1.4 Save both plots 
ggsave("plot1_series.png",           plot1_series,           width = 13, height = 4.5, dpi = 300)
ggsave("plot1_series_covidZoom.png", plot1_series_covidZoom, width = 13, height = 4.5, dpi = 300)

# 1.5 Combine both plots into one figure
plot_combined <- plot1_series / plot1_series_covidZoom
plot_combined
ggsave("plot1_combined.png", plot_combined, width = 13, height = 9, dpi = 300)


# =============================================================================
# (2). Estimate an MA model
# =============================================================================

# 2.1 First difference of inflation 
# dy_t = y_t - y_{t-1}: 

sample$dy_t <- sample$DPCCRV1Q225SBEA - dplyr::lag(sample$DPCCRV1Q225SBEA)

# 2.2 Identification: ACF and PACF
# For an MA(q) process the ACF cuts off after lag q and the PACF decays.
#Just for me to check if the MA(1) is actually justified.
par(mfrow = c(1, 2))
acf(sample$dy_t,  lag.max = 20, main = "ACF",  na.action = na.pass)
pacf(sample$dy_t, lag.max = 20, main = "PACF", na.action = na.pass)
par(mfrow = c(1, 1))

# 3.3 Estimate MA(1) by maximum likelihood 
# order = c(p, d, q) = c(0, 0, 1): no AR terms, no extra differencing, 1 MA lag.
# Model: dy_t = mu + e_t + theta_1 * e_{t-1}
# Note that R uses dy_t = mu + e_t + theta_1 * e_{t-1}, so in the docuyment, I will flip the sign
# I add include.mean = False, to not apply the estimation to dy_t - mu but dy_t as in Part 1. But it does not also affect the estimated value of theta
ma1 <- arima(sample$dy_t, order = c(0, 0, 1), include.mean = FALSE, method = "ML")
ma1   # theta_1 (ma1), intercept (= mu), sigma^2, log-likelihood, AIC

# 3.4 Inference ---------------------------------------------------------------
# Standard errors, z-statistics and p-values for theta_1 and mu.
coeftest(ma1)

# 3.5 Values of sigma_eta and sigma_nu (in percentage points of annualized quarterly inflation)

theta_in_document = - ma1[["coef"]][["ma1"]]
sigma_a   <- sqrt(ma1[["sigma2"]])
sigma_eta <- sqrt((theta_in_document)*sigma_a^2)
sigma_nu  <- sqrt((1+theta_in_document^2)*sigma_a^2 
                  -2*sigma_eta^2 )

cat(sprintf("Theta is : %.4f\n",     theta_in_document))
cat(sprintf("Sigma_a is : %.4f\n",   sigma_a))
cat(sprintf("Sigma_eta is : %.4f\n", sigma_eta))
cat(sprintf("Sigma_nu is : %.4f\n",  sigma_nu))

# Some extra interpreatation
#1. Signal-to-noise ratio
#sigma^2 = sigma_nu^2 / sigma_eta^2 = 0.363 / 0.2585 = 1.40. 
# Permanent (trend) shocks are somewhat larger than transitory ones, 
# so a lot of each quarter's inflation surprise persists.

# 2. Value of theta: We have theta = 0.3247, closer to 0 (RW) than to 1 (Inflation as Pure Noise around a mean)
# Because if theta were 0, there would be no transitory term (as the variance would be 0): Every shock would be permanent.
# Theta = 1 would deliver the opposite: the permanent variance would be zero and every sjhock would be transitory
# Theta = 0.32 neither variance is zero, so the conclusion is that permanent shocks dominate, not that inflation is a random walk.

Box.test(residuals(ma1), lag = 12, type = "Ljung-Box", fitdf = 1)

# p-value>0.05, so we can assume no leftover correlation: can't reject white noise, so the MA(1) is adequate.

# (3). Using the estimated parameters from (2), compute the filtered (τt|t) and smoothed (τt|T ) estimates of τ . 
# What do these say about the timing of the Covid inflation surge?

y     <- sample$DPCCRV1Q225SBEA          # level of inflation y_t
a_hat <- as.numeric(residuals(ma1))      # fitted innovations a_t (NA in first obs)
theta <- theta_in_document               # document sign convention

# 3.1 Filtered trend = Beveridge-Nelson trend: tau_{t|t} = y_t - theta * a_t (using results from Part1)
sample$tau_filt <- y - theta * a_hat

# 3.2 Smoothed trend: backward recursion tau_{t|T} = (1-theta) tau_{t|t} + theta tau_{t+1|T}
n <- length(y)
tau_smooth <- sample$tau_filt
for (t in (n - 1):1) {
  tau_smooth[t] <- (1 - theta) * sample$tau_filt[t] + theta * tau_smooth[t + 1]
}
sample$tau_smooth <- tau_smooth

# 3.3 Plot y_t with both trend estimates

filter_smoother_plot <- sample |>
  select(observation_date, y = DPCCRV1Q225SBEA, tau_filt, tau_smooth) |>
  pivot_longer(-observation_date) |>
  ggplot(aes(observation_date, value, colour = name, linewidth = name)) +
  geom_line(na.rm = TRUE) +
  scale_colour_manual(values = c(y = "grey70", tau_filt = "steelblue", tau_smooth = "firebrick"),
                      labels = c(y = expression(y[t]),
                                 tau_filt = expression(tau[t*"|"*t]),
                                 tau_smooth = expression(tau[t*"|"*T]))) +
  scale_linewidth_manual(values = c(y = 0.4, tau_filt = 0.6, tau_smooth = 0.8), guide = "none") +
  labs(title = "Core PCE inflation: filtered vs smoothed trend (MA(1) / local level)",
       x = NULL, y = "Percent", colour = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

filter_smoother_plot
filter_smoother_plot_covidZoom <-  filter_smoother_plot + 
  coord_cartesian(xlim = as.Date(c("2010-01-01", NA))) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")

filter_smoother_plot_covidZoom

# 3.4 Save both plots 
ggsave("filter_smoother_plot_covidZoom.png",           filter_smoother_plot_covidZoom,           width = 13, height = 4.5, dpi = 300)
ggsave("filter_smoother_plot.png", filter_smoother_plot, width = 13, height = 4.5, dpi = 300)

# lOOKING at exact numbers for interpretation:
sample |>
  filter(observation_date >= as.Date("2019-01-01"),
         observation_date <= as.Date("2023-12-31")) |>
  transmute(quarter = zoo::as.yearqtr(observation_date),
            y = round(DPCCRV1Q225SBEA, 2),
            tau_filt = round(tau_filt, 2),
            tau_smooth = round(tau_smooth, 2)) |>
  print(n = Inf)

# Volatility of quarterly changes: y_t vs filtered vs smoothed trend
sd_changes <- c(
  y          = sd(diff(sample$DPCCRV1Q225SBEA),  na.rm = TRUE),
  tau_filt   = sd(diff(sample$tau_filt),   na.rm = TRUE),
  tau_smooth = sd(diff(sample$tau_smooth), na.rm = TRUE)
)
round(sd_changes, 2)


# =============================================================================
# (4). Compute y^{BP}_{t|T} (using the result from Part 1, Question 3)
# -----------------------------------------------------------------------------
# y^{BP}_{t|T} = w(T-t) tau_{T|T} + sum_{j=t-T}^{t} w_{|j|} y_{t-j} + w(t) tau_{0|T}
# =============================================================================
Tn <- n - 1                                   # T

# 4.1 End-point trends
# tau_{T|T}: filter and smoother coincide at T
tau_TT <- sample$tau_smooth[n]

tau_0T <- (1 - theta) * y[1] + theta * sample$tau_smooth[2]

# 4.2 Low-pass weights
# w_0 = 1/16, w_j = sin(pi j / 16) / (pi j)
j <- 1:Tn
w <- c(1/16, sin(pi * j / 16) / (pi * j))     # w_0, w_1, ..., w_T

# 4.3 Tail sums: w(s) = sum_{j>s} w_j = (1 - w_0)/2 - sum_{j=1}^{s} w_j
w_tail <- (1 - w[1]) / 2 - c(0, cumsum(w[-1]))   # w(0), ..., w(T)

# 4.4 In-sample part: row t, column k holds w_{|t-k|}, so (Wmat %*% y)_t = sum_k w_{|t-k|} y_k
t_idx <- 0:Tn
Wmat  <- matrix(w[abs(outer(t_idx, t_idx, "-")) + 1], n, n)

# 4.5 y^{BP}_{t|T} for t = 0, ..., T
sample$y_bp <- w_tail[Tn - t_idx + 1] * tau_TT +   # w(T-t) tau_{T|T}
  as.numeric(Wmat %*% y) +            # observed y_0, ..., y_T
  w_tail[t_idx + 1] * tau_0T          # w(t)   tau_{0|T}

# 4.6 Plot y_t, the smoothed trend and y^{BP}_{t|T}
bp_plot <- sample |>
  select(observation_date, y = DPCCRV1Q225SBEA, tau_smooth, y_bp) |>
  pivot_longer(-observation_date) |>
  mutate(name = factor(name, levels = c("y", "tau_smooth", "y_bp"))) |>
  ggplot(aes(observation_date, value, colour = name, linewidth = name)) +
  geom_line(na.rm = TRUE) +
  scale_colour_manual(values = c(y = "grey70", tau_smooth = "firebrick", y_bp = "black"),
                      labels = c(y = expression(y[t]),
                                 tau_smooth = expression(tau[t*"|"*T]),
                                 y_bp = expression(y[t*"|"*T]^{BP}))) +
  scale_linewidth_manual(values = c(y = 0.4, tau_smooth = 0.6, y_bp = 0.8), guide = "none") +
  labs(title = "Core PCE inflation: smoothed trend vs low-pass component (periods > 8 years)",
       x = NULL, y = "Percent", colour = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

bp_plot

bp_plot_covidZoom <- bp_plot +
  coord_cartesian(xlim = as.Date(c("2010-01-01", NA))) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")

bp_plot_covidZoom

# 4.7 Save both plots
ggsave("bp_plot.png",           bp_plot,           width = 13, height = 4.5, dpi = 300)
ggsave("bp_plot_covidZoom.png", bp_plot_covidZoom, width = 13, height = 4.5, dpi = 300)

# 4.8 Numbers around COVID for interpretation
sample |>
  filter(observation_date >= as.Date("2019-01-01"),
         observation_date <= as.Date("2023-12-31")) |>
  transmute(quarter    = zoo::as.yearqtr(observation_date),
            y          = round(DPCCRV1Q225SBEA, 2),
            tau_smooth = round(tau_smooth, 2),
            y_bp       = round(y_bp, 2)) |>
  print(n = Inf)

# deviation quarterly changes of the BP
sd(diff(sample$y_bp))


# =============================================================================
# Estimate the MA(1) model for the first difference of inflation over different sample periods. Does it
# appear that the coefficients have changed? 
# =============================================================================

# Window1: 1960-1970
# Window2: 1970-1985
# Window3: 1986-to 2000
# Window: 2001 to 2020

# Estimate the MA(1) on dy_t restricted to [start, end] and return the key numbers
estimate_window <- function(start, end) {
  sel <- sample$observation_date >= as.Date(start) &
    sample$observation_date <= as.Date(end)
  m   <- arima(sample$dy_t[sel], order = c(0, 0, 1),
               include.mean = FALSE, method = "ML")
  theta <- -m$coef[["ma1"]]                 # document sign convention
  s2a   <- m$sigma2
  lb    <- Box.test(residuals(m), lag = 12, type = "Ljung-Box", fitdf = 1)
  data.frame(
    window     = paste(format(as.Date(start), "%Y"), format(as.Date(end), "%Y"), sep = "-"),
    n          = m$nobs,
    theta      = theta,
    se         = sqrt(m$var.coef[1, 1]),
    sigma2_a   = s2a,
    sigma2_eta = theta * s2a,
    sigma2_nu  = (1 - theta)^2 * s2a,
    snr        = (1 - theta)^2 / theta,     # sigma2_nu / sigma2_eta
    K          = 1 - theta,                 # steady-state gain
    LB_p       = lb$p.value
  )
}

windows <- list(
  c("1960-01-01", "1969-12-31"),   # Window 1
  c("1970-01-01", "1985-12-31"),   # Window 2
  c("1986-01-01", "2000-12-31"),   # Window 3
  c("2001-01-01", "2020-12-31")    # Window 4
)

ma_windows <- do.call(rbind, lapply(windows, function(w) estimate_window(w[1], w[2])))
print(ma_windows, digits = 3, row.names = FALSE)

# =============================================================================
# (6). Testing for breaks in theta
# -----------------------------------------------------------------------------
# Model: dy_t = a_t - theta_t a_{t-1}, with theta_t (and sigma2_a) constant within
# each regime. H0: theta equal across regimes; sigma2_a is left free in each regime.
# Requires: sample (with dy_t), windows, ma_windows from section (5).
# =============================================================================

# 6.0 Log-likelihood of an MA(1) on a vector x.
#     theta = NULL -> theta estimated; otherwise theta fixed (document sign convention).
#     sigma2_a is always concentrated out, i.e. free in each segment.
ma1_loglik <- function(x, theta = NULL) {
  if (is.null(theta)) {
    arima(x, order = c(0, 0, 1), include.mean = FALSE, method = "ML")$loglik
  } else {
    arima(x, order = c(0, 0, 1), include.mean = FALSE, method = "ML",
          fixed = -theta, transform.pars = FALSE)$loglik          # R's ma1 = -theta
  }
}

# 6.0b LR test of a common theta across a list of segments (df = #segments - 1)
lr_common_theta <- function(segments) {
  ll_u <- sum(sapply(segments, ma1_loglik))                       # separate thetas
  opt  <- optimize(function(t) sum(sapply(segments, ma1_loglik, theta = t)),
                   interval = c(0, 0.99), maximum = TRUE)         # common theta
  LR   <- 2 * (ll_u - opt$objective)
  df   <- length(segments) - 1
  c(LR = LR, df = df, p = pchisq(LR, df, lower.tail = FALSE), theta_common = opt$maximum)
}

# 6.1 Known break dates (Chow-type): my four windows ---------------------------
dy_windows <- lapply(windows, function(w) {
  sel <- sample$observation_date >= as.Date(w[1]) & sample$observation_date <= as.Date(w[2])
  sample$dy_t[sel]
})

# (a) Likelihood ratio, H0: theta_1 = ... = theta_4, chi2(3)
lr_windows <- lr_common_theta(dy_windows)
round(lr_windows, 4)

# (b) Wald, H0: R theta = 0 with R = adjacent differences, chi2(3)
#     V is diagonal because estimates from non-overlapping windows are asymptotically independent
th <- ma_windows$theta
V  <- diag(ma_windows$se^2)
R  <- cbind(diag(length(th) - 1), 0) - cbind(0, diag(length(th) - 1))
W  <- as.numeric(t(R %*% th) %*% solve(R %*% V %*% t(R)) %*% (R %*% th))
c(W = W, df = nrow(R), p = pchisq(W, nrow(R), lower.tail = FALSE))

# (c) Pairwise z-tests
pairs <- combn(length(th), 2)
pairwise <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(k) {
  i <- pairs[1, k]; j <- pairs[2, k]
  z <- (th[i] - th[j]) / sqrt(ma_windows$se[i]^2 + ma_windows$se[j]^2)
  data.frame(pair = paste(ma_windows$window[i], "vs", ma_windows$window[j]),
             diff = th[i] - th[j], z = z, p = 2 * pnorm(-abs(z)))
}))
print(pairwise, digits = 3, row.names = FALSE)

# # 6.2 Single known break on the full sample (Chow-type at 2000Q4) --------------
# dy   <- sample$dy_t[!is.na(sample$dy_t)]
# d_dy <- sample$observation_date[!is.na(sample$dy_t)]
# 
# lr_at <- function(k) lr_common_theta(list(dy[1:k], dy[(k + 1):length(dy)]))["LR"]
# 
# k_2000 <- max(which(d_dy <= as.Date("2000-12-31")))
# lr_2000 <- lr_at(k_2000)
# c(LR = unname(lr_2000), p = pchisq(lr_2000, 1, lower.tail = FALSE))
# 6.2 Full-sample series and LR for a single break after observation k ----------
# (shared by 6.2 and 6.3)
dy   <- sample$dy_t[!is.na(sample$dy_t)]                     # dy_1, ..., dy_T
d_dy <- sample$observation_date[!is.na(sample$dy_t)]         # their dates
n_dy <- length(dy)

lr_at <- function(k) lr_common_theta(list(dy[1:k], dy[(k + 1):n_dy]))[["LR"]]

# Chow-type test at a single known date (2000Q4), for comparison
k_2000  <- max(which(d_dy <= as.Date("2000-12-31")))
lr_2000 <- lr_at(k_2000)
c(LR = lr_2000, p = pchisq(lr_2000, 1, lower.tail = FALSE))

# 6.3 Unknown break date: QLR / sup-LR (Andrews 1993), 15% trimming -------------
k_grid <- floor(0.15 * n_dy):ceiling(0.85 * n_dy)
lr_seq <- sapply(k_grid, lr_at)

qlr       <- max(lr_seq)
k_star    <- k_grid[which.max(lr_seq)]
date_star <- d_dy[k_star]
cat(sprintf("QLR = %.2f, break after %s\n", qlr, zoo::as.yearqtr(date_star)))
# Andrews (1993) critical values, 1 restriction, 15% trimming: 10%: 7.12, 5%: 8.68, 1%: 12.16

# theta before / after the estimated break
theta_split <- rbind(
  before = estimate_window(min(d_dy), date_star)[, c("n", "theta", "se", "sigma2_a")],
  after  = estimate_window(d_dy[k_star + 1], max(d_dy))[, c("n", "theta", "se", "sigma2_a")]
)
print(theta_split, digits = 3)

# 6.4 Plot the LR sequence with the Andrews critical values
qlr_df <- data.frame(date = d_dy[k_grid], LR = lr_seq)
qlr_plot <- ggplot(qlr_df, aes(date, LR)) +
  geom_line(colour = "steelblue", linewidth = 0.7) +
  geom_hline(yintercept = 8.68,  linetype = "dashed", colour = "firebrick") +
  geom_hline(yintercept = 12.16, linetype = "dotted", colour = "firebrick") +
  annotate("point", x = date_star, y = qlr, colour = "firebrick", size = 2) +
  annotate("text",  x = date_star, y = qlr,
           label = sprintf("QLR = %.1f (%s)", qlr, zoo::as.yearqtr(date_star)),
           hjust = -0.1, vjust = -0.4, size = 3.2, colour = "firebrick") +
  annotate("text", x = min(qlr_df$date), y = 8.68,  label = "5% (Andrews)", hjust = 0, vjust = -0.4, size = 3.2) +
  annotate("text", x = min(qlr_df$date), y = 12.16, label = "1% (Andrews)", hjust = 0, vjust = -0.4, size = 3.2) +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +   # room for the label
  labs(x = NULL, y = "LR") +
  theme_minimal(base_size = 12)
qlr_plot
ggsave("qlr_plot.png", qlr_plot, width = 13, height = 4.5, dpi = 300)


# =============================================================================
# (7). What if: real-time trend tau_{t|t} with parameters from 2000-2019
# -----------------------------------------------------------------------------
# A researcher estimates the MA(1) on 2000-2019, fixes theta, and runs the
# steady-state filter tau_{t|t} = tau_{t-1|t-1} + (1 - theta)(y_t - tau_{t-1|t-1})
# through the COVID surge. Uses y, theta and sample$tau_filt (section 3)
# and estimate_window() (section 5).
# =============================================================================

# 7.1 Steady-state filter for a fixed theta, started at tau_{0|0} = y_0
filter_fixed_theta <- function(x, theta_fixed) {
  K   <- 1 - theta_fixed
  out <- numeric(length(x))
  out[1] <- x[1]
  for (s in 2:length(x)) out[s] <- out[s - 1] + K * (x[s] - out[s - 1])
  out
}

# 7.2 Parameters from 2000-2019
est_pre   <- estimate_window("2000-01-01", "2019-12-31")
theta_pre <- est_pre$theta
c(theta_2000_2019 = theta_pre, K_2000_2019 = 1 - theta_pre,
  theta_full = theta,          K_full = 1 - theta)

sample$tau_filt_pre <- filter_fixed_theta(y, theta_pre)
# full-sample comparison: sample$tau_filt from section 3

# 7.3 Numbers over the surge
sample |>
  filter(observation_date >= as.Date("2019-01-01"),
         observation_date <= as.Date("2024-12-31")) |>
  transmute(quarter  = zoo::as.yearqtr(observation_date),
            y        = DPCCRV1Q225SBEA,
            tau_pre  = round(tau_filt_pre, 2),
            tau_full = round(tau_filt, 2),
            gap      = round(tau_filt - tau_filt_pre, 2)) |>
  print(n = Inf)

# 7.4 Plot: y_t and the two real-time trends, 2015 onward
whatif_plot <- sample |>
  filter(observation_date >= as.Date("2015-01-01")) |>
  select(observation_date, y = DPCCRV1Q225SBEA, tau_filt, tau_filt_pre) |>
  pivot_longer(-observation_date) |>
  mutate(name = factor(name, levels = c("y", "tau_filt", "tau_filt_pre"))) |>
  ggplot(aes(observation_date, value, colour = name, linewidth = name)) +
  annotate("rect", xmin = covid_start, xmax = covid_end, ymin = -Inf, ymax = Inf,
           fill = "grey70", alpha = 0.35) +
  geom_hline(yintercept = 2, linetype = "dashed", colour = "grey40") +
  geom_line() +
  scale_colour_manual(values = c(y = "grey70", tau_filt = "steelblue", tau_filt_pre = "darkorange"),
                      labels = c(y = expression(y[t]),
                                 tau_filt     = expression(tau[t*"|"*t] ~ "(full-sample " * theta * ")"),
                                 tau_filt_pre = expression(tau[t*"|"*t] ~ "(2000-2019 " * theta * ")"))) +
  scale_linewidth_manual(values = c(y = 0.4, tau_filt = 0.7, tau_filt_pre = 0.9), guide = "none") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = NULL, y = "Percent", colour = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")
whatif_plot
ggsave("whatif_pre_covid_plot.png", whatif_plot, width = 13, height = 4.5, dpi = 300)

# 7.5 Other estimation windows: 2015-2025 and 1965-1985
alt_windows <- list(
  "2000-2019" = c("2000-01-01", "2019-12-31"),
  "2015-2025" = c("2015-01-01", "2025-12-31"),
  "1965-1985" = c("1965-01-01", "1985-12-31")
)
alt_est <- do.call(rbind, lapply(alt_windows, function(w) estimate_window(w[1], w[2])))
print(alt_est[, c("window", "n", "theta", "se", "K", "LB_p")], digits = 3, row.names = FALSE)

theta_2015 <- alt_est$theta[alt_est$window == "2015-2025"]
theta_1965 <- alt_est$theta[alt_est$window == "1965-1985"]
sample$tau_filt_2015 <- filter_fixed_theta(y, theta_2015)
sample$tau_filt_1965 <- filter_fixed_theta(y, theta_1965)

# 7.6 Numbers over the surge, all windows
sample |>
  filter(observation_date >= as.Date("2019-10-01"),
         observation_date <= as.Date("2024-12-31")) |>
  transmute(quarter  = zoo::as.yearqtr(observation_date),
            y        = DPCCRV1Q225SBEA,
            full     = round(tau_filt, 2),
            w2000_19 = round(tau_filt_pre, 2),
            w2015_25 = round(tau_filt_2015, 2),
            w1965_85 = round(tau_filt_1965, 2)) |>
  print(n = Inf)

# 7.7 Plot: y_t and the four real-time trends, 2015 onward
whatif_all_plot <- sample |>
  filter(observation_date >= as.Date("2015-01-01")) |>
  select(observation_date, y = DPCCRV1Q225SBEA, tau_filt, tau_filt_pre, tau_filt_2015, tau_filt_1965) |>
  pivot_longer(-observation_date) |>
  mutate(name = factor(name, levels = c("y", "tau_filt", "tau_filt_pre", "tau_filt_2015", "tau_filt_1965"))) |>
  ggplot(aes(observation_date, value, colour = name, linewidth = name)) +
  annotate("rect", xmin = covid_start, xmax = covid_end, ymin = -Inf, ymax = Inf,
           fill = "grey70", alpha = 0.35) +
  geom_hline(yintercept = 2, linetype = "dashed", colour = "grey40") +
  geom_line() +
  scale_colour_manual(
    values = c(y = "grey70", tau_filt = "steelblue", tau_filt_pre = "darkorange",
               tau_filt_2015 = "forestgreen", tau_filt_1965 = "purple"),
    labels = c(y             = expression(y[t]),
               tau_filt      = sprintf("Full sample (K = %.2f)", 1 - theta),
               tau_filt_pre  = sprintf("2000-2019 (K = %.2f)", 1 - theta_pre),
               tau_filt_2015 = sprintf("2015-2025 (K = %.2f)", 1 - theta_2015),
               tau_filt_1965 = sprintf("1965-1985 (K = %.2f)", 1 - theta_1965))) +
  scale_linewidth_manual(values = c(y = 0.4, tau_filt = 0.7, tau_filt_pre = 0.8,
                                    tau_filt_2015 = 0.8, tau_filt_1965 = 0.8), guide = "none") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = NULL, y = "Percent", colour = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")
whatif_all_plot
ggsave("whatif_all_windows_plot.png", whatif_all_plot, width = 13, height = 4.5, dpi = 300)

#Checking LB tests:
print(alt_est[, c("window", "n", "theta", "se", "K", "LB_p")], digits = 3, row.names = FALSE)
#Remark: wINDOW 2000-2019, the MA does not seem to fit well since LB test is rejected.

# window  n theta    se     K   LB_p
# 2000-2019 80 0.815 0.103 0.185 0.0344
# 2015-2025 44 0.376 0.120 0.624 0.5222
# 1965-1985 84 0.236 0.100 0.764 0.5778