# PART 1: Statistical analysis
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(stringr)
  library(purrr)
  library(tibble)
})

set.seed(123)

DATA_PATH <- "ScanRecords.csv"

raw <- read_csv(DATA_PATH, show_col_types = FALSE)

df <- raw %>%
  mutate(
    PatientType = factor(PatientType),
    Date        = as.Date(Date),
    Time        = as.numeric(Time),
    Duration    = as.numeric(Duration)
  ) %>%
  filter(!is.na(Date), !is.na(Time), !is.na(Duration), Duration > 0) %>%
  mutate(
    hrs  = floor(Time),
    mins = round((Time - hrs) * 60),
    hrs  = hrs + (mins %/% 60),
    mins = mins %% 60,
    CallDateTime = as.POSIXct(Date, tz = "UTC") + hours(hrs) + minutes(mins),
    Weekday = wday(CallDateTime, label = TRUE, week_start = 1)
  ) %>%
  select(-hrs, -mins)

# double check
message("Rows after cleaning: ", nrow(df))
print(df %>% count(PatientType))
message("Date range: ", min(as.Date(df$CallDateTime)), " to ", max(as.Date(df$CallDateTime)))
range(hour(df$CallDateTime) + minute(df$CallDateTime)/60) |> print()
table(wday(df$CallDateTime, week_start = 1) > 5) |> print()


# -------- SECTION 1 - Daily arrivals --------

daily_arrivals <- df %>%
  mutate(Date = as.Date(CallDateTime)) %>%
  count(PatientType, Date, name = "arrivals") %>%
  arrange(PatientType, Date)

all_workdays <- tibble(Date = seq(min(df$Date), max(df$Date), by = "day")) %>%
  filter(wday(Date, week_start = 1) <= 5)

daily_arrivals <- daily_arrivals %>%
  group_by(PatientType) %>%
  right_join(all_workdays, by = "Date") %>%
  mutate(arrivals = replace_na(arrivals, 0L)) %>%
  ungroup() %>%
  arrange(PatientType, Date)

arrivals_summary <- daily_arrivals %>%
  group_by(PatientType) %>%
  summarise(
    days = n(),
    mean_arrivals = mean(arrivals),
    sd_arrivals = sd(arrivals),
    q50 = quantile(arrivals, 0.50),
    q90 = quantile(arrivals, 0.90),
    q95 = quantile(arrivals, 0.95),
    min = min(arrivals),
    max = max(arrivals),
    .groups = "drop"
  )

print(arrivals_summary)

p_arrivals <- ggplot(daily_arrivals, aes(x = Date, y = arrivals)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.4) +
  facet_wrap(~ PatientType, scales = "free_y") +
  labs(
    title = "Daily arrivals by patient type",
    x = "Date", y = "Arrivals/day"
  ) +
  theme_minimal()

print(p_arrivals)

p_arrivals_hist <- ggplot(daily_arrivals, aes(x = arrivals)) +
  geom_histogram(bins = 15, color = "white") +
  facet_wrap(~ PatientType, scales = "free_y") +
  labs(
    title = "Distribution of daily arrivals",
    x = "Arrivals/day", y = "Number of days"
  ) +
  theme_minimal()

print(p_arrivals_hist)

# -------- SECTION 1.1 - Scan durations --------

dur_summary <- df %>%
  group_by(PatientType) %>%
  summarise(
    n = n(),
    mean_h = mean(Duration),
    sd_h = sd(Duration),
    q50_h = quantile(Duration, 0.50),
    q90_h = quantile(Duration, 0.90),
    q95_h = quantile(Duration, 0.95),
    min_h = min(Duration),
    max_h = max(Duration),
    .groups = "drop"
  ) %>%
  mutate(across(ends_with("_h"), ~ .x * 60, .names = "{.col}_min"))

print(dur_summary)

p_dur_hist <- ggplot(df, aes(x = Duration)) +
  geom_histogram(bins = 30, color = "white") +
  facet_wrap(~ PatientType, scales = "free_y") +
  labs(
    title = "Scan duration distribution by patient type",
    x = "Duration (hours)", y = "Number of scans"
  ) +
  theme_minimal()

print(p_dur_hist)


# -------- SECTION 2 - Type 1 arrivals: Poisson fit + parametric bootstrap CIs --------

# SECTION 2.0 — Extract Type 1 daily arrivals
type1_daily <- daily_arrivals %>%
  filter(PatientType == "Type 1") %>%
  arrange(Date)

stopifnot(nrow(type1_daily) > 0)

x <- type1_daily$arrivals
n_days <- length(x)

message("Type 1: number of workdays = ", n_days)
message("Type 1: total arrivals = ", sum(x))


lambda_hat <- mean(x)

message("Type 1: lambda_hat (mean arrivals/day) = ", round(lambda_hat, 3))


busy_thresholds <- c(18, 20, 23)
q_levels <- c(0.90, 0.95, 0.99)

type1_arrivals_point <- tibble(
  quantity = c(
    "Mean arrivals/day (lambda)",
    paste0("P(arrivals > ", busy_thresholds, ") under Poisson(lambda_hat)"),
    paste0(q_levels * 100, "% quantile of arrivals/day under Poisson(lambda_hat)")
  ),
  estimate = c(
    lambda_hat,
    sapply(busy_thresholds, function(k) 1 - ppois(k, lambda_hat)),
    sapply(q_levels, function(p) qpois(p, lambda_hat))
  )
)

print(type1_arrivals_point)

# -------- SECTION 2.1 - Parametric bootstrap for uncertainty --------

B <- 5000 

set.seed(123)

boot_mat <- replicate(B, {
  xb <- rpois(n_days, lambda_hat)
  lambda_b <- mean(xb)
  
  c(
    lambda_b,
    sapply(busy_thresholds, function(k) 1 - ppois(k, lambda_b)),
    sapply(q_levels, function(p) qpois(p, lambda_b))
  )
})

boot_mat <- t(boot_mat)

colnames(boot_mat) <- c(
  "lambda",
  paste0("P_gt_", busy_thresholds),
  paste0("q_", q_levels * 100)
)

boot_df <- as_tibble(boot_mat)

# -------- SECTION 2.2 - Build 95% bootstrap CI --------

ci_level <- 0.95
alpha <- 1 - ci_level

ci_tbl <- tibble(
  quantity = c(
    "Mean arrivals/day (lambda)",
    paste0("P(arrivals > ", busy_thresholds, ")"),
    paste0(q_levels * 100, "% quantile of arrivals/day")
  ),
  estimate = c(
    lambda_hat,
    sapply(busy_thresholds, function(k) 1 - ppois(k, lambda_hat)),
    sapply(q_levels, function(p) qpois(p, lambda_hat))
  ),
  ci_low = c(
    quantile(boot_df$lambda, alpha/2),
    sapply(busy_thresholds, function(k) quantile(boot_df[[paste0("P_gt_", k)]], alpha/2)),
    sapply(q_levels * 100, function(qp) quantile(boot_df[[paste0("q_", qp)]], alpha/2))
  ),
  ci_high = c(
    quantile(boot_df$lambda, 1 - alpha/2),
    sapply(busy_thresholds, function(k) quantile(boot_df[[paste0("P_gt_", k)]], 1 - alpha/2)),
    sapply(q_levels * 100, function(qp) quantile(boot_df[[paste0("q_", qp)]], 1 - alpha/2))
  )
)

print(ci_tbl)


# -------- SECTION 3 - Type 1 scan durations --------

# SECTION 3.0 — Extract Type 1 durations
type1_dur <- df %>%
  filter(PatientType == "Type 1") %>%
  pull(Duration)

stopifnot(length(type1_dur) > 0)

n_scans <- length(type1_dur)
message("Type 1: number of scans = ", n_scans)


mu_hat    <- mean(type1_dur)
sigma_hat <- sd(type1_dur)

message("Type 1: mu_hat (mean duration, hours) = ", round(mu_hat, 4))
message("Type 1: sigma_hat (sd duration, hours) = ", round(sigma_hat, 4))


slot_lengths <- c(0.45, 0.50, 0.55, 0.60) 

q_levels <- c(0.90, 0.95, 0.99)

type1_dur_point <- tibble(
  quantity = c(
    "Mean scan duration (hours)",
    paste0("P(duration > ", slot_lengths, " h)"),
    paste0(q_levels * 100, "% quantile of duration (hours)")
  ),
  estimate = c(
    mu_hat,
    sapply(slot_lengths, function(L) 1 - pnorm(L, mu_hat, sigma_hat)),
    sapply(q_levels, function(p) qnorm(p, mu_hat, sigma_hat))
  )
)

print(type1_dur_point)

# -------- SECTION 3.1 - Parametric bootstrap for uncertainty --------

B <- 5000
set.seed(123)

boot_mat <- replicate(B, {
  yb <- rnorm(n_scans, mu_hat, sigma_hat)
  
  mu_b    <- mean(yb)
  sigma_b <- sd(yb)
  
  c(
    mu_b,
    sigma_b,
    sapply(slot_lengths, function(L) 1 - pnorm(L, mu_b, sigma_b)),
    sapply(q_levels, function(p) qnorm(p, mu_b, sigma_b))
  )
})

boot_mat <- t(boot_mat)

colnames(boot_mat) <- c(
  "mu",
  "sigma",
  paste0("P_gt_", slot_lengths),
  paste0("q_", q_levels * 100)
)

boot_df <- as_tibble(boot_mat)

# -------- SECTION 3.2 - 95% bootstrap CI --------

ci_level <- 0.95
alpha <- 1 - ci_level

ci_tbl <- tibble(
  quantity = c(
    "Mean scan duration (hours)",
    "SD of scan duration (hours)",
    paste0("P(duration > ", slot_lengths, " h)"),
    paste0(q_levels * 100, "% quantile of duration (hours)")
  ),
  estimate = c(
    mu_hat,
    sigma_hat,
    sapply(slot_lengths, function(L) 1 - pnorm(L, mu_hat, sigma_hat)),
    sapply(q_levels, function(p) qnorm(p, mu_hat, sigma_hat))
  ),
  ci_low = c(
    quantile(boot_df$mu, alpha/2),
    quantile(boot_df$sigma, alpha/2),
    sapply(slot_lengths, function(L) quantile(boot_df[[paste0("P_gt_", L)]], alpha/2)),
    sapply(q_levels * 100, function(qp) quantile(boot_df[[paste0("q_", qp)]], alpha/2))
  ),
  ci_high = c(
    quantile(boot_df$mu, 1 - alpha/2),
    quantile(boot_df$sigma, 1 - alpha/2),
    sapply(slot_lengths, function(L) quantile(boot_df[[paste0("P_gt_", L)]], 1 - alpha/2)),
    sapply(q_levels * 100, function(qp) quantile(boot_df[[paste0("q_", qp)]], 1 - alpha/2))
  )
)

print(ci_tbl)

# -------- SECTION 3.3 — Convert to minutes --------

ci_tbl_minutes <- ci_tbl %>%
  mutate(
    estimate = if_else(
      str_detect(quantity, "hours"),
      estimate * 60,
      estimate
    ),
    ci_low = if_else(
      str_detect(quantity, "hours"),
      ci_low * 60,
      ci_low
    ),
    ci_high = if_else(
      str_detect(quantity, "hours"),
      ci_high * 60,
      ci_high
    ),
    quantity = str_replace(quantity, "\\(hours\\)", "(minutes)")
  ) %>%
  mutate(across(c(estimate, ci_low, ci_high), ~ round(.x, 2)))

print(ci_tbl_minutes)

# -------- SECTION 3.4 - Visual diagnostics --------

p_type1_dur_fit <- ggplot(tibble(Duration = type1_dur), aes(x = Duration)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, color = "white") +
  stat_function(
    fun = dnorm,
    args = list(mean = mu_hat, sd = sigma_hat),
    linewidth = 1
  ) +
  labs(
    title = "Type 1 scan durations: empirical vs Normal fit",
    subtitle = paste0(
      "mu_hat = ", round(mu_hat * 60, 1), " min, ",
      "sigma_hat = ", round(sigma_hat * 60, 1), " min"
    ),
    x = "Duration (hours)",
    y = "Density"
  ) +
  theme_minimal()

print(p_type1_dur_fit)

# Q–Q plot
p_type1_qq <- ggplot(tibble(Duration = type1_dur), aes(sample = Duration)) +
  stat_qq() +
  stat_qq_line() +
  labs(
    title = "Type 1 scan durations: Q–Q plot (Normal reference)",
    x = "Theoretical quantiles",
    y = "Sample quantiles"
  ) +
  theme_minimal()

print(p_type1_qq)

# -------- SECTION 3.5 - Slot-length recommendation --------

slot_summary <- tibble(
  slot_length_min = slot_lengths * 60,
  prob_overrun = sapply(slot_lengths, function(L) 1 - pnorm(L, mu_hat, sigma_hat))
) %>%
  mutate(prob_overrun_pct = round(100 * prob_overrun, 1))

print(slot_summary)


# -------- SECTION 4 - Type 2 arrivals --------

# SECTION 4.0 - Extract Type 2 daily arrivals

type2_daily <- daily_arrivals %>%
  filter(PatientType == "Type 2") %>%
  arrange(Date)

stopifnot(nrow(type2_daily) > 0)

x <- type2_daily$arrivals
n_days <- length(x)

message("Type 2: number of workdays = ", n_days)
message("Type 2: total arrivals = ", sum(x))

# -------- SECTION 4.1 - Plug-in point estimates (distribution-free) --------

mean_hat <- mean(x)

busy_thresholds <- c(9, 10, 11)

q_levels <- c(0.90, 0.95)

type2_arrivals_point <- tibble(
  quantity = c(
    "Mean arrivals/day",
    paste0("P(arrivals > ", busy_thresholds, ")"),
    paste0(q_levels * 100, "% quantile of arrivals/day")
  ),
  estimate = c(
    mean_hat,
    sapply(busy_thresholds, function(k) mean(x > k)),
    sapply(q_levels, function(p) as.numeric(quantile(x, p)))
  )
)

print(type2_arrivals_point)

# -------- SECTION 4.2 - Nonparametric bootstrap for uncertainty --------

B <- 5000
set.seed(123)

boot_mat <- replicate(B, {
  xb <- sample(x, size = n_days, replace = TRUE)
  
  c(
    mean(xb),
    sapply(busy_thresholds, function(k) mean(xb > k)),
    sapply(q_levels, function(p) as.numeric(quantile(xb, p)))
  )
})

boot_mat <- t(boot_mat)

colnames(boot_mat) <- c(
  "mean",
  paste0("P_gt_", busy_thresholds),
  paste0("q_", q_levels * 100)
)

boot_df <- as_tibble(boot_mat)

# -------- SECTION 4.3 - 95% bootstrap CI --------

ci_level <- 0.95
alpha <- 1 - ci_level

ci_tbl <- tibble(
  quantity = c(
    "Mean arrivals/day",
    paste0("P(arrivals > ", busy_thresholds, ")"),
    paste0(q_levels * 100, "% quantile of arrivals/day")
  ),
  estimate = c(
    mean_hat,
    sapply(busy_thresholds, function(k) mean(x > k)),
    sapply(q_levels, function(p) as.numeric(quantile(x, p)))
  ),
  ci_low = c(
    quantile(boot_df$mean, alpha/2),
    sapply(busy_thresholds, function(k) quantile(boot_df[[paste0("P_gt_", k)]], alpha/2)),
    sapply(q_levels * 100, function(qp) quantile(boot_df[[paste0("q_", qp)]], alpha/2))
  ),
  ci_high = c(
    quantile(boot_df$mean, 1 - alpha/2),
    sapply(busy_thresholds, function(k) quantile(boot_df[[paste0("P_gt_", k)]], 1 - alpha/2)),
    sapply(q_levels * 100, function(qp) quantile(boot_df[[paste0("q_", qp)]], 1 - alpha/2))
  )
)

print(ci_tbl)

# -------- SECTION 4.4 - Visual diagnostics (distribution-free) --------

# Empirical CDF
p_type2_ecdf <- ggplot(tibble(arrivals = x), aes(x = arrivals)) +
  stat_ecdf(geom = "step", linewidth = 1) +
  scale_x_continuous(breaks = seq(min(x), max(x), by = 1)) +
  labs(
    title = "Type 2 daily arrivals: empirical CDF",
    x = "Arrivals/day",
    y = "P(Arrivals ≤ x)"
  ) +
  theme_minimal()

print(p_type2_ecdf)


# -------- SECTION 5 - Type 2 scan durations --------

# SECTION 5.0 — Extract Type 2 durations

type2_dur <- df %>%
  filter(PatientType == "Type 2") %>%
  pull(Duration)

stopifnot(length(type2_dur) > 0)

n_scans <- length(type2_dur)
message("Type 2: number of scans = ", n_scans)

# -------- SECTION 5.1 - Plug-in point estimates (distribution-free) --------

mean_hat <- mean(type2_dur)

slot_lengths <- c(0.60, 0.70, 0.80, 0.90)  # ≈ 36, 42, 48, 54 minutes

q_levels <- c(0.90, 0.95)

type2_dur_point <- tibble(
  quantity = c(
    "Mean scan duration (hours)",
    paste0("P(duration > ", slot_lengths, " h)"),
    paste0(q_levels * 100, "% quantile of duration (hours)")
  ),
  estimate = c(
    mean_hat,
    sapply(slot_lengths, function(L) mean(type2_dur > L)),
    sapply(q_levels, function(p) as.numeric(quantile(type2_dur, p)))
  )
)

print(type2_dur_point)

# -------- SECTION 5.2 - Nonparametric bootstrap for uncertainty --------

B <- 5000
set.seed(123)

boot_mat <- replicate(B, {
  yb <- sample(type2_dur, size = n_scans, replace = TRUE)
  
  c(
    mean(yb),
    sapply(slot_lengths, function(L) mean(yb > L)),
    sapply(q_levels, function(p) as.numeric(quantile(yb, p)))
  )
})

boot_mat <- t(boot_mat)

colnames(boot_mat) <- c(
  "mean",
  paste0("P_gt_", slot_lengths),
  paste0("q_", q_levels * 100)
)

boot_df <- as_tibble(boot_mat)

# -------- SECTION 5.3 - 95% bootstrap CI --------

ci_level <- 0.95
alpha <- 1 - ci_level

ci_tbl <- tibble(
  quantity = c(
    "Mean scan duration (hours)",
    paste0("P(duration > ", slot_lengths, " h)"),
    paste0(q_levels * 100, "% quantile of duration (hours)")
  ),
  estimate = c(
    mean_hat,
    sapply(slot_lengths, function(L) mean(type2_dur > L)),
    sapply(q_levels, function(p) as.numeric(quantile(type2_dur, p)))
  ),
  ci_low = c(
    quantile(boot_df$mean, alpha/2),
    sapply(slot_lengths, function(L) quantile(boot_df[[paste0("P_gt_", L)]], alpha/2)),
    sapply(q_levels * 100, function(qp) quantile(boot_df[[paste0("q_", qp)]], alpha/2))
  ),
  ci_high = c(
    quantile(boot_df$mean, 1 - alpha/2),
    sapply(slot_lengths, function(L) quantile(boot_df[[paste0("P_gt_", L)]], 1 - alpha/2)),
    sapply(q_levels * 100, function(qp) quantile(boot_df[[paste0("q_", qp)]], 1 - alpha/2))
  )
)

print(ci_tbl)

# -------- SECTION 5.4 - Convert to minutes --------

ci_tbl_minutes <- ci_tbl %>%
  mutate(
    estimate = if_else(
      str_detect(quantity, "hours"),
      estimate * 60,
      estimate
    ),
    ci_low = if_else(
      str_detect(quantity, "hours"),
      ci_low * 60,
      ci_low
    ),
    ci_high = if_else(
      str_detect(quantity, "hours"),
      ci_high * 60,
      ci_high
    ),
    quantity = str_replace(quantity, "\\(hours\\)", "(minutes)")
  ) %>%
  mutate(across(c(estimate, ci_low, ci_high), ~ round(.x, 2)))

print(ci_tbl_minutes)

#  -------- SECTION 5.5 - Visual diagnostics (distribution-free) --------

p_type2_dur_ecdf <- ggplot(tibble(Duration = type2_dur), aes(x = Duration)) +
  stat_ecdf(geom = "step", linewidth = 1) +
  labs(
    title = "Type 2 scan durations: empirical CDF",
    x = "Duration (hours)",
    y = "P(Duration ≤ x)"
  ) +
  theme_minimal()

print(p_type2_dur_ecdf)

# -------- SECTION 5.6 - Slot-length --------

slot_summary <- tibble(
  slot_length_min = slot_lengths * 60,
  prob_overrun = sapply(slot_lengths, function(L) mean(type2_dur > L))
) %>%
  mutate(prob_overrun_pct = round(100 * prob_overrun, 1))

print(slot_summary)


# -------- MONTE CARLO --------

set.seed(123)

type2_dur <- df %>% filter(PatientType == "Type 2") %>% pull(Duration)
n <- length(type2_dur)

m_data <- mean(type2_dur)
v_data <- var(type2_dur)

B_boot <- 400   
R <- 400        

plugin_estimates <- function(y) {
  c(
    mean = mean(y),
    q90  = as.numeric(quantile(y, 0.90)),
    q95  = as.numeric(quantile(y, 0.95)),
    P_gt_0.8 = mean(y > 0.8),
    P_gt_0.9 = mean(y > 0.9)
  )
}

bootstrap_ci <- function(y, B = 400, alpha = 0.05) {
  n <- length(y)
  boot_stats <- replicate(B, {
    yb <- sample(y, size = n, replace = TRUE)
    plugin_estimates(yb)
  })
  boot_stats <- t(boot_stats)
  
  est <- plugin_estimates(y)
  low <- apply(boot_stats, 2, quantile, probs = alpha/2)
  high <- apply(boot_stats, 2, quantile, probs = 1 - alpha/2)
  
  list(est = est, low = low, high = high)
}

# 1) Lognormal
s2_logn <- log(1 + v_data / (m_data^2))
s_logn <- sqrt(s2_logn)
mu_logn <- log(m_data) - 0.5 * s2_logn

r_lognormal <- function(n) rlnorm(n, meanlog = mu_logn, sdlog = s_logn)

# 2) Gamma
shape_g <- (m_data^2) / v_data
scale_g <- v_data / m_data

r_gamma <- function(n) rgamma(n, shape = shape_g, scale = scale_g)

# 3) Weibull
weibull_mom_fit <- function(target_mean, target_var) {
  obj <- function(par) {
    k <- par[1]; lam <- par[2]
    if (k <= 0.1 || lam <= 0) return(1e9)
    m <- lam * gamma(1 + 1/k)
    v <- lam^2 * (gamma(1 + 2/k) - gamma(1 + 1/k)^2)
    (m - target_mean)^2 + (v - target_var)^2
  }
  fit <- optim(par = c(2, target_mean), fn = obj, method = "Nelder-Mead")
  list(shape = fit$par[1], scale = fit$par[2])
}

wb <- weibull_mom_fit(m_data, v_data)
r_weibull <- function(n) rweibull(n, shape = wb$shape, scale = wb$scale)

truth_logn <- function() c(
  mean = m_data,
  q90 = qlnorm(0.90, mu_logn, s_logn),
  q95 = qlnorm(0.95, mu_logn, s_logn),
  P_gt_0.8 = 1 - plnorm(0.8, mu_logn, s_logn),
  P_gt_0.9 = 1 - plnorm(0.9, mu_logn, s_logn)
)

truth_gamma <- function() c(
  mean = m_data,
  q90 = qgamma(0.90, shape = shape_g, scale = scale_g),
  q95 = qgamma(0.95, shape = shape_g, scale = scale_g),
  P_gt_0.8 = 1 - pgamma(0.8, shape = shape_g, scale = scale_g),
  P_gt_0.9 = 1 - pgamma(0.9, shape = shape_g, scale = scale_g)
)

truth_weib <- function() c(
  mean = wb$scale * gamma(1 + 1/wb$shape),
  q90 = qweibull(0.90, shape = wb$shape, scale = wb$scale),
  q95 = qweibull(0.95, shape = wb$shape, scale = wb$scale),
  P_gt_0.8 = 1 - pweibull(0.8, shape = wb$shape, scale = wb$scale),
  P_gt_0.9 = 1 - pweibull(0.9, shape = wb$shape, scale = wb$scale)
)

scenarios <- list(
  Lognormal = list(r = r_lognormal, truth = truth_logn),
  Gamma     = list(r = r_gamma,     truth = truth_gamma),
  Weibull   = list(r = r_weibull,   truth = truth_weib)
)

run_mc <- function(rfun, truth_fun, R = 400) {
  qnames <- names(plugin_estimates(type2_dur))
  
  truth <- truth_fun()
  truth <- truth[qnames]
  
  if (any(is.na(truth))) {
    missing <- qnames[is.na(truth)]
    stop("Truth function is missing these quantities: ", paste(missing, collapse = ", "))
  }
  
  results <- replicate(R, {
    y <- rfun(n)
    out <- bootstrap_ci(y, B = B_boot, alpha = 0.05)
    
    est  <- out$est[qnames]
    low  <- out$low[qnames]
    high <- out$high[qnames]
    
    c(
      est,
      setNames(low,  paste0("low.",  qnames)),
      setNames(high, paste0("high.", qnames))
    )
  })
  
  results <- as_tibble(t(results))
  
  metrics <- map_dfr(qnames, function(q) {
    est  <- results[[q]]
    low  <- results[[paste0("low.", q)]]
    high <- results[[paste0("high.", q)]]
    tru  <- as.numeric(truth[q])
    
    tibble(
      quantity = q,
      truth = tru,
      bias = mean(est - tru),
      rmse = sqrt(mean((est - tru)^2)),
      coverage_95 = mean(low <= tru & tru <= high)
    )
  })
  
  metrics
}

mc_summary <- imap_dfr(scenarios, function(sc, sc_name) {
  m <- run_mc(sc$r, sc$truth, R = R)
  m$scenario <- sc_name
  m
}) %>%
  select(scenario, quantity, truth, bias, rmse, coverage_95) %>%
  arrange(quantity, scenario)

print(mc_summary)