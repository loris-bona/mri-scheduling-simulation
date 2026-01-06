# ============================================================
# PART 1 — Statistical analysis of ScanRecords.csv
# -----------------------------
# SECTION 0 — Setup
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(stringr)
  library(purrr)
})

set.seed(123)

DATA_PATH <- "ScanRecords.csv"

WORK_START <- 8
WORK_END   <- 17

# -----------------------------
# SECTION 0.1 — Helper functions

decimal_hours_to_hms <- function(x) {
  x <- as.numeric(x)
  hrs <- floor(x)
  mins <- round((x - hrs) * 60)
  # handle rounding edge case: 9.999 -> mins 60
  carry <- mins %/% 60
  hrs <- hrs + carry
  mins <- mins %% 60
  list(h = hrs, m = mins, s = 0)
}

clamp_to_working_hours <- function(dt, work_start = 8, work_end = 17) {
  d <- as.Date(dt)
  h <- hour(dt) + minute(dt)/60 + second(dt)/3600
  h2 <- pmin(pmax(h, work_start), work_end)
  hrs <- floor(h2)
  mins <- round((h2 - hrs) * 60)
  carry <- mins %/% 60
  hrs <- hrs + carry
  mins <- mins %% 60
  as.POSIXct(d) + hours(hrs) + minutes(mins)
}


is_work_time <- function(dt, work_start = 8, work_end = 17) {
  wd <- wday(dt, week_start = 1)
  h  <- hour(dt) + minute(dt)/60 + second(dt)/3600
  (wd <= 5) & (h >= work_start) & (h <= work_end)
}


next_work_time <- function(dt, work_start = 8, work_end = 17) {
  dt <- as.POSIXct(dt)
  wd <- wday(dt, week_start = 1)
  h  <- hour(dt) + minute(dt)/60 + second(dt)/3600
  
  if (wd >= 6) {
    days_to_mon <- 8 - wd
    d <- as.Date(dt) + days(days_to_mon)
    return(as.POSIXct(d) + hours(work_start))
  }
  
  if (h < work_start) {
    d <- as.Date(dt)
    return(as.POSIXct(d) + hours(work_start))
  }
  
  if (h > work_end) {
    d <- as.Date(dt) + days(1)
    return(next_work_time(as.POSIXct(d) + hours(work_start), work_start, work_end))
  }
  
  dt
}

biz_diff_hours <- function(start_dt, end_dt, work_start = 8, work_end = 17) {
  start_dt <- as.POSIXct(start_dt)
  end_dt   <- as.POSIXct(end_dt)
  
  if (is.na(start_dt) || is.na(end_dt)) return(NA_real_)
  if (end_dt < start_dt) return(NA_real_)
  
  s <- next_work_time(start_dt, work_start, work_end)
  e <- end_dt
  
  total <- 0
  
  while (s < e) {
    wd <- wday(s, week_start = 1)
    if (wd >= 6) {
      s <- next_work_time(s, work_start, work_end)
      next
    }
    
    day_date <- as.Date(s)
    day_start <- as.POSIXct(day_date) + hours(work_start)
    day_end   <- as.POSIXct(day_date) + hours(work_end)
    
    s <- max(s, day_start)
    
    if (s >= day_end) {
      s <- next_work_time(as.POSIXct(day_date + days(1)) + hours(work_start),
                          work_start, work_end)
      next
    }
    
    chunk_end <- min(day_end, e)
    total <- total + as.numeric(difftime(chunk_end, s, units = "hours"))
    
    s <- next_work_time(as.POSIXct(day_date + days(1)) + hours(work_start),
                        work_start, work_end)
  }
  
  total
}

# -----------------------------
# SECTION 0.2 — Load + basic cleaning

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

# sanity checks (optional)
message("Rows after cleaning: ", nrow(df))
print(df %>% count(PatientType))
message("Date range: ", min(as.Date(df$CallDateTime)), " to ", max(as.Date(df$CallDateTime)))
range(hour(df$CallDateTime) + minute(df$CallDateTime)/60) |> print()
table(wday(df$CallDateTime, week_start = 1) > 5) |> print()

# -----------------------------
# SECTION 1 — Daily arrivals (EDA for both types)

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

# Plot: daily arrivals over time
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

# Plot: histogram of arrivals/day
p_arrivals_hist <- ggplot(daily_arrivals, aes(x = arrivals)) +
  geom_histogram(bins = 15, color = "white") +
  facet_wrap(~ PatientType, scales = "free_y") +
  labs(
    title = "Distribution of daily arrivals",
    x = "Arrivals/day", y = "Number of days"
  ) +
  theme_minimal()

print(p_arrivals_hist)

# -----------------------------
# SECTION 1.1 — Scan durations (EDA for both types)

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

# Optional: Q-Q plots (especially relevant for Type 1 normality check)
p_qq <- df %>%
  ggplot(aes(sample = Duration)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ PatientType, scales = "free") +
  labs(title = "Q-Q plots of scan durations (Normal reference)") +
  theme_minimal()

print(p_qq)

# -----------------------------
# SECTION 1.2 — OPTIONAL diagnostic:
# Interarrival times (working-time only)
# Why optional?
#  - Part 1 does not require interarrival-time analysis.
#  - I include it only as a sanity check consistent with:
#       Type 1 Poisson arrivals/day => approximately exponential interarrivals.
# Requires helper function biz_diff_hours() defined earlier.

interarrival_df <- df %>%
  arrange(PatientType, CallDateTime) %>%
  group_by(PatientType) %>%
  mutate(
    prev_dt = lag(CallDateTime),
    interarrival_hours = purrr::map2_dbl(
      prev_dt,
      CallDateTime,
      ~ if (is.na(.x)) NA_real_ else biz_diff_hours(.x, .y, WORK_START, WORK_END)
    )
  ) %>%
  ungroup()

interarrival_non_na <- interarrival_df %>%
  filter(!is.na(interarrival_hours))

# Quick numeric summaries (optional)
interarrival_summary <- interarrival_non_na %>%
  group_by(PatientType) %>%
  summarise(
    n = n(),
    mean_hours = mean(interarrival_hours),
    median_hours = median(interarrival_hours),
    q90_hours = quantile(interarrival_hours, 0.90),
    q95_hours = quantile(interarrival_hours, 0.95),
    .groups = "drop"
  )

print(interarrival_summary)

p_interarrival <- ggplot(interarrival_non_na, aes(x = interarrival_hours)) +
  geom_histogram(bins = 30, color = "white") +
  facet_wrap(~ PatientType, scales = "free_y") +
  labs(
    title = "OPTIONAL: Interarrival times (hours, working-time only)",
    subtitle = "Sanity check (Type 1 Poisson/day implies exponential-like interarrivals)",
    x = "Interarrival time (hours)",
    y = "Count"
  ) +
  theme_minimal()

print(p_interarrival)

# -----------------------------
# SECTION 2 — Type 1 arrivals: Poisson fit + parametric bootstrap CIs

# -----------------------------
# SECTION 2.0 — Extract Type 1 daily arrivals

type1_daily <- daily_arrivals %>%
  filter(PatientType == "Type 1") %>%
  arrange(Date)

stopifnot(nrow(type1_daily) > 0)

x <- type1_daily$arrivals
n_days <- length(x)

message("Type 1: number of workdays = ", n_days)
message("Type 1: total arrivals = ", sum(x))

# -----------------------------
# SECTION 2.1 — Point estimate for lambda (MLE)
lambda_hat <- mean(x)

message("Type 1: lambda_hat (mean arrivals/day) = ", round(lambda_hat, 3))

# -----------------------------
# SECTION 2.2 — Management-friendly quantities (point estimates)
busy_thresholds <- c(18, 20, 23)  # e.g., "more than 20 patients/day"
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

# -----------------------------
# SECTION 2.3 — Parametric bootstrap for uncertainty
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

# -----------------------------
# SECTION 2.4 — Build 95% bootstrap confidence intervals
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

# -----------------------------
# SECTION 2.5 — Better formatting (optional but recommended for report)

ci_tbl_pretty <- ci_tbl %>%
  mutate(
    estimate = if_else(str_detect(quantity, "^P\\("), 100 * estimate, estimate),
    ci_low   = if_else(str_detect(quantity, "^P\\("), 100 * ci_low, ci_low),
    ci_high  = if_else(str_detect(quantity, "^P\\("), 100 * ci_high, ci_high)
  ) %>%
  mutate(
    estimate = if_else(str_detect(quantity, "quantile"), round(estimate), estimate),
    ci_low   = if_else(str_detect(quantity, "quantile"), round(ci_low), ci_low),
    ci_high  = if_else(str_detect(quantity, "quantile"), round(ci_high), ci_high)
  ) %>%
  mutate(across(c(estimate, ci_low, ci_high), ~ round(.x, 3)))

print(ci_tbl_pretty)

# -----------------------------
# SECTION 2.6 — Visual check: fitted Poisson vs empirical histogram

emp_counts <- tibble(arrivals = x)

p_type1_poisson_fit <- ggplot(emp_counts, aes(x = arrivals)) +
  geom_histogram(binwidth = 1, boundary = -0.5, color = "white") +
  stat_function(
    fun = function(k) dpois(round(k), lambda_hat) * n_days,
    linewidth = 1
  ) +
  scale_x_continuous(breaks = seq(min(x), max(x), by = 1)) +
  labs(
    title = "Type 1 daily arrivals: empirical vs Poisson fit",
    subtitle = paste0("Fitted Poisson mean (lambda_hat) = ", round(lambda_hat, 2)),
    x = "Arrivals/day",
    y = "Number of days (empirical); Poisson curve scaled to days"
  ) +
  theme_minimal()

print(p_type1_poisson_fit)

# -----------------------------
# SECTION 2.7 — (Optional) Quick dispersion check (Poisson sanity)

dispersion_ratio <- var(x) / mean(x)
message("Type 1: Var/Mean (dispersion ratio) = ", round(dispersion_ratio, 3))

# -----------------------------
# OUTPUTS you will reuse later (Part 2 inputs)
type1_lambda_hat <- lambda_hat
type1_arrivals_ci <- ci_tbl_pretty



# -----------------------------
# SECTION 3 — Type 1 scan durations:

# -----------------------------
# SECTION 3.0 — Extract Type 1 durations
type1_dur <- df %>%
  filter(PatientType == "Type 1") %>%
  pull(Duration)

stopifnot(length(type1_dur) > 0)

n_scans <- length(type1_dur)
message("Type 1: number of scans = ", n_scans)

# -----------------------------
# SECTION 3.1 — Parameter estimates (MLEs)
mu_hat    <- mean(type1_dur)
sigma_hat <- sd(type1_dur)

message("Type 1: mu_hat (mean duration, hours) = ", round(mu_hat, 4))
message("Type 1: sigma_hat (sd duration, hours) = ", round(sigma_hat, 4))

# -----------------------------
# SECTION 3.2 — Management-friendly quantities (point estimates)
slot_lengths <- c(0.45, 0.50, 0.55, 0.60)  # ≈ 27, 30, 33, 36 minutes

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

# -----------------------------
# SECTION 3.3 — Parametric bootstrap for uncertainty

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

# -----------------------------
# SECTION 3.4 — 95% bootstrap confidence intervals
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

# -----------------------------
# SECTION 3.5 — Convert to minutes (report-ready table)
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
    )
  ) %>%
  mutate(across(c(estimate, ci_low, ci_high), ~ round(.x, 2)))

print(ci_tbl_minutes)

# -----------------------------
# SECTION 3.6 — Visual diagnostics

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

# -----------------------------
# SECTION 3.7 — Slot-length recommendation (interpretation helper)
slot_summary <- tibble(
  slot_length_min = slot_lengths * 60,
  prob_overrun = sapply(slot_lengths, function(L) 1 - pnorm(L, mu_hat, sigma_hat))
) %>%
  mutate(prob_overrun_pct = round(100 * prob_overrun, 1))

print(slot_summary)

# -----------------------------
# OUTPUTS you will reuse later (Part 2 inputs)
type1_mu_hat    <- mu_hat
type1_sigma_hat <- sigma_hat
type1_duration_ci <- ci_tbl_minutes
type1_slot_summary <- slot_summary


# -----------------------------
# SECTION 4 — Type 2 arrivals:
# -----------------------------
# SECTION 4.0 — Extract Type 2 daily arrivals

type2_daily <- daily_arrivals %>%
  filter(PatientType == "Type 2") %>%
  arrange(Date)

stopifnot(nrow(type2_daily) > 0)

x <- type2_daily$arrivals
n_days <- length(x)

message("Type 2: number of workdays = ", n_days)
message("Type 2: total arrivals = ", sum(x))

# -----------------------------
# SECTION 4.1 — Plug-in point estimates (distribution-free)
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

# -----------------------------
# SECTION 4.2 — Nonparametric bootstrap for uncertainty

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

# -----------------------------
# SECTION 4.3 — 95% bootstrap confidence intervals

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

# -----------------------------
# SECTION 4.4 — Visual diagnostics (distribution-free)

p_type2_arrivals_hist <- ggplot(tibble(arrivals = x), aes(x = arrivals)) +
  geom_histogram(binwidth = 1, boundary = -0.5, color = "white") +
  scale_x_continuous(breaks = seq(min(x), max(x), by = 1)) +
  labs(
    title = "Type 2 daily arrivals (empirical distribution)",
    x = "Arrivals/day",
    y = "Number of days"
  ) +
  theme_minimal()

print(p_type2_arrivals_hist)

# Empirical CDF (very intuitive for management)
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

# -----------------------------
# SECTION 4.5 — Interpretation helper

summary_tbl <- tibble(
  arrivals_day = sort(unique(x)),
  prob_exceed  = sapply(sort(unique(x)), function(k) mean(x > k))
)

print(summary_tbl)

# -----------------------------
# OUTPUTS you will reuse later (Part 2 inputs)
type2_mean_arrivals <- mean_hat
type2_arrivals_ci <- ci_tbl
type2_arrivals_ecdf <- summary_tbl


# -----------------------------
# SECTION 5 — Type 2 scan durations:

# -----------------------------
# SECTION 5.0 — Extract Type 2 durations

type2_dur <- df %>%
  filter(PatientType == "Type 2") %>%
  pull(Duration)

stopifnot(length(type2_dur) > 0)

n_scans <- length(type2_dur)
message("Type 2: number of scans = ", n_scans)

# -----------------------------
# SECTION 5.1 — Plug-in point estimates (distribution-free)

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

# -----------------------------
# SECTION 5.2 — Nonparametric bootstrap for uncertainty

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

# -----------------------------
# SECTION 5.3 — 95% bootstrap confidence intervals

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

# -----------------------------
# SECTION 5.4 — Convert to minutes (report-ready)

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
    )
  ) %>%
  mutate(across(c(estimate, ci_low, ci_high), ~ round(.x, 2)))

print(ci_tbl_minutes)

# -----------------------------
# SECTION 5.5 — Visual diagnostics (distribution-free)

p_type2_dur_hist <- ggplot(tibble(Duration = type2_dur), aes(x = Duration)) +
  geom_histogram(bins = 30, color = "white") +
  labs(
    title = "Type 2 scan durations (empirical distribution)",
    x = "Duration (hours)",
    y = "Number of scans"
  ) +
  theme_minimal()

print(p_type2_dur_hist)

p_type2_dur_ecdf <- ggplot(tibble(Duration = type2_dur), aes(x = Duration)) +
  stat_ecdf(geom = "step", linewidth = 1) +
  labs(
    title = "Type 2 scan durations: empirical CDF",
    x = "Duration (hours)",
    y = "P(Duration ≤ x)"
  ) +
  theme_minimal()

print(p_type2_dur_ecdf)

# -----------------------------
# SECTION 5.6 — Slot-length interpretation helper

slot_summary <- tibble(
  slot_length_min = slot_lengths * 60,
  prob_overrun = sapply(slot_lengths, function(L) mean(type2_dur > L))
) %>%
  mutate(prob_overrun_pct = round(100 * prob_overrun, 1))

print(slot_summary)

# -----------------------------
# OUTPUTS you will reuse later (Part 2 inputs)
type2_mean_duration <- mean_hat
type2_duration_ci <- ci_tbl_minutes
type2_slot_summary <- slot_summary



# ============================================================
# MONTE CARLO STUDY (Part 1 add-on)
# Robustness of Type 2 duration estimators + bootstrap CIs
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
})

set.seed(123)

# --- Inputs from your data (Type 2 durations in HOURS) ---
type2_dur <- df %>% filter(PatientType == "Type 2") %>% pull(Duration)
n <- length(type2_dur)

m_data <- mean(type2_dur)
v_data <- var(type2_dur)

# Quantities you care about (match Part 1)
q_levels <- c(0.90, 0.95)
slot_lengths <- c(0.8, 0.9)   # hours (48 and 54 minutes)

# Bootstrap settings
B_boot <- 400   # per Monte Carlo replicate (keep moderate for speed)
R <- 400        # Monte Carlo replicates (increase if you want)

# --- Helper: compute plug-in estimates on a sample y ---
plugin_estimates <- function(y) {
  c(
    mean = mean(y),
    q90  = as.numeric(quantile(y, 0.90)),
    q95  = as.numeric(quantile(y, 0.95)),
    P_gt_0.8 = mean(y > 0.8),
    P_gt_0.9 = mean(y > 0.9)
  )
}

# --- Helper: bootstrap percentile CI for plug-in quantities ---
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

# ============================================================
# Scenario parameter calibration (match mean/variance of data)
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
