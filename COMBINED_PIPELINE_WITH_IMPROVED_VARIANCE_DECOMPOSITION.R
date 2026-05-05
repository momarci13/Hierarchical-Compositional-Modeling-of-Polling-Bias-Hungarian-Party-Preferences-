# ============================================================
# COMBINED PIPELINE (paste-and-run) -
#   1) Frequentist "old methodology": mblogit month FE (+ optional covars) + institute RE
#   2) Bayesian "new methodology": brms Dirichlet–Multinomial with spline time + method + institute RE
#   3) Plots for BOTH:
#        - Frequentist: house effects (robust), month-fitted shares (safe fallback)
#        - Bayesian: time trends (with N_eff in newdata), house effects, method effects, PPC (manual; pp_check not supported)
#
# Dataset: valasztasok_vegso.xlsx
# Folder:  C:/Users/molna/Documents/KUTATÁSOK/2_Választás_bayesi_minta/Új bayesi
# Output:  ./plots
# ============================================================
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(lubridate)
  library(stringr)
  library(tibble)
  library(rlang)
  library(mclogit)
  library(brms)
  library(ggplot2)
  library(tidyr)
  library(patchwork)
})

# ---- set your working directory (keep yours) ----
setwd("C:/Users/molna/Documents/KUTATÁSOK/2_Választás_bayesi_minta/Új bayesi")

# ---- output folder ----
dir.create("plots", showWarnings = FALSE)

theme_pub <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

# ============================================================
# 0) Shared helpers
# ============================================================
hun_date <- function(x){
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_replace_all(x, intToUtf8(160), " ")
  x <- str_replace_all(x, "\\s+", " ")
  x <- str_replace_all(x, "\\.$", "")
  
  m <- c("január"="01","február"="02","március"="03","április"="04","május"="05","június"="06",
         "július"="07","augusztus"="08","szeptember"="09","október"="10","november"="11","december"="12")
  for (nm in names(m)) x <- str_replace_all(x, paste0("\\b", nm, "\\b"), m[[nm]])
  
  x <- str_replace_all(x, "\\.", "")
  x <- str_replace_all(x, " ", "-")
  suppressWarnings(ymd(x))
}

clean_share <- function(x){
  x <- as.character(x)
  x <- str_replace_all(x, intToUtf8(160), " ")
  x <- str_trim(x)
  x[x %in% c("", "-", "–", "—", "NA", "N/A", "n/a", "null", "NULL")] <- NA_character_
  x <- str_replace_all(x, "%", "")
  x <- str_replace_all(x, "[^0-9,\\.\\-]", "")
  x <- str_replace_all(x, ",", ".")
  suppressWarnings(as.numeric(x))
}

clean_N <- function(x){
  x <- as.character(x)
  x <- str_replace_all(x, intToUtf8(160), " ")
  x <- str_replace_all(x, "\\s+", " ")
  x <- str_trim(x)
  x <- str_replace_all(x, "[^0-9]", "")
  suppressWarnings(as.numeric(x))
}

to_safe <- function(s){
  s %>%
    iconv(from = "", to = "ASCII//TRANSLIT") %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "")
}

safe_levels <- function(f){
  if (!is.factor(f)) return(f)
  lv <- levels(f)
  lv2 <- to_safe(lv)
  lv2[lv2 == "" | is.na(lv2)] <- "lvl"
  levels(f) <- make.unique(lv2, sep = "_")
  f
}

lump_rare_levels <- function(f, min_n = 5){
  if (!is.factor(f)) return(f)
  tab <- table(f, useNA = "no")
  rare <- names(tab)[tab < min_n]
  f <- as.character(f)
  f[f %in% rare] <- "other"
  factor(f)
}

# Deterministic largest remainder rounding -> counts (sums exactly to N)
largest_remainder_counts <- function(p, N){
  p <- as.numeric(p)
  p[p < 0 | is.na(p)] <- 0
  s <- sum(p)
  if (s <= 0) p <- rep(1/length(p), length(p)) else p <- p/s
  
  raw <- p * N
  base <- floor(raw)
  rem <- as.integer(N - sum(base))
  
  if (rem > 0) {
    frac <- raw - base
    ord <- order(frac, seq_along(frac), decreasing = TRUE)  # deterministic tie-break
    base[ord[seq_len(rem)]] <- base[ord[seq_len(rem)]] + 1L
  }
  
  # final safeguard
  diff <- as.integer(N - sum(base))
  if (diff != 0) base[which.max(base)] <- base[which.max(base)] + diff
  
  as.integer(base)
}

# ============================================================
# 1) Load + common cleaning (done ONCE)
# ============================================================
df <- read_excel("valasztasok_vegso_tisztitott.xlsx") %>%
  rename(institute = `Intézet`, date_end = `Felmérés vége`) %>%
  mutate(date = hun_date(date_end), institute = as.factor(institute)) %>%
  arrange(date)

# Party column renames to ASCII-safe
if ("Egyéb(MM+LMP+Jobbik+Momentum+Párbeszéd+MSZP+DK)" %in% names(df)) {
  df <- df %>% rename(Egyeb = `Egyéb(MM+LMP+Jobbik+Momentum+Párbeszéd+MSZP+DK)`)
}

mi_candidates <- c("Mi_Hazánk", "Mi Hazánk", "MiHazánk", "Mi_Hazank", "Mi Hazank", "MiHazank")
mi_found <- mi_candidates[mi_candidates %in% names(df)]
if (length(mi_found) == 1) {
  df <- df %>% rename(MiHazank = !!sym(mi_found))
} else if (!("MiHazank" %in% names(df))) {
  stop("Could not find a Mi Hazánk column. Check the exact column name in the Excel.")
}

party_cols <- c("Fidesz","Tisza","MiHazank","MKKP","Egyeb")
missing_parties <- setdiff(party_cols, names(df))
if (length(missing_parties) > 0) stop("Missing party columns: ", paste(missing_parties, collapse = ", "))

# shares -> numeric -> proportions -> row-normalize
df[party_cols] <- lapply(df[party_cols], clean_share)
df[party_cols] <- lapply(df[party_cols], function(x){ x[is.na(x)] <- 0; x })

if (median(unlist(df[party_cols]), na.rm = TRUE) > 1) {
  df[party_cols] <- lapply(df[party_cols], function(x) x / 100)
}

df <- df %>%
  mutate(row_sum = rowSums(across(all_of(party_cols)))) %>%
  mutate(across(all_of(party_cols),
                ~ ifelse(row_sum > 0, .x / row_sum, 1/length(party_cols)))) %>%
  select(-row_sum)

# N_eff
deff <- 1.7
df <- df %>%
  mutate(
    N_raw = clean_N(N),
    N_raw = ifelse(is.na(N_raw) | N_raw <= 0, 1000, N_raw),
    N_eff = pmax(50L, round(N_raw / deff)),
    N_eff = pmin(N_eff, 2000L),
    N_eff = as.integer(N_eff)
  )

# safe institute levels once
df$institute <- safe_levels(as.factor(df$institute))

# ============================================================
# 2) FREQUENTIST BLOCK (mblogit "old methodology")
# ============================================================
message("\n=============================\nFREQUENTIST (mblogit) FIT\n=============================")

df_f <- df %>% mutate(month = factor(format(date, "%Y-%m")))

# covariates for frequentist (method + scenarios), sanitized & lumped
covars_raw <- c("Módszer",
                "MSZP nem indul","Párbeszéd nem indul","Humanisták nem indul",
                "Második reformkor nem indul","LMP nem indul","Megoldás Mozgalom nem indul",
                "Humanisták megjelennek","Momentum nem indul","Minenki Magyaroszága nem indul","Feri off")
covars_raw <- covars_raw[covars_raw %in% names(df_f)]

covars_safe <- setNames(to_safe(covars_raw), covars_raw)
rename_map <- setNames(as.list(covars_raw), unname(covars_safe[covars_raw]))
if (length(rename_map) > 0) df_f <- df_f %>% rename(!!!rename_map)
covars_safe_names <- unname(covars_safe[covars_raw])

for (v in covars_safe_names) {
  if (v %in% names(df_f)) {
    if (is.character(df_f[[v]])) df_f[[v]] <- as.factor(df_f[[v]])
    df_f[[v]] <- safe_levels(df_f[[v]])
    df_f[[v]] <- lump_rare_levels(df_f[[v]], min_n = 5)
  }
}

# drop non-varying
covars_ok_f <- covars_safe_names[covars_safe_names %in% names(df_f)]
covars_ok_f <- covars_ok_f[sapply(covars_ok_f, function(v){
  x <- df_f[[v]]
  if (is.factor(x) || is.character(x)) dplyr::n_distinct(x, na.rm = TRUE) > 1
  else stats::var(as.numeric(x), na.rm = TRUE) > 0
})]

# counts from proportions
counts_f <- t(vapply(seq_len(nrow(df_f)), function(i){
  largest_remainder_counts(as.numeric(df_f[i, party_cols]), df_f$N_eff[i])
}, integer(length(party_cols))))
colnames(counts_f) <- party_cols

df_counts_f <- bind_cols(
  df_f %>% select(date, institute, month, N_eff, all_of(covars_ok_f)),
  as.data.frame(counts_f)
)

stopifnot(all(rowSums(df_counts_f[, party_cols, drop=FALSE]) == df_counts_f$N_eff))

# pseudocount (frequentist stabilization)
USE_PSEUDOCOUNT <- TRUE
PSEUDOCOUNT <- 1L
if (USE_PSEUDOCOUNT) {
  df_counts_f[, party_cols] <- lapply(df_counts_f[, party_cols], function(x) as.integer(x + PSEUDOCOUNT))
  df_counts_f$N_eff <- rowSums(df_counts_f[, party_cols, drop=FALSE])
}

fit_try <- function(use_random = TRUE, use_covars = TRUE){
  rhs2 <- "0 + month"
  if (use_covars && length(covars_ok_f) > 0) rhs2 <- paste(rhs2, paste(covars_ok_f, collapse = " + "), sep = " + ")
  f2 <- as.formula(paste0("cbind(", paste(party_cols, collapse = ", "), ") ~ ", rhs2))
  if (use_random) {
    mblogit(f2, random = ~ 1 | institute, catCov = "diagonal", data = df_counts_f)
  } else {
    mblogit(f2, catCov = "diagonal", data = df_counts_f)
  }
}

fit_freq <- tryCatch(fit_try(use_random = TRUE,  use_covars = TRUE),  error = function(e) NULL)
fit_label <- "A: month + covars + random institute"
if (is.null(fit_freq)) {
  fit_freq <- tryCatch(fit_try(use_random = FALSE, use_covars = TRUE), error = function(e) NULL)
  fit_label <- "B: month + covars (no random institute)"
}
if (is.null(fit_freq)) {
  fit_freq <- tryCatch(fit_try(use_random = TRUE,  use_covars = FALSE), error = function(e) NULL)
  fit_label <- "C: month only + random institute"
}
if (is.null(fit_freq)) {
  fit_freq <- tryCatch(fit_try(use_random = FALSE, use_covars = FALSE), error = function(e) NULL)
  fit_label <- "D: month only (no random institute)"
}
if (is.null(fit_freq)) stop("mblogit failed in all configurations.")

message("mblogit converged using: ", fit_label)
print(summary(fit_freq))

# ============================================================
# 2A) Frequentist plots
# ============================================================

# (F1) Frequentist house effects: robust extraction (no fragile party filtering)
if (grepl("random institute", fit_label)) {
  
  reF_obj <- ranef(fit_freq)$institute
  reF_mat <- tryCatch(as.matrix(reF_obj), error = function(e) NULL)
  if (is.null(reF_mat)) reF_mat <- tryCatch(as.matrix(unclass(reF_obj)), error = function(e) NULL)
  
  if (is.null(reF_mat) || nrow(reF_mat) == 0) {
    message("Could not coerce ranef(fit_freq)$institute to a matrix; skipping frequentist house-effects plot.")
  } else {
    
    message("ranef(fit_freq)$institute columns: ", paste(colnames(reF_mat), collapse = ", "))
    
    reF_long <- as.data.frame(reF_mat)
    reF_long$institute <- rownames(reF_long)
    
    reF_long <- tidyr::pivot_longer(
      reF_long,
      cols = -institute,
      names_to = "coef",
      values_to = "Estimate"
    )
    
    # party label from coefficient label (keep as-is if parsing fails)
    reF_long$party <- reF_long$coef
    reF_long$party <- sub(":\\(Intercept\\)$", "", reF_long$party)
    reF_long$party <- sub("\\(Intercept\\)$", "", reF_long$party)
    reF_long$party <- sub("~1$", "", reF_long$party)
    reF_long$party <- sub("^eq", "", reF_long$party)
    reF_long$party <- sub("\\..*$", "", reF_long$party)
    
    reF_long <- reF_long %>%
      group_by(party) %>%
      arrange(Estimate, .by_group = TRUE) %>%
      mutate(institute_ord = factor(institute, levels = unique(institute))) %>%
      ungroup()
    
    pF_house <- ggplot(reF_long, aes(x = Estimate, y = institute_ord)) +
      geom_vline(xintercept = 0, linetype = 2, alpha = 0.6) +
      geom_point(size = 1.8) +
      facet_wrap(~ party, scales = "free_y", ncol = 2) +
      labs(
        title = "Frequentist house effects (mblogit random institute; point estimates)",
        x = "Random intercept deviation (log-odds scale)",
        y = "Institute"
      ) +
      theme_pub
    
    print(pF_house)
    ggsave("plots/F_01_house_effects.png", pF_house, width = 10, height = 8, dpi = 300)
    ggsave("plots/F_01_house_effects.pdf", pF_house, width = 10, height = 8)
  }
  
} else {
  message("Frequentist model has no random institute; skipping frequentist house-effects plot.")
}

# (F2) Frequentist fitted month shares:
months <- levels(df_counts_f$month)
newdata_month <- data.frame(month = factor(months, levels = months), stringsAsFactors = FALSE)

for (v in covars_ok_f) {
  if (is.factor(df_counts_f[[v]])) newdata_month[[v]] <- factor(levels(df_counts_f[[v]])[1], levels = levels(df_counts_f[[v]]))
  else newdata_month[[v]] <- df_counts_f[[v]][1]
}
newdata_month$institute <- df_counts_f$institute[1]

p_hat <- tryCatch(predict(fit_freq, newdata = newdata_month, type = "prob"), error = function(e) NULL)

if (!is.null(p_hat)) {
  p_hat_df <- as.data.frame(p_hat)
  # try to keep only party columns if they exist
  keep_cols <- intersect(colnames(p_hat_df), party_cols)
  if (length(keep_cols) == 0) keep_cols <- colnames(p_hat_df)
  
  p_hat_df$month <- months
  p_hat_long <- pivot_longer(p_hat_df, cols = all_of(keep_cols), names_to = "party", values_to = "prob") %>%
    mutate(date = as.Date(paste0(month, "-15")))
  
  pF_time <- ggplot(p_hat_long, aes(x = date, y = prob)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ party, scales = "free_y", ncol = 2) +
    labs(title = "Frequentist fitted month effects (mblogit)", x = "Month", y = "Fitted vote share") +
    theme_pub
  
  print(pF_time)
  ggsave("plots/F_02_time_trends.png", pF_time, width = 10, height = 7, dpi = 300)
  ggsave("plots/F_02_time_trends.pdf", pF_time, width = 10, height = 7)
} else {
  
  # empirical monthly means from observed shares (weighted by N_eff)
  df_emp <- df_counts_f %>%
    mutate(across(all_of(party_cols), ~ .x / N_eff)) %>%
    group_by(month) %>%
    summarise(
      Nw = sum(N_eff),
      across(all_of(party_cols), ~ weighted.mean(.x, w = N_eff)),
      .groups = "drop"
    ) %>%
    mutate(date = as.Date(paste0(as.character(month), "-15")))
  
  emp_long <- pivot_longer(df_emp, cols = all_of(party_cols), names_to = "party", values_to = "share")
  
  pF_time <- ggplot(emp_long, aes(x = date, y = share)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ party, scales = "free_y", ncol = 2) +
    labs(title = "Empirical monthly mean shares (fallback; mblogit predict() failed)", x = "Month", y = "Weighted mean share") +
    theme_pub
  
  print(pF_time)
  ggsave("plots/F_02_time_trends_empirical.png", pF_time, width = 10, height = 7, dpi = 300)
  ggsave("plots/F_02_time_trends_empirical.pdf", pF_time, width = 10, height = 7)
}

# ============================================================
# 3) BAYESIAN BLOCK (brms Dirichlet–Multinomial)
# ============================================================
message("\n=============================\nBAYESIAN (brms DM) FIT\n=============================")

df_b <- df %>%
  mutate(
    t_day  = as.numeric(date - min(date, na.rm = TRUE)),
    t_week = t_day / 7
  )

# Keep ONLY method (modszer) for Bayesian
if ("Módszer" %in% names(df_b)) {
  df_b <- df_b %>% rename(modszer = `Módszer`)
} else if ("Modszer" %in% names(df_b)) {
  df_b <- df_b %>% rename(modszer = Modszer)
}
if ("modszer" %in% names(df_b)) {
  df_b$modszer <- as.factor(df_b$modszer)
  df_b$modszer <- safe_levels(df_b$modszer)
}

df_b$institute <- safe_levels(as.factor(df_b$institute))

# Deterministic counts (no pseudocount for Bayes)
counts_b <- t(vapply(seq_len(nrow(df_b)), function(i){
  largest_remainder_counts(as.numeric(df_b[i, party_cols]), df_b$N_eff[i])
}, integer(length(party_cols))))
colnames(counts_b) <- party_cols

base_cols_b <- c("date", "institute", "t_week", "N_eff")
if ("modszer" %in% names(df_b)) base_cols_b <- c(base_cols_b, "modszer")

df_counts_b <- bind_cols(
  df_b %>% select(all_of(base_cols_b)),
  as.data.frame(counts_b)
)

stopifnot(all(rowSums(df_counts_b[, party_cols, drop=FALSE]) == df_counts_b$N_eff))

k_spline <- 6
rhs <- c(paste0("s(t_week, k = ", k_spline, ")"))
if ("modszer" %in% names(df_counts_b)) rhs <- c(rhs, "modszer")
rhs <- c(rhs, "(1 | institute)")

form_bayes <- as.formula(paste0(
  "cbind(", paste(party_cols, collapse=", "), ") | trials(N_eff) ~ ",
  paste(rhs, collapse=" + ")
))

dp <- default_prior(form_bayes, data = df_counts_b, family = dirichlet_multinomial())
dp$prior[dp$class == "phi"] <- "gamma(2, 0.01)"
dp$prior[(dp$class == "b") & (dp$prior == "(flat)")] <- "normal(0, 0.7)"
dp$prior[dp$class == "Intercept"] <- "student_t(3, 0, 5)"
dp$prior[(dp$class == "sd") & !is.na(dp$group) & dp$group == "institute"] <- "exponential(3)"
dp$prior[dp$class == "sds"] <- "exponential(1)"
priors_dm <- as.brmsprior(dp)

fit_bayes_dm <- brm(
  form_bayes,
  data = df_counts_b,
  family = dirichlet_multinomial(),
  prior = priors_dm,
  chains = 4, cores = 4, iter = 3000, warmup = 1500,
  control = list(adapt_delta = 0.98, max_treedepth = 12)
)

print(summary(fit_bayes_dm))

# ============================================================
# 3A) Bayesian plots
# ============================================================

# (B1) Time trends with 95% credible bands on SHARE scale
if ("modszer" %in% names(df_counts_b)) {
  ref_modszer <- levels(df_counts_b$modszer)[1]
} else {
  ref_modszer <- NULL
}

t_grid <- seq(
  min(df_counts_b$t_week, na.rm = TRUE),
  max(df_counts_b$t_week, na.rm = TRUE),
  length.out = 200
)
# IMPORTANT: include N_eff in newdata because model uses trials(N_eff)
newdata_time <- data.frame(
  t_week = t_grid,
  institute = df_counts_b$institute[1],   # ignored when re_formula=NA
  N_eff = median(df_counts_b$N_eff)
)
if (!is.null(ref_modszer)) newdata_time$modszer <- factor(ref_modszer, levels = levels(df_counts_b$modszer))

ep <- posterior_epred(fit_bayes_dm, newdata = newdata_time, re_formula = NA)
# ep dims: draws x time x categories

summ_party <- function(mat){
  data.frame(
    mean = apply(mat, 2, mean),
    lo   = apply(mat, 2, quantile, probs = 0.025),
    hi   = apply(mat, 2, quantile, probs = 0.975)
  )
}

party_names <- dimnames(ep)[[3]]
plot_df <- do.call(rbind, lapply(seq_along(party_names), function(k){
  s <- summ_party(ep[,,k])
  s$t_week <- t_grid
  s$party <- party_names[k]
  s
}))

min_date <- min(df_counts_b$date, na.rm = TRUE)
plot_df <- plot_df %>% mutate(date = min_date + round(t_week * 7))

pB_time <- ggplot(plot_df, aes(x = date, y = mean)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ party, scales = "free_y", ncol = 2) +
  labs(title = "Bayesian time trends (brms DM; population-level)", x = "Date", y = "Vote share") +
  theme_pub

print(pB_time)
ggsave("plots/B_01_time_trends.png", pB_time, width = 10, height = 7, dpi = 300)
ggsave("plots/B_01_time_trends.pdf", pB_time, width = 10, height = 7)

# (B2) Bayesian house effects (ranef with 95% CrI)
reB <- ranef(fit_bayes_dm)$institute

reB_df <- lapply(dimnames(reB)[[3]], function(par){
  x <- reB[,,par]
  data.frame(
    institute = rownames(x),
    Estimate = x[,"Estimate"],
    Q2.5     = x[,"Q2.5"],
    Q97.5    = x[,"Q97.5"],
    param    = par,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

reB_df$party <- sub("^mu", "", reB_df$param)
reB_df$party <- sub("_Intercept$", "", reB_df$party)

reB_df <- reB_df %>%
  group_by(party) %>%
  arrange(Estimate, .by_group = TRUE) %>%
  mutate(institute_ord = factor(institute, levels = unique(institute))) %>%
  ungroup()

pB_house <- ggplot(reB_df, aes(x = Estimate, y = institute_ord)) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.6) +
  geom_errorbar(aes(xmin = Q2.5, xmax = Q97.5), orientation = "y", linewidth = 0.6) +
  geom_point(size = 1.8) +
  facet_wrap(~ party, scales = "free_y", ncol = 2) +
  labs(
    title = "Bayesian house effects (brms DM; 95% credible intervals)",
    x = "Log-odds deviation (vs Fidesz baseline)",
    y = "Institute"
  ) +
  theme_pub

print(pB_house)
ggsave("plots/B_02_house_effects.png", pB_house, width = 10, height = 8, dpi = 300)
ggsave("plots/B_02_house_effects.pdf", pB_house, width = 10, height = 8)

# (B3) Bayesian method effects (modszer) if present
if ("modszer" %in% names(df_counts_b)) {
  
  fe <- as.data.frame(fixef(fit_bayes_dm))
  fe$term <- rownames(fe)
  
  fe_m <- fe %>% filter(grepl("_modszer", term))
  if (nrow(fe_m) > 0) {
    fe_m$party <- sub("^mu", "", fe_m$term)
    fe_m$party <- sub("_modszer.*$", "", fe_m$party)
    fe_m$method <- sub("^.*_modszer", "", fe_m$term)
    
    pB_method <- ggplot(fe_m, aes(x = Estimate, y = reorder(method, Estimate))) +
      geom_vline(xintercept = 0, linetype = 2, alpha = 0.6) +
      geom_errorbar(aes(xmin = Q2.5, xmax = Q97.5), orientation = "y", linewidth = 0.6) +
      geom_point(size = 2) +
      facet_wrap(~ party, scales = "free_y", ncol = 2) +
      labs(
        title = "Bayesian method effects (brms DM; vs baseline method)",
        x = "Effect on log-odds ratio (vs Fidesz)",
        y = "Method contrast"
      ) +
      theme_pub
    
    print(pB_method)
    ggsave("plots/B_03_method_effects.png", pB_method, width = 10, height = 6, dpi = 300)
    ggsave("plots/B_03_method_effects.pdf", pB_method, width = 10, height = 6)
  } else {
    message("No modszer coefficients found in fixef(); skipping Bayesian method plot.")
  }
  
} else {
  message("No 'modszer' in Bayesian data; skipping Bayesian method plot.")
}

# (B4) Manual PPC for Dirichlet-multinomial (pp_check not implemented)
y_rep <- posterior_predict(fit_bayes_dm, ndraws = 200)  # draws x N x K
y_obs <- as.matrix(df_counts_b[, party_cols, drop = FALSE])

rep_shares <- sweep(y_rep, c(1,2), rowSums(y_rep), "/")  # draws x N x K
obs_shares <- y_obs / rowSums(y_obs)

ppc_df <- do.call(rbind, lapply(seq_along(party_cols), function(k){
  data.frame(
    party = party_cols[k],
    rep_mean = apply(rep_shares[,,k], 1, mean),
    obs_mean = mean(obs_shares[,k])
  )
}))

pB_ppc <- ggplot(ppc_df, aes(x = rep_mean)) +
  geom_histogram(bins = 30) +
  geom_vline(aes(xintercept = obs_mean), linewidth = 0.9) +
  facet_wrap(~ party, scales = "free", ncol = 2) +
  labs(
    title = "Posterior predictive check (means of party shares)",
    x = "Replicated mean share across polls",
    y = "Count of posterior draws"
  ) +
  theme_pub

print(pB_ppc)
ggsave("plots/B_04_ppc.png", pB_ppc, width = 10, height = 6, dpi = 300)
ggsave("plots/B_04_ppc.pdf", pB_ppc, width = 10, height = 6)

# ============================================================
# (B5) Shrinkage plot: raw (no-pooling) vs shrunk (partial-pooling) house effects
# ============================================================
# Observed proportions per poll
obs_b_shares <- df_counts_b %>%
  mutate(across(all_of(party_cols), ~ .x / N_eff))

# Grand mean proportion per party
grand_mean_b <- colMeans(obs_b_shares[, party_cols, drop = FALSE])

# Which parties appear in the house-effects table?
# (brms fits K-1 equations; the remaining party is the reference)
parties_in_reB <- unique(reB_df$party)
ref_p <- setdiff(party_cols, parties_in_reB)
if (length(ref_p) != 1) ref_p <- party_cols[length(party_cols)]   # fallback: last

# Per-institute mean proportions (simple unweighted mean across polls)
inst_means_b <- obs_b_shares %>%
  select(institute, all_of(party_cols)) %>%
  group_by(institute) %>%
  summarise(across(all_of(party_cols), mean), .groups = "drop")

# Reference-category grand mean and per-institute mean, used to form log-odds
grand_ref  <- grand_mean_b[[ref_p]]

raw_logodds_df <- inst_means_b %>%
  pivot_longer(-institute, names_to = "party", values_to = "inst_share") %>%
  filter(party != ref_p) %>%
  left_join(
    data.frame(party = setdiff(party_cols, ref_p),
               grand_k = grand_mean_b[setdiff(party_cols, ref_p)]),
    by = "party"
  ) %>%
  left_join(
    inst_means_b %>% transmute(institute, inst_ref = .data[[ref_p]]),
    by = "institute"
  ) %>%
  mutate(
    # raw log-odds of inst share vs grand-mean log-odds (deviation from overall)
    raw = log(inst_share / inst_ref) - log(grand_k / grand_ref)
  ) %>%
  select(institute, party, raw)

shrinkage_df <- raw_logodds_df %>%
  left_join(reB_df %>% select(institute, party, shrunk = Estimate),
            by = c("institute", "party")) %>%
  filter(is.finite(raw) & !is.na(shrunk))

pB_shrinkage <- ggplot(shrinkage_df, aes(x = raw, y = shrunk)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55",
              linewidth = 0.6) +                          # no-shrinkage diagonal
  geom_hline(yintercept = 0, linetype = 3, colour = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey40", linewidth = 0.5) +
  geom_segment(aes(xend = raw, yend = 0),
               alpha = 0.25, colour = "grey60", linewidth = 0.5) + # shrinkage arrow
  geom_point(size = 2.2, alpha = 0.85) +
  facet_wrap(~ party, scales = "free", ncol = 2) +
  labs(
    title    = "Hierarchical shrinkage of pollster house effects",
    subtitle = paste0(
      "Each point = one institute. Dashed line: no shrinkage (raw = shrunk).\n",
      "Vertical segments show pull toward zero (complete-pooling baseline)."
    ),
    x = "Raw log-odds deviation (no-pooling estimate)",
    y = "Posterior mean house effect (partial pooling)"
  ) +
  theme_pub

print(pB_shrinkage)
ggsave("plots/B_05_pollster_bias_shrinkage.png", pB_shrinkage,
       width = 10, height = 8, dpi = 300)
ggsave("plots/B_05_pollster_bias_shrinkage.pdf", pB_shrinkage,
       width = 10, height = 8)

# ============================================================
# (B6) Party-wise posterior SD of pollster house effects (sigma_u)
# ============================================================
# as_draws_df gives one column per sd_institute__ parameter
sigma_draws_raw <- as_draws_df(fit_bayes_dm) %>%
  select(starts_with("sd_institute__"))

if (ncol(sigma_draws_raw) > 0) {

  sigma_long <- sigma_draws_raw %>%
    pivot_longer(everything(), names_to = "param", values_to = "sd_draw")

  # Strip brms naming convention: sd_institute__muTisza_Intercept -> Tisza
  sigma_long$party <- sigma_long$param
  sigma_long$party <- sub("^sd_institute__mu", "",    sigma_long$party)
  sigma_long$party <- sub("_Intercept$",       "",    sigma_long$party)
  sigma_long$party <- sub("^sd_institute__",   "",    sigma_long$party)  # fallback

  sigma_summ <- sigma_long %>%
    group_by(party) %>%
    summarise(
      mean_sd = mean(sd_draw),
      lo       = quantile(sd_draw, 0.025),
      hi       = quantile(sd_draw, 0.975),
      .groups  = "drop"
    ) %>%
    arrange(desc(mean_sd))

  pB_sigma <- ggplot(sigma_summ,
                     aes(x = reorder(party, mean_sd), y = mean_sd)) +
    geom_col(fill = "#619CFF", alpha = 0.75, width = 0.6) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25, linewidth = 0.7) +
    coord_flip() +
    labs(
      title    = "Posterior SD of pollster house effects by party",
      subtitle = "Bar = posterior mean; whiskers = 95% credible interval",
      x        = "Party",
      y        = expression(hat(sigma)[u] ~ "(log-odds scale)")
    ) +
    theme_pub

  print(pB_sigma)
  ggsave("plots/B_06_sigma_house_effect_sd.png", pB_sigma,
         width = 8, height = 5, dpi = 300)
  ggsave("plots/B_06_sigma_house_effect_sd.pdf", pB_sigma,
         width = 8, height = 5)

} else {
  message("No sd_institute__ parameters found in posterior draws; skipping B_06.")
}

# ============================================================
# (B7) Posterior distribution of the concentration parameter phi
# ============================================================
phi_draws_df <- as_draws_df(fit_bayes_dm) %>%
  select(starts_with("phi"))
plot_data_dir <- file.path("plots", "plot_data")
dir.create(plot_data_dir, showWarnings = FALSE, recursive = TRUE)
if (ncol(phi_draws_df) > 0) {
  
  phi_vec  <- phi_draws_df[[1]]
  phi_mean <- mean(phi_vec)
  phi_lo   <- quantile(phi_vec, 0.025)
  phi_hi   <- quantile(phi_vec, 0.975)
  
  pB_phi <- ggplot(data.frame(phi = phi_vec), aes(x = phi)) +
    geom_histogram(bins = 60, fill = "#619CFF", alpha = 0.70, colour = NA) +
    geom_vline(xintercept = phi_mean,
               linewidth = 0.9, colour = "black") +
    geom_vline(xintercept = c(phi_lo, phi_hi),
               linetype = 2, linewidth = 0.7, colour = "grey35") +
    annotate("text", x = phi_mean, y = Inf,
             label = sprintf("Mean = %.1f", phi_mean),
             hjust = -0.12, vjust = 1.6, size = 4) +
    labs(
      title    = expression("Posterior distribution of concentration parameter " * phi),
      subtitle = sprintf("Posterior mean = %.1f;  95%% CrI [%.1f,\u2009%.1f]",
                         phi_mean, phi_lo, phi_hi),
      x        = expression(phi),
      y        = "Number of posterior draws"
    ) +
    theme_pub
  
  print(pB_phi)
  ggsave("plots/B_07_phi_posterior.png", pB_phi, width = 8, height = 5, dpi = 300)
  ggsave("plots/B_07_phi_posterior.pdf", pB_phi, width = 8, height = 5)
  
  # ---- Save B7 plot data ----
  # Full posterior draw vector (allows replotting with any bin width or quantile)
  write.csv(
    data.frame(phi = phi_vec),
    file.path(plot_data_dir, "B_07_phi_draws.csv"),
    row.names = FALSE
  )
  # Scalar summary (mean + 95% CrI) for quick reference
  write.csv(
    data.frame(
      phi_mean = phi_mean,
      phi_lo   = phi_lo,
      phi_hi   = phi_hi
    ),
    file.path(plot_data_dir, "B_07_phi_summary.csv"),
    row.names = FALSE
  )
  message("Saved: B_07_phi_draws.csv and B_07_phi_summary.csv")
  
} else {
  message("No phi parameter found in posterior draws; skipping B_07.")
}

# ============================================================
# 4) SAVE ALL PLOT DATA
# ============================================================
message("\n=============================\nSAVING PLOT DATA\n=============================")

plot_data_dir <- file.path("plots", "plot_data")
dir.create(plot_data_dir, showWarnings = FALSE, recursive = TRUE)

# ---- F1: Frequentist house effects ----
if (exists("reF_long")) {
  write.csv(reF_long, file.path(plot_data_dir, "F_01_house_effects.csv"), row.names = FALSE)
  message("Saved: F_01_house_effects.csv")
}

# ---- F2: Frequentist time trends ----
if (exists("p_hat_long")) {
  write.csv(p_hat_long, file.path(plot_data_dir, "F_02_time_trends_fitted.csv"), row.names = FALSE)
  message("Saved: F_02_time_trends_fitted.csv")
} else if (exists("emp_long")) {
  write.csv(emp_long, file.path(plot_data_dir, "F_02_time_trends_empirical.csv"), row.names = FALSE)
  message("Saved: F_02_time_trends_empirical.csv")
}

# ---- B1: Bayesian time trends ----
write.csv(plot_df, file.path(plot_data_dir, "B_01_time_trends.csv"), row.names = FALSE)
message("Saved: B_01_time_trends.csv")

# ---- B2: Bayesian house effects ----
write.csv(reB_df, file.path(plot_data_dir, "B_02_house_effects.csv"), row.names = FALSE)
message("Saved: B_02_house_effects.csv")

# ---- B3: Bayesian method effects ----
if (exists("fe_m") && nrow(fe_m) > 0) {
  write.csv(fe_m, file.path(plot_data_dir, "B_03_method_effects.csv"), row.names = FALSE)
  message("Saved: B_03_method_effects.csv")
}

# ---- B4: PPC data ----
write.csv(ppc_df, file.path(plot_data_dir, "B_04_ppc.csv"), row.names = FALSE)
message("Saved: B_04_ppc.csv")

ppc_rep_matrix <- do.call(cbind, lapply(seq_along(party_cols), function(k) {
  apply(rep_shares[,,k], 1, mean)
}))
colnames(ppc_rep_matrix) <- party_cols
saveRDS(ppc_rep_matrix, file.path(plot_data_dir, "B_04_ppc_rep_matrix.rds"))
message("Saved: B_04_ppc_rep_matrix.rds")

# ---- B5: Shrinkage plot data (partial-pooling estimates) ----
if (exists("shrinkage_df")) {
  write.csv(shrinkage_df, file.path(plot_data_dir, "B_05_shrinkage.csv"), row.names = FALSE)
  message("Saved: B_05_shrinkage.csv")
}

# ---- B5 supplement: No-pooling (raw) log-odds estimates underlying B_05 ----
if (exists("raw_logodds_df")) {
  write.csv(raw_logodds_df, file.path(plot_data_dir, "B_05_raw_logodds_nopooling.csv"), row.names = FALSE)
  message("Saved: B_05_raw_logodds_nopooling.csv")
}

# ---- B6: Posterior SD of house effects by party ----
if (exists("sigma_summ")) {
  write.csv(sigma_summ, file.path(plot_data_dir, "B_06_sigma_house_effect_sd.csv"), row.names = FALSE)
  message("Saved: B_06_sigma_house_effect_sd.csv")
}

# ---- Raw modelling inputs ----
write.csv(df_counts_f, file.path(plot_data_dir, "input_frequentist_counts.csv"), row.names = FALSE)
message("Saved: input_frequentist_counts.csv")

write.csv(df_counts_b, file.path(plot_data_dir, "input_bayesian_counts.csv"), row.names = FALSE)
message("Saved: input_bayesian_counts.csv")

# ---- Posterior epred arrays (full fidelity) ----
saveRDS(ep,           file.path(plot_data_dir, "B_epred_time_trends.rds"))
message("Saved: posterior epred RDS arrays")

message("\n✓ All plot data saved into: ", normalizePath(plot_data_dir))
message("\n✓ DONE. Plots saved into: ", normalizePath("plots"))
message("\n=============================")
message("SUMMARY")
message("=============================")
message("Frequentist: 2 plots")
message("Bayesian: 7 plots + 6 improved variance decomposition plots = 13 plots")
message("TOTAL: 15 publication-ready plots")
message("=============================\n")

# ============================================================
# 5) EXTENSION: PRIOR GRID SEARCH + MODEL COMPARISON
# ============================================================

message("\n=============================")
message("EXTENSION: PRIOR GRID SEARCH")
message("=============================")

# ---- ensure loo is available ----
if (!requireNamespace("loo", quietly = TRUE)) {
  stop("Package 'loo' is required for model comparison. Install it first.")
}
library(loo)

# ============================================================
# 5A) PRIOR GRID
# ============================================================

prior_grid <- expand.grid(
  phi_prior = c("gamma(0.01,0.01)", "gamma(2,0.01)", "gamma(1,1)"),
  beta_sd   = c(0.2, 0.35, 0.7),
  intercept_scale = c(2.5, 5),
  sd_institute = c("exponential(2)", "exponential(3)"),
  sd_spline = c("exponential(1)", "exponential(2)"),
  stringsAsFactors = FALSE
)

prior_grid$model_id <- paste0("M", seq_len(nrow(prior_grid)))

message("Total models to fit: ", nrow(prior_grid))

# ============================================================
# 5B) PRIOR BUILDER
# ============================================================

build_priors <- function(cfg, form, data) {
  dp <- default_prior(form, data = data, family = dirichlet_multinomial())
  
  dp$prior[dp$class == "phi"] <- cfg$phi_prior
  dp$prior[(dp$class == "b") & (dp$prior == "(flat)")] <- paste0("normal(0,", cfg$beta_sd, ")")
  dp$prior[dp$class == "Intercept"] <- paste0("student_t(3,0,", cfg$intercept_scale, ")")
  
  dp$prior[(dp$class == "sd") & dp$group == "institute"] <- cfg$sd_institute
  dp$prior[dp$class == "sds"] <- cfg$sd_spline
  
  as.brmsprior(dp)
}

# ============================================================
# 5C) MODEL LOOP
# ============================================================

results_list <- list()

for (i in seq_len(nrow(prior_grid))) {
  
  cfg <- prior_grid[i, ]
  message("\n--- [", cfg$model_id, "] fitting ---")
  
  priors_i <- build_priors(cfg, form_bayes, df_counts_b)
  
  fit_i <- tryCatch({
    brm(
      form_bayes,
      data = df_counts_b,
      family = dirichlet_multinomial(),
      prior = priors_i,
      chains = 2, cores = 2,
      iter = 2000, warmup = 1000,
      refresh = 0,
      control = list(adapt_delta = 0.95)
    )
  }, error = function(e) {
    message("FAILED: ", cfg$model_id)
    return(NULL)
  })
  
  if (is.null(fit_i)) next
  
  # ---- diagnostics ----
  summ <- summary(fit_i)
  
  rhat_max <- suppressWarnings(max(summ$fixed[,"Rhat"], na.rm = TRUE))
  ess_min  <- suppressWarnings(min(summ$fixed[,"Bulk_ESS"], na.rm = TRUE))
  
  loo_i <- tryCatch(loo(fit_i), error = function(e) NULL)
  elpd <- if (!is.null(loo_i)) loo_i$estimates["elpd_loo","Estimate"] else NA
  
  results_list[[cfg$model_id]] <- list(
    config = cfg,
    fit = fit_i,
    rhat_max = rhat_max,
    ess_min = ess_min,
    elpd = elpd
  )
}

# ============================================================
# 5D) BUILD COMPARISON TABLE
# ============================================================

if (length(results_list) == 0) {
  stop("All models failed. No comparison possible.")
}

results_df <- do.call(rbind, lapply(results_list, function(x){
  data.frame(
    model_id = x$config$model_id,
    phi_prior = x$config$phi_prior,
    beta_sd = x$config$beta_sd,
    intercept_scale = x$config$intercept_scale,
    sd_institute = x$config$sd_institute,
    sd_spline = x$config$sd_spline,
    rhat_max = x$rhat_max,
    ess_min = x$ess_min,
    elpd = x$elpd
  )
}))

results_df <- results_df %>% arrange(desc(elpd))

write.csv(results_df, "plots/C_01_model_comparison.csv", row.names = FALSE)

message("\nTop models by ELPD:")
print(head(results_df, 10))

# ============================================================
# 5E) DIAGNOSTIC PLOTS
# ============================================================

# ---- (1) ELPD ranking ----
p_elpd <- ggplot(results_df, aes(x = reorder(model_id, elpd), y = elpd)) +
  geom_col(fill = "#619CFF", alpha = 0.8) +
  coord_flip() +
  labs(
    title = "Model comparison (LOO-ELPD)",
    x = "Model",
    y = "ELPD"
  ) +
  theme_pub

ggsave("plots/C_02_elpd_comparison.png", p_elpd, width = 8, height = 6, dpi = 300)

# ---- (2) Rhat distribution ----
p_rhat <- ggplot(results_df, aes(x = rhat_max)) +
  geom_histogram(bins = 30, fill = "#A23B72", alpha = 0.7) +
  geom_vline(xintercept = 1.01, linetype = 2) +
  labs(
    title = "Max R-hat across models",
    x = "R-hat",
    y = "Count"
  ) +
  theme_pub

ggsave("plots/C_03_rhat_distribution.png", p_rhat, width = 7, height = 5, dpi = 300)

# ---- (3) Phi sensitivity ----
p_phi <- ggplot(results_df, aes(x = phi_prior, y = elpd)) +
  geom_boxplot(fill = "#F18F01", alpha = 0.7) +
  labs(
    title = "Sensitivity to phi prior",
    x = "Phi prior",
    y = "ELPD"
  ) +
  theme_pub

ggsave("plots/C_04_phi_sensitivity.png", p_phi, width = 7, height = 5, dpi = 300)

# ---- (4) Heatmap ----
p_heat <- ggplot(results_df, aes(x = factor(beta_sd), y = sd_institute, fill = elpd)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(
    title = "Performance landscape",
    x = "Beta SD",
    y = "Institute prior",
    fill = "ELPD"
  ) +
  theme_pub

ggsave("plots/C_05_performance_heatmap.png", p_heat, width = 7, height = 5, dpi = 300)

# ============================================================
# 5F) SAVE BEST MODEL
# ============================================================

best_model_id <- results_df$model_id[1]
best_fit <- results_list[[best_model_id]]$fit

saveRDS(best_fit, "plots/C_best_model.rds")

message("\n✓ EXTENSION COMPLETE")
message("Best model: ", best_model_id)
message("=============================\n")

