# ============================================================================
#  IMF INFLATION TIME SERIES ANALYSIS — PRODUCTION SCRIPT
#  Source : IMF World Economic Outlook (WEO) Database, April 2025
#  Dataset: IMF_data_cleaned.csv  (wide format, year columns 1980-2024)
#  Combines: Base ARIMA script + Extended GARCH & Forecast Evaluation script
#
#  Sections:
#    1.  Setup & Libraries
#    2.  Data Loading & Validation
#    3.  Exploratory Data Analysis (EDA)
#    4.  Time Series Diagnostics (ACF/PACF · ADF · KPSS · PP)
#    5.  ARIMA Modelling (auto + grid search)
#    6.  Forecasting 2024-2026
#    7.  Residual Diagnostics
#    8.  Correlation Analysis
#    9.  Volatility Diagnostics (Rolling SD · ARCH-LM · EWMA)
#    10. GARCH Modelling (GARCH 1,1 + 1,2; normal & Student-t)
#    11. Forecast Evaluation (RMSE · MAE · MAPE · Skill Score)
#    12. Final Visual Dashboard
# ============================================================================


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# SECTION 1 · SETUP & LIBRARIES
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(scales)
  library(tseries)
  library(forecast)
  library(urca)
  library(lmtest)
  library(gridExtra)
  library(corrplot)
  library(broom)
  library(strucchange)
  library(zoo)
  library(grid)
  library(cowplot)
  library(fGarch)
  library(e1071)
})

# Output folder for all plots and CSVs
OUT <- "plots/"
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

set.seed(42)

# ── Colour palettes ──────────────────────────────────────────────────────────
PALETTE <- c(
  USA = "#1f77b4",
  GBR = "#e377c2",
  DEU = "#2ca02c",
  FRA = "#d62728",
  JPN = "#ff7f0e",
  CHN = "#9467bd",
  IND = "#8c564b",
  BRA = "#17becf",
  RUS = "#bcbd22",
  TUR = "#7f7f7f",
  ARG = "#e6ab02",
  MEX = "#66a61e",
  ZAF = "#a6761d",
  NGA = "#e7298a",
  EGY = "#1b9e77"
)

# ── Custom ggplot2 theme ─────────────────────────────────────────────────────
theme_ts <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(
      plot.title    = element_text(face = "bold", size = base + 2, colour = "#1a1a2e"),
      plot.subtitle = element_text(size = base - 2, colour = "#555555", margin = margin(b = 6)),
      plot.caption  = element_text(size = base - 3.5, colour = "#888888"),
      axis.title    = element_text(size = base - 1),
      legend.title  = element_text(size = base - 2, face = "bold"),
      legend.text   = element_text(size = base - 2.5),
      panel.grid.minor = element_blank(),
      strip.text    = element_text(face = "bold", size = base - 1),
      plot.margin   = margin(10, 14, 8, 10),
      panel.border  = element_rect(colour = "#e0e0e0", fill = NA, linewidth = 0.4)
    )
}

# Helper: save ggplot to OUT folder
sp <- function(p, name, w = 14, h = 8) {
  ggsave(paste0(OUT, name), plot = p, width = w, height = h, dpi = 150, bg = "white")
  cat("  Saved:", name, "\n")
}

# Event shading rectangles (Asian Crisis, GFC, COVID, Post-COVID)
events <- data.frame(
  xmin = c(1997,   2008,   2020,   2021.3),
  xmax = c(1998.5, 2009.5, 2021,   2023.2),
  lbl  = c("Asian Crisis", "GFC", "COVID-19", "Post-COVID Surge"),
  fill = c("#f39c12", "#c0392b", "#2980b9", "#e67e22"),
  stringsAsFactors = FALSE
)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# SECTION 2 · DATA LOADING & VALIDATION
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

cat("\n══ SECTION 2: DATA LOADING & VALIDATION ══\n")

# ── 2.1  Read raw wide CSV ───────────────────────────────────────────────────
raw <- read.csv("C:/Users/hp/Desktop/IMF data cleaned.csv", stringsAsFactors = FALSE, check.names = FALSE)
cat("Raw dimensions:", nrow(raw), "rows x", ncol(raw), "columns\n")
cat("Column names:\n"); print(names(raw))

# ── 2.2  Auto-detect key column names ────────────────────────────────────────
country_col <- names(raw)[sapply(names(raw), function(n)
  grepl("^country$", n, ignore.case = TRUE))][1]
iso_col     <- names(raw)[grepl("^iso$", names(raw), ignore.case = TRUE)][1]
subj_col    <- names(raw)[grepl("descriptor", names(raw), ignore.case = TRUE)][1]
units_col   <- names(raw)[grepl("^units$", names(raw), ignore.case = TRUE)][1]
est_col     <- names(raw)[grepl("estimates", names(raw), ignore.case = TRUE)][1]
year_cols   <- names(raw)[grepl("^[12][0-9]{3}$", names(raw))]
year_nums   <- as.integer(year_cols)

cat("\nDetected columns:\n")
cat("  Country  :", country_col, "\n")
cat("  ISO      :", iso_col, "\n")
cat("  Subject  :", subj_col, "\n")
cat("  Units    :", units_col, "\n")
cat("  Estimates:", est_col, "\n")
cat("  Years    :", min(year_nums), "to", max(year_nums), "\n")
cat("  Subjects :", paste(unique(raw[[subj_col]]), collapse = " | "), "\n")
cat("  Units    :", paste(unique(raw[[units_col]]), collapse = " | "), "\n")

# ── 2.2b  Strip thousand-separator commas from all year columns ─────────────
# The IMF CSV formats hyperinflation values with commas (e.g. "2,947.73").
# These must be stripped BEFORE pivot_longer so they parse correctly as numeric.
for (yc in year_cols) {
  raw[[yc]] <- gsub(",", "", raw[[yc]])
}

# ── 2.3  Filter: CPI average, percent change (the key analytic series) ───────
cpi_filter <- raw[
  grepl("average consumer", raw[[subj_col]], ignore.case = TRUE) &
  grepl("percent change",   raw[[units_col]], ignore.case = TRUE),
]
cat("\nRows after CPI % change filter:", nrow(cpi_filter), "\n")
cat("Countries                      :", length(unique(cpi_filter[[country_col]])), "\n")

# ── 2.4  Pivot wide -> long ──────────────────────────────────────────────────
df_long_raw <- tidyr::pivot_longer(
  cpi_filter,
  cols      = all_of(year_cols),
  names_to  = "Year",
  values_to = "Inflation_Rate"
)
df_long_raw$Year           <- as.integer(df_long_raw$Year)
# Strip thousand-separator commas before parsing (e.g. "2,947.73" -> 2947.73)
# R as.numeric() silently returns NA for comma-formatted numbers, causing hyperinflation
# values (BRA, RUS, TUR etc.) to be lost and incorrectly interpolated later.
df_long_raw$Inflation_Rate <- suppressWarnings(as.numeric(gsub(",", "", df_long_raw$Inflation_Rate)))

# Flag years at or after Estimates Start After as projected
df_long_raw$Estimate <- !is.na(df_long_raw[[est_col]]) & df_long_raw$Year >= df_long_raw[[est_col]]

# Rename to standard columns
df_long_raw <- df_long_raw[, c(country_col, iso_col, "Year", "Inflation_Rate", "Estimate")]
names(df_long_raw)[names(df_long_raw) == country_col] <- "Country"
names(df_long_raw)[names(df_long_raw) == iso_col]     <- "ISO"

# ── 2.5  Clean: remove duplicates ───────────────────────────────────────────
n_before <- nrow(df_long_raw)
df_long_raw <- df_long_raw[!duplicated(df_long_raw[, c("Country","ISO","Year")]), ]
cat("Duplicates removed:", n_before - nrow(df_long_raw), "\n")

# ── 2.6  Filter to analytic window ──────────────────────────────────────────
df_long <- df_long_raw[df_long_raw$Year >= 1980 & df_long_raw$Year <= 2024, ]
cat("Final rows:", nrow(df_long), "| Years:", min(df_long$Year), "-", max(df_long$Year), "\n")
cat("Countries :", length(unique(df_long$Country)), "\n")
print(head(df_long, 6))

# ── 2.7  Country lookup & ISO sets ──────────────────────────────────────────
nm <- unique(df_long[, c("ISO","Country")])

all_iso <- c("USA","GBR","DEU","FRA","JPN","CHN","IND","BRA","RUS","TUR","ARG","MEX","ZAF","NGA","EGY")
missing_iso <- setdiff(all_iso, unique(df_long$ISO))
if (length(missing_iso) > 0) {
  cat("WARNING - ISOs not found:", paste(missing_iso, collapse = ", "), "\n")
  all_iso <- intersect(all_iso, unique(df_long$ISO))
}

model_iso <- intersect(c("USA","DEU","JPN","GBR","CHN","IND"), all_iso)
model_nm <- setNames(nm$Country[match(model_iso, nm$ISO)], model_iso)

cat("\nAnalysis ISOs :", paste(all_iso,   collapse = " "), "\n")
cat("Modelling ISOs:", paste(model_iso, collapse = " "), "\n")

# Economy group splits for aggregates and EDA plots
adv_isos      <- intersect(c("USA","DEU","JPN","GBR","FRA"), all_iso)
em_isos       <- intersect(c("CHN","IND","BRA","MEX","ZAF","TUR","ARG","NGA","EGY","RUS"), all_iso)

# EDA sub-groups by inflation magnitude — prevents scale-domination in charts:
# moderate_isos: max inflation <35%, fully comparable on a linear axis
# moderate_isos: max inflation <35% — all comparable on a shared linear axis
# NGA max=29%, EGY max=24%, CHN max=24%, ZAF max=15%, IND max=13%, MEX max=35%
moderate_isos <- intersect(c("CHN","IND","MEX","ZAF","EGY","NGA"), all_iso)
# high_isos    : 50-135%, y-axis capped at 150%
# TUR max=104%, ARG max=133%  (NGA moved to moderate — its max is only 29%)
high_isos     <- intersect(c("TUR","ARG"), all_iso)
# hyper_isos   : >150% peak, require log10 scale (BRA peak 2948%, RUS peak 874%)
hyper_isos    <- intersect(c("BRA","RUS"), all_iso)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# REUSABLE HELPER FUNCTIONS
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# Extract annual ts() for a single ISO
get_ts <- function(iso, start_yr = 1990, end_yr = 2023) {
  tryCatch({
    sub_df <- df_long[df_long$ISO == iso & df_long$Year >= start_yr & df_long$Year <= end_yr, ]
    sub_df <- sub_df[order(sub_df$Year), ]
    vals   <- sub_df$Inflation_Rate
    # Linearly interpolate internal NAs; leave leading/trailing NAs (trimmed below)
    if (any(is.na(vals)) && !all(is.na(vals))) {
      vals <- zoo::na.approx(vals, na.rm = FALSE)
    }
    # Trim leading/trailing NAs
    first_ok <- which(!is.na(vals))[1]
    last_ok  <- tail(which(!is.na(vals)), 1)
    if (is.na(first_ok)) return(NULL)
    actual_start <- sub_df$Year[first_ok]
    vals <- vals[first_ok:last_ok]
    ts(vals, start = actual_start, frequency = 1)
  }, error = function(e) { warning("get_ts: ", iso, " - ", e$message); NULL })
}

# Country name from ISO
iso_to_name <- function(iso) {
  v <- nm$Country[nm$ISO == iso]
  if (length(v) == 0 || is.na(v[1])) iso else v[1]
}

# Safely fit ARIMA
# Safely fit ARIMA
# CRITICAL: Use Arima() from {forecast} + as.numeric() — NOT base arima().
# Base arima() stores a ts reference that becomes NULL inside lapply(), making
# forecast() fail with "attempt to set attribute on NULL".
# This is the root cause of the "no valid grid models" / "None" in fig06.
safe_arima <- function(series, order) {
  tryCatch(Arima(as.numeric(series), order = order), error = function(e) NULL)
}

# Compute RMSE / MAE / sMAPE / Bias
# sMAPE replaces MAPE: avoids division-by-zero for near-zero actuals (e.g. Japan).
eval_metrics <- function(actual, fcst) {
  e     <- fcst - actual
  smape <- mean(2 * abs(e) / (abs(actual) + abs(fcst) + 1e-10), na.rm = TRUE) * 100
  data.frame(
    RMSE  = round(sqrt(mean(e^2, na.rm = TRUE)), 3),
    MAE   = round(mean(abs(e), na.rm = TRUE),    3),
    sMAPE = round(smape, 1),
    Bias  = round(mean(e, na.rm = TRUE), 3)
  )
}

# Event rectangles layer for ggplot
event_rects <- function() {
  geom_rect(
    data = events,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
    inherit.aes = FALSE,
    alpha = 0.11
  )
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# SECTION 3 · EXPLORATORY DATA ANALYSIS
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

cat("\n══ SECTION 3: EDA ══\n")

df_eda <- merge(
  df_long[df_long$ISO %in% all_iso & df_long$Year >= 1990 & df_long$Year <= 2023, ],
  setNames(nm, c("ISO","CountryName")),
  by = "ISO"
)

# ── FIG 01: Faceted country panels ──────────────────────────────────────────
cat("[FIG 01] Country panels\n")

p01 <- ggplot(df_eda, aes(Year, Inflation_Rate, colour = ISO,
                           linetype = as.character(Estimate))) +
  event_rects() +
  scale_fill_identity() +
  annotate("rect", xmin = 2023.5, xmax = 2024.5, ymin = -Inf, ymax = Inf,
           fill = "grey70", alpha = 0.07) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = PALETTE, guide = "none") +
  scale_linetype_manual(values = c("FALSE" = "solid", "TRUE" = "dashed"),
                        labels = c("Historical", "IMF Projection"), name = "") +
  scale_x_continuous(breaks = seq(1990, 2024, 4)) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  facet_wrap(~CountryName, scales = "free_y", ncol = 5) +
  labs(title    = "IMF WEO — Country Inflation Rate Panels (1990-2023)",
       subtitle = "Indicator: CPI % change (average consumer prices) | Source: IMF WEO April 2025",
       x = "Year", y = "Inflation Rate (%)") +
  theme_ts() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
sp(p01, "fig01_country_panels.png", 18, 14)

# ── FIG 02: Grouped line charts by economy type ──────────────────────────────
cat("[FIG 02] Grouped line charts\n")

make_group_plot <- function(isos, title_str) {
  df_sub  <- df_long[df_long$ISO %in% isos & df_long$Year >= 1990 & df_long$Year <= 2024, ]
  cnames  <- nm$Country[match(isos, nm$ISO)]
  ggplot(df_sub, aes(Year, Inflation_Rate, colour = ISO, linetype = as.character(Estimate))) +
    event_rects() +
    scale_fill_identity() +
    annotate("rect", xmin = 2023.5, xmax = 2024.5, ymin = -Inf, ymax = Inf, fill = "grey70", alpha = 0.07) +
    geom_line(linewidth = 1.0) +
    geom_point(data = df_sub[df_sub$Year == 2022, ], size = 2.5, show.legend = FALSE) +
    scale_colour_manual(values = PALETTE[isos], labels = cnames, name = "Country") +
    scale_linetype_manual(values = c("FALSE" = "solid", "TRUE" = "dashed"),
                          labels = c("Historical", "IMF Projection"), name = "") +
    scale_x_continuous(breaks = seq(1990, 2024, 5)) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(title = title_str, subtitle = "Source: IMF WEO April 2025 | CPI % change",
         x = "Year", y = "Inflation Rate (%)") +
    theme_ts() +
    theme(legend.position = "bottom")
}
# ── fig02a: Advanced economies (all max <10%) ────────────────────────────────
sp(make_group_plot(adv_isos,
  "Advanced Economies - Annual Inflation (1990-2024)"), "fig02a_advanced.png", 14, 7)

# ── fig02b: Moderate emerging (max <35%, comparable linear scale) ─────────────
# BRA (2948%), RUS (874%), TUR (104%), ARG (133%) are EXCLUDED from this panel —
# their extreme values would compress all other lines to near zero on a shared axis.
sp(make_group_plot(moderate_isos,
  "Moderate Emerging Markets — Annual Inflation (1990-2024)\n(CHN, IND, MEX, ZAF, EGY, NGA  |  max inflation <35%  |  comparable linear scale)"),
  "fig02b_moderate.png", 16, 7)

# ── fig02c: High-inflation (50-135%) capped at 150pp ─────────────────────────
make_high_plot <- function(isos, title_str, cap = 150) {
  df_sub        <- df_long[df_long$ISO %in% isos & df_long$Year >= 1990 & df_long$Year <= 2024, ]
  cnames        <- nm$Country[match(isos, nm$ISO)]
  df_sub$PlotY  <- pmin(df_sub$Inflation_Rate, cap)
  df_sub$AtCap  <- !is.na(df_sub$Inflation_Rate) & df_sub$Inflation_Rate > cap
  ggplot(df_sub, aes(Year, PlotY, colour = ISO, linetype = as.character(Estimate))) +
    event_rects() +
    scale_fill_identity() +
    geom_line(linewidth = 1.0) +
    geom_point(data = df_sub[df_sub$AtCap == TRUE, ], aes(Year, PlotY),
               shape = 24, size = 3, fill = "white", show.legend = FALSE) +
    scale_colour_manual(values = PALETTE[isos], labels = cnames, name = "Country") +
    scale_linetype_manual(values = c("FALSE" = "solid", "TRUE" = "dashed"),
                          labels = c("Historical", "IMF Projection"), name = "") +
    scale_x_continuous(breaks = seq(1990, 2024, 5)) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(NA, cap * 1.06)) +
    labs(title    = title_str,
         subtitle = paste0("Y-axis capped at ", cap, "pp | White triangles = values exceed cap | Source: IMF WEO April 2025"),
         x = "Year", y = paste0("Inflation Rate (%, capped at ", cap, "%)"),
         caption  = "Note: Argentina missing data 1990-1997 and 2014-2016 in IMF WEO dataset") +
    theme_ts() + theme(legend.position = "bottom")
}
sp(make_high_plot(high_isos, "High-Inflation Economies — Annual Inflation (1990-2024)\n(TUR & ARG only  |  max 50-135%  |  y-axis capped at 150%)", cap = 150),
  "fig02c_high_inflation.png", 14, 7)

# ── fig02d: Hyperinflation on log10 scale (BRA, RUS) ─────────────────────────
# BRA peaked at 2,948% (1990), RUS at 874% (1993). A linear scale makes both
# unreadable after stabilisation. Log10 scale shows the full trajectory clearly.
make_hyper_plot <- function(isos, title_str) {
  df_sub <- df_long[df_long$ISO %in% isos & df_long$Year >= 1988 & df_long$Year <= 2010, ]
  df_sub <- df_sub[!is.na(df_sub$Inflation_Rate) & df_sub$Inflation_Rate > 0, ]
  cnames <- nm$Country[match(isos, nm$ISO)]
  ggplot(df_sub, aes(Year, Inflation_Rate, colour = ISO)) +
    event_rects() + scale_fill_identity() +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.2) +
    scale_colour_manual(values = PALETTE[isos], labels = cnames, name = "Country") +
    scale_x_continuous(breaks = seq(1988, 2010, 2)) +
    scale_y_log10(labels = function(x) paste0(format(x, big.mark = ",", scientific = FALSE), "%"),
                  breaks  = c(1, 5, 10, 50, 100, 500, 1000, 3000)) +
    annotation_logticks(sides = "l", colour = "grey70") +
    labs(title    = title_str,
         subtitle = "Log\u2081\u2080 scale | BRA peak: 2,948% (1990) | RUS peak: 874% (1993)",
         x = "Year", y = "Inflation Rate (%, log\u2081\u2080 scale)",
         caption  = "Note: Plotted 1988-2010 to capture stabilisation. Linear scale impossible — both lines would be unreadable post-2000.") +
    theme_ts() + theme(legend.position = "bottom")
}
sp(make_hyper_plot(hyper_isos, "Hyperinflation Economies 1988-2010 — BRA & RUS (Log\u2081\u2080 Scale)"),
  "fig02d_hyperinflation.png", 14, 7)

# ── FIG 03a: Box plots ───────────────────────────────────────────────────────
cat("[FIG 03] Distribution plots\n")

df_eda$Capped <- pmin(df_eda$Inflation_Rate, 80)
p03a <- ggplot(df_eda, aes(reorder(CountryName, Inflation_Rate, FUN = median), Capped, fill = ISO)) +
  geom_boxplot(alpha = 0.8, outlier.size = 1.2, outlier.alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_fill_manual(values = PALETTE, guide = "none") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  coord_flip() +
  labs(title    = "Distribution of Inflation Rates by Country (1990-2023)",
       subtitle = "Values capped at 80% for scale | Source: IMF WEO April 2025",
       x = NULL, y = "Inflation Rate (%)") +
  theme_ts()
sp(p03a, "fig03a_boxplots.png", 12, 8)

# ── FIG 03b: Violin plots ────────────────────────────────────────────────────
p03b <- ggplot(df_eda, aes(ISO, Capped, fill = ISO)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_jitter(width = 0.07, size = 1.2, alpha = 0.55, colour = "#333333") +
  geom_boxplot(width = 0.12, fill = "white", alpha = 0.8, outlier.shape = NA) +
  scale_fill_manual(values = PALETTE, guide = "none") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(title    = "Inflation Distribution - Violin Plots (1990-2023)",
       subtitle = "Each dot = one year | Source: IMF WEO April 2025",
       x = NULL, y = "Inflation Rate (%)") +
  theme_ts() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sp(p03b, "fig03b_violins.png", 14, 7)

# ── FIG 03c: Histograms ──────────────────────────────────────────────────────
hist_isos <- c("USA","DEU","JPN","CHN","IND","BRA")
df_hist   <- df_eda[df_eda$ISO %in% hist_isos, ]
mean_df   <- aggregate(Inflation_Rate ~ ISO + CountryName, df_hist, mean)
names(mean_df)[3] <- "m"

p03c <- ggplot(df_hist, aes(Inflation_Rate, fill = ISO)) +
  geom_histogram(bins = 12, colour = "white", alpha = 0.85) +
  geom_vline(data = mean_df, aes(xintercept = m), colour = "red", linetype = "dashed", linewidth = 0.9) +
  facet_wrap(~CountryName, scales = "free", ncol = 3) +
  scale_fill_manual(values = PALETTE, guide = "none") +
  labs(title    = "Inflation Rate Histograms - Selected Countries (1990-2023)",
       subtitle = "Red dashed = sample mean | Source: IMF WEO April 2025",
       x = "Inflation Rate (%)", y = "Frequency") +
  theme_ts()
sp(p03c, "fig03c_histograms.png", 14, 9)

# ── FIG 03d: Heatmap ─────────────────────────────────────────────────────────
df_heat <- df_eda[df_eda$Year >= 2000 & df_eda$Year <= 2023, ]
df_heat$Cap50 <- pmin(df_heat$Inflation_Rate, 50)

p03d <- ggplot(df_heat, aes(Year, reorder(CountryName, Inflation_Rate, FUN = mean), fill = Cap50)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradientn(
    colours = c("#2166ac","#f7f7f7","#f4a582","#d6604d","#b2182b"),
    values  = c(0, 0.03, 0.12, 0.25, 1),
    limits  = c(0, 50), oob = squish, name = "Inflation %\n(cap 50%)"
  ) +
  scale_x_continuous(breaks = seq(2000, 2023, 2)) +
  labs(title    = "Inflation Rate Heatmap - Countries x Years (2000-2023)",
       subtitle = "Colour intensity = inflation level | Source: IMF WEO April 2025",
       x = "Year", y = NULL) +
  theme_ts() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sp(p03d, "fig03d_heatmap.png", 15, 8)

# ── Descriptive statistics CSV ────────────────────────────────────────────────
cat("\nDescriptive statistics:\n")
desc_df <- do.call(rbind, lapply(all_iso, function(iso) {
  s <- na.omit(as.numeric(get_ts(iso)))
  if (length(s) < 5) return(NULL)
  data.frame(
    Country  = iso_to_name(iso), N = length(s),
    Mean     = round(mean(s), 2), Median = round(median(s), 2),
    SD       = round(sd(s), 2),   Min    = round(min(s), 2),
    Max      = round(max(s), 2),
    Skewness = round(e1071::skewness(s, type = 2), 3),
    CV       = round(sd(s) / abs(mean(s)) * 100, 1)
  )
}))
print(desc_df)
write.csv(desc_df, paste0(OUT, "descriptive_stats.csv"), row.names = FALSE)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# SECTION 4 · TIME SERIES DIAGNOSTICS
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

cat("\n══ SECTION 4: TIME SERIES DIAGNOSTICS ══\n")

# ── FIG 04: ACF / PACF panels ────────────────────────────────────────────────
cat("[FIG 04] ACF / PACF\n")

png(paste0(OUT, "fig04_acf_pacf.png"), width = 1800, height = 2200, res = 130, bg = "white")
par(mfrow = c(length(model_iso), 2), mar = c(3.5, 3.5, 2.5, 1),
    oma = c(1, 1, 3, 1), bg = "white")
for (iso in model_iso) {
  s <- get_ts(iso)
  if (is.null(s)) next
  acf(s,  lag.max = 20, main = paste("ACF - ",  iso_to_name(iso)),
      col = "#1f77b4", lwd = 2, ci.col = "#d62728")
  pacf(s, lag.max = 15, main = paste("PACF - ", iso_to_name(iso)),
       col = "#d62728", lwd = 2, ci.col = "#1f77b4")
}
mtext("ACF & PACF - IMF WEO Inflation Rates (1990-2023)",
      outer = TRUE, cex = 1.3, font = 2, col = "#1a1a2e")
dev.off()
cat("  Saved: fig04_acf_pacf.png\n")

# ── Stationarity tests: ADF, KPSS, PP ────────────────────────────────────────
cat("\nStationarity tests:\n")
stat_results <- lapply(all_iso, function(iso) {
  s    <- get_ts(iso)
  name <- iso_to_name(iso)
  sv <- as.numeric(na.omit(s))
  if (is.null(s) || length(sv) < 10) return(data.frame(Country = name, ISO = iso))
  tryCatch({
    adf_r  <- adf.test(sv)
    pp_r   <- pp.test(sv)
    kpss_r <- kpss.test(sv, null = "Level")
    data.frame(
      Country = name, ISO = iso,
      ADF_stat  = round(adf_r$statistic,  3), ADF_pval  = round(adf_r$p.value,  4),
      PP_stat   = round(pp_r$statistic,   3), PP_pval   = round(pp_r$p.value,   4),
      KPSS_stat = round(kpss_r$statistic, 3), KPSS_pval = round(kpss_r$p.value, 4),
      Stationary = adf_r$p.value < 0.05
    )
  }, error = function(e) {
    warning("Stationarity failed: ", iso, " - ", e$message)
    data.frame(Country = name, ISO = iso)
  })
})
stat_df <- do.call(rbind, stat_results)
cols_to_show <- intersect(c("Country","ADF_pval","PP_pval","KPSS_pval","Stationary"), names(stat_df))
print(stat_df[, cols_to_show])
write.csv(stat_df, paste0(OUT, "stationarity_tests.csv"), row.names = FALSE)

# ── FIG 05a: ADF p-value bar chart ──────────────────────────────────────────
stat_plot_df <- stat_df[!is.na(stat_df$ADF_pval), ]
p05a <- ggplot(stat_plot_df,
               aes(reorder(Country, ADF_pval), ADF_pval, fill = Stationary)) +
  geom_col(alpha = 0.85, colour = "white") +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "black", linewidth = 1) +
  geom_hline(yintercept = 0.01, linetype = "dotted", colour = "navy", linewidth = 0.8) +
  geom_text(aes(label = paste0("p=", ADF_pval)), hjust = -0.1, size = 2.8) +
  scale_fill_manual(
    values = c("TRUE" = "#27ae60", "FALSE" = "#e74c3c"),
    labels = c("TRUE" = "Stationary (p<0.05)", "FALSE" = "Unit Root (p>=0.05)"),
    name   = "ADF Result"
  ) +
  scale_y_continuous(limits = c(0, 0.75)) +
  coord_flip() +
  labs(title    = "Augmented Dickey-Fuller (ADF) Unit Root Test",
       subtitle = "H0: Unit root exists | Reject H0 at p<0.05 = stationary",
       x = NULL, y = "ADF p-value", caption = "IMF WEO April 2025 | 1990-2023") +
  theme_ts()
sp(p05a, "fig05a_adf_pvalues.png", 12, 8)

# ── FIG 05b: Original vs first-differenced (USA) ────────────────────────────
cat("[FIG 05b] Differencing example\n")
usa_raw  <- get_ts("USA")
usa_diff <- diff(usa_raw)
df_usa_d <- rbind(
  data.frame(Year = as.numeric(time(usa_raw)),  Value = as.numeric(usa_raw),  Series = "Original Series"),
  data.frame(Year = as.numeric(time(usa_diff)), Value = as.numeric(usa_diff), Series = "First Difference (Delta)")
)
p05b <- ggplot(df_usa_d, aes(Year, Value, colour = Series)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  facet_wrap(~Series, scales = "free_y", nrow = 2) +
  scale_colour_manual(
    values = c("Original Series" = "#1f77b4", "First Difference (Delta)" = "#d62728"),
    guide  = "none"
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(title    = "USA Inflation - Original vs. First-Differenced Series",
       subtitle = "First differencing removes trend/unit root",
       x = "Year", y = "Inflation Rate (%)") +
  theme_ts()
sp(p05b, "fig05b_differencing.png", 12, 7)

# ── FIG 05c: CUSUM structural stability ─────────────────────────────────────
cat("[FIG 05c] CUSUM stability test\n")
tryCatch({
  cusum_obj <- efp(usa_raw ~ 1, type = "Rec-CUSUM")
  bp_obj    <- breakpoints(usa_raw ~ 1)
  png(paste0(OUT, "fig05c_cusum.png"), width = 1400, height = 600, res = 130, bg = "white")
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 2), oma = c(0, 0, 2.5, 0), bg = "white")
  plot(cusum_obj, main = "CUSUM Test - USA Inflation", col = "#1f77b4", lwd = 2)
  plot(bp_obj,    main = "Structural Breakpoints - USA Inflation")
  if (!is.na(bp_obj$breakpoints[1]))
    abline(v = time(usa_raw)[bp_obj$breakpoints], col = "#e67e22", lwd = 2, lty = 2)
  mtext("Structural Stability Tests - USA Inflation (1990-2023)",
        outer = TRUE, cex = 1.2, font = 2)
  dev.off()
  cat("  Saved: fig05c_cusum.png\n")
}, error = function(e) warning("CUSUM failed: ", e$message))


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# SECTION 5 · ARIMA MODELLING
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

cat("\n══ SECTION 5: ARIMA MODELLING ══\n")
cat("Primary evaluation  : Train 1990-2014 | Pre-COVID test 2015-2019\n")
cat("Stress-test         : Train 1990-2019 | COVID-era  test 2020-2023\n\n")

# DESIGN NOTES:
#  (a) PRIMARY   train=1990-2014, test=2015-2019: clean pre-COVID window.
#      This is the headline evaluation of ARIMA forecast skill.
#  (b) STRESS    train=1990-2019, test=2020-2023: COVID structural break.
#      Large errors reflect the shock, not model failure. Reported separately.
#  Model SELECTION: AICc on training set (penalises complexity on small n).
#  Model EVALUATION: RMSE on held-out test set (true out-of-sample accuracy).
#  Full-sample model (1990-2023): used only for the 2024-2026 point forecast.

model_results <- lapply(model_iso, function(iso) {
  cat("  Fitting:", iso_to_name(iso), "\n")

  # Primary split (pre-COVID)
  s_primary   <- get_ts(iso, 1990, 2019)
  s_train_pre <- window(s_primary, end = 2014)
  s_test_pre  <- window(s_primary, start = 2015)

  # Stress split (COVID era)
  s_full_to23 <- get_ts(iso, 1990, 2023)
  s_train_cov <- window(s_full_to23, end = 2019)
  s_test_cov  <- window(s_full_to23, start = 2020)

  # auto.arima on primary training set
  auto_m <- tryCatch(
    auto.arima(as.numeric(s_train_pre), seasonal = FALSE, stepwise = TRUE, ic = "aicc"),
    error = function(e) auto.arima(as.numeric(s_train_pre))
  )
  cat("    auto.arima:", capture.output(auto_m)[1],
      "| AICc:", round(auto_m$aicc, 1), "\n")

  # Grid search p in 0:2, d in 0:1, q in 0:2 — ranked by AICc
  cands <- list()
  gc()
  for (p in 0:2) {
    for (d in 0:1) {
      for (q in 0:2) {
        key <- paste0("ARIMA(", p, ",", d, ",", q, ")")
        m   <- safe_arima(s_train_pre, c(p, d, q))
        if (!is.null(m)) {
          fc_v <- tryCatch(
            as.numeric(forecast(m, h = length(s_test_pre))$mean),
            error = function(e) NULL
          )
          if (!is.null(fc_v)) {
            em       <- eval_metrics(as.numeric(s_test_pre), fc_v)
            k        <- length(m$coef) + 1
            n        <- length(s_train_pre)
            aicc_val <- AIC(m) + 2 * k * (k + 1) / max(n - k - 1, 1)
            cands[[key]] <- data.frame(
              Model = key, AICc = round(aicc_val, 2),
              AIC = round(AIC(m), 2), BIC = round(BIC(m), 2),
              RMSE_preCOVID = em$RMSE, MAE_preCOVID = em$MAE
            )
          }
        }
      }
    }
  }
  gc()
  cand_df <- do.call(rbind, cands)
  if (is.null(cand_df) || nrow(cand_df) == 0) {
    cat("    WARNING: No valid grid models for", iso, "\n")
    cand_df <- data.frame(Model = "None", AICc = NA, AIC = NA, BIC = NA,
                          RMSE_preCOVID = NA, MAE_preCOVID = NA)
  } else {
    cand_df <- cand_df[order(cand_df$AICc), ]
    cat("    Best by AICc:", cand_df$Model[1], "| AICc =", cand_df$AICc[1], "\n")
  }

  # Refit on full series 1990-2023 for forecasting 2024-2026
  m_full <- tryCatch(
    Arima(as.numeric(s_full_to23), order = auto_m$arma[c(1, 6, 2)]),
    error = function(e) auto.arima(as.numeric(s_full_to23), seasonal = FALSE, stepwise = TRUE)
  )
  cat("    Full-sample ARIMA(", paste(m_full$arma[c(1,6,2)], collapse = ","), ")",
      "| AIC:", round(AIC(m_full), 1), "\n")

  list(iso = iso, name = iso_to_name(iso),
       model = auto_m, full_model = m_full,
       train_pre = s_train_pre, test_pre  = s_test_pre,
       train_cov = s_train_cov, test_cov  = s_test_cov,
       full_ts   = s_full_to23, cand_table = cand_df)
})
names(model_results) <- model_iso

# Model comparison CSV (grid, ranked by AICc)
all_cands_df <- do.call(rbind, lapply(model_iso, function(iso) {
  df <- head(model_results[[iso]]$cand_table, 10)
  df$Country <- iso_to_name(iso)
  df
}))
write.csv(all_cands_df, paste0(OUT, "arima_model_comparison.csv"), row.names = FALSE)

# ── FIG 06: AICc comparison bar ───────────────────────────────────────────────
cat("[FIG 06] ARIMA AICc comparison\n")
plot_cands <- all_cands_df[!is.na(all_cands_df$AICc), ]
if (nrow(plot_cands) > 0) {
  p06 <- ggplot(plot_cands, aes(reorder(Model, AICc), AICc, fill = Country)) +
    geom_col(alpha = 0.85, colour = "white") +
    facet_wrap(~Country, scales = "free_x", nrow = 2) +
    scale_fill_manual(
      values = setNames(unname(PALETTE[model_iso]), sapply(model_iso, iso_to_name)),
      guide  = "none"
    ) +
    coord_flip() +
    labs(title    = "ARIMA Model Selection — AICc by Country",
         subtitle = "Lower AICc preferred | Train: 1990-2014 (pre-COVID primary window)",
         x = "Model", y = "AICc", caption = "Source: IMF WEO April 2025") +
    theme_ts() +
    theme(axis.text.y = element_text(size = 8))
  sp(p06, "fig06_model_aic_comparison.png", 16, 10)
}

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# SECTION 6 · FORECASTING 2024-2026
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

cat("\n══ SECTION 6: FORECASTING 2024-2026 ══\n")
cat("NOTE: This dataset contains only ONE year of IMF figures beyond 2023:\n")
cat("  2024 = IMF near-final estimate (Estimates Start After = 2024).\n")
cat("  No 2025 or 2026 IMF data exists in this CSV. The green square marks\n")
cat("  the 2024 IMF estimate only — NOT a 3-year IMF projection line.\n\n")

fc_list <- lapply(model_iso, function(iso) {
  res   <- model_results[[iso]]
  m_f   <- res$full_model
  fc    <- forecast(m_f, h = 3, level = c(80, 95))
  fc_df <- data.frame(
    Year  = 2024:2026,
    Point = as.numeric(fc$mean),
    Lo80  = as.numeric(fc$lower[, 1]), Hi80 = as.numeric(fc$upper[, 1]),
    Lo95  = as.numeric(fc$lower[, 2]), Hi95 = as.numeric(fc$upper[, 2])
  )
  # Only 2024 available as IMF estimate — do NOT plot a 2025-2026 IMF line
  imf_2024 <- df_long[df_long$ISO == iso & df_long$Year == 2024,
                       c("Year","Inflation_Rate")]
  names(imf_2024)[2] <- "IMF"
  hist_df <- data.frame(Year  = as.numeric(time(res$full_ts)),
                        Value = as.numeric(res$full_ts))
  ord_str <- paste0("ARIMA(", paste(m_f$arma[c(1, 6, 2)], collapse = ","), ")")

  p <- ggplot() +
    geom_ribbon(data = fc_df, aes(Year, ymin = Lo95, ymax = Hi95),
                fill = "#1f77b4", alpha = 0.12) +
    geom_ribbon(data = fc_df, aes(Year, ymin = Lo80, ymax = Hi80),
                fill = "#1f77b4", alpha = 0.22) +
    geom_line(data  = hist_df, aes(Year, Value),
              colour = "#1f77b4", linewidth = 1) +
    geom_point(data = hist_df, aes(Year, Value),
               colour = "#1f77b4", size = 1.2) +
    geom_line(data  = fc_df, aes(Year, Point),
              colour = "#e74c3c", linewidth = 1.2, linetype = "dashed") +
    geom_point(data = fc_df, aes(Year, Point),
               colour = "#e74c3c", size = 3, shape = 17) +
    geom_vline(xintercept = 2023.5, linetype = "dotted", colour = "grey50") +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    scale_x_continuous(breaks = seq(1990, 2026, 5)) +
    labs(title    = res$name,
         subtitle = paste0(ord_str, " | AIC=", round(AIC(m_f), 1),
                           " | trained 1990-2023"),
         x = NULL, y = "Inflation (%)") +
    theme_ts() +
    theme(plot.title = element_text(size = 11))

  # Add 2024 IMF estimate as a single reference point (not a trend line)
  if (nrow(imf_2024) > 0 && !is.na(imf_2024$IMF[1])) {
    p <- p +
      geom_point(data = imf_2024, aes(Year, IMF),
                 colour = "#2ca02c", size = 4.5, shape = 15) +
      annotate("text", x = 2024.15, y = imf_2024$IMF[1],
               label = paste0("IMF\nest.\n", round(imf_2024$IMF[1], 1), "%"),
               hjust = 0, size = 2.6, colour = "#2ca02c", fontface = "bold")
  }
  p
})

cat("[FIG 07] ARIMA forecasts\n")
png(paste0(OUT, "fig07_arima_forecasts.png"),
    width = 2400, height = 1600, res = 100, bg = "white")
do.call(grid.arrange, c(fc_list, list(
  ncol = 3,
  top  = textGrob(
    "ARIMA Forecasts 2024-2026  |  Green square = 2024 IMF near-final estimate only",
    gp = gpar(fontsize = 12, fontface = "bold", col = "#1a1a2e")
  )
)))
dev.off()
cat("  Saved: fig07_arima_forecasts.png\n")

# SECTION 7 · RESIDUAL DIAGNOSTICS
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

cat("\n══ SECTION 7: RESIDUAL DIAGNOSTICS ══\n")

diag_list <- lapply(model_iso, function(iso) {
  m_f  <- model_results[[iso]]$full_model
  r    <- residuals(m_f)
  lb   <- Box.test(r, lag = 10, type = "Ljung-Box")
  sw   <- shapiro.test(r)
  cat(sprintf("  %-18s | LB p=%.4f | SW p=%.4f\n", iso_to_name(iso), lb$p.value, sw$p.value))
  list(iso = iso, name = iso_to_name(iso), model = m_f,
       resid = r, lb_p = lb$p.value, sw_p = sw$p.value)
})
names(diag_list) <- model_iso

cat("[FIG 08] Residual panels\n")
png(paste0(OUT, "fig08_residual_diagnostics.png"),
    width = 2000, height = 2200, res = 130, bg = "white")
par(mfrow = c(length(model_iso), 4), mar = c(3, 3, 2.5, 1),
    oma = c(0, 0, 4, 0), bg = "white")
for (iso in model_iso) {
  d  <- diag_list[[iso]]; r <- d$resid
  plot(r, main = paste0(d$name, "\nResiduals"), type = "l",
       col = "#1f77b4", lwd = 1.5, xlab = "Year", ylab = "Resid")
  abline(h = 0, col = "red", lty = 2)
  acf(r,  lag.max = 15, main = "ACF of Residuals",
      col = "#2ca02c", lwd = 2, ci.col = "#e74c3c")
  qqnorm(r, main = "Normal Q-Q", pch = 20, col = "#9467bd"); qqline(r, col = "red", lwd = 2)
  hist(r, breaks = 10, main = "Residual Histogram", probability = TRUE,
       col = "#aec7e8", border = "white", xlab = "Residual")
  curve(dnorm(x, mean(r), sd(r)), add = TRUE, col = "red", lwd = 2)
}
mtext("Residual Diagnostics - ARIMA Models | IMF WEO Inflation",
      outer = TRUE, cex = 1.4, font = 2, col = "#1a1a2e")
dev.off()
cat("  Saved: fig08_residual_diagnostics.png\n")

hyp_df <- do.call(rbind, lapply(model_iso, function(iso) {
  d           <- diag_list[[iso]]
  adf_series  <- adf.test(as.numeric(na.omit(get_ts(iso))))
  adf_resid   <- adf.test(as.numeric(na.omit(d$resid)))
  data.frame(
    Country                        = d$name,
    Series_Stationary_ADF          = ifelse(adf_series$p.value < 0.05, "Yes", "No"),
    Residuals_Stationary           = ifelse(adf_resid$p.value  < 0.05, "Yes", "No"),
    No_Autocorr_LB                 = ifelse(d$lb_p > 0.05, "Yes (good)", "No"),
    Normal_Residuals_SW            = ifelse(d$sw_p > 0.05, "Yes", "No"),
    LB_pval = round(d$lb_p, 4), SW_pval = round(d$sw_p, 4)
  )
}))
print(hyp_df)
write.csv(hyp_df, paste0(OUT, "hypothesis_tests.csv"), row.names = FALSE)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# SECTION 8 · CORRELATION ANALYSIS
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

cat("\n══ SECTION 8: CORRELATION ANALYSIS ══\n")

df_cor_sub <- df_long[df_long$ISO %in% all_iso & df_long$Year >= 1990 & df_long$Year <= 2023, ]
df_cor_sub <- merge(df_cor_sub, setNames(nm, c("ISO","CountryName")), by = "ISO")
df_wide_cor <- tidyr::pivot_wider(
  df_cor_sub[, c("CountryName","Year","Inflation_Rate")],
  names_from  = CountryName,
  values_from = Inflation_Rate
)
df_wide_cor$Year <- NULL
corr_mat <- cor(df_wide_cor, use = "pairwise.complete.obs")

cat("[FIG 09a] Correlation matrix\n")
png(paste0(OUT, "fig09a_corrplot.png"), width = 1600, height = 1400, res = 130, bg = "white")
par(mar = c(2, 2, 3, 2), bg = "white")
corrplot(corr_mat,
         method = "ellipse", type = "upper", tl.cex = 0.85,
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         number.cex = 0.65,
         col = colorRampPalette(c("#2166ac","#f7f7f7","#b2182b"))(200),
         title = "Pairwise Inflation Correlations (1990-2023) | IMF WEO",
         mar   = c(0, 0, 2, 0))
dev.off()
cat("  Saved: fig09a_corrplot.png\n")

# Advanced vs Emerging scatter
df_adv_agg <- aggregate(Inflation_Rate ~ Year,
  df_long[df_long$ISO %in% adv_isos & df_long$Year >= 1990 & df_long$Year <= 2023, ],
  mean, na.rm = TRUE)
names(df_adv_agg)[2] <- "Advanced"
df_em_agg  <- aggregate(Inflation_Rate ~ Year,
  df_long[df_long$ISO %in% em_isos  & df_long$Year >= 1990 & df_long$Year <= 2023, ],
  mean, na.rm = TRUE)
names(df_em_agg)[2] <- "Emerging"
df_sc  <- merge(df_adv_agg, df_em_agg, by = "Year")
r_val  <- cor(df_sc$Advanced, df_sc$Emerging, use = "complete.obs")

p09b <- ggplot(df_sc, aes(Advanced, Emerging)) +
  geom_point(aes(colour = Year), size = 3.5, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, colour = "#d62728", fill = "#f4c2c2", linewidth = 1.2) +
  scale_colour_viridis_c(option = "plasma", name = "Year") +
  geom_text(data = df_sc[df_sc$Year %in% c(2000,2008,2009,2020,2021,2022,2023), ],
            aes(label = Year), vjust = -0.8, size = 3, colour = "#333333") +
  annotate("text", x = min(df_sc$Advanced) + 0.1, y = max(df_sc$Emerging) - 1,
           label = paste0("Pearson r = ", round(r_val, 3)),
           size = 4.5, hjust = 0, fontface = "bold", colour = "#2c3e50") +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(title    = "Advanced vs. Emerging Market Average Inflation (1990-2023)",
       subtitle = "Positive co-movement with heterogeneous magnitude",
       x = "Advanced avg (%)", y = "Emerging avg (%)",
       caption = "Source: IMF WEO April 2025") +
  theme_ts()
sp(p09b, "fig09b_adv_em_scatter.png", 10, 7)

corr_m <- reshape2::melt(corr_mat)
p09c <- ggplot(corr_m, aes(Var1, Var2, fill = value)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(value, 2)), size = 2.2) +
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  labs(title = "Inflation Correlation Matrix (1990-2023)", x = NULL, y = NULL) +
  theme_ts() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7))
sp(p09c, "fig09c_corr_heatmap.png", 13, 11)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# SECTION 9 · VOLATILITY DIAGNOSTICS
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

cat("\n══ SECTION 9: VOLATILITY DIAGNOSTICS ══\n")

# ── Rolling 5-year SD ────────────────────────────────────────────────────────
cat("[FIG 10] Rolling 5-year volatility\n")
df_roll <- do.call(rbind, lapply(all_iso, function(iso) {
  ts_obj <- get_ts(iso)
  if (is.null(ts_obj)) return(NULL)
  s  <- as.numeric(ts_obj)
  yr <- as.integer(time(ts_obj))
  if (length(s) < 5) return(NULL)
  roll_sd <- sapply(5:length(s), function(i) sd(s[max(1, i-4):i], na.rm = TRUE))
  data.frame(Year = yr[5:length(yr)], RollSD = roll_sd,
             Country = iso_to_name(iso), ISO = iso)
}))

p10 <- ggplot(df_roll, aes(Year, RollSD, colour = ISO)) +
  event_rects() +
  scale_fill_identity() +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = PALETTE,
                      labels = sapply(names(PALETTE), iso_to_name), name = "Country") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_x_continuous(breaks = seq(1994, 2023, 4)) +
  labs(title    = "Rolling 5-Year Volatility of Inflation - All 15 Countries",
       subtitle = "sigma(t) = SD over rolling 5-year window | Source: IMF WEO April 2025",
       x = "Year", y = "Rolling Std. Deviation (%)") +
  theme_ts() +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 3))
sp(p10, "fig10_rolling_volatility.png", 14, 8)

# ── ARCH-LM tests ────────────────────────────────────────────────────────────
cat("\nARCH-LM tests on ARIMA residuals:\n")
arch_results <- lapply(model_iso, function(iso) {
  r        <- residuals(model_results[[iso]]$full_model)
  r2       <- r^2
  lb5      <- Box.test(r2, lag = 5,  type = "Ljung-Box")
  lb10     <- Box.test(r2, lag = 10, type = "Ljung-Box")
  sk       <- e1071::skewness(as.numeric(r))
  ku       <- e1071::kurtosis(as.numeric(r))
  cat(sprintf("  %-18s | Kurt=%.2f | Skew=%.3f | ARCH-LM(5) p=%.4f | ARCH-LM(10) p=%.4f\n",
              iso_to_name(iso), ku, sk, lb5$p.value, lb10$p.value))
  list(iso = iso, name = iso_to_name(iso), resid = r, resid_sq = r2,
       skew = sk, kurt = ku, arch5_p = lb5$p.value, arch10_p = lb10$p.value,
       arch_present = lb5$p.value < 0.10)
})
names(arch_results) <- model_iso

# ── FIG 11: Squared residuals panels ────────────────────────────────────────
cat("[FIG 11] Squared residuals\n")
png(paste0(OUT, "fig11_squared_residuals.png"),
    width = 2000, height = 1600, res = 130, bg = "white")
par(mfrow = c(length(model_iso), 3), mar = c(3, 3.5, 2.5, 1),
    oma = c(0, 0, 4, 0), bg = "white")
for (iso in model_iso) {
  a  <- arch_results[[iso]]
  r  <- as.numeric(a$resid); r2 <- as.numeric(a$resid_sq)
  yr <- seq(start(model_results[[iso]]$full_model$x)[1],
            by = 1, length.out = length(r))
  plot(yr, r, type = "l", col = "#1f77b4", lwd = 1.8,
       main = paste0(a$name, "\nARIMA Residuals"), xlab = "", ylab = "Residual", las = 1)
  abline(h = 0, col = "grey50", lty = 2)
  abline(h = c(-2*sd(r), 2*sd(r)), col = "#e74c3c", lty = 3, lwd = 1.2)
  text(2020, max(r, na.rm = TRUE)*0.85, paste0("Kurt=", round(a$kurt, 2)), cex = 0.85, col = "#9467bd", font = 2)
  plot(yr, r2, type = "l", col = "#e74c3c", lwd = 1.8,
       main = "Squared Residuals (Volatility Proxy)", xlab = "", ylab = "e^2", las = 1)
  abline(h = mean(r2, na.rm = TRUE), col = "#2ca02c", lty = 2, lwd = 1.5)
  text(2020, max(r2, na.rm = TRUE)*0.88,
       paste0("ARCH-LM\np=", round(a$arch5_p, 3)), cex = 0.8,
       col = if (a$arch_present) "#e74c3c" else "#27ae60")
  acf(r2, lag.max = 12, main = "ACF of Squared Residuals",
      col = "#9467bd", lwd = 2, ci.col = "#e74c3c", xlab = "Lag", las = 1)
}
mtext("Volatility Diagnostics - Squared ARIMA Residuals & ARCH Effect Testing",
      outer = TRUE, cex = 1.3, font = 2, col = "#1a1a2e")
dev.off()
cat("  Saved: fig11_squared_residuals.png\n")

arch_sum_df <- do.call(rbind, lapply(arch_results, function(a) {
  data.frame(Country = a$name, Skewness = round(a$skew, 3),
             ExcessKurt = round(a$kurt, 3), ARCH5_p = round(a$arch5_p, 4),
             ARCH10_p = round(a$arch10_p, 4), ARCH_Present = a$arch_present)
}))
print(arch_sum_df)
write.csv(arch_sum_df, paste0(OUT, "arch_test_results.csv"), row.names = FALSE)

p12 <- ggplot(arch_sum_df, aes(reorder(Country, -ExcessKurt), ExcessKurt, fill = ARCH_Present)) +
  geom_col(alpha = 0.85, colour = "white", width = 0.65) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_text(aes(label = paste0("Kurt=", ExcessKurt, "\np=", ARCH5_p)),
            vjust = -0.4, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = c("TRUE" = "#e74c3c", "FALSE" = "#27ae60"),
                    labels = c("TRUE" = "ARCH Effects (p<0.10)", "FALSE" = "No ARCH"), name = "") +
  scale_y_continuous(limits = c(-3, 9)) +
  labs(title = "Excess Kurtosis & ARCH-LM Test - ARIMA Residuals",
       subtitle = "Excess kurtosis > 0 = heavy tails | ARCH-LM(5) p-value annotated",
       x = NULL, y = "Excess Kurtosis") +
  theme_ts() +
  theme(legend.position = "bottom")
sp(p12, "fig12_arch_summary.png", 12, 7)

# ── EWMA conditional variance ────────────────────────────────────────────────
cat("[FIG 13] EWMA volatility\n")
ewma_vol <- do.call(rbind, lapply(model_iso, function(iso) {
  r  <- as.numeric(residuals(model_results[[iso]]$full_model))
  lambda <- 0.94
  h2 <- numeric(length(r)); h2[1] <- r[1]^2
  for (t in 2:length(r)) h2[t] <- lambda*h2[t-1] + (1-lambda)*r[t-1]^2
  yr_ewma <- seq(start(model_results[[iso]]$full_model$x)[1],
               by = 1, length.out = length(r))
  data.frame(Year = yr_ewma, EWMA_Vol = sqrt(h2),
             Country = iso_to_name(iso), ISO = iso)
}))

p13 <- ggplot(ewma_vol, aes(Year, EWMA_Vol, colour = ISO)) +
  event_rects() +
  scale_fill_identity() +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = PALETTE[model_iso], labels = model_nm, name = "Country") +
  scale_y_continuous(labels = function(x) paste0(round(x, 1), "%")) +
  labs(title    = "EWMA Conditional Volatility of ARIMA Residuals (lambda=0.94)",
       subtitle = "h^2(t) = lambda*h^2(t-1) + (1-lambda)*epsilon^2(t-1) | Source: IMF WEO April 2025",
       x = "Year", y = "Conditional Volatility (%)") +
  theme_ts() +
  theme(legend.position = "bottom")
sp(p13, "fig13_ewma_volatility.png", 14, 7)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# SECTION 10 · GARCH MODELLING
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

cat("\n══ SECTION 10: GARCH MODELLING ══\n")
cat("WARNING: GARCH fitted on ~33 annual observations (100+ recommended).\n")
cat("Results are ILLUSTRATIVE. Treat parameter estimates with caution.\n\n")

garch_results <- lapply(model_iso, function(iso) {
  cat("\n  Country:", iso_to_name(iso), "\n")
  res      <- model_results[[iso]]
  s        <- as.numeric(get_ts(iso))
  d_ord    <- res$full_model$arma[6]
  s_diff   <- if (d_ord > 0) diff(s, differences = d_ord) else s

  g11_n <- tryCatch(fGarch::garchFit(~garch(1,1), data = s_diff, trace = FALSE, cond.dist = "norm"),
                    error = function(e) { cat("   GARCH(1,1)-norm failed\n"); NULL })
  g11_t <- tryCatch(fGarch::garchFit(~garch(1,1), data = s_diff, trace = FALSE, cond.dist = "std"),
                    error = function(e) { cat("   GARCH(1,1)-std failed\n"); NULL })
  g12_t <- tryCatch(fGarch::garchFit(~garch(1,2), data = s_diff, trace = FALSE, cond.dist = "std"),
                    error = function(e) NULL)

  models_fit <- Filter(Negate(is.null),
                       list("GARCH(1,1)-N" = g11_n, "GARCH(1,1)-t" = g11_t, "GARCH(1,2)-t" = g12_t))

  best_name <- if (length(models_fit) > 0) {
    aics  <- sapply(models_fit, function(m) tryCatch(m@fit$ics["AIC"], error = function(e) Inf))
    sub("[.]AIC$", "", names(which.min(aics)))
  } else "None"
  best_mod <- models_fit[[best_name]]

  cat("   Best:", best_name, "\n")
  if (!is.null(best_mod)) {
    cf <- coef(best_mod)
    persist <- tryCatch(cf["alpha1"] + cf["beta1"], error = function(e) NA)
    cat(sprintf("   omega=%.5f  alpha1=%.4f  beta1=%.4f  alpha+beta=%.4f\n",
                tryCatch(cf["omega"],  error = function(e) NA),
                tryCatch(cf["alpha1"], error = function(e) NA),
                tryCatch(cf["beta1"],  error = function(e) NA),
                tryCatch(persist, error = function(e) NA)))
  }

  list(iso = iso, name = iso_to_name(iso), s_diff = s_diff,
       g11_n = g11_n, g11_t = g11_t, g12_t = g12_t,
       best_name = best_name, best_mod = best_mod, models = models_fit)
})
names(garch_results) <- model_iso

# GARCH model comparison
garch_comp <- do.call(rbind, lapply(model_iso, function(iso) {
  g <- garch_results[[iso]]
  do.call(rbind, Filter(Negate(is.null), lapply(names(g$models), function(mn) {
    m <- g$models[[mn]]
    tryCatch({
      ics <- m@fit$ics; cf <- coef(m)
      data.frame(Country = g$name, Model = mn,
                 AIC = round(ics["AIC"], 3), BIC = round(ics["BIC"], 3),
                 omega  = round(tryCatch(cf["omega"],  error = function(e) NA), 6),
                 alpha1 = round(tryCatch(cf["alpha1"], error = function(e) NA), 4),
                 beta1  = round(tryCatch(cf["beta1"],  error = function(e) NA), 4),
                 Persistence = round(tryCatch(cf["alpha1"] + cf["beta1"], error = function(e) NA), 4))
    }, error = function(e) NULL)
  })))
}))
print(garch_comp)
write.csv(garch_comp, paste0(OUT, "garch_model_comparison.csv"), row.names = FALSE)

# ── FIG 14: GARCH conditional volatility ─────────────────────────────────────
cat("[FIG 14] GARCH conditional volatility\n")
png(paste0(OUT, "fig14_garch_volatility.png"),
    width = 2000, height = 2000, res = 130, bg = "white")
par(mfrow = c(length(model_iso), 2), mar = c(3.5, 4, 3, 1.5),
    oma = c(0, 0, 4, 0), bg = "white")
for (iso in model_iso) {
  g   <- garch_results[[iso]]; bm <- g$best_mod
  d_o <- model_results[[iso]]$full_model$arma[6]
  yr  <- (1990 + ifelse(d_o > 0, 1, 0)):2023
  if (!is.null(bm) && length(bm@sigma.t) == length(yr)) {
    h_t <- bm@sigma.t
    plot(yr, g$s_diff, type = "l", col = "#1f77b4", lwd = 1.5,
         main = paste0(g$name, " - Differenced Series & +/-2sigma"),
         xlab = "", ylab = "Delta Inflation (%)", las = 1)
    polygon(c(yr, rev(yr)),
            c(as.numeric(bm@fitted) + 2*h_t, rev(as.numeric(bm@fitted) - 2*h_t)),
            col = adjustcolor("#e74c3c", 0.15), border = NA)
    lines(yr, as.numeric(bm@fitted) + 2*h_t, col = "#e74c3c", lwd = 1.2, lty = 2)
    lines(yr, as.numeric(bm@fitted) - 2*h_t, col = "#e74c3c", lwd = 1.2, lty = 2)
    legend("topleft", legend = c("Delta Inflation", "+/-2 Cond. SD"),
           col = c("#1f77b4","#e74c3c"), lty = c(1, 2), lwd = 1.5, cex = 0.75, bty = "n")
    plot(yr, h_t, type = "l", col = "#e74c3c", lwd = 2,
         main = paste0("GARCH Cond. Volatility - ", g$best_name),
         xlab = "", ylab = "sigma(t) (%)", las = 1)
    polygon(c(yr[1], yr, yr[length(yr)]), c(0, h_t, 0),
            col = adjustcolor("#e74c3c", 0.15), border = NA)
    abline(h = mean(h_t), col = "#2ca02c", lwd = 1.5, lty = 2)
    cf_g    <- coef(bm)
    persist <- tryCatch(round(cf_g["alpha1"] + cf_g["beta1"], 3), error = function(e) "?")
    text(yr[4], max(h_t)*0.92, paste0("alpha+beta=", persist), cex = 0.9, font = 2, col = "#9467bd")
  } else {
    plot.new(); text(0.5, 0.5, paste("GARCH convergence\nfailed:", g$name), cex = 1.1, col = "red")
    plot.new()
  }
}
mtext("ARIMA-GARCH(1,1) Conditional Volatility - Inflation Residuals",
      outer = TRUE, cex = 1.3, font = 2, col = "#1a1a2e")
dev.off()
cat("  Saved: fig14_garch_volatility.png\n")

# ── FIG 15: GARCH standardised residuals ─────────────────────────────────────
cat("[FIG 15] GARCH standardised residuals\n")
png(paste0(OUT, "fig15_garch_std_resid.png"),
    width = 2000, height = 1800, res = 130, bg = "white")
par(mfrow = c(length(model_iso), 3), mar = c(3, 3.5, 2.5, 1),
    oma = c(0, 0, 4, 0), bg = "white")
for (iso in model_iso) {
  g   <- garch_results[[iso]]; bm <- g$best_mod
  d_o <- model_results[[iso]]$full_model$arma[6]
  yr  <- (1990 + ifelse(d_o > 0, 1, 0)):2023
  if (!is.null(bm) && length(bm@sigma.t) == length(yr)) {
    sr <- bm@residuals / bm@sigma.t
    plot(yr, sr, type = "l", col = "#1f77b4", lwd = 1.5,
         main = paste0(g$name, "\nStd. Residuals"), xlab = "", ylab = "z(t)", las = 1)
    abline(h = c(-2, 0, 2), col = c("#e74c3c","grey50","#e74c3c"), lty = c(2, 1, 2))
    acf(sr^2, lag.max = 12, main = "ACF of z^2(t)",
        col = "#9467bd", lwd = 2, ci.col = "#e74c3c", xlab = "Lag")
    qqnorm(sr, main = "Normal Q-Q (Std Resid)", pch = 20, col = "#2ca02c", cex = 0.9)
    qqline(sr, col = "red", lwd = 2)
  } else {
    for (k in 1:3) { plot.new(); text(0.5, 0.5, paste("N/A -", g$name), col = "red") }
  }
}
mtext("GARCH Standardised Residual Diagnostics",
      outer = TRUE, cex = 1.3, font = 2, col = "#1a1a2e")
dev.off()
cat("  Saved: fig15_garch_std_resid.png\n")

# GARCH persistence table
persist_df <- do.call(rbind, Filter(Negate(is.null), lapply(model_iso, function(iso) {
  g <- garch_results[[iso]]; bm <- g$best_mod
  if (is.null(bm)) return(NULL)
  cf      <- coef(bm)
  omega   <- tryCatch(cf["omega"],  error = function(e) NA)
  alpha1  <- tryCatch(cf["alpha1"], error = function(e) 0)
  beta1   <- tryCatch(cf["beta1"],  error = function(e) 0)
  persist <- alpha1 + beta1
  lr_var  <- if (!is.na(persist) && persist < 1) omega / (1 - persist) else NA
  data.frame(Country = g$name, Model = g$best_name,
             omega = round(omega, 5), alpha1 = round(alpha1, 4), beta1 = round(beta1, 4),
             Persistence = round(persist, 4), LR_Variance = round(lr_var, 4),
             LR_Vol = round(sqrt(lr_var), 4))
})))
print(persist_df)
write.csv(persist_df, paste0(OUT, "garch_persistence.csv"), row.names = FALSE)

if (nrow(persist_df) > 0) {
  p16 <- ggplot(persist_df, aes(reorder(Country, -Persistence), Persistence, fill = Persistence)) +
    geom_col(alpha = 0.85, colour = "white", width = 0.6) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "#e74c3c", linewidth = 1.1) +
    geom_text(aes(label = paste0("a+b=", Persistence, "\nLR s=", LR_Vol, "%")),
              vjust = -0.3, size = 3.2, fontface = "bold", colour = "#1a1a2e") +
    scale_fill_gradient2(low = "#27ae60", mid = "#f39c12", high = "#e74c3c",
                         midpoint = 0.85, name = "Persistence") +
    scale_y_continuous(limits = c(0, 1.4)) +
    labs(title    = "GARCH Volatility Persistence (alpha + beta) by Country",
         subtitle = "Values close to 1 = long-memory volatility | >1 = explosive",
         x = NULL, y = "Persistence (alpha + beta)",
         caption = "Long-run var = omega / (1 - alpha - beta) | Source: IMF WEO April 2025") +
    theme_ts() +
    theme(legend.position = "right")
  sp(p16, "fig16_garch_persistence.png", 12, 7)
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# SECTION 11 · FORECAST EVALUATION
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

gc()
cat("\n══ SECTION 11: FORECAST EVALUATION ══\n")
cat("PRIMARY   evaluation: pre-COVID test 2015-2019 (headline forecast skill)\n")
cat("STRESS    evaluation: COVID-era  test 2020-2023 (structural-break context)\n\n")
cat("Naive benchmark: random walk (I(1)) or historical mean (I(0)) per ADF result.\n\n")

# ── Adaptive naive benchmark ──────────────────────────────────────────────────
naive_forecast <- function(train_series, h, is_stationary) {
  if (is_stationary)
    rep(mean(as.numeric(train_series), na.rm = TRUE), h)  # mean for I(0)
  else
    rep(as.numeric(tail(train_series, 1)), h)              # random walk for I(1)
}

# Determine stationarity per ISO from ADF on 1990-2023
stationary_flag <- setNames(
  sapply(model_iso, function(iso) {
    sv <- as.numeric(na.omit(get_ts(iso)))
    if (length(sv) < 10) return(FALSE)
    p  <- tryCatch(adf.test(sv)$p.value, error = function(e) 1)
    p < 0.05
  }), model_iso
)

# ── Evaluation helper (runs on any train/test window) ─────────────────────────
eval_on_window <- function(iso, s_train, s_test, label) {
  s_actual   <- as.numeric(s_test)
  h_test     <- length(s_actual)
  is_stat    <- stationary_flag[iso]
  naive_fc   <- naive_forecast(s_train, h_test, is_stat)
  naive_rmse <- sqrt(mean((naive_fc - s_actual)^2, na.rm = TRUE))
  naive_type <- if (is_stat) "Mean (I(0))" else "Random Walk (I(1))"

  cand_models <- as.character(model_results[[iso]]$cand_table$Model)
  eval_models <- Filter(Negate(is.null), lapply(cand_models, function(mn) {
    ord <- suppressWarnings(
      as.integer(unlist(strsplit(gsub("ARIMA\\(|\\)", "", mn), ","))))
    if (any(is.na(ord))) return(NULL)
    m  <- safe_arima(s_train, ord)
    if (is.null(m)) return(NULL)
    fc <- tryCatch(as.numeric(forecast(m, h = h_test)$mean),
                   error = function(e) NULL)
    if (is.null(fc)) return(NULL)
    cbind(data.frame(Model = mn), eval_metrics(s_actual, fc))
  }))

  eval_df <- do.call(rbind, eval_models)
  if (is.null(eval_df) || nrow(eval_df) == 0) {
    auto_fc <- tryCatch(
      as.numeric(forecast(model_results[[iso]]$model, h = h_test)$mean),
      error = function(e) rep(NA, h_test))
    eval_df <- cbind(data.frame(Model = "auto"), eval_metrics(s_actual, auto_fc))
  } else {
    eval_df <- eval_df[order(eval_df$RMSE), ]
  }

  best_ord <- suppressWarnings(
    as.integer(unlist(strsplit(gsub("ARIMA\\(|\\)", "", eval_df$Model[1]), ","))))
  best_m  <- safe_arima(s_train, best_ord)
  best_fc <- if (!is.null(best_m))
    as.numeric(forecast(best_m, h = h_test)$mean) else rep(NA, h_test)

  skill <- round(1 - eval_df$RMSE[1] / naive_rmse, 3)
  cat(sprintf("  %-18s [%s] RMSE=%.3f | Naive(%s)=%.3f | Skill=%+.3f\n",
              iso_to_name(iso), label,
              eval_df$RMSE[1], naive_type, naive_rmse, skill))

  list(iso = iso, name = iso_to_name(iso), label = label,
       actuals = s_actual, years = as.integer(time(s_test)),
       best_arima_fc = best_fc, best_arima_nm = eval_df$Model[1],
       naive_fc = naive_fc, naive_type = naive_type,
       eval_table = eval_df, naive_rmse = naive_rmse,
       best_rmse = eval_df$RMSE[1], skill = skill)
}

# ── PRIMARY evaluation (pre-COVID 2015-2019) ──────────────────────────────────
cat("--- PRIMARY evaluation (pre-COVID 2015-2019) ---\n")
eval_pre <- setNames(lapply(model_iso, function(iso) {
  res <- model_results[[iso]]
  eval_on_window(iso, res$train_pre, res$test_pre, "pre-COVID")
}), model_iso)

# ── STRESS evaluation (COVID-era 2020-2023) ───────────────────────────────────
cat("\n--- STRESS evaluation (COVID-era 2020-2023) ---\n")
eval_cov <- setNames(lapply(model_iso, function(iso) {
  res <- model_results[[iso]]
  eval_on_window(iso, res$train_cov, res$test_cov, "COVID-era")
}), model_iso)

# Alias used by dashboard
eval_results <- eval_pre

# ── FIG 17: Combined evaluation chart ────────────────────────────────────────
cat("[FIG 17] Forecast evaluation chart\n")

hist_rows <- do.call(rbind, lapply(model_iso, function(iso) {
  s <- as.numeric(get_ts(iso, 1990, 2014))
  yr <- 1990:(1990 + length(s) - 1)
  data.frame(Year = yr, Value = s, Series = "Historical (train)",
             ISO = iso, Country = iso_to_name(iso))
}))

comp_rows <- do.call(rbind, lapply(model_iso, function(iso) {
  ep  <- eval_pre[[iso]];  ec <- eval_cov[[iso]]
  m_f <- model_results[[iso]]$full_model
  fc3 <- as.numeric(forecast(m_f, h = 3)$mean)
  imf_val <- df_long$Inflation_Rate[df_long$ISO == iso & df_long$Year == 2024]
  imf_row <- if (length(imf_val) > 0 && !is.na(imf_val[1]))
    data.frame(Year = 2024, Value = imf_val[1], Series = "2024 IMF estimate",
               ISO = iso, Country = ep$name) else NULL
  rbind(
    data.frame(Year = ep$years, Value = ep$actuals,
               Series = "Actual 2015-2019", ISO = iso, Country = ep$name),
    data.frame(Year = ep$years, Value = ep$best_arima_fc,
               Series = "ARIMA (primary)",   ISO = iso, Country = ep$name),
    data.frame(Year = ep$years, Value = ep$naive_fc,
               Series = "Naive (primary)",   ISO = iso, Country = ep$name),
    data.frame(Year = ec$years, Value = ec$actuals,
               Series = "Actual 2020-2023",  ISO = iso, Country = ep$name),
    data.frame(Year = ec$years, Value = ec$best_arima_fc,
               Series = "ARIMA (stress)",    ISO = iso, Country = ep$name),
    data.frame(Year = 2024:2026, Value = fc3,
               Series = "ARIMA 2024-26",     ISO = iso, Country = ep$name),
    imf_row
  )
}))
plot_df <- rbind(hist_rows, comp_rows)

pal17 <- c("Historical (train)" = "#aaaaaa", "Actual 2015-2019" = "#1a1a2e",
           "ARIMA (primary)"    = "#e74c3c", "Naive (primary)"  = "#ff7f0e",
           "Actual 2020-2023"   = "#555577", "ARIMA (stress)"   = "#c0392b",
           "ARIMA 2024-26"      = "#e74c3c", "2024 IMF estimate"= "#2ca02c")
lty17  <- c("Historical (train)" = "solid",   "Actual 2015-2019" = "solid",
            "ARIMA (primary)"    = "dashed",  "Naive (primary)"  = "dotdash",
            "Actual 2020-2023"   = "solid",   "ARIMA (stress)"   = "longdash",
            "ARIMA 2024-26"      = "dashed",  "2024 IMF estimate"= "blank")

p17 <- ggplot(plot_df[plot_df$Year >= 2010, ],
              aes(Year, Value, colour = Series, linetype = Series)) +
  annotate("rect", xmin=2014.5, xmax=2019.5, ymin=-Inf, ymax=Inf,
           fill="#e8f4fd", alpha=0.55) +
  annotate("rect", xmin=2019.5, xmax=2023.5, ymin=-Inf, ymax=Inf,
           fill="#fdf0e8", alpha=0.55) +
  annotate("rect", xmin=2023.5, xmax=2026.5, ymin=-Inf, ymax=Inf,
           fill="#e8f8f5", alpha=0.55) +
  annotate("text", x=2017,   y=Inf, label="Primary\n2015-19",
           vjust=1.5, size=2.4, colour="#2980b9") +
  annotate("text", x=2021.5, y=Inf, label="COVID\nstress",
           vjust=1.5, size=2.4, colour="#e67e22") +
  annotate("text", x=2025,   y=Inf, label="Forecast\nhorizon",
           vjust=1.5, size=2.4, colour="#27ae60") +
  geom_line(linewidth = 0.9) +
  geom_point(data = plot_df[plot_df$Series == "2024 IMF estimate", ],
             size = 3.5, shape = 15) +
  facet_wrap(~Country, scales = "free_y", nrow = 2) +
  scale_colour_manual(values = pal17, name = "") +
  scale_linetype_manual(values = lty17, name = "") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_x_continuous(breaks = seq(2010, 2026, 3)) +
  labs(title    = "ARIMA Forecast Evaluation: Pre-COVID (primary) vs COVID-era (stress)",
       subtitle = "Blue = primary eval 2015-19 | Orange = COVID stress 2020-23 | Green pt = 2024 IMF estimate only",
       x = "Year", y = "Inflation Rate (%)",
       caption = "Naive: random walk for I(1) series, mean for I(0) | Source: IMF WEO April 2025") +
  theme_ts() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))
sp(p17, "fig17_forecast_evaluation.png", 18, 11)

# ── Accuracy metrics tables ───────────────────────────────────────────────────
make_metrics_row <- function(e, window_label) {
  et <- if (nrow(e$eval_table) > 0) e$eval_table[1, ] else
    data.frame(RMSE = NA, MAE = NA, sMAPE = NA, Bias = NA)
  data.frame(Country = e$name, ISO = e$iso, Window = window_label,
             ARIMA_RMSE = et$RMSE, ARIMA_MAE = et$MAE,
             ARIMA_sMAPE = if ("sMAPE" %in% names(et)) et$sMAPE else NA,
             Naive_RMSE = round(e$naive_rmse, 3), Naive_type = e$naive_type,
             Skill_Score = e$skill, ARIMA_Bias = et$Bias)
}
metrics_pre <- do.call(rbind, lapply(model_iso, function(iso)
  make_metrics_row(eval_pre[[iso]], "Pre-COVID (2015-2019)")))
metrics_cov <- do.call(rbind, lapply(model_iso, function(iso)
  make_metrics_row(eval_cov[[iso]], "COVID-era (2020-2023)")))
metrics_df <- rbind(metrics_pre, metrics_cov)
cat("\nForecast accuracy:\n"); print(metrics_df)
write.csv(metrics_df, paste0(OUT, "forecast_accuracy.csv"), row.names = FALSE)

# ── FIG 18: RMSE bar — primary window ────────────────────────────────────────
p18 <- tidyr::pivot_longer(metrics_pre, cols = c(ARIMA_RMSE, Naive_RMSE),
                            names_to = "Method", values_to = "RMSE") %>%
  dplyr::mutate(Method = gsub("_RMSE","", Method)) %>%
  ggplot(aes(reorder(Country, -RMSE), RMSE, fill = Method)) +
  geom_col(position = "dodge", alpha = 0.85, colour = "white", width = 0.7) +
  geom_text(aes(label = round(RMSE, 2)), position = position_dodge(0.7),
            vjust = -0.4, size = 3.1) +
  scale_fill_manual(values = c("ARIMA" = "#1f77b4", "Naive" = "#e74c3c"),
                    name = "Model") +
  labs(title    = "Forecast RMSE: ARIMA vs. Naive — Pre-COVID Window (2015-2019)",
       subtitle = "Primary evaluation: avoids COVID structural break | lower = better",
       x = NULL, y = "RMSE (pp)") +
  theme_ts() + theme(legend.position = "top")
sp(p18, "fig18_rmse_comparison.png", 12, 7)

# ── FIG 19: Skill scores — both windows ──────────────────────────────────────
metrics_df$Window <- factor(metrics_df$Window,
  levels = c("Pre-COVID (2015-2019)", "COVID-era (2020-2023)"))

p19 <- ggplot(metrics_df, aes(reorder(Country, -Skill_Score), Skill_Score,
                               fill = Skill_Score > 0)) +
  geom_col(alpha = 0.82, colour = "white", width = 0.65) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
  geom_text(aes(label = paste0(round(Skill_Score*100, 1), "%")),
            vjust = ifelse(metrics_df$Skill_Score >= 0, -0.4, 1.2),
            size = 3.0, fontface = "bold") +
  facet_wrap(~Window, nrow = 2) +
  scale_fill_manual(values = c("TRUE" = "#27ae60", "FALSE" = "#e74c3c"),
                    labels = c("TRUE" = "Beats naive", "FALSE" = "Naive wins"),
                    name = "") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(-2.5, 1.0)) +
  labs(title    = "ARIMA Skill Score vs. Adaptive Naive Benchmark",
       subtitle = "Skill = 1 - RMSE(ARIMA)/RMSE(Naive) | >0 means ARIMA adds value",
       x = NULL, y = "Skill Score") +
  theme_ts() + theme(legend.position = "bottom")
sp(p19, "fig19_skill_scores.png", 14, 10)

# ── ARIMA vs 2024 IMF estimate (single year available) ───────────────────────
div_df <- do.call(rbind, Filter(Negate(is.null), lapply(model_iso, function(iso) {
  m_f     <- model_results[[iso]]$full_model
  fc3     <- as.numeric(forecast(m_f, h = 3)$mean)
  imf_val <- df_long$Inflation_Rate[df_long$ISO == iso & df_long$Year == 2024]
  if (length(imf_val) == 0 || is.na(imf_val[1])) return(NULL)
  data.frame(Year = 2024, ARIMA = round(fc3[1], 2),
             IMF_estimate = round(imf_val[1], 2),
             Diff = round(fc3[1] - imf_val[1], 2),
             Country = iso_to_name(iso), ISO = iso)
})))
print(div_df)
write.csv(div_df, paste0(OUT, "arima_vs_imf_2024.csv"), row.names = FALSE)

if (nrow(div_df) > 0) {
  p20 <- ggplot(div_df, aes(Country, Diff, fill = Diff > 0)) +
    geom_col(alpha = 0.80, colour = "white", width = 0.6) +
    geom_hline(yintercept = 0, linewidth = 1, colour = "#333333") +
    geom_text(aes(label = paste0(ifelse(Diff > 0, "+", ""), round(Diff, 2), " pp")),
              vjust = ifelse(div_df$Diff >= 0, -0.5, 1.3),
              size = 3.5, fontface = "bold") +
    scale_fill_manual(values = c("TRUE" = "#e74c3c", "FALSE" = "#2980b9"),
                      labels = c("TRUE" = "ARIMA > IMF", "FALSE" = "ARIMA < IMF"),
                      name = "") +
    labs(title    = "ARIMA 2024 Forecast vs IMF Near-Final Estimate",
         subtitle = "Only 2024 is available in this dataset — 2025-2026 IMF data not present",
         x = NULL, y = "Difference (ARIMA - IMF, pp)") +
    theme_ts() + theme(legend.position = "bottom")
  sp(p20, "fig20_arima_vs_imf_divergence.png", 12, 6)
}

# ── FIG 21: Forecast error bars — primary window ─────────────────────────────
cat("[FIG 21] Forecast error bars\n")
png(paste0(OUT, "fig21_forecast_errors.png"),
    width = 1800, height = 1600, res = 100, bg = "white")
par(mfrow = c(ceiling(length(model_iso)/2), 2),
    mar = c(4, 4, 3, 2), oma = c(0, 0, 3.5, 0), bg = "white")
for (iso in model_iso) {
  ep      <- eval_pre[[iso]]
  arima_e <- ep$best_arima_fc - ep$actuals
  naive_e <- ep$naive_fc      - ep$actuals
  barplot(rbind(arima_e, naive_e), beside = TRUE,
          names.arg = ep$years, col = c("#1f77b4","#e74c3c"), border = NA,
          main = paste0(ep$name, "\nForecast Errors (Primary: 2015-2019)"),
          ylab = "Error (pp)", xlab = "Year", las = 1)
  abline(h = 0, lty = 1, col = "grey40")
  legend("topright",
         legend = c(paste0("ARIMA RMSE=", round(ep$best_rmse, 2)),
                    paste0(ep$naive_type, " RMSE=", round(ep$naive_rmse, 2))),
         fill = c("#1f77b4","#e74c3c"), bty = "n", cex = 0.8)
}
mtext("ARIMA vs Naive Forecast Errors — Primary Evaluation (2015-2019)",
      outer = TRUE, cex = 1.2, font = 2, col = "#1a1a2e")
dev.off()
cat("  Saved: fig21_forecast_errors.png\n")

# SECTION 12 · FINAL VISUAL DASHBOARD
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

cat("\n══ SECTION 12: FINAL DASHBOARD ══\n")

# Panel A: Avg inflation trends — Advanced vs Emerging
df_adv_d <- aggregate(Inflation_Rate ~ Year,
  df_long[df_long$ISO %in% adv_isos & df_long$Year >= 1990 & df_long$Year <= 2024, ], mean)
df_em_d  <- aggregate(Inflation_Rate ~ Year,
  df_long[df_long$ISO %in% em_isos  & df_long$Year >= 1990 & df_long$Year <= 2024, ], mean)
df_all_d <- aggregate(Inflation_Rate ~ Year,
  df_long[df_long$ISO %in% all_iso  & df_long$Year >= 1990 & df_long$Year <= 2024, ], mean)

df_lines_d <- rbind(
  data.frame(Year = df_adv_d$Year, Value = df_adv_d$Inflation_Rate, Group = "Advanced Economies"),
  data.frame(Year = df_em_d$Year,  Value = df_em_d$Inflation_Rate,  Group = "Emerging Markets"),
  data.frame(Year = df_all_d$Year, Value = df_all_d$Inflation_Rate, Group = "All 15 Countries")
)

pDA <- ggplot(df_lines_d, aes(Year, Value, colour = Group)) +
  event_rects() +
  scale_fill_identity() +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("Advanced Economies" = "#1f77b4",
                                 "Emerging Markets" = "#d62728",
                                 "All 15 Countries" = "#333333"), name = "Group") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_x_continuous(breaks = seq(1990, 2024, 6)) +
  labs(title = "Global Inflation Timeline (1990-2024)", x = NULL, y = "Rate (%)") +
  theme_ts() +
  theme(legend.position = "bottom", legend.text = element_text(size = 7.5))

# Panel B: 2022 peak ranking
p2022 <- df_long[df_long$Year == 2022 & df_long$ISO %in% all_iso, ]
p2022 <- merge(p2022, setNames(nm, c("ISO","CountryName")), by = "ISO")
p2022 <- p2022[!is.na(p2022$Inflation_Rate), ]

pDB <- ggplot(p2022, aes(reorder(CountryName, Inflation_Rate), Inflation_Rate, fill = ISO)) +
  geom_col(alpha = 0.85, colour = "white") +
  geom_text(aes(label = paste0(round(Inflation_Rate, 0), "%")), hjust = -0.1, size = 2.8) +
  scale_fill_manual(values = PALETTE, guide = "none") +
  coord_flip() +
  scale_y_continuous(limits = c(0, max(p2022$Inflation_Rate, na.rm = TRUE) * 1.2)) +
  labs(title = "2022 Peak Inflation Ranking", x = NULL, y = "%") +
  theme_ts()

# Panel C: ADF p-values
pDC <- p05a + labs(title = "ADF Stationarity Test", subtitle = NULL) +
  theme(legend.position = "bottom", legend.text = element_text(size = 7.5))

# Panel D: Correlation heatmap
pDD <- p09c + labs(title = "Inflation Correlations (2000-2023)") +
  theme(axis.text = element_text(size = 6))

cat("[FIG 22] Final dashboard\n")
dash <- grid.arrange(pDA, pDB, pDC, pDD, ncol = 2, nrow = 2,
  top = textGrob(
    "IMF WEO - Inflation Analysis Dashboard | April 2025 (CPI % Change)",
    gp = gpar(fontsize = 14, fontface = "bold", col = "#1a1a2e")
  )
)
ggsave(paste0(OUT, "fig22_dashboard.png"), plot = dash,
       width = 20, height = 16, dpi = 150, bg = "white")
cat("  Saved: fig22_dashboard.png\n")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# FINAL SUMMARY
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("ANALYSIS COMPLETE - IMF WEO Inflation Time Series\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

cat("-- Descriptive Statistics --\n")
print(desc_df[, c("Country","Mean","Median","SD","Max")])

cat("\n-- Stationarity (ADF) --\n")
if ("ADF_pval" %in% names(stat_df))
  print(stat_df[, c("Country","ADF_pval","Stationary")])

cat("\n-- ARIMA Residual Diagnostics --\n")
print(hyp_df[, c("Country","No_Autocorr_LB","Normal_Residuals_SW")])

cat("\n-- GARCH Persistence --\n")
if (nrow(persist_df) > 0)
  print(persist_df[, c("Country","Model","alpha1","beta1","Persistence","LR_Vol")])

cat("\n-- Forecast Accuracy (Primary: Pre-COVID 2015-2019) --\n")
print(metrics_pre[, c("Country","ARIMA_RMSE","Naive_RMSE","Skill_Score","Naive_type")])

cat("\n-- ARIMA vs 2024 IMF Estimate --\n")
if (exists("div_df") && nrow(div_df) > 0)
  print(div_df[, c("Country","ARIMA","IMF_estimate","Diff")])

cat("\n-- Output files --\n")
output_files <- list.files(OUT, full.names = FALSE)
cat("Folder    :", OUT, "\n")
cat("File count:", length(output_files), "\n")
for (f in sort(output_files)) cat("  *", f, "\n")
cat("\nDone.\n")
