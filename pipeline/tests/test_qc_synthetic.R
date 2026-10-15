# Synthetic check for qc-outlier-fixes.
# Run from repo root:  Rscript pipeline/tests/test_qc_synthetic.R
# Builds a clean seasonal EVI series, injects a single-date dip, a 3-obs dip
# and a run of exact-zero observations, runs DoPhenologyPlanet with the new
# defaults and with legacy settings, and checks that (a) the bad points are
# removed from the raw series and (b) the summer spline is not pulled down.
suppressWarnings(suppressMessages(library(jsonlite)))
source("pipeline/wetlsp-functions.R")
params <- fromJSON("pipeline/wetlsp-parameters.json")
params$phenology_parameters$dormStart <- as.Date(params$phenology_parameters$dormStart)
params$phenology_parameters$dormEnd   <- as.Date(params$phenology_parameters$dormEnd)

set.seed(1)
dates <- sort(as.Date("2023-01-01") + sort(sample(0:(3*365-1), 330)))
doy   <- as.numeric(format(dates, "%j"))
evi_true <- 0.15 + 0.45 * exp(-((doy - 200) / 45)^2)
evi <- evi_true + rnorm(length(dates), 0, 0.015)

# Invert EVI to bands with fixed blue/red so CheckSpike sees plausible values
blue <- rep(0.04, length(dates)); red <- rep(0.05, length(dates))
# evi = 2.5(n-r)/(n+6r-7.5b+1)  ->  n = (evi*(6r-7.5b+1) + 2.5r)/(2.5-evi)
nir_from <- function(e, i = seq_along(e)) (e*(6*red[i] - 7.5*blue[i] + 1) + 2.5*red[i]) / (2.5 - e)
nir <- nir_from(evi)

y24 <- format(dates, "%Y") == "2024"
summer <- which(y24 & doy > 180 & doy < 230)
single <- summer[3]
multi  <- summer[8:10]
zeros  <- which(dates >= as.Date("2024-07-21") & dates <= as.Date("2024-08-07"))
zeros  <- setdiff(zeros, c(single, multi))
nir[single] <- nir_from(0.10, single)
nir[multi]  <- nir_from(c(0.12, 0.10, 0.13), multi)
# exact zeros as delivered by a fill/zero band (all bands 0)
blue[zeros] <- 0; red[zeros] <- 0; nir[zeros] <- 0

run <- function(p) {
  acc <- new.env(); acc$raw <- NULL; acc$qc <- NULL; acc$sd <- NULL; acc$se <- NULL
  sink_fn <- function(pix_meta, dates_raw, evi_raw, pred_dates, evi_spline, qc_raw = NULL) {
    if (!is.null(dates_raw)) { acc$raw <- evi_raw; acc$qc <- qc_raw }
    if (!is.null(pred_dates)) { acc$sd <- c(acc$sd, as.Date(pred_dates)); acc$se <- c(acc$se, evi_spline) }
  }
  DoPhenologyPlanet(blue*1e4, rep(500, length(dates)), red*1e4, nir*1e4, dates,
                    2023:2025, p, 0, ts_sink = sink_fn, pix_meta = list(cell = 1))
  acc
}

new <- run(params)
legacy_p <- params
legacy_p$phenology_parameters[c("secondSplineQC","envelopeReweight")] <- list(FALSE, FALSE)
legacy_p$phenology_parameters$dormantFloorMode <- "all"
legacy_p$phenology_parameters$maskNegativeValues <- FALSE
legacy_p$phenology_parameters[c("eviValidMin","eviValidMax")] <- list(-Inf, Inf)
legacy_p$phenology_parameters$MADspikeThresh <- 1e6   # approximates the old, effectively-disabled madTest
legacy <- run(legacy_p)

bad <- c(single, multi, zeros)
cat("Injected bad obs:", length(bad), "\n")
cat("New:    bad obs remaining non-NA in raw =", sum(!is.na(new$raw[bad])), "\n")
cat("        qc codes on bad obs:", paste(table(new$qc[bad]), names(table(new$qc[bad])), sep="x", collapse=" "), "\n")
cat("        good obs removed =", sum(is.na(new$raw[-bad])), "of", length(dates) - length(bad), "\n")
cat("Legacy: bad obs remaining non-NA in raw =", sum(!is.na(legacy$raw[bad])), " values:", paste(round(legacy$raw[bad],3), collapse=" "), "\n")

win <- which(as.Date(new$sd, origin="1970-01-01") >= as.Date("2024-07-10") &
             as.Date(new$sd, origin="1970-01-01") <= as.Date("2024-08-20"))
truth <- function(d) { j <- as.numeric(format(d, "%j")); 0.15 + 0.45*exp(-((j-200)/45)^2) }
sd_new <- as.Date(new$sd, origin="1970-01-01")
sd_leg <- as.Date(legacy$sd, origin="1970-01-01")
wl <- which(sd_leg >= as.Date("2024-07-10") & sd_leg <= as.Date("2024-08-20"))
bias_new <- mean(new$se[win] - truth(sd_new[win]))
bias_leg <- mean(legacy$se[wl] - truth(sd_leg[wl]))
cat(sprintf("Mean spline bias Jul10-Aug20 2024: new = %+.3f, legacy = %+.3f\n", bias_new, bias_leg))

stopifnot(sum(!is.na(new$raw[bad])) == 0)
stopifnot(abs(bias_new) < 0.03, abs(bias_new) < abs(bias_leg))
cat("PASS\n")

# Case 2: second-spline QC alone (madTest effectively off) must still catch the dips
pp <- params; pp$phenology_parameters$MADspikeThresh <- 1e6
ss <- run(pp)
cat("SS-only: bad obs remaining =", sum(!is.na(ss$raw[bad])),
    "| good obs removed =", sum(is.na(ss$raw[-bad])), "\n")
stopifnot(sum(!is.na(ss$raw[bad])) == 0)
cat("PASS (second-spline only)\n")
