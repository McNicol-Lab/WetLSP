# Unit test for UDM2 scene matching (step 01).
# Run from repo root:  Rscript pipeline/tests/test_udm2_matching.R
# Same-date, multi-scene PlanetScope filenames: each SR scene must get its own
# UDM2, and a scene without one must get NA (not another scene's mask).
suppressWarnings(suppressMessages(source("pipeline/wetlsp-functions.R")))

ok <- TRUE
check <- function(cond, msg) {
  cat(if (cond) "PASS " else "FAIL ", msg, "\n", sep = "")
  if (!cond) ok <<- FALSE
}

d <- "ord1/PSScene/"
sr <- paste0(d, c(
  "20201022_191920_74_227a_3B_AnalyticMS_SR_8b_clip.tif",
  "20201022_191923_12_227a_3B_AnalyticMS_SR_8b_clip.tif",
  "20201022_181500_1049_3B_AnalyticMS_SR_harmonized_clip.tif",   # no UDM2
  "20201205_191349_38_240f_3B_AnalyticMS_SR_8b_clip.tif"))
udm <- paste0(d, c(
  "20201205_191349_38_240f_3B_udm2_clip.tif",
  "20201022_191923_12_227a_3B_udm2_clip.tif",
  "20201022_191920_74_227a_3B_udm2_clip.tif"))

m <- MatchUDM2(sr, udm)
check(identical(m, c(3L, 2L, NA_integer_, 1L)), "each same-date scene gets its own UDM2; missing -> NA")
check(identical(SceneIdFromFile(sr[1]), "20201022_191920_74_227a"), "scene ID = stem before _3B_")

# old date-only rule gives scene 1 the first 20201022 UDM2 listed (scene 2's)
old <- sapply(sr[1:2], function(f) grep(substr(basename(f), 1, 8), udm)[1])
check(old[[1]] != m[1], "date-only rule picks the wrong mask (regression guard)")

# same scene ID in two orders: prefer UDM2 in same directory
sr2  <- c("a/PSScene/20210101_100000_00_1111_3B_AnalyticMS_SR_clip.tif")
udm2 <- c("b/PSScene/20210101_100000_00_1111_3B_udm2_clip.tif",
          "a/PSScene/20210101_100000_00_1111_3B_udm2_clip.tif")
check(MatchUDM2(sr2, udm2) == 2L, "same-directory UDM2 preferred")
check(all(is.na(MatchUDM2(sr, character(0)))), "no UDM2 files -> all NA")

if (!ok) quit(status = 1)
cat("All UDM2 matching tests passed.\n")
