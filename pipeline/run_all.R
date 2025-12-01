#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript run_all.R <siteIndex1> [<siteIndex2> ...]")
}

site_indices <- as.numeric(args)
if (any(is.na(site_indices))) stop("❌ Site indices must be numeric")

scripts <- c(
  "01-wetlsp-process.R",
  "02-wetlsp-chunks.R",
  "03-wetlsp-phenometrics.R",
  "04-wetlsp-generate.R",
  "05-wetlsp-format.R"
)

for (numSite in site_indices) {
  cat("\n🌎============================================\n")
  cat("▶ Starting full pipeline for site", numSite, "\n")
  cat("============================================\n\n")
  
  # Optional logging directory
  log_dir <- file.path("logs", paste0("site_", numSite))
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (s in scripts) {
    step_name <- gsub("^\\d{2}-|\\.R$", "", s)
    log_file <- file.path(log_dir, paste0(step_name, ".log"))
    cmd <- sprintf("Rscript pipeline/code/%s %d > %s 2>&1", s, numSite, log_file)
    
    cat("🚀 Running", s, "for site", numSite, "\n")
    status <- system(cmd)
    
    if (status != 0) {
      cat("❌ Error in", s, "for site", numSite, "\n")
      cat("📄 See log:", log_file, "\n")
      quit(status = status)
    } else {
      cat("✅ Completed", s, "for site", numSite, "\n\n")
    }
  }
  
  cat("🎉 All stages complete for site", numSite, "\n")
}

cat("\n🏁 All requested sites processed successfully.\n")