Understood — if the phenometrics live in the per-year WetLSP_<year>.nc files, the clean approach is:
  1.	extract the phenometrics at each pixel location in sf_all (for each year you want)
2.	assemble a feature matrix (metrics × years)
3.	cluster
4.	join class back to sf_all
5.	optionally write a class raster and split sf_all by class

Below is a working template you can paste in.

⸻

1) Extract phenometrics from NetCDF at sf_all pixel points

library(terra)
library(sf)
library(dplyr)
library(purrr)

# ---- inputs ----
years <- 2021:2024
nc_dir <- file.path(site_dir, "netcdf")  # adjust if needed
nc_paths <- setNames(
  file.path(nc_dir, paste0("WetLSP_", years, ".nc")),
  years
)

# Choose a sane subset first (you can add more later)
vars <- c("NumCycles","EVImax","EVIamp","EVIarea","Peak","OGI","OGD","OGMn","numObs","QA")

# Ensure you have a stable id
if (!("pixel_id" %in% names(sf_all))) sf_all$pixel_id <- seq_len(nrow(sf_all))

# Convert sf points to SpatVector for terra::extract
v <- terra::vect(sf_all)

extract_year <- function(nc_path, vars, year) {
  r <- terra::rast(nc_path)
  
  missing <- setdiff(vars, names(r))
  if (length(missing) > 0) stop("Missing vars in ", nc_path, ": ", paste(missing, collapse=", "))
  
  ex <- terra::extract(r[[vars]], v)  # returns data.frame with ID + vars
  ex <- ex |>
    as_tibble() |>
    rename(extract_id = ID) |>
    mutate(pixel_id = sf_all$pixel_id[extract_id]) |>
    select(-extract_id)
  
  # suffix with year so features are unique
  names(ex)[names(ex) %in% vars] <- paste0(vars, "_", year)
  ex
}

feat_list <- imap(nc_paths, ~extract_year(.x, vars, .y))
feat <- reduce(feat_list, left_join, by = "pixel_id")

Quick sanity checks

dplyr::glimpse(feat)
colSums(is.na(feat))

If you see tons of NAs, that’s usually because some pixels are outside the valid footprint for that year (or masked). That’s not fatal; you just need an NA strategy.

⸻

2) Build feature matrix + cluster

Minimal NA handling (recommended first pass)
•	Drop pixels with “too many” missing values.
•	Impute remaining NAs with column medians.

X <- feat |> select(-pixel_id)

# Drop rows with >25% missing (tune this)
keep <- rowMeans(is.na(X)) <= 0.25
feat2 <- feat[keep, ]
X2 <- X[keep, ]

# Median impute
med <- apply(X2, 2, \(z) median(z, na.rm = TRUE))
for (j in seq_len(ncol(X2))) {
  idx <- is.na(X2[[j]])
  if (any(idx)) X2[[j]][idx] <- med[j]
}

# Standardize
X2s <- scale(as.matrix(X2))

set.seed(1)
k <- 8  # tune later
km <- kmeans(X2s, centers = k, nstart = 20)

feat2$pheno_class <- factor(km$cluster)

Join back to sf_all (pixels that were dropped remain NA class):
  
  sf_all$pheno_class <- NA
sf_all$pheno_class[match(feat2$pixel_id, sf_all$pixel_id)] <- as.character(feat2$pheno_class)
sf_all$pheno_class <- factor(sf_all$pheno_class)


⸻

3) Make a class map raster (very useful)

# Use any year's nc as the template grid
tmpl <- terra::rast(nc_paths[[1]])[[1]]  # first variable as template

# Rasterize the point classes (points must align to pixel centers; your sf_all should)
pts <- terra::vect(sf_all[!is.na(sf_all$pheno_class), c("pheno_class")])
pts$pheno_class_num <- as.integer(factor(pts$pheno_class))

class_r <- terra::rasterize(pts, tmpl, field = "pheno_class_num", touches = FALSE)

plot(class_r, main = paste0(site_id, " phenology classes (k=", k, ")"))

Save it:
  
  terra::writeRaster(class_r, file.path(site_dir, paste0(site_id, "_pheno_class_k", k, ".tif")),
                     overwrite = TRUE)


⸻

4) Re-chunk sf_all into class files

out_dir <- file.path(site_dir, "tables_sf_pheno_classes")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

by_class <- split(sf_all[!is.na(sf_all$pheno_class), ], sf_all$pheno_class[!is.na(sf_all$pheno_class)])

purrr::iwalk(by_class, \(obj, cls) {
  saveRDS(obj, file.path(out_dir, paste0(site_id, "_class_", cls, ".rds")))
})


⸻

Two important performance notes
1.	terra::extract() is doing the heavy lifting. If sf_all is huge and this is slow, do it in batches:
  
  batch_ids <- split(seq_len(nrow(sf_all)), ceiling(seq_len(nrow(sf_all))/50000))
# loop batches: subset sf_all, subset v, extract, rbind

2.	Don’t start with all metrics. Start with a compact set (the vars above), validate the class map, then expand.

⸻

If you paste the names(rast(nc_path)) for one of your NetCDFs (just the variable names), I can recommend the best minimal feature set for clustering (often fewer, less collinear features gives cleaner classes).