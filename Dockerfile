# Use a Rocker image with base R and common packages
FROM harbor.cyverse.org/vice/rstudio/geospatial:latest

# Install system dependencies for spatial packages
RUN apt-get update && apt-get install -y \
    libproj-dev \
    libgdal-dev \
    libgeos-dev \
    libudunits2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install additional R packages your code needs
RUN Rscript -e "install.packages(c('raster', 'terra', 'ncdf4', 'foreach', 'doParallel', 'geojsonR', 'gdalUtilities', 'rjson'))"

# Create and set a working directory
# WORKDIR /app
RUN cd /app

# Copy all local files to the container
COPY . .

# Set default command (can be changed at run time)
CMD ["Rscript", "Code_LSP/01_img_process_gm_v1.R"]