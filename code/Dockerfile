# Dockerfile for reproducible environment
FROM rocker/r-ver:4.5.3

# system dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev libssl-dev libcurl4-openssl-dev gdal-bin libgdal-dev \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /home/rstudio/project

# copy files
COPY . /home/rstudio/project

# install renv and restore
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org')"
RUN R -e "install.packages('renv', repos='https://cloud.r-project.org')"
RUN R -e "if (file.exists('renv.lock')) renv::restore()"

ENTRYPOINT ["Rscript","run_pipeline.R"]

