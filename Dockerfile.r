# R environment for the manuscript's figure-generation scripts
# (notebooks/figures/**/*.R, notebooks/figures/**/*.Rmd).
#
# Runs on Linux, where grDevices::cairo_pdf() links against the system
# libcairo2 package directly -- unlike the CRAN macOS build of R, which links
# cairo_pdf() against XQuartz's X11 libs. That's what makes this image able to
# produce cairo_pdf() output without installing XQuartz on the host.
#
# Base is rocker/r-ver (not rocker/tidyverse) because rocker/tidyverse's
# versioned tags (e.g. 4.4.2) are amd64-only -- no arm64 manifest -- which
# breaks `docker build` on Apple Silicon. rocker/r-ver publishes multi-arch
# manifests for the same version tags, so we install tidyverse ourselves.
FROM rocker/r-ver:4.4.2

# System libraries needed to build tidyverse/ggplot2's compiled deps
# (ragg, textshaping, systemfonts, xml2, curl, ssl) plus cairo (for
# grDevices::cairo_pdf) and pandoc (for rmarkdown::render() on the .Rmd
# figure scripts). cabextract + ttf-mscorefonts-installer provide real Arial,
# matching the scripts' `family = 'Arial'` theme settings (ggplot2 theme
# text= / element_geom()) -- the figures were designed/reviewed on macOS,
# where Arial ships as a system font, rather than a metric-compatible
# substitute like Liberation Sans.
RUN echo "ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true" | debconf-set-selections \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
       cabextract \
       ttf-mscorefonts-installer \
       fontconfig \
       pandoc \
       libcairo2-dev \
       libfontconfig1-dev \
       libfreetype6-dev \
       libharfbuzz-dev \
       libfribidi-dev \
       libpng-dev \
       libtiff5-dev \
       libjpeg-dev \
       libxml2-dev \
       libcurl4-openssl-dev \
       libssl-dev \
    && fc-cache -f \
    && rm -rf /var/lib/apt/lists/*

# ttf-mscorefonts-installer's postinst only queues a request for
# update-notifier's package-data-downloader hook to fetch and cabextract the
# actual font archives; that hook is normally fired later by a running
# system, which a `docker build` never has, so the package "installs"
# successfully but /usr/share/fonts/truetype/msttcorefonts is left empty.
# extrafont::font_import() below then crashes (0 ttf files found) instead of
# reporting a missing-font error. Run the downloader synchronously so the
# fonts are actually in place before that step.
RUN /usr/lib/update-notifier/package-data-downloader \
    && fc-cache -f

# CRAN packages used across the figure scripts (see library()/require() calls
# under notebooks/figures/). rocker/r-ver:4.4.2 pins its
# default repo to a Posit Package Manager snapshot dated 2025-02-27 (predates
# ggplot2 4.0.0, released September 2025) -- override to a later snapshot so
# we get ggplot2 >= 4.0, required by notebooks/figures/figure_5_6/Figure_6b.R's
# theme(geom = element_geom(...)), which doesn't exist before 4.0. Deliberately
# still a *dated* p3m.dev snapshot (not the rolling cloud.r-project.org
# mirror): that live mirror only serves source packages and choked on
# transitive deps (googledrive/googlesheets4/reprex for tidyverse,
# htmlwidgets for plotly) when tested; p3m.dev's per-date snapshots serve
# prebuilt binaries and installed cleanly. Bump this date if a script starts
# needing a CRAN package/version newer than what's on 2025-10-15.
RUN install2.r --error --skipinstalled --repos https://p3m.dev/cran/2025-10-15 \
    tidyverse \
    patchwork \
    ggh4x \
    extrafont \
    ggforce \
    ggrepel \
    ggsci \
    plotly \
    readxl \
    eulerr \
    reticulate \
    svglite \
    remotes

# ggsankey has no CRAN release; install from GitHub (per
# notebooks/figures/extended_data_figure_4_6_7_8_9/README_Extended_data_figures.md).
RUN Rscript -e 'remotes::install_github("davidsjoberg/ggsankey")'

# Populate extrafont's font database at build time (font_import() scans
# system TTFs and builds metrics -- a slow, one-time step) so containers start
# with Arial etc. already registered instead of paying this cost, or hitting
# "font family 'Arial' not found" render failures, on every run.
RUN Rscript -e 'extrafont::font_import(prompt = FALSE)'

# The figure scripts call print(<plot>) to preview the plot before ggsave()'ing
# it. With no display attached, R's default device for that print() falls
# back to base grDevices::pdf() -- which, unlike cairo_pdf()/ggsave(), resolves
# fonts against R's built-in Type1/AFM table rather than via fontconfig, so it
# doesn't see whatever extrafont registered and dies with "invalid font type"
# on the first non-Type1 glyph (e.g. Arial, or the ≤/≥ characters in the
# figure's classification labels). Point the default device at cairo_pdf
# instead, which -- like ggsave() -- resolves fonts via fontconfig/systemfonts,
# so print() renders the same way the saved files will.
RUN echo 'options(device = function(...) grDevices::cairo_pdf(filename = tempfile(fileext = ".pdf"), ...))' \
    >> /usr/local/lib/R/etc/Rprofile.site

WORKDIR /usr/src/app

CMD ["bash", "-c", "tail -f /dev/null"]
