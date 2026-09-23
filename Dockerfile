FROM python:3.12-slim

# cabextract + ttf-mscorefonts-installer provide real Arial, matching the
# Python figure scripts' matplotlib styling (e.g.
# src/make_extended_data_figure_7alt.py's Arial-based rcParams) -- the
# figures were designed/reviewed on macOS, where Arial ships as a system
# font, rather than a metric-compatible substitute. Same reasoning as
# Dockerfile.r's ttf-mscorefonts-installer step. ttf-mscorefonts-installer
# lives in Debian's contrib component (non-free EULA'd fonts), not enabled
# by default in this slim base's sources, so add it first. Unlike
# Dockerfile.r's Ubuntu-based rocker/r-ver image, this Debian base's
# ttf-mscorefonts-installer extracts the fonts synchronously during `apt-get
# install` itself -- no separate update-notifier hook to trigger.
RUN sed -i 's/^Components: main/Components: main contrib/' /etc/apt/sources.list.d/debian.sources \
    && echo "ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true" | debconf-set-selections \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
       cabextract \
       ttf-mscorefonts-installer \
       fontconfig \
    && fc-cache -f \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir poetry==2.1.2

WORKDIR /usr/src/app

COPY pyproject.toml poetry.lock ./
# --extras notebooks: pulls in ipykernel/jupyterlab so `jupyter nbconvert
# --execute` (Stage 2 analysis notebooks, Stage 3 figure notebooks) works.
# Not --all-extras: dev/test tooling (ruff, pytest, pre-commit) doesn't
# belong in the runtime image.
RUN poetry config virtualenvs.create false \
    && poetry install --no-root --no-interaction --extras notebooks

COPY . .

ENV PYTHONUNBUFFERED=1

# Entrypoint is set per-service in compose.yaml (mirrors variant-annotation's
# Dockerfile, which stays generic and lets each compose service pick its own
# `python -m src.<module>` entrypoint).
CMD ["bash", "-c", "tail -f /dev/null"]
