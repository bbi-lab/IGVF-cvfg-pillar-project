FROM python:3.12-slim

RUN pip install --no-cache-dir poetry==2.1.2

WORKDIR /usr/src/app

COPY pyproject.toml poetry.lock ./
RUN poetry config virtualenvs.create false \
    && poetry install --no-root --no-interaction

COPY . .

ENV PYTHONUNBUFFERED=1

# Entrypoint is set per-service in compose.yaml (mirrors variant-annotation's
# Dockerfile, which stays generic and lets each compose service pick its own
# `python -m src.<module>` entrypoint).
CMD ["bash", "-c", "tail -f /dev/null"]
