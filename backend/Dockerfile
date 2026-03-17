FROM mambaorg/micromamba:2.5.0

ENV MAMBA_DOCKERFILE_ACTIVATE=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /app

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl ca-certificates tini \
  && rm -rf /var/lib/apt/lists/*
USER mambauser

COPY --chown=mambauser:mambauser backend/environment.yml /app/environment.yml
COPY --chown=mambauser:mambauser backend/requirements.txt /app/requirements.txt

RUN if [ -f /app/environment.yml ]; then \
      micromamba create -y -n app -f /app/environment.yml; \
    else \
      micromamba create -y -n app python=3.11; \
    fi && \
    micromamba run -n app python -V && \
    if [ -f /app/requirements.txt ]; then \
      micromamba run -n app pip install --no-cache-dir -r /app/requirements.txt; \
    fi && \
    micromamba clean -a -y

COPY --chown=mambauser:mambauser backend /app
COPY --chown=mambauser:mambauser backend/worker.sh /app/worker.sh
RUN chmod +x /app/worker.sh

EXPOSE 8000

ENTRYPOINT ["/usr/bin/tini","--"]

ENV PORT=8000
CMD ["sh","-c","exec micromamba run -n app uvicorn app.main:app --host 0.0.0.0 --port ${PORT}"]
