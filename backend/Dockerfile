FROM mambaorg/micromamba:1.5.10

WORKDIR /app

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /app/environment.yml

RUN micromamba create -y -n chembind -f /app/environment.yml && \
    micromamba clean -a -y

COPY --chown=$MAMBA_USER:$MAMBA_USER . /app

EXPOSE 8000

CMD ["micromamba", "run", "-n", "chembind", "uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
