FROM mambaorg/micromamba:1.5.10

WORKDIR /app

# Copy env first for caching
COPY environment.yml /app/environment.yml

# Create env and install runtime deps
RUN micromamba create -y -n chembind -f /app/environment.yml && \
    micromamba run -n chembind pip install --no-cache-dir fastapi uvicorn && \
    micromamba clean -a -y

# Copy app code
COPY . /app

EXPOSE 8000

# Run uvicorn inside the conda env
CMD ["micromamba", "run", "-n", "chembind", "uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
