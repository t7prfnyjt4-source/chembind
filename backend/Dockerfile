FROM mambaorg/micromamba:1.5.10

WORKDIR /app

COPY environment.yml /app/environment.yml
COPY requirements-dev.txt /app/requirements-dev.txt

RUN micromamba create -y -n chembind -f /app/environment.yml && \
    micromamba run -n chembind pip install --no-cache-dir fastapi uvicorn && \
    micromamba run -n chembind pip install --no-cache-dir -r /app/requirements-dev.txt && \
    micromamba clean -a -y

COPY . /app

EXPOSE 8000

CMD ["micromamba", "run", "-n", "chembind", "uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
