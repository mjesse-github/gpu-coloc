FROM pytorch/pytorch:2.6.0-cuda12.4-cudnn9-runtime

RUN apt-get update && apt-get install -y --no-install-recommends git && \
    rm -rf /var/lib/apt/lists/*

RUN pip install -U pip && pip install --no-cache-dir gpu-coloc

WORKDIR /workspace

ENTRYPOINT ["gpu-coloc"]