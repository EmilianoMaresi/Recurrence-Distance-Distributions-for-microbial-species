FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV CONDA_DIR=/opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH

# ---- System dependencies ----
RUN apt-get update && apt-get install -y \
    wget \
    git \
    ca-certificates \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# ---- Install Miniconda ----
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p $CONDA_DIR && \
    rm /tmp/miniconda.sh

# ---- Accept Conda Terms of Service ----
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# ---- Create conda environment ----
COPY environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml && conda clean -a -y


# ---- Activate conda environment by default ----
ENV CONDA_DEFAULT_ENV=RDD_env
ENV PATH=$CONDA_DIR/envs/RDD_env/bin:$PATH

# ---- Set working directory ----
WORKDIR /app

# ---- Copy project code ----
COPY . /app

# ---- Default command ----
ENTRYPOINT ["python"]
CMD ["main.py"]

