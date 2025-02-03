FROM mambaorg/micromamba:1.4.2-bullseye

# We want all RUN commands to use Bash.
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Copy over your environment file
COPY environment.yml /tmp/environment.yml

# Create the environment
RUN micromamba create -n hf39 -f /tmp/environment.yml --yes && \
    micromamba clean --all --yes

ARG WORKDIR=/usr/local/roadmap2030
ENV WORKDIR=${WORKDIR}

RUN micromamba shell init -s bash -p /opt/conda

# If needed, ensure the file exists and append your activation line
RUN touch /home/mambauser/.bashrc
RUN echo 'micromamba activate hf39' >> /home/mambauser/.bashrc

USER root
RUN apt-get update -y
RUN apt install git -y
RUN apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    make \
    build-essential \
    zlib1g-dev \
    libpython3-dev \
    libssl-dev \
    libffi-dev \
    libsqlite3-dev \
    libgdal-dev

RUN usermod -aG sudo mambauser && \
    echo "mambauser ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

RUN apt-get update -y && apt-get install htop gpg curl -y
RUN curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list
RUN apt-get update -y && apt-get install google-cloud-cli -y

ENV GEE_KEY_PATH=/usr/local/secrets/service-account-key.json
ENV GOOGLE_APPLICATION_CREDENTIALS=/usr/local/secrets/service-account-key.json

ARG CACHEBUST
RUN echo $(date +%s) > /tmp/cachebuster

RUN git clone https://github.com/springinnovate/ecoshard.git /usr/local/ecoshard && \
    cd /usr/local/ecoshard && \
    micromamba run -n hf39 pip install . && \
    git log -1 --format="%h on %ci" > /usr/local/ecoshard.gitversion

RUN git clone https://github.com/springinnovate/inspring.git /usr/local/inspring && \
    cd /usr/local/inspring && \
    micromamba run -n hf39 pip install . && \
    git log -1 --format="%h on %ci" > /usr/local/inspring.gitversion

# This shows the timedate of the ecoshard and inspring repos
RUN echo 'if [ -f "/usr/local/ecoshard.gitversion" ]; then' >> /home/mambauser/.bashrc && \
    echo '  echo "ecoshard: commit on $(cat /usr/local/ecoshard.gitversion)"' >> /home/mambauser/.bashrc && \
    echo 'fi' >> /home/mambauser/.bashrc && \
    echo 'if [ -f "/usr/local/inspring.gitversion" ]; then' >> /home/mambauser/.bashrc && \
    echo '  echo "inspring: commit on $(cat /usr/local/inspring.gitversion)"' >> /home/mambauser/.bashrc && \
    echo 'fi' >> /home/mambauser/.bashrc

USER mambauser
WORKDIR ${WORKDIR}
CMD ["/bin/bash"]
#CMD ["micromamba", "run", "-n", "hf39", "python", "-m", "ecoshard.geosharding.geosharding", "swy_amazon_1992.ini", "--debug_n_aoi", "1"]
