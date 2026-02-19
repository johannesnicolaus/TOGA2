FROM continuumio/anaconda3

LABEL Author="JN Wibisana"
LABEL Version="20260219_1.0.0"

# Set environment variables
ENV PATH=/opt/TOGA2/toga2/bin:$HOME/.cargo/bin:/opt/nextflow:/opt/conda/bin:/opt/TOGA2:/opt/TOGA2/bin:/opt/TOGA2/bin/prank:/opt/TOGA2/bin/iqtree2:/opt/TOGA2/bin/intronIC:$PATH
ENV SHELL=/bin/bash
ENV VIRTUAL_ENV=/opt/TOGA2/toga2

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    openjdk-17-jdk \
    git \
    bsd-mailx \
    mailutils \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install Rust/Cargo
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Install Nextflow
RUN mkdir -p /opt/nextflow && \
    curl -sSL https://get.nextflow.io | bash && \
    mv nextflow /opt/nextflow/nextflow && \
    chmod +x /opt/nextflow/nextflow

# Facilitate installation by removing conflicting packages
RUN python3 -m pip uninstall -y gensim numba datashader && \
    python3 -m pip install streamlit --upgrade --root-user-action=ignore && \
    python3 -m pip install pyarrow --upgrade --root-user-action=ignore && \
    python3 -m pip install bottleneck --upgrade --root-user-action=ignore

# Make bash the default shell
RUN sed -i 's/required/sufficient/g' /etc/pam.d/chsh && \
    rm /bin/sh && ln -f /bin/bash /bin/sh

# Clone TOGA2 repository (main branch for v2.0.7a)
RUN git clone --recurse-submodules -b main https://github.com/hillerlab/TOGA2 /opt/TOGA2

# Build TOGA2 (create virtual environment for maturin)
WORKDIR /opt/TOGA2
RUN python3 -m venv toga2 && \
    . toga2/bin/activate && \
    make

# Create writable directory for Nextflow temporary files (fallback)
RUN mkdir -p /.nextflow && chmod 777 /.nextflow

# Create dummy user for Slurm
RUN adduser --disabled-password --gecos "" slurm && \
    echo "slurmadmin:x:300:300::/opt/slurm/slurm:/bin/false" >> /etc/passwd && \
    echo "slurmadmin:x:300:" >> /etc/group

# Set working directory
WORKDIR /opt/TOGA2

# Create startup script that sets up .nextflow directory
RUN echo '#!/bin/bash\n\
source /opt/TOGA2/toga2/bin/activate\n\
\n\
# Set NXF_HOME to user home .nextflow if it exists, otherwise use container fallback\n\
if [ -n "$HOME" ] && [ -d "$HOME" ]; then\n\
    export NXF_HOME="$HOME/.nextflow"\n\
    mkdir -p "$NXF_HOME"\n\
else\n\
    export NXF_HOME="/.nextflow"\n\
fi\n\
\n\
exec "$@"\n' > /opt/entrypoint.sh && \
    chmod +x /opt/entrypoint.sh

ENTRYPOINT ["/opt/entrypoint.sh"]

CMD ["/bin/bash"]
