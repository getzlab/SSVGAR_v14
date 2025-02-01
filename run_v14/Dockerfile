# Use cicirello/pyaction:4.14.0 as base image for docker 
FROM cicirello/pyaction:4.14.0

# Install/set up from root directory
WORKDIR /

### INSTALLATIONS:
# Clone SSVGAR Github repo
RUN git clone -b v14 https://ssvgar-dummy-user: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@github.com/getzlab/SSVGAR.git
# Install requirements.txt for this version of SSVGAR
RUN python -m pip install --upgrade pip && \
    pip install --no-cache-dir -r SSVGAR/run_v14/requirements.txt
# Install bedtools 2.27.1 for SSVGAR
RUN curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg | dd of=/usr/share/keyrings/githubcli-archive-keyring.gpg && \
    chmod go+r /usr/share/keyrings/githubcli-archive-keyring.gpg && \
    echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" | tee /etc/apt/sources.list.d/github-cli.list > /dev/null && \
    apt update && \ 
    apt-get install -y wget && apt-get install -y build-essential && apt-get install -y libz-dev
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz && \
    tar -zxvf bedtools-2.27.1.tar.gz
WORKDIR /bedtools2
RUN make
ENV PATH="/bedtools2/bin:$PATH"
# Install text editors for convenience (vim and emacs)
RUN apt install -y vim && apt install -y emacs

# Set the working directory to be the SSVGAR Github repo
WORKDIR /SSVGAR
