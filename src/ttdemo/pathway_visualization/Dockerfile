FROM openjdk:11-jre-slim

# Set environment
ENV PATHVISIO_VERSION=3.3.0

# Install wget and unzip
RUN apt-get update && apt-get install -y wget unzip xvfb && rm -rf /var/lib/apt/lists/*

# Download PathVisio
RUN wget https://github.com/PathVisio/pathvisio/releases/download/v${PATHVISIO_VERSION}/PathVisio-${PATHVISIO_VERSION}-linux.tar.gz \
    && tar -xzf PathVisio-${PATHVISIO_VERSION}-linux.tar.gz -C /opt/ \
    && rm PathVisio-${PATHVISIO_VERSION}-linux.tar.gz

# Copy any plugins, if needed
# COPY plugins/*.jar /opt/PathVisio/plugins/

# Set working directory
WORKDIR /opt/PathVisio

# Start PathVisio with RPC server (headless with Xvfb)
CMD ["xvfb-run", "-a", "./pathvisio.sh", "--rpc"]
