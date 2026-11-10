FROM python:3.11-slim

WORKDIR /app

# Prevent Python from writing .pyc files & buffering output
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Install basic system build essentials
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    && rm -rf /lib/apt/lists/*

# Append to Dockerfile
RUN apt-get update && apt-get install -y --no-install-recommends git \
    && rm -rf /lib/apt/lists/* \
    && git clone https://github.com/DanYev/Martini-3-DNA-RNA.git /opt/martini_dna_rna \
    && chmod +x /opt/martini_dna_rna/*.py

# Copy and install lightweight frontend requirements
COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Copy application source code
COPY . /app

# Expose Streamlit port
EXPOSE 8501

# Command to launch Streamlit frontend
CMD ["streamlit", "run", "frontend/app.py", "--server.port=8501", "--server.address=0.0.0.0"]
