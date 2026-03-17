FROM python:3.11-slim
 
LABEL description="Telos-S MCP Server — Genomic intelligence tools for LLMs"
 
RUN apt-get update && apt-get install -y --no-install-recommends curl \
    && rm -rf /var/lib/apt/lists/*
 
WORKDIR /app
 
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
 
COPY mcp_server.py .
 
ENV PYTHONUNBUFFERED=1
 
EXPOSE 8001
 
HEALTHCHECK --interval=15s --timeout=5s --retries=3 \
    CMD curl -f http://localhost:8001/health || exit 1
 
CMD ["python", "mcp_server.py"]