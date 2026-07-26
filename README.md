# Telos-S · Deploy

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21381037.svg)](https://doi.org/10.5281/zenodo.21381037)

### Powered by Telos Genomics

Docker Compose orchestration for the full **Telos-S** stack: the genomic analysis backend, the frontend, an MCP (Model Context Protocol) server, and, optionally, a local LLM stack (Ollama + OpenWebUI) for exploring results conversationally.

This repository **contains no application code**. It orchestrates containers built from two other repositories (`telos-s-backend`, `telos-s-frontend`) plus the MCP server included here.

---

## Architecture

```
                        ┌─────────────────┐
   User ──HTTP─────────▶│  telos-frontend  │  (Nginx + Svelte 5 build)
                        └────────┬─────────┘
                                 │ VITE_API_URL
                                 ▼
                        ┌─────────────────┐
                        │  telos-backend   │  (FastAPI + ESM-2 650M)
                        └────────┬─────────┘
                                 │ output/ (shared volume, read-only)
                                 ▼
                        ┌─────────────────┐
                        │   telos-mcp      │  (MCP server over results)
                        └────────┬─────────┘
                                 │
                     ┌───────────┴────────────┐
                     ▼                        ▼
              ┌─────────────┐         ┌───────────────┐
              │   Ollama     │◀───────▶│   OpenWebUI    │  (optional stack)
              └─────────────┘         └───────────────┘
```

- **telos-backend**: exposes the REST API (`/api/v1/analysis/...`), runs the analysis pipeline (Spike extraction → alignment → imputation → Telos Prophet → ESM-2 → report), and automatically downloads the Wuhan reference (`NC_045512.2`) on first boot.
- **telos-frontend**: SPA served by Nginx, consumes the backend API.
- **telos-mcp**: exposes analysis results (JSON files under `output/`) as MCP _tools_ so a local or remote LLM can query, summarize, and compare them in natural language.
- **ollama / openwebui**: optional local conversational stack for interacting with `telos-mcp` without depending on an external provider.

---

## Prerequisites

- Docker and Docker Compose v2
- Clone the application repos as sibling folders of this repo:

```
telos-s/                        ← working root folder (any name)
├── telos-s-deploy/              ← this repository
│   ├── docker-compose.yml
│   ├── Dockerfile               (MCP server image)
│   ├── mcp_server.py
│   ├── requirements.txt
│   ├── mcpo-config.json
│   └── config-claude.example.json
├── telos-s-backend/              ← git clone of the backend repo
│   ├── backend.py
│   ├── requirements.txt
│   ├── modules/
│   └── data/
└── telos-s-frontend/              ← git clone of the frontend repo
    ├── src/
    └── package.json
```

> `docker-compose.yml` references `../backend` and `../frontend` as build `context` — adjust these paths if you use different folder names.

---

## Configuration

Create a `.env` file next to `docker-compose.yml`:

```bash
# Public backend URL, reachable from the end user's BROWSER
# (not Docker's internal URL)
PUBLIC_URL=http://localhost:6002

# URL the frontend uses to call the backend (build-time, Vite)
VITE_API_URL=http://localhost:6002

# Mapbox token for the frontend's simulation tab
VITE_MAPBOX_TOKEN=pk.eyJ1Ijoi...
```

| Variable            | Service              | Description                                                                               |
| ------------------- | -------------------- | ----------------------------------------------------------------------------------------- |
| `PUBLIC_URL`        | backend              | Base URL used to build the download links (CSV, heatmap, report) returned to the frontend |
| `VITE_API_URL`      | frontend (build arg) | Backend URL that the user's browser should use                                            |
| `VITE_MAPBOX_TOKEN` | frontend (build arg) | Public Mapbox token                                                                       |
| `ESM_2_SIZE`        | backend              | ESM-2 model to use (default `facebook/esm2_t33_650M_UR50D`)                               |

---

## Bringing up the stack

```bash
docker compose up -d
```

This builds and starts:

| Service          | Host port       | Description                  |
| ---------------- | --------------- | ---------------------------- |
| `telos-frontend` | `6001` → `80`   | Web interface                |
| `telos-backend`  | `6002` → `8000` | REST API + pipeline          |
| `telos-mcp`      | `8001`          | MCP server (Streamable HTTP) |

On first boot, `telos-backend` automatically downloads the Wuhan reference from NCBI (`entrypoint.sh`) and extracts the reference Spike protein — no manual steps required beyond having internet access from the container.

### Checking service health

```bash
curl http://localhost:6002/health   # backend
curl http://localhost:8001/health   # MCP server
```

---

## Extended stack (Ollama + OpenWebUI + pentesting MCP)

`docker-compose-example.yml` additionally includes:

- **`ollama`**: local LLM (`llama3.1:8b` or higher recommended for interpreting results). On Apple Silicon with OrbStack it uses native Metal; on Linux with an NVIDIA GPU, uncomment the `deploy.resources` block to pass the GPU through to the container.
- **`openwebui`**: conversational interface over Ollama, preconfigured to consume `telos-mcp` as a tool source.
  To use this extended stack:

```bash
docker compose -f docker-compose-example.yml up -d ollama
docker compose -f docker-compose-example.yml exec ollama ollama pull llama3.1:8b
docker compose -f docker-compose-example.yml up -d
```

---

## The MCP server (`telos-mcp`)

Exposes pipeline results to any MCP client (OpenWebUI, Claude Desktop, etc.) via **Streamable HTTP** (spec `2025-03-26`), with additional legacy SSE support for compatibility.

### Available tools

| Tool                      | Description                                                                                                         |
| ------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| `get_analysis_results`    | Full results of an analysis by `job_id`: score, lineage, mutations, epidemiological parameters, Prophet predictions |
| `list_recent_analyses`    | Lists the N most recent analyses, with status filtering                                                             |
| `get_variant_summary`     | Structured narrative summary, optimized for an LLM to explain to non-experts                                        |
| `compare_variants`        | Side-by-side comparison of two analyses (score, shared/unique mutations, R0)                                        |
| `get_prophet_predictions` | Only Telos Prophet's evolutionary predictions (positions 452, 484, 501, 681)                                        |

### Connecting an MCP client

**OpenWebUI** (via `mcpo`, included in the extended stack): uses `mcpo-config.json`, which points to `http://telos-mcp:8001/mcp`.

**Claude Desktop** (or another external MCP client): copy `config-claude.example.json` into your MCP server configuration, replacing `TU_USUARIO` with your system user:

```json
{
  "mcpServers": {
    "telos-s-mcp": {
      "command": "/etc/profiles/per-user/TU_USUARIO/bin/npx",
      "args": ["-y", "mcp-remote", "http://localhost:8001/mcp"]
    }
  }
}
```

---

## Volumes

| Volume              | Contents                                                                                                           | Persistence |
| ------------------- | ------------------------------------------------------------------------------------------------------------------ | ----------- |
| `telos_output`      | Pipeline results (CSVs, heatmaps, reports, JSON) — shared between backend and MCP server (MCP mounts it read-only) | Persistent  |
| `telos_reference`   | Automatically downloaded Wuhan reference                                                                           | Persistent  |
| `huggingface_cache` | ESM-2 model cache (~2.5 GB) — avoids re-downloading on every rebuild                                               | Persistent  |
| `ollama_models`     | Local LLM models (extended stack)                                                                                  | Persistent  |
| `openwebui_data`    | OpenWebUI configuration and conversations (extended stack)                                                         | Persistent  |

---

## Troubleshooting

**Backend fails to start / Wuhan reference download fails**
Check that the `telos-backend` container has outbound internet access (`entrypoint.sh` runs `curl` against `eutils.ncbi.nlm.nih.gov`). Alternatively, manually mount `spike_wuhan.txt` into `/app/data/` via a volume.

**Frontend can't reach the backend**
`VITE_API_URL` is injected at **build time**, not runtime — if you change the backend URL, you must rebuild the frontend image (`docker compose build telos-frontend`).

**`telos-mcp` not responding / inconsistent session in OpenWebUI**
The server re-registers unknown `Mcp-Session-Id` sessions instead of rejecting them with a 404 (to tolerate container restarts). If the issue persists, check `docker logs telos-mcp`.

---

## Related repositories

- `telos-s-backend` — FastAPI API + analysis pipeline (ESM-2, Telos Prophet)
- `telos-s-frontend` — Svelte 5 + Vite SPA

---
