"""
Telos-S MCP Server
------------------
A Model Context Protocol server that exposes the results of the Telos-S pipeline to local language models (Ollama) via OpenWebUI.

It implements the MCP protocol over HTTP/SSE, compatible with OpenWebUI >= 0.5 and any standard MCP client.

Available tools:
    - get_analysis_results → Complete results for a specific analysis by job_id
    - list_recent_analyses → Lists the N most recent analyses
    - get_variant_summary → An executive summary in natural language, ready for LLMs
    - compare_variants → Comparison between two analyses
    - get_prophet_predictions → Only future evolution predictions

Standalone usage (without Docker):
    pip install fastapi uvicorn sse-starlette
    python mcp_server.py
 
Default port: 8001
"""
 
import json
import os
import glob
import uuid
from datetime import datetime
from pathlib import Path
from typing import Optional, Any
 
from fastapi import FastAPI, HTTPException, Request, Response
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from sse_starlette.sse import EventSourceResponse
import uvicorn
 
# =============================================================================
# CONFIGURATION
# =============================================================================
 
OUTPUT_DIR = Path(os.getenv("TELOS_OUTPUT_DIR", "/app/output"))
JOBS_DIR = OUTPUT_DIR / "jobs"
MCP_PORT = int(os.getenv("MCP_PORT", "8001"))
 
# Información del servidor MCP
SERVER_INFO = {
    "name": "telos-s-mcp",
    "version": "0.1.1",
    "description": "Telos-S genomic intelligence — SARS-CoV-2 variant analysis tools",
    "vendor": "Telos Genomics",
}
 
# Definición de tools MCP
TOOLS = [
    {
        "name": "get_analysis_results",
        "description": (
            "Retrieves the complete results of a Telos-S variant analysis by job ID. "
            "Returns aggression score, lineage classification, mutations list, "
            "epidemiological parameters (R0, incubation period), and Prophet predictions "
            "for future evolution at key positions (452, 484, 501, 681). "
            "Use this when the user asks about a specific analysis or variant."
        ),
        "inputSchema": {
            "type": "object",
            "properties": {
                "job_id": {
                    "type": "string",
                    "description": "The job ID returned when the analysis was submitted (e.g. 'job_abc123def456')"
                }
            },
            "required": ["job_id"]
        }
    },
    {
        "name": "list_recent_analyses",
        "description": (
            "Lists the most recent Telos-S analyses with their key results. "
            "Returns variant name, aggression score, lineage, and completion time. "
            "Use this when the user asks what analyses have been done, or wants to "
            "explore available results."
        ),
        "inputSchema": {
            "type": "object",
            "properties": {
                "limit": {
                    "type": "integer",
                    "description": "Maximum number of analyses to return (default: 10, max: 50)",
                    "default": 10
                },
                "status_filter": {
                    "type": "string",
                    "description": "Filter by status: 'completed', 'processing', 'failed', or 'all'",
                    "default": "completed"
                }
            }
        }
    },
    {
        "name": "get_variant_summary",
        "description": (
            "Returns a structured narrative summary of a variant analysis, "
            "optimized for language model interpretation and explanation to non-experts. "
            "Includes risk assessment, key mutations in plain language, evolutionary "
            "pressure predictions, and epidemiological implications. "
            "Use this when the user asks you to explain, interpret, or describe a variant."
        ),
        "inputSchema": {
            "type": "object",
            "properties": {
                "job_id": {
                    "type": "string",
                    "description": "The job ID of the completed analysis"
                }
            },
            "required": ["job_id"]
        }
    },
    {
        "name": "compare_variants",
        "description": (
            "Compares two Telos-S variant analyses side by side. "
            "Returns differences in aggression score, lineage, key mutations, "
            "R0 estimates, and Prophet predictions. "
            "Use this when the user asks to compare two variants or asks which is more dangerous."
        ),
        "inputSchema": {
            "type": "object",
            "properties": {
                "job_id_a": {
                    "type": "string",
                    "description": "Job ID of the first variant analysis"
                },
                "job_id_b": {
                    "type": "string",
                    "description": "Job ID of the second variant analysis"
                }
            },
            "required": ["job_id_a", "job_id_b"]
        }
    },
    {
        "name": "get_prophet_predictions",
        "description": (
            "Returns only the Telos Prophet evolutionary predictions for a variant. "
            "These are ESM-2-based structural stability predictions for the 4 critical "
            "positions: RBM-452, RBM-484, RBM-501, and Furin-681. "
            "Use this when discussing viral evolution, immune escape potential, "
            "or future mutation paths."
        ),
        "inputSchema": {
            "type": "object",
            "properties": {
                "job_id": {
                    "type": "string",
                    "description": "The job ID of the completed analysis"
                }
            },
            "required": ["job_id"]
        }
    }
]
 
# =============================================================================
# FASTAPI APP
# =============================================================================
 
app = FastAPI(
    title="Telos-S MCP Server",
    description="MCP server for Telos-S genomic intelligence",
    version="0.1.1"
)
 
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)
 
# =============================================================================
# Load the JSON status of a job.
# =============================================================================
 
def load_job(job_id: str) -> dict:
    """Load the JSON status of a job."""
    job_file = JOBS_DIR / f"{job_id}.json"
    if not job_file.exists():
        raise HTTPException(status_code=404, detail=f"Job '{job_id}' not found")
    with open(job_file) as f:
        return json.load(f)
 
 
def load_prophet_data(job_id: str, variant_name: str) -> Optional[list]:
    """Load the Prophet predictions for a specific variant."""
    prophet_path = OUTPUT_DIR / "prophet" / f"mutation_predictions_spike_{job_id}_{variant_name}.json"
    if prophet_path.exists():
        with open(prophet_path) as f:
            return json.load(f)
    return None
 
 
def score_to_risk_label(score: float) -> str:
    if score > 1200:
        return "MAXIMUM ALERT — High immune evasion capacity"
    elif score > 600:
        return "ACTIVE MONITORING — Variant of interest"
    else:
        return "OBSERVATION — Moderate mutations"
 
 
def get_aggression_score(results: dict) -> float:
    """
    Read the `aggression_score` from wherever it is available.
    New jobs: results["aggression_score"] (promoted by the backend)
    Old jobs: results["epi_params"]["aggression_score"]
    """
    score = results.get("aggression_score")
    if score is None or score == 0:
        score = results.get("epi_params", {}).get("aggression_score", 0.0)
    return float(score or 0.0)
 
 
def format_prophet_for_llm(prophet_data: list) -> str:
    """Convert Prophet predictions into structured text for use with large language models (LLMs)."""
    if not prophet_data:
        return "No Prophet predictions available."
 
    lines = []
    for target in prophet_data:
        pos = target.get("wuhan_position")
        original = target.get("original_aa", "?")
        name = target.get("target_site", f"Position {pos}")
        predictions = target.get("predictions", [])
 
        top_candidates = [p for p in predictions if p.get("amino") != original]
        if top_candidates:
            top = top_candidates[0]
            alert = "⚠️ HIGH" if top["confidence"] > 20 else "LOW"
            lines.append(
                f"  {name} (pos {pos}): currently {original} → "
                f"most likely evolution to {top['amino']} "
                f"({top['confidence']:.1f}% structural probability) — "
                f"evolutionary pressure: {alert}"
            )
        else:
            lines.append(f"  {name} (pos {pos}): structurally stable, no dominant mutation path")
 
    return "\n".join(lines)
 
 
# =============================================================================
# TOOLS IMPLEMENTATION
# =============================================================================
 
def tool_get_analysis_results(job_id: str) -> dict:
    job = load_job(job_id)
 
    if job["status"] != "completed":
        return {
            "status": job["status"],
            "message": f"Analysis is {job['status']}. Current step: {job.get('current_step', 'unknown')}",
            "progress": job.get("progress", 0)
        }
 
    results = job.get("results", {})
    variant_name = results.get("variant_name", "Unknown")
    prophet_data = load_prophet_data(job_id, variant_name)
 
    return {
        "job_id": job_id,
        "status": "completed",
        "variant_name": variant_name,
        "aggression_score": get_aggression_score(results),
        "risk_level": score_to_risk_label(get_aggression_score(results)),
        "lineage": results.get("lineage", "Unknown"),
        "lineage_confidence": results.get("lineage_confidence", 0),
        "sequence_quality": results.get("sequence_quality", 0),
        "mutations": results.get("mutations", []),
        "epi_params": results.get("epi_params", {}),
        "prophet_predictions": prophet_data,
        "completed_at": job.get("completed_at"),
        "files": results.get("files", {})
    }
 
 
def tool_list_recent_analyses(limit: int = 10, status_filter: str = "completed") -> list:
    job_files = sorted(JOBS_DIR.glob("*.json"), key=lambda f: f.stat().st_mtime, reverse=True)
    results = []
 
    for job_file in job_files[:100]:  # Read up to 100 characters to filter
        if len(results) >= limit:
            break
        try:
            with open(job_file) as f:
                job = json.load(f)
 
            if status_filter != "all" and job.get("status") != status_filter:
                continue
 
            job_results = job.get("results", {})
            score = get_aggression_score(job_results)
 
            results.append({
                "job_id": job_file.stem,
                "variant_name": job_results.get("variant_name", "Unknown"),
                "status": job.get("status"),
                "aggression_score": score,
                "risk_level": score_to_risk_label(score) if job.get("status") == "completed" else None,
                "lineage": job_results.get("lineage"),
                "completed_at": job.get("completed_at"),
                "started_at": job.get("started_at")
            })
        except Exception:
            continue
 
    return results
 
 
def tool_get_variant_summary(job_id: str) -> dict:
    job = load_job(job_id)
 
    if job["status"] != "completed":
        return {"error": f"Analysis not completed yet. Status: {job['status']}"}
 
    results = job.get("results", {})
    variant_name = results.get("variant_name", "Unknown")
    score = get_aggression_score(results)
    lineage = results.get("lineage", "Unknown")
    lineage_conf = results.get("lineage_confidence", 0)
    quality = results.get("sequence_quality", 0)
    epi = results.get("epi_params", {})
    mutations = results.get("mutations", [])
    prophet_data = load_prophet_data(job_id, variant_name)
 
    # Top mutations by score
    reliable_muts = [m for m in mutations if m.get("Reliability") == "RELIABLE"]
    top_muts = sorted(reliable_muts, key=lambda m: abs(m.get("Score", 0)), reverse=True)[:5]
 
    top_muts_text = "\n".join([
        f"  {m['Mutation']} — zone: {m['Context']}, score: {m['Score']}"
        for m in top_muts
    ]) or "  No reliable mutations detected"

    # lineage_confidence: Convert to float with a safe fallback
    try:
        lineage_conf = float(lineage_conf or 0)
    except (TypeError, ValueError):
        lineage_conf = 0.0
 
    prophet_text = format_prophet_for_llm(prophet_data) if prophet_data else "  Not available"
 
    # Creating a summary as structured text
    summary_text = f"""
    TELOS-S VARIANT INTELLIGENCE REPORT
    =====================================
    Variant: {variant_name}
    Probable Lineage: {lineage} ({lineage_conf}% match)
    Sequencing Quality: {quality:.1f}%
 
    RISK ASSESSMENT
    ---------------
    Aggression Score: {score:.1f}
    Risk Level: {score_to_risk_label(score)}
 
    Interpretation:
    - Scores > 1200 indicate high immune evasion potential (comparable to Omicron BA.1/BA.2 emergence)
    - Scores 600–1200 indicate an active variant of interest requiring monitoring
    - Scores < 600 indicate moderate mutational burden
 
    EPIDEMIOLOGICAL PARAMETERS (ESM-2 derived)
    ------------------------------------------
    Estimated R0: {epi.get('r0_estimated', 'N/A')}
    Incubation period: {epi.get('incubation_period_days', 'N/A')} days
    Base transmissibility: {epi.get('transmissibility_base', 'N/A')}
 
    Note: These parameters are derived computationally from protein structural 
    stability analysis and should be validated with epidemiological field data.
 
    TOP 5 CRITICAL MUTATIONS (reliable positions only)
    ---------------------------------------------------
    {top_muts_text}
 
    TELOS PROPHET — FUTURE EVOLUTION PREDICTIONS
    (ESM-2 structural stability at 4 key positions)
    ---------------------------------------------------
    {prophet_text}
 
    METHODOLOGY NOTE
    ----------------
    Analysis performed with ESM-2 650M (Meta AI) protein language model.
    Mutations near sequencing gaps (X positions ±5 residues) are excluded
    from risk scoring to prevent artifact-driven false positives.
    Data source: Confirmed sequencing only. Imputed positions flagged.
    """.strip()
 
    return {
        "job_id": job_id,
        "variant_name": variant_name,
        "summary": summary_text,
        "key_numbers": {
            "aggression_score": score,
            "r0_estimated": epi.get("r0_estimated"),
            "lineage": lineage,
            "reliable_mutations_count": len(reliable_muts)
        }
    }
 
 
def tool_compare_variants(job_id_a: str, job_id_b: str) -> dict:
    job_a = load_job(job_id_a)
    job_b = load_job(job_id_b)
 
    if job_a["status"] != "completed" or job_b["status"] != "completed":
        return {"error": "Both analyses must be completed to compare"}
 
    res_a = job_a.get("results", {})
    res_b = job_b.get("results", {})
 
    score_a = get_aggression_score(res_a)
    score_b = get_aggression_score(res_b)
    score_diff = score_b - score_a
    score_diff_pct = ((score_b - score_a) / score_a * 100) if score_a > 0 else 0
 
    epi_a = res_a.get("epi_params", {})
    epi_b = res_b.get("epi_params", {})
 
    # Unique mutations in each variant
    muts_a = {m["Mutation"] for m in res_a.get("mutations", []) if m.get("Reliability") == "RELIABLE"}
    muts_b = {m["Mutation"] for m in res_b.get("mutations", []) if m.get("Reliability") == "RELIABLE"}
    shared = muts_a & muts_b
    unique_to_a = muts_a - muts_b
    unique_to_b = muts_b - muts_a
 
    return {
        "variant_a": {
            "job_id": job_id_a,
            "name": res_a.get("variant_name", "Unknown"),
            "aggression_score": score_a,
            "risk_level": score_to_risk_label(score_a),
            "lineage": res_a.get("lineage"),
            "r0_estimated": epi_a.get("r0_estimated"),
        },
        "variant_b": {
            "job_id": job_id_b,
            "name": res_b.get("variant_name", "Unknown"),
            "aggression_score": score_b,
            "risk_level": score_to_risk_label(score_b),
            "lineage": res_b.get("lineage"),
            "r0_estimated": epi_b.get("r0_estimated"),
        },
        "comparison": {
            "score_difference": score_diff,
            "score_difference_pct": round(score_diff_pct, 1),
            "more_aggressive": res_b.get("variant_name") if score_diff > 0 else res_a.get("variant_name"),
            "r0_difference": (epi_b.get("r0_estimated", 0) - epi_a.get("r0_estimated", 0))
                if epi_a.get("r0_estimated") and epi_b.get("r0_estimated") else None,
            "shared_mutations": sorted(shared),
            "unique_to_a": sorted(unique_to_a),
            "unique_to_b": sorted(unique_to_b),
        }
    }
 
 
def tool_get_prophet_predictions(job_id: str) -> dict:
    job = load_job(job_id)
    results = job.get("results", {})
    variant_name = results.get("variant_name", "Unknown")
    prophet_data = load_prophet_data(job_id, variant_name)
 
    if not prophet_data:
        return {
            "job_id": job_id,
            "variant_name": variant_name,
            "predictions": None,
            "message": "Prophet predictions not found. Run analysis with oraculo_mutaciones enabled."
        }
 
    return {
        "job_id": job_id,
        "variant_name": variant_name,
        "aggression_score": get_aggression_score(results),
        "predictions": prophet_data,
        "summary": format_prophet_for_llm(prophet_data)
    }
 
 
# =============================================================================
# TOOLS DISPATCHER
# =============================================================================
 
TOOL_HANDLERS = {
    "get_analysis_results": lambda args: tool_get_analysis_results(**args),
    "list_recent_analyses": lambda args: tool_list_recent_analyses(**{k: v for k, v in args.items() if v is not None}),
    "get_variant_summary": lambda args: tool_get_variant_summary(**args),
    "compare_variants": lambda args: tool_compare_variants(**args),
    "get_prophet_predictions": lambda args: tool_get_prophet_predictions(**args),
}
 
# =============================================================================
# MCP Session Management (Streamable HTTP – spec 2025-03-26)
# OpenWebUI sends the Mcp-Session-Id in each request after initialization.
# We maintain a set of active sessions for validation purposes.
# For local use in memory, this is sufficient; for production, use Redis.
# =============================================================================
_active_sessions: set = set()
 
 
# =============================================================================
# MCP ENDPOINTS (JSON-RPC protocol over HTTP)
# =============================================================================
 
@app.get("/health")
async def health():
    jobs_count = len(list(JOBS_DIR.glob("*.json"))) if JOBS_DIR.exists() else 0
    return {
        "status": "healthy",
        "server": SERVER_INFO["name"],
        "output_dir": str(OUTPUT_DIR),
        "total_analyses": jobs_count,
        "timestamp": datetime.now().isoformat()
    }
 
 
@app.get("/")
async def root():
    return {
        "mcp_server": SERVER_INFO,
        "tools_available": [t["name"] for t in TOOLS],
        "protocol": "MCP over HTTP",
        "docs": "/docs"
    }
 
 
@app.post("/mcp")
async def mcp_endpoint(request: Request):
    """
    Main endpoint MCP – Streamable HTTP (spec 2025-03-26).

    Session management:
      - initialize: creates a session, returns the Mcp-Session-Id in the response header
      - subsequent requests: validates the Mcp-Session-Id from the incoming header
      - OpenWebUI sends the session ID with each request after the initial handshake
    """
    from fastapi.responses import JSONResponse as _JSONResponse
 
    try:
        body = await request.json()
    except Exception:
        raise HTTPException(status_code=400, detail="Invalid JSON body")
 
    method = body.get("method")
    params = body.get("params", {})
    request_id = body.get("id", 1)
 
    # Guard: Request without a 'method' parameter
    if not method:
        return _JSONResponse({
            "jsonrpc": "2.0",
            "id": request_id,
            "error": {
                "code": -32600,
                "message": "Invalid request: missing 'method' field. Expected: initialize, tools/list, tools/call"
            }
        })
 
    # --- initialize — create a new session ---
    if method == "initialize":
        session_id = uuid.uuid4().hex
        _active_sessions.add(session_id)
 
        response_body = {
            "jsonrpc": "2.0",
            "id": request_id,
            "result": {
                "protocolVersion": "2024-11-05",
                "capabilities": {"tools": {}},
                "serverInfo": SERVER_INFO
            }
        }
        # Return the session ID in the header – OpenWebUI reuses it in subsequent requests.
        return _JSONResponse(
            content=response_body,
            headers={"Mcp-Session-Id": session_id}
        )
 
    # --- For all other methods: Optional session ---
    # If "Mcp-Session-Id" is provided, but it's not in memory (e.g., container restarted),
    # we accept it anyway instead of rejecting with a 404 error.
    # OpenWebUI doesn't handle the 404 session error correctly and gets stuck loading.
    incoming_session = request.headers.get("Mcp-Session-Id")
    if incoming_session and incoming_session not in _active_sessions:
        # Instead of rejecting the session, re-register it
        _active_sessions.add(incoming_session)
 
    # --- notifications/initialized — Client acknowledgment, does not require a response ---
    if method == "notifications/initialized":
        return _JSONResponse(content={}, status_code=200)
 
    # --- tools/list ---
    elif method == "tools/list":
        return _JSONResponse({
            "jsonrpc": "2.0",
            "id": request_id,
            "result": {"tools": TOOLS}
        })
 
    # --- tools/call ---
    elif method == "tools/call":
        tool_name = params.get("name")
        tool_args = params.get("arguments", {})
 
        if tool_name not in TOOL_HANDLERS:
            return _JSONResponse({
                "jsonrpc": "2.0",
                "id": request_id,
                "error": {"code": -32601, "message": f"Tool '{tool_name}' not found"}
            })
 
        try:
            result = TOOL_HANDLERS[tool_name](tool_args)
            return _JSONResponse({
                "jsonrpc": "2.0",
                "id": request_id,
                "result": {
                    "content": [{"type": "text", "text": json.dumps(result, indent=2, default=str)}]
                }
            })
        except HTTPException as e:
            return _JSONResponse({
                "jsonrpc": "2.0",
                "id": request_id,
                "error": {"code": -32000, "message": e.detail}
            })
        except Exception as e:
            return _JSONResponse({
                "jsonrpc": "2.0",
                "id": request_id,
                "error": {"code": -32000, "message": str(e)}
            })
 
    # --- Unrecognized method ---
    else:
        return _JSONResponse({
            "jsonrpc": "2.0",
            "id": request_id,
            "error": {"code": -32601, "message": f"Method '{method}' not found"}
        })
 
 
@app.delete("/mcp")
async def mcp_session_terminate(request: Request):
    """
    Explicitly ends an MCP session.
    When the chat is closed, OpenWebUI sends a DELETE /mcp request with the Mcp-Session-Id.
    """
    session_id = request.headers.get("Mcp-Session-Id")
    if session_id and session_id in _active_sessions:
        _active_sessions.discard(session_id)
    return Response(status_code=200)
 
 
@app.get("/mcp")
async def mcp_get(request: Request):
    """
    GET /mcp — MCP streamable HTTP specification (server→client channel).
    OpenWebUI opens this after initialization to receive proactive notifications from the server. 
    We maintain it active with periodic keep-alive requests.
    Without this endpoint, the client receives a 405 error and the session becomes inconsistent.
    """
    import asyncio
 
    incoming_session = request.headers.get("Mcp-Session-Id")
    if incoming_session and incoming_session not in _active_sessions:
        _active_sessions.add(incoming_session)
 
    async def server_events():
        try:
            while True:
                # Send a "keep-alive" message every 15 seconds to prevent proxies and 
                # load balancers from closing the connection due to inactivity.
                await asyncio.sleep(15)
                yield {
                    "event": "ping",
                    "data": "{}"
                }
        except asyncio.CancelledError:
            pass
 
    return EventSourceResponse(server_events())
 
 
@app.get("/mcp/sse")
async def mcp_sse(request: Request):
    """
    Bidirectional SSE endpoint for legacy MCP clients (pre-2025 protocol).
    Compatibility with Claude Desktop and other clients that use SSE transport.

    Flow:
        1. Client GET /mcp/sse → Server sends "endpoint" event with the message URL
        2. Client POST to that URL with JSON-RPC requests
        3. Server responds via SSE with the results
    """
    import asyncio
 
    session_id = uuid.uuid4().hex[:12]
    messages_url = f"/mcp/sse/messages/{session_id}"
 
    # Temporary storage for responses for this session
    # (In production, use Redis; for local use, a dictionary in memory is sufficient)
    if not hasattr(app.state, "sse_sessions"):
        app.state.sse_sessions = {}
    
    queue: asyncio.Queue = asyncio.Queue()
    app.state.sse_sessions[session_id] = queue
 
    async def event_generator():
        # 1. Announce the URL where the client should send their messages.
        yield {
            "event": "endpoint",
            "data": messages_url
        }
 
        # 2. Maintain an open line of communication and relay responses.
        try:
            while True:
                try:
                    message = await asyncio.wait_for(queue.get(), timeout=30.0)
                    yield {
                        "event": "message",
                        "data": json.dumps(message)
                    }
                except asyncio.TimeoutError:
                    # A mechanism to prevent the proxy from disconnecting the connection.
                    yield {
                        "event": "ping",
                        "data": "{}"
                    }
        except asyncio.CancelledError:
            pass
        finally:
            app.state.sse_sessions.pop(session_id, None)
 
    return EventSourceResponse(event_generator())
 
 
@app.post("/mcp/sse/messages/{session_id}")
async def mcp_sse_message(session_id: str, request: Request):
    """
    Receive JSON-RPC messages from the client using SSE and process them, 
    sending the response back over the open SSE connection.
    """
    if not hasattr(app.state, "sse_sessions") or session_id not in app.state.sse_sessions:
        raise HTTPException(status_code=404, detail=f"SSE session '{session_id}' not found")
 
    try:
        body = await request.json()
    except Exception:
        raise HTTPException(status_code=400, detail="Invalid JSON body")
 
    # Reusing the main dispatcher
    method = body.get("method")
    params = body.get("params", {})
    request_id = body.get("id", 1)
 
    # Process it in the same way as the POST /mcp endpoint.
    fake_request = type("R", (), {"json": lambda self: body})()
    
    if method == "initialize":
        response = {"jsonrpc": "2.0", "id": request_id, "result": {
            "protocolVersion": "2024-11-05",
            "capabilities": {"tools": {}},
            "serverInfo": SERVER_INFO
        }}
    elif method == "tools/list":
        response = {"jsonrpc": "2.0", "id": request_id, "result": {"tools": TOOLS}}
    elif method == "tools/call":
        tool_name = params.get("name")
        tool_args = params.get("arguments", {})
        if tool_name not in TOOL_HANDLERS:
            response = {"jsonrpc": "2.0", "id": request_id,
                       "error": {"code": -32601, "message": f"Tool '{tool_name}' not found"}}
        else:
            try:
                result = TOOL_HANDLERS[tool_name](tool_args)
                response = {"jsonrpc": "2.0", "id": request_id,
                           "result": {"content": [{"type": "text", "text": json.dumps(result, indent=2, default=str)}]}}
            except Exception as e:
                response = {"jsonrpc": "2.0", "id": request_id,
                           "error": {"code": -32000, "message": str(e)}}
    else:
        response = {"jsonrpc": "2.0", "id": request_id,
                   "error": {"code": -32601, "message": f"Method '{method}' not found"}}
 
    # Send response via the SSE queue
    queue = app.state.sse_sessions[session_id]
    await queue.put(response)
 
    return {"ok": True}
 
 
# =============================================================================
# MAIN
# =============================================================================
 
if __name__ == "__main__":
    print(f"🧬 Telos-S MCP Server v{SERVER_INFO['version']}")
    print(f"   Output dir: {OUTPUT_DIR}")
    print(f"   Port: {MCP_PORT}")
    print(f"   Tools: {', '.join(t['name'] for t in TOOLS)}")
    print()
 
    uvicorn.run(
        "mcp_server:app",
        host="0.0.0.0",
        port=MCP_PORT,
        reload=False
    )
 