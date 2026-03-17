"""
Telos-S MCP Server
------------------
Model Context Protocol server que expone los resultados del pipeline de
Telos-S a modelos de lenguaje locales (Ollama) vía OpenWebUI.
 
Implementa el protocolo MCP sobre HTTP/SSE, compatible con OpenWebUI >= 0.5
y con cualquier cliente MCP estándar.
 
Tools expuestas:
    - get_analysis_results      → Resultados completos de un análisis por job_id
    - list_recent_analyses      → Lista los N análisis más recientes
    - get_variant_summary       → Resumen ejecutivo en lenguaje natural listo para LLM
    - compare_variants          → Comparación entre dos análisis
    - get_prophet_predictions   → Solo las predicciones de evolución futura
 
Uso standalone (sin Docker):
    pip install fastapi uvicorn sse-starlette
    python mcp_server.py
 
Puerto por defecto: 8001
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
# CONFIGURACIÓN
# =============================================================================
 
OUTPUT_DIR = Path(os.getenv("TELOS_OUTPUT_DIR", "/app/output"))
JOBS_DIR = OUTPUT_DIR / "jobs"
MCP_PORT = int(os.getenv("MCP_PORT", "8001"))
 
# Información del servidor MCP
SERVER_INFO = {
    "name": "telos-s-mcp",
    "version": "1.0.0",
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
    version="1.0.0"
)
 
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)
 
# =============================================================================
# HELPERS — Lectura de datos del pipeline
# =============================================================================
 
def load_job(job_id: str) -> dict:
    """Carga el JSON de estado de un job."""
    job_file = JOBS_DIR / f"{job_id}.json"
    if not job_file.exists():
        raise HTTPException(status_code=404, detail=f"Job '{job_id}' not found")
    with open(job_file) as f:
        return json.load(f)
 
 
def load_prophet_data(variant_name: str) -> Optional[list]:
    """Carga las predicciones Prophet para una variante."""
    prophet_path = OUTPUT_DIR / "prophet" / f"mutation_predictions_spike_{variant_name}.json"
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
    Lee aggression_score de donde esté disponible.
    Jobs nuevos: results["aggression_score"] (promovido a raíz por el backend)
    Jobs viejos: results["epi_params"]["aggression_score"]
    """
    score = results.get("aggression_score")
    if score is None or score == 0:
        score = results.get("epi_params", {}).get("aggression_score", 0.0)
    return float(score or 0.0)
 
 
def format_prophet_for_llm(prophet_data: list) -> str:
    """Convierte predicciones Prophet a texto estructurado para el LLM."""
    if not prophet_data:
        return "No Prophet predictions available."
 
    lines = []
    for target in prophet_data:
        pos = target.get("detected_position")
        original = target.get("original", "?")
        name = target.get("target", f"Position {pos}")
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
# IMPLEMENTACIÓN DE TOOLS
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
    prophet_data = load_prophet_data(variant_name)
 
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
 
    for job_file in job_files[:100]:  # Leer hasta 100 para filtrar
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
    prophet_data = load_prophet_data(variant_name)
 
    # Top mutations por score
    reliable_muts = [m for m in mutations if m.get("confidence") == "CONFIABLE"]
    top_muts = sorted(reliable_muts, key=lambda m: abs(m.get("score", 0)), reverse=True)[:5]
 
    top_muts_text = "\n".join([
        f"  {m['mutation']} — zone: {m['zone']}, score: {m['score']:.1f}"
        for m in top_muts
    ]) or "  No reliable mutations detected"
 
    prophet_text = format_prophet_for_llm(prophet_data) if prophet_data else "  Not available"
 
    # Construir el summary como texto estructurado
    summary_text = f"""
TELOS-S VARIANT INTELLIGENCE REPORT
=====================================
Variant: {variant_name}
Probable Lineage: {lineage} ({lineage_conf:.0f}% match)
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
 
    # Mutaciones únicas en cada variante
    muts_a = {m["mutation"] for m in res_a.get("mutations", []) if m.get("confidence") == "CONFIABLE"}
    muts_b = {m["mutation"] for m in res_b.get("mutations", []) if m.get("confidence") == "CONFIABLE"}
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
    prophet_data = load_prophet_data(variant_name)
 
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
# DISPATCHER DE TOOLS
# =============================================================================
 
TOOL_HANDLERS = {
    "get_analysis_results": lambda args: tool_get_analysis_results(**args),
    "list_recent_analyses": lambda args: tool_list_recent_analyses(**{k: v for k, v in args.items() if v is not None}),
    "get_variant_summary": lambda args: tool_get_variant_summary(**args),
    "compare_variants": lambda args: tool_compare_variants(**args),
    "get_prophet_predictions": lambda args: tool_get_prophet_predictions(**args),
}
 
# =============================================================================
# GESTIÓN DE SESIONES MCP (Streamable HTTP — spec 2025-03-26)
# OpenWebUI envía Mcp-Session-Id en cada request después del initialize.
# Mantenemos un set de sesiones activas para validarlas.
# Para uso local en memoria es suficiente; para producción usar Redis.
# =============================================================================
_active_sessions: set = set()
 
 
# =============================================================================
# ENDPOINTS MCP (protocolo JSON-RPC sobre HTTP)
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
    Endpoint principal MCP — Streamable HTTP (spec 2025-03-26).
 
    Gestión de sesiones:
      - initialize:  crea sesión, devuelve Mcp-Session-Id en el header de respuesta
      - resto:       valida Mcp-Session-Id del header entrante
      - OpenWebUI envía el session ID en cada request tras el handshake inicial
    """
    from fastapi.responses import JSONResponse as _JSONResponse
 
    try:
        body = await request.json()
    except Exception:
        raise HTTPException(status_code=400, detail="Invalid JSON body")
 
    method = body.get("method")
    params = body.get("params", {})
    request_id = body.get("id", 1)
 
    # Guard: request sin method
    if not method:
        return _JSONResponse({
            "jsonrpc": "2.0",
            "id": request_id,
            "error": {
                "code": -32600,
                "message": "Invalid request: missing 'method' field. Expected: initialize, tools/list, tools/call"
            }
        })
 
    # --- initialize — crear sesión nueva ---
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
        # Devolver session ID en header — OpenWebUI lo reutiliza en requests siguientes
        return _JSONResponse(
            content=response_body,
            headers={"Mcp-Session-Id": session_id}
        )
 
    # --- Para todos los demás métodos: sesión opcional ---
    # Si viene Mcp-Session-Id pero no está en memoria (ej: contenedor reiniciado),
    # lo aceptamos igual en lugar de rechazar con 404.
    # OpenWebUI no maneja el 404 de sesión correctamente y se queda cargando.
    incoming_session = request.headers.get("Mcp-Session-Id")
    if incoming_session and incoming_session not in _active_sessions:
        # Re-registrar la sesión en lugar de rechazarla
        _active_sessions.add(incoming_session)
 
    # --- notifications/initialized — ACK del cliente, no requiere respuesta ---
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
 
    # --- Método no reconocido ---
    else:
        return _JSONResponse({
            "jsonrpc": "2.0",
            "id": request_id,
            "error": {"code": -32601, "message": f"Method '{method}' not found"}
        })
 
 
@app.delete("/mcp")
async def mcp_session_terminate(request: Request):
    """
    Termina una sesión MCP explícitamente.
    OpenWebUI envía DELETE /mcp con Mcp-Session-Id al cerrar el chat.
    """
    session_id = request.headers.get("Mcp-Session-Id")
    if session_id and session_id in _active_sessions:
        _active_sessions.discard(session_id)
    return Response(status_code=200)
 
 
@app.get("/mcp")
async def mcp_get(request: Request):
    """
    GET /mcp — canal SSE servidor→cliente (MCP Streamable HTTP spec).
    OpenWebUI lo abre después del initialize para recibir notificaciones
    proactivas del servidor. Lo mantenemos vivo con keepalives periódicos.
    Sin este endpoint el cliente recibe 405 y la sesión queda inconsistente.
    """
    import asyncio
 
    incoming_session = request.headers.get("Mcp-Session-Id")
    if incoming_session and incoming_session not in _active_sessions:
        _active_sessions.add(incoming_session)
 
    async def server_events():
        try:
            while True:
                # Keepalive cada 15s — evita que proxies y load balancers
                # cierren la conexión por inactividad
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
    Endpoint SSE bidireccional para clientes MCP legacy (protocolo pre-2025).
    Compatibilidad con Claude Desktop y otros clientes que usen SSE transport.
 
    Flujo:
      1. Cliente GET /mcp/sse → servidor envía evento "endpoint" con la URL de mensajes
      2. Cliente POST a esa URL con requests JSON-RPC
      3. Servidor responde via SSE con los resultados
    """
    import asyncio
 
    session_id = uuid.uuid4().hex[:12]
    messages_url = f"/mcp/sse/messages/{session_id}"
 
    # Almacén temporal de respuestas para esta sesión
    # (en producción usar Redis; para uso local un dict en memoria es suficiente)
    if not hasattr(app.state, "sse_sessions"):
        app.state.sse_sessions = {}
    
    queue: asyncio.Queue = asyncio.Queue()
    app.state.sse_sessions[session_id] = queue
 
    async def event_generator():
        # 1. Anunciar la URL donde el cliente debe enviar sus mensajes
        yield {
            "event": "endpoint",
            "data": messages_url
        }
 
        # 2. Mantener la conexión abierta y retransmitir respuestas
        try:
            while True:
                try:
                    message = await asyncio.wait_for(queue.get(), timeout=30.0)
                    yield {
                        "event": "message",
                        "data": json.dumps(message)
                    }
                except asyncio.TimeoutError:
                    # Keepalive para evitar que el proxy cierre la conexión
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
    Recibe mensajes JSON-RPC del cliente SSE y los procesa,
    enviando la respuesta de vuelta por la conexión SSE abierta.
    """
    if not hasattr(app.state, "sse_sessions") or session_id not in app.state.sse_sessions:
        raise HTTPException(status_code=404, detail=f"SSE session '{session_id}' not found")
 
    try:
        body = await request.json()
    except Exception:
        raise HTTPException(status_code=400, detail="Invalid JSON body")
 
    # Reusar el dispatcher principal
    method = body.get("method")
    params = body.get("params", {})
    request_id = body.get("id", 1)
 
    # Procesar igual que el endpoint POST /mcp
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
 
    # Enviar respuesta por la cola SSE
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
 