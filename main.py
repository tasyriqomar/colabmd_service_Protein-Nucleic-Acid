from fastapi import FastAPI, UploadFile, File, Form, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from celery.result import AsyncResult
from pathlib import Path
import uuid
import os

from tasks.worker import (
    celery_app,
    run_aptamer_task,
    run_haddock_docking_task,
    run_backmapping_task,
    run_trajectory_analysis_task,
    run_plip_analysis_task,
    run_mmpbsa_task,
    run_pca_fel_task
)

app = FastAPI(title="ColabMD Protein-Nucleic Acid Service", version="2.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

SHARED_DATA = Path("/app/shared_data")
SHARED_DATA.mkdir(parents=True, exist_ok=True)

@app.get("/api/v1/status/{task_id}")
def get_status(task_id: str):
    res = AsyncResult(task_id, app=celery_app)
    response = {"task_id": task_id, "status": res.status}
    if res.ready():
        response["result"] = res.result
    return response

@app.get("/api/v1/download-prep-result/{task_id}")
def download_result(task_id: str):
    res = AsyncResult(task_id, app=celery_app)
    if not res.ready():
        raise HTTPException(status_code=400, detail="Task not complete")
    result = res.result
    if isinstance(result, dict) and "result_file" in result:
        file_path = result["result_file"]
        if os.path.exists(file_path):
            return FileResponse(file_path, filename=Path(file_path).name)
    raise HTTPException(status_code=404, detail="File not found")

@app.post("/api/v1/run-aptamer")
async def api_aptamer(aptamer_pdb: UploadFile = File(...), convert_rna: bool = Form(True)):
    job_id = f"apt_{uuid.uuid4().hex[:10]}"
    work_dir = SHARED_DATA / job_id
    work_dir.mkdir(parents=True, exist_ok=True)
    content = await aptamer_pdb.read()
    
    # Save both casing variants directly before Celery picks it up
    with open(work_dir / "Aptamer.pdb", "wb") as f:
        f.write(content)
    with open(work_dir / "aptamer.pdb", "wb") as f:
        f.write(content)

    task = run_aptamer_task.delay(str(work_dir), content, convert_rna)
    return {"task_id": task.id, "job_id": job_id}

@app.post("/api/v1/run-haddock-docking")
async def api_docking(protein_pdb: UploadFile = File(...), dna_pdb: UploadFile = File(...), sampling: int = Form(20), cmrest: bool = Form(True)):
    job_id = f"dock_{uuid.uuid4().hex[:10]}"
    work_dir = SHARED_DATA / job_id
    work_dir.mkdir(parents=True, exist_ok=True)
    p_bytes = await protein_pdb.read()
    d_bytes = await dna_pdb.read()
    task = run_haddock_docking_task.delay(str(work_dir), p_bytes, d_bytes, sampling, cmrest)
    return {"task_id": task.id, "job_id": job_id}

@app.post("/api/v1/run-backmapping")
async def api_backmapping(
    xtc: UploadFile = File(...), tpr: UploadFile = File(...),
    mdp: UploadFile = File(...), ndx: UploadFile = File(...),
    fixed_pdb: UploadFile = File(...), aptamer_pdb: UploadFile = File(...)
):
    job_id = f"bm_{uuid.uuid4().hex[:10]}"
    work_dir = SHARED_DATA / job_id
    work_dir.mkdir(parents=True, exist_ok=True)
    task = run_backmapping_task.delay(
        str(work_dir), await xtc.read(), await tpr.read(),
        await mdp.read(), await ndx.read(), await fixed_pdb.read(), await aptamer_pdb.read()
    )
    return {"task_id": task.id, "job_id": job_id}

@app.post("/api/v1/analysis/trajectory")
def api_trajectory(payload: dict):
    task = run_trajectory_analysis_task.delay(
        payload["job_id"],
        payload.get("rmsd_group", 4),
        payload.get("rmsf_group", 1),
        payload.get("rg_group", 1)
    )
    return {"task_id": task.id}

@app.post("/api/v1/analysis/plip")
def api_plip(payload: dict):
    task = run_plip_analysis_task.delay(payload["job_id"])
    return {"task_id": task.id}

@app.post("/api/v1/analysis/mmpbsa")
def api_mmpbsa(payload: dict):
    task = run_mmpbsa_task.delay(
        payload["job_id"],
        payload.get("protein_idx", 1),
        payload.get("ligand_idx", 12)
    )
    return {"task_id": task.id}

@app.post("/api/v1/analysis/pca-fel")
def api_pca_fel(payload: dict):
    task = run_pca_fel_task.delay(payload["job_id"], payload.get("backbone_group", 4))
    return {"task_id": task.id}
