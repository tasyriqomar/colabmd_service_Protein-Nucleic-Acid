from fastapi import FastAPI, UploadFile, File, Form, HTTPException, APIRouter
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from pydantic import BaseModel
from celery.result import AsyncResult
from pathlib import Path
from typing import Optional
import shutil
import uuid
import os
import time

from tasks.worker import (
    celery_app,
    run_aptamer_task,
    run_hdock_docking_task,
    run_haddock_docking_task,
    run_complex_cg_prep_task,
    run_backmapping_workflow,
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

SHARED_DATA = Path("/content/drive/MyDrive/ColabMD_Data/shared_data")
SHARED_DATA.mkdir(parents=True, exist_ok=True)

router = APIRouter()

class CleanupRequest(BaseModel):
    username: str
    tabs: list[int] = [1, 3, 4]
    task_ids: list[str] = []

@router.post("/cleanup-workspace")
def cleanup_workspace(req: CleanupRequest):
    # 1. Terminate any lingering Celery tasks so worker slots (--concurrency=2) are freed
    for tid in req.task_ids:
        celery_app.control.revoke(tid, terminate=True, signal="SIGKILL")

    # 2. Remove Tab 1, 3, 4 artifacts for this user in shared_data
    base_dir = os.path.join("/app/shared_data", req.username)
    tab_folders = {
        1: ["tab1", "preparation"],
        3: ["tab3", "cgmd_setup"],
        4: ["tab4", "simulation", "cgmd_run"],
    }

    deleted = []
    for tab in req.tabs:
        for folder in tab_folders.get(tab, []):
            target = os.path.join(base_dir, folder)
            if os.path.exists(target):
                shutil.rmtree(target, ignore_errors=True)
                deleted.append(target)

    return {"status": "cleaned", "deleted": deleted}

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

@app.post("/api/v1/run-hdock-docking")
async def run_hdock_docking(
    protein_pdb: UploadFile = File(...),
    dna_pdb: UploadFile = File(...),
    nmax: int = Form(100),
    interaction_cutoff: float = Form(4.5),
):
    job_id = f"hdock_{int(time.time())}"
    work_dir = f"/app/shared_data/{job_id}"
    os.makedirs(work_dir, exist_ok=True)

    prot_bytes = await protein_pdb.read()
    dna_bytes = await dna_pdb.read()

    task = run_hdock_docking_task.delay(work_dir, prot_bytes, dna_bytes, nmax, interaction_cutoff)
    return {"message": "HDOCK docking queued.", "task_id": task.id}

@app.post("/api/v1/run-haddock-docking")
async def api_docking(protein_pdb: UploadFile = File(...), dna_pdb: UploadFile = File(...), sampling: int = Form(20), cmrest: bool = Form(True)):
    job_id = f"dock_{uuid.uuid4().hex[:10]}"
    work_dir = SHARED_DATA / job_id
    work_dir.mkdir(parents=True, exist_ok=True)
    p_bytes = await protein_pdb.read()
    d_bytes = await dna_pdb.read()
    task = run_haddock_docking_task.delay(str(work_dir), p_bytes, d_bytes, sampling, cmrest)
    return {"task_id": task.id, "job_id": job_id}

@app.post("/api/v1/run-complex-cg-prep")
async def api_complex_cg_prep(
    complex_pdb: UploadFile = File(...),
    dnatype: str = Form("ss"),
    elastic: bool = Form(True)
):
    job_id = f"cgprep_{uuid.uuid4().hex[:10]}"
    work_dir = SHARED_DATA / job_id
    work_dir.mkdir(parents=True, exist_ok=True)
    content = await complex_pdb.read()
    task = run_complex_cg_prep_task.delay(str(work_dir), content, dnatype, elastic)
    return {"task_id": task.id, "job_id": job_id}

@app.post("/api/v1/run-backmapping")
async def api_backmapping(
    protein_pdb: UploadFile = File(...),
    aptamer_pdb: UploadFile = File(...),
    xtc_file: UploadFile = File(...),
    tpr_file: UploadFile = File(...),
    mdp_file: UploadFile = File(...),
    ndx_file: UploadFile = File(...),
    system_gro: UploadFile = File(...),
    aptamer_map: Optional[UploadFile] = File(None),
    topol_top: Optional[UploadFile] = File(None),
    loop_from: int = Form(0),
    loop_until: int = Form(10),
    solute_group: int = Form(1)
):
    job_id = f"bm_{uuid.uuid4().hex[:10]}"
    work_dir = SHARED_DATA / job_id
    work_dir.mkdir(parents=True, exist_ok=True)

    # Save uploaded files into the job directory
    file_map = {
        protein_pdb.filename: await protein_pdb.read(),
        aptamer_pdb.filename: await aptamer_pdb.read(),
        xtc_file.filename: await xtc_file.read(),
        tpr_file.filename: await tpr_file.read(),
        mdp_file.filename: await mdp_file.read(),
        ndx_file.filename: await ndx_file.read(),
        system_gro.filename: await system_gro.read()
    }
    for filename, content in file_map.items():
        with open(work_dir / filename, "wb") as f:
            f.write(content)

    map_name = None
    if aptamer_map:
        map_name = aptamer_map.filename
        with open(work_dir / map_name, "wb") as f:
            f.write(await aptamer_map.read())

    top_name = None
    if topol_top:
        top_name = topol_top.filename
        with open(work_dir / top_name, "wb") as f:
            f.write(await topol_top.read())

    task = run_backmapping_workflow.delay(
        str(work_dir),
        protein_pdb.filename,
        aptamer_pdb.filename,
        xtc_file.filename,
        tpr_file.filename,
        mdp_file.filename,
        ndx_file.filename,
        system_gro.filename,
        map_name,
        top_name,
        loop_from,
        loop_until,
        solute_group
    )
    return {"message": "Backmapping queued.", "task_id": task.id, "job_id": job_id}

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
