from fastapi import APIRouter, UploadFile, File, Form, HTTPException
from fastapi.responses import FileResponse
from pydantic import BaseModel
from typing import Optional
import time
import os
from pathlib import Path
from tasks.worker import (
    generate_mol2_task,
    run_aa_prep_task,
    run_cg_prep_task,
    run_aptamer_task,
    run_cgmd_prep_workflow,
    run_backmapping_workflow,
    celery_app,
)

router = APIRouter()
SHARED_DIR = Path("/app/shared_data")


def clean_aptamer_pdb(input_pdb: str, output_pdb: str, aptamer_type: str = "DNA"):
    """
    Standardizes residue and atom naming for DNA or RNA aptamers.
    Ensures prime notation (') is used instead of asterisks (*).
    """
    dna_map = {
        "A": "DA", "T": "DT", "U": "DT", "G": "DG", "C": "DC",
        "DA": "DA", "DT": "DT", "DG": "DG", "DC": "DC",
        "RA": "DA", "RU": "DT", "RG": "DG", "RC": "DC",
    }
    rna_map = {
        "A": "A", "U": "U", "G": "G", "C": "C",
        "RA": "A", "RU": "U", "RG": "G", "RC": "C",
    }

    with open(input_pdb, "r") as infile, open(output_pdb, "w") as outfile:
        for line in infile:
            if not line.startswith(("ATOM", "HETATM")):
                continue

            res_name = line[17:20].strip()
            atom_name = line[12:16].strip()

            # Standardize sugar atom names (replace * with ')
            if "*" in atom_name:
                atom_name = atom_name.replace("*", "'")

            if aptamer_type.upper() == "DNA":
                # Strip 2'-hydroxyl if an RNA structure was loaded as DNA
                if atom_name in ["O2'", "HO2'"]:
                    continue
                # Map to DA, DT, DC, DG
                new_res = dna_map.get(res_name, res_name)
            else:
                # RNA processing: preserve O2' and map to standard RNA naming
                new_res = rna_map.get(res_name, res_name)

            # Reconstruct PDB record line
            atom_formatted = f" {atom_name:<3}" if len(atom_name) < 4 else f"{atom_name:<4}"
            res_formatted = f"{new_res:>3}"
            new_line = line[:12] + atom_formatted + " " + res_formatted + line[20:]
            outfile.write(new_line)


@router.post("/run-aptamer")
async def run_aptamer(
    aptamer_pdb: UploadFile = File(...),
    aptamer_type: str = Form("DNA"),
    convert_rna: bool = Form(True),
    session_id: Optional[str] = Form(None),
):
    job_id = f"aptamer_{int(time.time())}"
    wd = SHARED_DIR / job_id
    wd.mkdir(parents=True, exist_ok=True)

    input_path = wd / aptamer_pdb.filename
    output_path = wd / "aptamer_fixed.pdb"

    with open(input_path, "wb") as f:
        f.write(await aptamer_pdb.read())

    # Enforce False for RNA mode so RNA structures are never converted to DNA
    should_convert_rna = convert_rna if aptamer_type.upper() == "DNA" else False

    task = run_aptamer_task.delay(
        str(wd),
        str(input_path),
        str(output_path),
        aptamer_type,
        should_convert_rna,
    )
    return {"message": f"{aptamer_type.upper()} aptamer prep queued.", "task_id": task.id}


# 1. Add ligand_name to the Mol2 Request
class Mol2Request(BaseModel):
    task_id: str
    ligand_name: str


@router.post("/generate-mol2")
def generate_mol2(req: Mol2Request):
    docking_task = celery_app.AsyncResult(req.task_id)
    if not docking_task.ready() or docking_task.result.get("status") != "success":
        raise HTTPException(status_code=400, detail="Docking incomplete.")

    ligaa_pdb_path = docking_task.result.get("result_file")

    job_id = f"mol2_{int(time.time())}"
    wd = SHARED_DIR / job_id
    wd.mkdir(parents=True, exist_ok=True)

    # Pass the dynamic ligand_name here
    task = generate_mol2_task.delay(str(wd), ligaa_pdb_path, req.ligand_name)
    return {"message": "MOL2 queued.", "task_id": task.id}


@router.post("/run-aa-prep")
async def run_aa_prep(
    ligand_name: str = Form(...),  # 2. Accept ligand_name as a Form field
    mol2_file: UploadFile = File(...),
    str_file: UploadFile = File(...),
):
    job_id = f"aa_{int(time.time())}"
    wd = SHARED_DIR / job_id
    wd.mkdir(parents=True, exist_ok=True)

    mol2_path = wd / mol2_file.filename
    str_path = wd / str_file.filename

    with open(mol2_path, "wb") as f:
        f.write(await mol2_file.read())
    with open(str_path, "wb") as f:
        f.write(await str_file.read())

    # Pass the dynamic ligand_name here
    task = run_aa_prep_task.delay(str(wd), str(mol2_path), str(str_path), ligand_name)
    return {"message": "AA prep queued.", "task_id": task.id}


@router.post("/run-cg-prep")
async def run_cg_prep(
    ligand_name: str = Form(...),
    smiles: str = Form(...),
):
    job_id = f"cg_{int(time.time())}"
    wd = SHARED_DIR / job_id
    wd.mkdir(parents=True, exist_ok=True)

    task = run_cg_prep_task.delay(str(wd), smiles, ligand_name)
    return {"message": "CG prep queued.", "task_id": task.id}


@router.get("/download-prep-result/{task_id}")
def download_prep_result(task_id: str):
    task_result = celery_app.AsyncResult(task_id)
    if task_result.ready() and task_result.result.get("status") == "success":
        file_path = task_result.result.get("result_file")
        if file_path and os.path.exists(file_path):
            if file_path.endswith(".zip"):
                mt = "application/zip"
            elif file_path.endswith(".pdb"):
                mt = "chemical/x-pdb"
            else:
                mt = "chemical/x-mol2"
            return FileResponse(path=file_path, filename=Path(file_path).name, media_type=mt)
    return {"error": "File not ready or task failed."}


class CGMDPrepRequest(BaseModel):
    job_id: str
    protein_pdb: str
    ligand_itp: str


@router.post("/run-cgmd-prep")
async def start_cgmd_prep(
    protein_pdb: UploadFile = File(...),
    ligand_itp: UploadFile = File(...),
    ligand_gro: UploadFile = File(...),
    ligand_code: str = Form(...),
):
    # 1. Create a fresh, unique folder for this CGMD job
    timestamp = str(int(time.time()))
    work_dir = f"/app/shared_data/cgmd_{timestamp}"
    os.makedirs(work_dir, exist_ok=True)

    # 2. Save the uploaded files
    for upload in [protein_pdb, ligand_itp, ligand_gro]:
        file_path = os.path.join(work_dir, upload.filename)
        with open(file_path, "wb") as f:
            f.write(await upload.read())

    # 3. Dispatch to the Celery worker
    task = run_cgmd_prep_workflow.delay(
        work_dir,
        protein_pdb.filename,
        ligand_itp.filename,
        ligand_gro.filename,
        ligand_code,
    )

    return {"message": "CGMD Prep queued", "task_id": task.id}


@router.post("/run-backmapping")
async def start_backmapping(
    protein_pdb: UploadFile = File(...),
    ligand_map: UploadFile = File(...),
    xtc_file: UploadFile = File(...),
    tpr_file: UploadFile = File(...),
    mdp_file: UploadFile = File(...),
    ndx_file: UploadFile = File(...),
    lig_itp_file: UploadFile = File(...),
    lig_prm_file: UploadFile = File(...),
    lig_cg_itp_file: UploadFile = File(...),
    system_gro: UploadFile = File(...),
    ligand_code: str = Form(...),
    loop_from: int = Form(...),
    loop_until: int = Form(...),
):
    timestamp = str(int(time.time()))
    work_dir = f"/app/shared_data/backmap_{timestamp}"
    os.makedirs(work_dir, exist_ok=True)

    # Save all files to the working directory
    file_list = [
        protein_pdb, ligand_map, xtc_file, tpr_file, mdp_file,
        ndx_file, lig_itp_file, lig_prm_file, lig_cg_itp_file, system_gro,
    ]
    for upload in file_list:
        file_path = os.path.join(work_dir, upload.filename)
        with open(file_path, "wb") as f:
            f.write(await upload.read())

    task = run_backmapping_workflow.delay(
        work_dir, protein_pdb.filename, ligand_map.filename,
        xtc_file.filename, tpr_file.filename, mdp_file.filename, ndx_file.filename,
        lig_itp_file.filename, lig_prm_file.filename, lig_cg_itp_file.filename, system_gro.filename,
        ligand_code, loop_from, loop_until,
    )

    return {"message": "Backmapping queued", "task_id": task.id}
