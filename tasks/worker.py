from celery import Celery
import subprocess
from pathlib import Path
import os
import zipfile
import re
import shutil
import glob
import sys

# Configure Redis connection with dynamic fallback for single-container hosting
REDIS_URL = os.getenv("REDIS_URL", "redis://127.0.0.1:6379/0")

celery_app = Celery(
    "colabmd_tasks",
    broker=REDIS_URL,
    backend=REDIS_URL
)

# ==============================================================================
# 1. APTAMER & NUCLEIC ACID PROCESSING
# ==============================================================================

@celery_app.task
def run_aptamer_processing_task(work_dir: str, aptamer_pdb_content: bytes, convert_rna: bool):
    """Converts RNA to DNA and renames residues using local repository scripts."""
    wd = Path(work_dir)
    try:
        os.chdir("/app")
    except Exception:
        pass

    try:
        with open(wd / "aptamer.pdb", "wb") as f:
            f.write(aptamer_pdb_content)
            
        shutil.copy("/app/tools/Aptamer-Processing/rna_to_dna_with_methyl.py", wd / "rna_to_dna_with_methyl.py")
        shutil.copy("/app/tools/Aptamer-Processing/rename_dna.py", wd / "rename_dna.py")
        
        if convert_rna:
            subprocess.run([sys.executable, "rna_to_dna_with_methyl.py"], cwd=wd, check=True, capture_output=True, text=True)
            
        subprocess.run([sys.executable, "rename_dna.py"], cwd=wd, check=True, capture_output=True, text=True)
        
        fixed_pdb = wd / "aptamer_fixed.pdb"
        if not fixed_pdb.exists():
            return {"status": "failed", "message": "Failed to generate aptamer_fixed.pdb"}
            
        return {"status": "success", "result_file": str(fixed_pdb)}
    except subprocess.CalledProcessError as e:
        err = e.stderr if e.stderr else e.stdout
        return {"status": "failed", "message": f"Aptamer Processing Error: {err}"}
    except Exception as e:
        return {"status": "failed", "message": str(e)}

# ==============================================================================
# 2. CGMD PREPARATION (PROTEIN + APTAMER)
# ==============================================================================

@celery_app.task
def run_cgmd_prep_task(work_dir: str, protein_pdb_content: bytes, aptamer_pdb_content: bytes):
    """Coarse-grains Protein and DNA Aptamer, solvates with insane.py, and creates CG system."""
    wd = Path(work_dir)
    try:
        os.chdir("/app")
    except Exception:
        pass

    try:
        with open(wd / "protein.pdb", "wb") as f: f.write(protein_pdb_content)
        with open(wd / "aptamer_fixed.pdb", "wb") as f: f.write(aptamer_pdb_content)

        shutil.copytree("/app/tools/CGMD", wd / "CGMD", dirs_exist_ok=True)

        subprocess.run(["gmx", "editconf", "-f", "aptamer_fixed.pdb", "-o", "aptamer_aa.gro"], cwd=wd, check=True, capture_output=True)
        subprocess.run(["python2", "CGMD/martinize-dna.py", "-dnatype", "ss", "-f", "aptamer_aa.gro", "-o", "dna.top", "-x", "dna_CG.pdb"], cwd=wd, check=True, capture_output=True)
        subprocess.run(["gmx", "editconf", "-f", "dna_CG.pdb", "-o", "aptamer_CG.gro"], cwd=wd, check=True, capture_output=True)

        subprocess.run(["python2", "CGMD/martinize.py", "-f", "protein.pdb", "-o", "system.top", "-x", "CG.pdb", "-p", "backbone", "-ff", "martini22", "-elastic"], cwd=wd, check=True, capture_output=True)
        subprocess.run(["gmx", "editconf", "-f", "CG.pdb", "-o", "CG.gro"], cwd=wd, check=True, capture_output=True)

        with open(wd / "system.top", "r") as f:
            top = f.read().replace('#include "martini.itp"', '#include "CGMD/martini_v2.1P-dna.itp"\n#include "CGMD/martini_v2.0_ions.itp"')
        with open(wd / "system.top", "w") as f:
            f.write(top)

        subprocess.run(["gmx", "insert-molecules", "-f", "CG.gro", "-ci", "aptamer_CG.gro", "-box", "10", "10", "10", "-nmol", "1", "-o", "CG_aptamer.gro"], cwd=wd, check=True, capture_output=True)
        subprocess.run(["python2", "CGMD/insane.py", "-f", "CG_aptamer.gro", "-o", "system.gro", "-pbc", "cubic", "-box", "10,10,10", "-salt", "0.15", "-charge", "auto", "-sol", "W"], cwd=wd, check=True, capture_output=True)

        with open(wd / "system.top", "r") as f:
            top = f.read().replace("NA+", "NA").replace("CL-", "CL")
        with open(wd / "system.top", "w") as f:
            f.write(top)

        zip_path = wd / "cgmd_prep.zip"
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for fname in ["system.top", "system.gro", "CG.pdb", "aptamer_CG.gro"]:
                if (wd / fname).exists(): zipf.write(wd / fname, arcname=fname)

        return {"status": "success", "result_file": str(zip_path)}
    except subprocess.CalledProcessError as e:
        err = e.stderr if e.stderr else e.stdout
        return {"status": "failed", "message": f"CGMD Prep Error: {err}"}
    except Exception as e:
        return {"status": "failed", "message": str(e)}

# ==============================================================================
# 3. BACKMAPPING TASK
# ==============================================================================

@celery_app.task
def run_backmapping_task(work_dir: str, xtc: bytes, tpr: bytes, mdp: bytes, ndx: bytes, fixed_pdb: bytes, apt_pdb: bytes):
    """Executes the automated Backmapping loop using Dr_Tasyriq_V4_M2.sh."""
    wd = Path(work_dir)
    try:
        os.chdir("/app")
    except Exception:
        pass

    try:
        shutil.copytree("/app/tools/Backmapping", wd / "Backmapping", dirs_exist_ok=True)
        
        with open(wd / "Backmapping" / "dynamic_fitted.xtc", "wb") as f: f.write(xtc)
        with open(wd / "Backmapping" / "dynamic.tpr", "wb") as f: f.write(tpr)
        with open(wd / "Backmapping" / "dynamic.mdp", "wb") as f: f.write(mdp)
        with open(wd / "Backmapping" / "index.ndx", "wb") as f: f.write(ndx)
        with open(wd / "Backmapping" / "fixed.pdb", "wb") as f: f.write(fixed_pdb)
        with open(wd / "Backmapping" / "aptamer_relaxed.pdb", "wb") as f: f.write(apt_pdb)

        bm_dir = wd / "Backmapping"

        # Combine PDB structures
        subprocess.run('cat fixed.pdb aptamer_relaxed.pdb | grep -v -E "TITLE|MODEL|ENDMDL|CRYST1|REMARK|END" > combined.pdb', shell=True, cwd=bm_dir, check=True)
        subprocess.run('echo "END" >> combined.pdb', shell=True, cwd=bm_dir, check=True)
        subprocess.run(["gmx", "pdb2gmx", "-f", "combined.pdb", "-o", "aa.gro", "-ignh"], cwd=bm_dir, input="1\n", check=True, text=True, capture_output=True)

        # Set permissions and execute looping script
        subprocess.run(["chmod", "+x", "Dr_Tasyriq_V4_M2.sh", "Dr_Tasyriq_V4_M2_Redo.sh", "initram-v5.sh", "backward.py"], cwd=bm_dir, check=True)
        subprocess.run(["bash", "Dr_Tasyriq_V4_M2.sh"], cwd=bm_dir, check=True, capture_output=True, text=True)
        
        # Concatenate XTC trajectories
        xtc_dir = bm_dir / "all_xtc"
        subprocess.run(["gmx", "trjcat", "-f", "*.xtc", "-o", "all.xtc"], shell=True, cwd=xtc_dir, check=True)
        
        if (bm_dir / "gro" / "aa_0ns.gro").exists():
            shutil.copy(bm_dir / "gro" / "aa_0ns.gro", bm_dir / "aa_0ns.gro")
        
        subprocess.run(["gmx", "grompp", "-f", "dynamic.mdp", "-c", "aa_0ns.gro", "-p", "topol.top", "-o", "dynamic.tpr", "-maxwarn", "5"], cwd=bm_dir, check=True, capture_output=True)
        shutil.copy(bm_dir / "dynamic.tpr", xtc_dir / "dynamic.tpr")
        
        subprocess.run(["gmx", "make_ndx", "-f", "dynamic.tpr", "-o", "index.ndx"], input="q\n", cwd=xtc_dir, check=True, text=True, capture_output=True)
        subprocess.run(["gmx", "trjconv", "-s", "dynamic.tpr", "-n", "index.ndx", "-f", "all.xtc", "-dt", "1000", "-o", "concatenated.pdb"], input="0\n", cwd=xtc_dir, check=True, text=True, capture_output=True)

        zip_path = wd / "backmap_results.zip"
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as z:
            for fname in ["all.xtc", "dynamic.tpr", "index.ndx", "concatenated.pdb"]:
                if (xtc_dir / fname).exists(): 
                    z.write(xtc_dir / fname, arcname=fname)

        return {"status": "success", "result_file": str(zip_path)}
    except subprocess.CalledProcessError as e:
        err = e.stderr if e.stderr else e.stdout
        return {"status": "failed", "message": f"Backmapping Error: {err}"}
    except Exception as e:
        return {"status": "failed", "message": str(e)}

# ==============================================================================
# 4. GROMACS TRAJECTORY ANALYSIS (RMSD, RMSF, Rg)
# ==============================================================================

@celery_app.task
def run_trajectory_analysis_task(job_id: str, rmsd_group: int = 4, rmsf_group: int = 1, rg_group: int = 1):
    """Run RMSD, RMSF, and Rg calculations on backmapped trajectories."""
    import matplotlib.pyplot as plt

    try:
        os.chdir("/app")
    except Exception:
        pass

    job_dir = Path(f"/app/shared_data/{job_id}")
    xtc_dir = job_dir / "all_xtc"
    analysis_dir = job_dir / "analysis" / "trajectory"
    analysis_dir.mkdir(parents=True, exist_ok=True)

    tpr_file = xtc_dir / "dynamic.tpr"
    xtc_file = xtc_dir / "all.xtc"
    ndx_file = xtc_dir / "index.ndx"

    if not (tpr_file.exists() and xtc_file.exists() and ndx_file.exists()):
        return {"status": "error", "message": f"Missing required files in {xtc_dir}."}

    try:
        subprocess.run(["gmx", "rms", "-s", str(tpr_file), "-f", str(xtc_file), "-n", str(ndx_file), "-o", str(analysis_dir / "rmsd.xvg"), "-tu", "ns"], input=f"{rmsd_group}\n{rmsd_group}\n", text=True, capture_output=True, cwd=str(analysis_dir), check=True)
        subprocess.run(["gmx", "rmsf", "-s", str(tpr_file), "-f", str(xtc_file), "-n", str(ndx_file), "-o", str(analysis_dir / "rmsf.xvg"), "-res"], input=f"{rmsf_group}\n", text=True, capture_output=True, cwd=str(analysis_dir), check=True)
        subprocess.run(["gmx", "gyrate", "-s", str(tpr_file), "-f", str(xtc_file), "-n", str(ndx_file), "-o", str(analysis_dir / "gyration.xvg")], input=f"{rg_group}\n", text=True, capture_output=True, cwd=str(analysis_dir), check=True)

        def plot_xvg(xvg_path, title, xlabel, ylabel, out_png, x_scale=1.0):
            x, y = [], []
            if Path(xvg_path).exists():
                with open(xvg_path, 'r') as f:
                    for line in f:
                        if not line.startswith(('@', '#')):
                            parts = line.split()
                            if len(parts) >= 2:
                                x.append(float(parts[0]) * x_scale)
                                y.append(float(parts[1]))
                plt.figure(figsize=(8, 5))
                plt.plot(x, y, color='blue', linewidth=1.5)
                plt.title(title, fontsize=14, fontweight='bold')
                plt.xlabel(xlabel, fontsize=12)
                plt.ylabel(ylabel, fontsize=12)
                plt.grid(True, linestyle='--', alpha=0.7)
                plt.tight_layout()
                plt.savefig(out_png, dpi=300)
                plt.close()

        plot_xvg(analysis_dir / "rmsd.xvg", "RMSD Over Time", "Time (ns)", "RMSD (nm)", analysis_dir / "RMSD.png")
        plot_xvg(analysis_dir / "rmsf.xvg", "Per-Residue RMSF", "Residue Number", "RMSF (nm)", analysis_dir / "RMSF.png")
        plot_xvg(analysis_dir / "gyration.xvg", "Radius of Gyration", "Time (ns)", "Rg (nm)", analysis_dir / "Gyration.png", x_scale=0.001)

        zip_path = analysis_dir / "trajectory_results.zip"
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for file_path in analysis_dir.glob("*"):
                if file_path.name != "trajectory_results.zip":
                    zipf.write(file_path, arcname=file_path.name)

        return {"status": "success", "message": "Trajectory Analysis completed successfully.", "result_file": str(zip_path)}
    except subprocess.CalledProcessError as e:
        err = e.stderr if e.stderr else e.stdout
        return {"status": "error", "message": f"GROMACS Error: {err}"}
    except Exception as e:
        return {"status": "error", "message": str(e)}

# ==============================================================================
# 5. TTCLUST CLUSTERING
# ==============================================================================

@celery_app.task
def run_ttclust_task(job_id: str, n_clusters: int = 5):
    """Run TTClust trajectory clustering."""
    try:
        os.chdir("/app")
    except Exception:
        pass

    job_dir = Path(f"/app/shared_data/{job_id}")
    xtc_file = job_dir / "all_xtc" / "all.xtc"
    analysis_dir = job_dir / "analysis" / "ttclust"
    analysis_dir.mkdir(parents=True, exist_ok=True)

    gro_files = list((job_dir / "gro").glob("aa_*ns.gro"))
    if not xtc_file.exists() or not gro_files:
        return {"status": "error", "message": "Missing all.xtc or .gro files."}

    try:
        subprocess.run([sys.executable, "-m", "pip", "install", "-q", "ttclust"], check=True)
        
        env = os.environ.copy()
        env["MPLBACKEND"] = "Agg"

        process = subprocess.run(
            ["ttclust", "-f", str(xtc_file), "-t", str(gro_files[0]), "-n", str(n_clusters)],
            cwd=analysis_dir, env=env, text=True, capture_output=True
        )

        zip_path = analysis_dir / "ttclust_results.zip"
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root, _, files in os.walk(analysis_dir):
                for file in files:
                    if file != "ttclust_results.zip":
                        fp = os.path.join(root, file)
                        zipf.write(fp, arcname=os.path.relpath(fp, analysis_dir))

        return {"status": "success", "message": "TTClust Complete", "result_file": str(zip_path)}
    except Exception as e:
        return {"status": "error", "message": str(e)}

# ==============================================================================
# 6. PLIP INTERACTION PROFILER
# ==============================================================================

@celery_app.task
def run_plip_analysis_task(job_id: str):
    """Run the complete PLIP analysis pipeline using local /app/tools/Analysis tools."""
    try:
        os.chdir("/app")
    except Exception:
        pass

    os.environ["MPLBACKEND"] = "Agg"
    env = os.environ.copy()
    env["MPLBACKEND"] = "Agg"

    work_dir = Path(f"/app/shared_data/{job_id}/analysis/plip")
    work_dir.mkdir(parents=True, exist_ok=True)
    backmap_dir = Path(f"/app/shared_data/{job_id}/all_xtc")

    tpr_file = backmap_dir / "dynamic.tpr"
    ndx_file = backmap_dir / "index.ndx"
    xtc_file = backmap_dir / "all.xtc"

    if not (tpr_file.exists() and xtc_file.exists()):
        return {"status": "error", "message": f"Missing required trajectory files in {backmap_dir}."}

    try:
        # Step 1: Convert trajectory to concatenated PDB
        subprocess.run(
            ["gmx", "trjconv", "-s", str(tpr_file), "-n", str(ndx_file), "-f", str(xtc_file), "-dt", "1000", "-o", "concatenated.pdb"],
            input="0\n", cwd=work_dir, check=True, text=True, capture_output=True, env=env
        )

        # Copy Analysis scripts from tools
        analysis_dir = work_dir / "Analysis"
        if analysis_dir.exists():
            shutil.rmtree(analysis_dir)
        shutil.copytree("/app/tools/Analysis", analysis_dir)
        
        shutil.move(work_dir / "concatenated.pdb", analysis_dir / "concatenated.pdb")

        # Execute analysis pipeline scripts
        subprocess.run([sys.executable, "split_pdb.py"], cwd=analysis_dir, check=True, capture_output=True, text=True, env=env)
        subprocess.run(["bash", "plip.sh"], cwd=analysis_dir, check=True, capture_output=True, text=True, env=env)
        subprocess.run([sys.executable, "convert_xml_to_json.py"], cwd=analysis_dir, check=True, capture_output=True, text=True, env=env)
        
        if (analysis_dir / "split_pdbs" / "extracted_data.json").exists():
            shutil.copy(analysis_dir / "split_pdbs" / "extracted_data.json", analysis_dir / "extracted_data.json")

        subprocess.run([sys.executable, "process_json.py"], cwd=analysis_dir, check=True, capture_output=True, text=True, env=env)

        bond_csv = analysis_dir / "bond_csv_files"
        split_bond_csv = analysis_dir / "split_pdbs" / "bond_csv_files"
        if bond_csv.exists():
            if split_bond_csv.exists():
                shutil.rmtree(split_bond_csv)
            shutil.copytree(bond_csv, split_bond_csv)

        for script in ["percentage_interactions.py", "type_interactions_color.py", "residue_interactions_color_csv_a.py", "timeline_interaction_color.py"]:
            if (analysis_dir / script).exists():
                subprocess.run([sys.executable, script], cwd=analysis_dir, check=True, capture_output=True, text=True, env=env)

        zip_path = work_dir / "plip_results.zip"
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root, _, files in os.walk(analysis_dir):
                for file in files:
                    fp = os.path.join(root, file)
                    zipf.write(fp, arcname=os.path.relpath(fp, analysis_dir))

        return {"status": "success", "result_file": str(zip_path)}
    except subprocess.CalledProcessError as e:
        err = e.stderr if e.stderr else e.stdout
        cmd_str = " ".join(e.cmd) if isinstance(e.cmd, list) else str(e.cmd)
        return {"status": "error", "message": f"PLIP Pipeline failed at command: {cmd_str}\nError Details:\n{err}"}
    except Exception as e:
        return {"status": "error", "message": str(e)}

# ==============================================================================
# 7. MM-GBSA BINDING FREE ENERGY
# ==============================================================================

@celery_app.task
def run_mmpbsa_task(job_id: str, protein_idx: int = 1, ligand_idx: int = 12, start_frame: int = 1, end_frame: int = 0):
    """Runs gmx_MMPBSA calculation."""
    try:
        os.chdir("/app")
    except Exception:
        pass

    extra_paths = ["/root/micromamba/bin", "/usr/local/bin", os.path.expanduser("~/.local/bin"), "/opt/conda/bin"]
    os.environ["PATH"] = ":".join(extra_paths) + ":" + os.environ.get("PATH", "")

    mmpbsa_exec = shutil.which("gmx_MMPBSA") or "/root/micromamba/bin/gmx_MMPBSA"
    if not os.path.exists(mmpbsa_exec):
        return {"status": "error", "message": "gmx_MMPBSA executable not found."}

    job_dir = Path(f"/app/shared_data/{job_id}")
    mmgbsa_dir = job_dir / "analysis" / "MMGBSA"
    mmgbsa_dir.mkdir(parents=True, exist_ok=True)
    xtc_dir = job_dir / "all_xtc"

    try:
        # Copy input files
        for src, dst in [
            (xtc_dir / "all.xtc", mmgbsa_dir / "all.xtc"),
            (xtc_dir / "dynamic.tpr", mmgbsa_dir / "dynamic.tpr"),
            (xtc_dir / "index.ndx", mmgbsa_dir / "index.ndx"),
            (job_dir / "Backmapping" / "aa.gro", mmgbsa_dir / "aa.gro"),
            (job_dir / "Backmapping" / "topol.top", mmgbsa_dir / "topol.top")
        ]:
            if src.exists():
                shutil.copy2(src, dst)

        # Copy mmpbsa.in script from /app/tools/Analysis or download fallback
        if (Path("/app/tools/Analysis/mmpbsa.in")).exists():
            shutil.copy("/app/tools/Analysis/mmpbsa.in", mmgbsa_dir / "mmpbsa.in")

        mmpbsa_cmd = [
            mmpbsa_exec, "-O", "-i", "mmpbsa.in", "-cs", "dynamic.tpr",
            "-ct", "all.xtc", "-ci", "index.ndx", "-cg", str(protein_idx), str(ligand_idx),
            "-cp", "topol.top", "-o", "FINAL_RESULTS_MMPBSA.dat",
            "-eo", "FINAL_RESULTS_MMPBSA.csv", "-do", "FINAL_DECOMP_MMPBSA.dat",
            "-deo", "FINAL_DECOMP_MMPBSA.csv", "-nogui"
        ]

        subprocess.run(mmpbsa_cmd, cwd=mmgbsa_dir, check=True, capture_output=True, text=True, env=os.environ)

        summary_data = ""
        res_dat = mmgbsa_dir / "FINAL_RESULTS_MMPBSA.dat"
        if res_dat.exists():
            with open(res_dat, "r") as f:
                summary_data = "".join(f.readlines()[-20:])

        zip_path = job_dir / "analysis" / "mmgbsa_results.zip"
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root, _, files in os.walk(mmgbsa_dir):
                for file in files:
                    fp = os.path.join(root, file)
                    zipf.write(fp, arcname=os.path.relpath(fp, mmgbsa_dir))

        return {"status": "success", "result_file": str(zip_path), "summary_data": summary_data}
    except subprocess.CalledProcessError as e:
        err = e.stderr if e.stderr else e.stdout
        return {"status": "error", "message": f"MM-GBSA Error: {err}"}
    except Exception as e:
        return {"status": "error", "message": str(e)}

# ==============================================================================
# 8. PCA & FREE ENERGY LANDSCAPE (FEL) (NO PYMOL)
# ==============================================================================

@celery_app.task
def run_pca_fel_task(job_id: str, backbone_group: int = 4):
    """Executes Covariance, Eigenvector projections, and Free Energy Landscape (FEL)."""
    try:
        os.chdir("/app")
    except Exception:
        pass

    os.environ["MPLBACKEND"] = "Agg"
    env = os.environ.copy()
    env["MPLBACKEND"] = "Agg"

    job_dir = Path(f"/app/shared_data/{job_id}")
    xtc_dir = job_dir / "all_xtc"
    pca_dir = job_dir / "analysis" / "PCA-FEL"
    pca_dir.mkdir(parents=True, exist_ok=True)

    tpr_file = xtc_dir / "dynamic.tpr"
    xtc_file = xtc_dir / "all.xtc"
    ndx_file = xtc_dir / "index.ndx"

    if not (tpr_file.exists() and xtc_file.exists() and ndx_file.exists()):
        return {"status": "error", "message": "Missing required trajectory files."}

    try:
        group_input = f"{backbone_group}\n{backbone_group}\n"

        # 1. Covariance Analysis
        covar_cmd = ["gmx", "covar", "-f", str(xtc_file), "-s", str(tpr_file), "-n", str(ndx_file), "-o", "eigenval.xvg", "-v", "eigenvec.trr", "-av", "average.pdb", "-l", "covar.log", "-b", "0", "-tu", "ns"]
        subprocess.run(covar_cmd, input=group_input, text=True, capture_output=True, cwd=pca_dir, check=True)

        # 2. Eigenvalue Projections
        projections = [("1", "2"), ("1", "3"), ("2", "3")]
        for first, last in projections:
            anaeig_cmd = ["gmx", "anaeig", "-v", "eigenvec.trr", "-f", str(xtc_file), "-s", str(tpr_file), "-n", str(ndx_file), "-comp", f"eigcomp{first}-{last}.xvg", "-rmsf", f"eigrmsf{first}-{last}.xvg", "-2d", f"2dproj{first}-{last}.xvg", "-b", "0", "-tu", "ns", "-first", first, "-last", last]
            subprocess.run(anaeig_cmd, input=group_input, text=True, capture_output=True, cwd=pca_dir, check=True)

        # 3. FEL Calculation
        subprocess.run(["gmx", "sham", "-f", "2dproj1-2.xvg", "-notime", "-ls", "gibb1-2.xpm"], cwd=pca_dir, check=True, capture_output=True)

        # Copy Analysis scripts from /app/tools/Analysis
        for script in ["xpm2txt.py", "PCA.py", "FEL.py"]:
            if (Path(f"/app/tools/Analysis/{script}")).exists():
                shutil.copy(f"/app/tools/Analysis/{script}", pca_dir / script)

        subprocess.run(["python2", "xpm2txt.py", "-f", "gibb1-2.xpm", "-o", "FEL.dat"], cwd=pca_dir, check=True, capture_output=True)
        subprocess.run([sys.executable, "PCA.py"], cwd=pca_dir, check=True, capture_output=True, env=env)
        subprocess.run([sys.executable, "FEL.py"], cwd=pca_dir, check=True, capture_output=True, env=env)

        zip_path = pca_dir / "pca_fel_results.zip"
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root, _, files in os.walk(pca_dir):
                for file in files:
                    if file != "pca_fel_results.zip":
                        fp = os.path.join(root, file)
                        zipf.write(fp, arcname=file)

        return {"status": "success", "message": "PCA & FEL Complete", "result_file": str(zip_path)}
    except subprocess.CalledProcessError as e:
        err = e.stderr if e.stderr else e.stdout
        return {"status": "error", "message": f"PCA/FEL Error: {err}"}
    except Exception as e:
        return {"status": "error", "message": str(e)}
