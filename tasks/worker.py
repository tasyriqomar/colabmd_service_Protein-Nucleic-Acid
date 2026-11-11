from celery import Celery
import subprocess
from pathlib import Path
import os
import io
import zipfile
import shutil
import glob
import sys
import gzip
import numpy as np
import pandas as pd

REDIS_URL = os.getenv("REDIS_URL", "redis://127.0.0.1:6379/0")

celery_app = Celery(
    "colabmd_tasks",
    broker=REDIS_URL,
    backend=REDIS_URL
)

# ==============================================================================
# TAB 1: APTAMER PREP
# ==============================================================================
@celery_app.task
def run_aptamer_task(work_dir: str, aptamer_pdb_content: bytes, convert_rna: bool):
    wd = Path(work_dir)
    wd.mkdir(parents=True, exist_ok=True)
    try:
        # Guarantee both filenames exist in the task working directory
        with open(wd / "Aptamer.pdb", "wb") as f:
            f.write(aptamer_pdb_content)
        with open(wd / "aptamer.pdb", "wb") as f:
            f.write(aptamer_pdb_content)

        shutil.copy("/app/tools/Aptamer-Processing/rna_to_dna_with_methyl.py", wd / "rna_to_dna_with_methyl.py")
        shutil.copy("/app/tools/Aptamer-Processing/rename_dna.py", wd / "rename_dna.py")

        if convert_rna:
            # Converts Aptamer.pdb -> aptamer_DNA.pdb
            subprocess.run([sys.executable, "rna_to_dna_with_methyl.py"], cwd=str(wd), check=True, capture_output=True, text=True)
        else:
            # If already DNA, create aptamer_DNA.pdb for rename_dna.py
            shutil.copy(wd / "Aptamer.pdb", wd / "aptamer_DNA.pdb")

        # Converts aptamer_DNA.pdb -> aptamer_fixed.pdb
        subprocess.run([sys.executable, "rename_dna.py"], cwd=str(wd), check=True, capture_output=True, text=True)

        fixed_pdb = wd / "aptamer_fixed.pdb"
        if not fixed_pdb.exists():
            return {"status": "failed", "message": "Failed to generate aptamer_fixed.pdb"}

        return {"status": "success", "result_file": str(fixed_pdb)}
    except subprocess.CalledProcessError as e:
        err = e.stderr if e.stderr else e.stdout
        return {"status": "failed", "message": f"Aptamer Processing Error: {err}"}
    except Exception as e:
        return {"status": "failed", "message": str(e)}

run_aptamer_processing_task = run_aptamer_task

# ==============================================================================
# TAB 3: HADDOCK3 DOCKING & CONTACT ANALYSIS (NON-STANDARD BASE COMPATIBLE)
# ==============================================================================
import numpy as np


def prepare_canonical_aptamer_for_haddock(
    input_pdb_path: Path, output_pdb_path: Path
):
    """Converts modified nucleobases (CFZ, UFT, UMS, ADS) into canonical RNA

    (C, U, U, A) and suppresses non-standard atoms so HADDOCK3 [topoaa] generates
    CNS topologies without stripping residues or crashing.
    """
    res_map = {
        'CFZ': 'C',
        'UFT': 'U',
        'UMS': 'U',
        'ADS': 'A',
    }

    # Atom types not recognized by CNS standard dna-rna.top
    skip_atoms = {"F2'", "SE2'", "CA'", "HA'1", "HA'2", "HA'3"}
    atom_rename = {'S2P': 'O1P', 'S3P': 'O2P'}

    cleaned_lines = []
    with open(input_pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()

                if atom_name in skip_atoms:
                    continue

                new_res = res_map.get(res_name, res_name)
                new_atom = atom_rename.get(atom_name, atom_name)

                # Format atom name adhering strictly to PDB column specification
                if len(new_atom) == 4 or (
                    len(new_atom) > 0 and new_atom[0].isdigit()
                ):
                    formatted_atom = f'{new_atom:<4}'
                else:
                    formatted_atom = f' {new_atom:<3}'

                formatted_res = f'{new_res:>3}'
                # Enforce chain B and canonical RNA representation
                new_line = f'ATOM  {line[6:12]}{formatted_atom} {formatted_res} B{line[22:]}'
                cleaned_lines.append(new_line)
            elif line.startswith(('TER', 'END')):
                cleaned_lines.append(line)

    with open(output_pdb_path, 'w') as f:
        f.writelines(cleaned_lines)


def restore_modified_aptamer_pose(
    docked_pose_pdb: Path,
    original_modified_pdb: Path,
    output_complex_pdb: Path,
):
    """Calculates the rigid-body transformation (R, t) between the docked canonical pose

    and the original relaxed aptamer, then restores all original non-standard nucleobases
    (CFZ, UFT, UMS, ADS, F2', SE2', etc.) into the final docked complex.
    """
    fit_atoms = ["C1'", "C3'", "C4'", "O4'", "C5'"]

    # 1. Extract docked backbone atoms (Chain B) and docked protein lines (Chain A)
    docked_coords = {}
    chain_a_lines = []
    with open(docked_pose_pdb, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                chain = line[21]
                if chain == 'A':
                    chain_a_lines.append(line)
                elif chain == 'B':
                    res_num = int(line[22:26])
                    atom_name = line[12:16].strip()
                    if atom_name in fit_atoms:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        docked_coords[(res_num, atom_name)] = np.array(
                            [x, y, z]
                        )

    # 2. Extract corresponding backbone atoms from original relaxed aptamer
    orig_coords = {}
    orig_all_lines = []
    with open(original_modified_pdb, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                orig_all_lines.append(line)
                res_num = int(line[22:26])
                atom_name = line[12:16].strip()
                if atom_name in fit_atoms:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    orig_coords[(res_num, atom_name)] = np.array([x, y, z])

    # 3. Compute optimal Kabsch rotation and translation
    common_keys = sorted(
        list(set(docked_coords.keys()) & set(orig_coords.keys()))
    )
    if len(common_keys) < 3:
        raise ValueError(
            f'Insufficient matching backbone atoms ({len(common_keys)}) to align modified aptamer.'
        )

    p_docked = np.array([docked_coords[k] for k in common_keys])
    p_orig = np.array([orig_coords[k] for k in common_keys])

    c_docked = np.mean(p_docked, axis=0)
    c_orig = np.mean(p_orig, axis=0)

    p_docked_c = p_docked - c_docked
    p_orig_c = p_orig - c_orig

    H = p_orig_c.T @ p_docked_c
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    t = c_docked - c_orig @ R.T

    # 4. Transform all atoms of the modified aptamer to the docked orientation
    chain_b_restored_lines = []
    for line in orig_all_lines:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        new_xyz = np.array([x, y, z]) @ R.T + t
        # Ensure Chain ID is explicitly set to 'B'
        new_line = f'{line[:21]}B{line[22:30]}{new_xyz[0]:8.3f}{new_xyz[1]:8.3f}{new_xyz[2]:8.3f}{line[54:]}'
        chain_b_restored_lines.append(new_line)

    # 5. Assemble final complex: Chain A (Protein) + Chain B (Full Modified Aptamer)
    with open(output_complex_pdb, 'w') as f:
        f.writelines(chain_a_lines)
        if chain_a_lines and not chain_a_lines[-1].startswith('TER'):
            f.write('TER\n')
        f.writelines(chain_b_restored_lines)
        f.write('TER\nEND\n')


@celery_app.task
def run_haddock_docking_task(
    work_dir: str,
    protein_pdb: bytes,
    dna_pdb: bytes,
    sampling: int = 20,
    cmrest: bool = True,
):
    wd = Path(work_dir)
    try:
        raw_rbd = wd / 'RBD.pdb'
        raw_dna = wd / 'DNA.pdb'
        with open(raw_rbd, 'wb') as f:
            f.write(protein_pdb)
        with open(raw_dna, 'wb') as f:
            f.write(dna_pdb)

        # 1. Format Chain IDs
        subprocess.run(
            'pdb_chain -A RBD.pdb > RBD_fixed.pdb',
            shell=True,
            cwd=wd,
            check=True,
        )
        subprocess.run(
            'pdb_chain -B DNA.pdb > DNA_modified.pdb',
            shell=True,
            cwd=wd,
            check=True,
        )

        # 2. Prepare Canonical Scaffold for HADDOCK3 CNS topology
        prepare_canonical_aptamer_for_haddock(
            wd / 'DNA_modified.pdb', wd / 'DNA_fixed.pdb'
        )

        # 3. HADDOCK3 Configuration
        docking_cfg = f"""
run_dir = "rigid_body_docking"
mode = "local"
molecules = ["RBD_fixed.pdb", "DNA_fixed.pdb"]

[topoaa]

[rigidbody]
sampling = {sampling}
cmrest = {"true" if cmrest else "false"}
"""
        with open(wd / 'docking.cfg', 'w') as f:
            f.write(docking_cfg)

        subprocess.run(
            ['haddock3', 'docking.cfg'],
            cwd=wd,
            check=True,
            capture_output=True,
            text=True,
        )

        # 4. Parse Docking Scores
        out_folder = wd / 'rigid_body_docking' / '1_rigidbody'
        docked_poses = sorted(glob.glob(str(out_folder / '*.pdb.gz')))
        if not docked_poses:
            return {
                'status': 'failed',
                'message': 'No docked poses generated by HADDOCK3.',
            }

        pose_data = []
        for p_path in docked_poses:
            p_name = Path(p_path).name.replace('.pdb.gz', '')
            metrics = {
                'Pose': p_name,
                'HADDOCK Score': None,
                'E_vdw (kcal/mol)': None,
                'E_elec (kcal/mol)': None,
                'E_desolv (kcal/mol)': None,
                'E_air (kcal/mol)': None,
                'BSA (Å²)': None,
            }
            with gzip.open(p_path, 'rt') as f:
                for line in f:
                    if line.startswith('REMARK HADDOCK score:'):
                        metrics['HADDOCK Score'] = float(
                            line.split(':')[1].strip()
                        )
                    elif line.startswith('REMARK vdw energy:'):
                        metrics['E_vdw (kcal/mol)'] = float(
                            line.split(':')[1].strip()
                        )
                    elif line.startswith('REMARK elec energy:'):
                        metrics['E_elec (kcal/mol)'] = float(
                            line.split(':')[1].strip()
                        )
                    elif line.startswith('REMARK desolv energy:'):
                        metrics['E_desolv (kcal/mol)'] = float(
                            line.split(':')[1].strip()
                        )
                    elif line.startswith('REMARK air energy:'):
                        metrics['E_air (kcal/mol)'] = float(
                            line.split(':')[1].strip()
                        )
                    elif line.startswith('REMARK buried surface area:'):
                        metrics['BSA (Å²)'] = float(line.split(':')[1].strip())
                    elif line.startswith('ATOM'):
                        break
            pose_data.append(metrics)

        df_scores = pd.DataFrame(pose_data)
        if df_scores['HADDOCK Score'].notnull().any():
            df_scores = df_scores.sort_values(
                by='HADDOCK Score', ascending=True
            ).reset_index(drop=True)
        df_scores.to_csv(wd / 'haddock3_docking_scores.csv', index=False)

        # 5. Decompress Best Pose and Restore Full Non-Standard Nucleobases
        best_pose_gz = docked_poses[0]
        temp_canonical_pdb = wd / 'docked_canonical_temp.pdb'
        with (
            gzip.open(best_pose_gz, 'rt') as f_in,
            open(temp_canonical_pdb, 'w') as f_out,
        ):
            f_out.write(f_in.read())

        best_pose_pdb = wd / 'docked_complex_best.pdb'
        restore_modified_aptamer_pose(
            docked_pose_pdb=temp_canonical_pdb,
            original_modified_pdb=wd / 'DNA_modified.pdb',
            output_complex_pdb=best_pose_pdb,
        )

        # 6. Extract Contacts & H-Bonds (Preserving F, SE, S for 5DO4 interface)
        from Bio.PDB import PDBParser

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('complex', str(best_pose_pdb))
        chain_a, chain_b = structure[0]['A'], structure[0]['B']
        contacts, hbonds = [], []
        # Added SE and F for halogen/selenium coordination with Thrombin Exosite I
        heteroatoms = {'N', 'O', 'S', 'F', 'P', 'SE'}

        for atom_a in chain_a.get_atoms():
            if atom_a.element == 'H':
                continue
            for atom_b in chain_b.get_atoms():
                if atom_b.element == 'H':
                    continue
                dist = atom_a - atom_b
                if dist <= 4.0:
                    r_a = f'{atom_a.get_parent().get_resname()}{atom_a.get_parent().get_id()[1]}'
                    r_b = f'{atom_b.get_parent().get_resname()}{atom_b.get_parent().get_id()[1]}'
                    contacts.append({
                        'Protein Residue': r_a,
                        'DNA Residue': r_b,
                        'Distance (Å)': round(float(dist), 2),
                    })
                    if (
                        dist <= 3.5
                        and atom_a.element in heteroatoms
                        and atom_b.element in heteroatoms
                    ):
                        hbonds.append({
                            'Protein Residue': r_a,
                            'Protein Atom': atom_a.get_name(),
                            'DNA Residue': r_b,
                            'DNA Atom': atom_b.get_name(),
                            'Distance (Å)': round(float(dist), 2),
                        })

        pd.DataFrame(contacts).drop_duplicates().to_csv(
            wd / 'interfacial_contacts.csv', index=False
        )
        pd.DataFrame(hbonds).drop_duplicates().to_csv(
            wd / 'hydrogen_bonds.csv', index=False
        )

        # 7. Zip Output
        zip_path = wd / 'haddock3_docking_results.zip'
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as z:
            for fname in [
                'docked_complex_best.pdb',
                'haddock3_docking_scores.csv',
                'interfacial_contacts.csv',
                'hydrogen_bonds.csv',
            ]:
                if (wd / fname).exists():
                    z.write(wd / fname, arcname=fname)

        return {
            'status': 'success',
            'result_file': str(zip_path),
            'scores': df_scores.to_dict(orient='records'),
        }
    except Exception as e:
        return {
            'status': 'failed',
            'message': f'HADDOCK3 Docking Error: {str(e)}',
        }

# ==============================================================================
# TAB 4: DUAL DNA/RNA COMPLEX COARSE-GRAINED PREPARATION
# ==============================================================================
@celery_app.task
def run_complex_cg_prep_task(work_dir: str, complex_pdb_content: bytes, dnatype: str = "ss", elastic: bool = True):
    wd = Path(work_dir)
    wd.mkdir(parents=True, exist_ok=True)
    try:
        cg_tools_dir = Path("/app/tools/CGMD")
        complex_pdb = wd / "docked_complex_best.pdb"
        with open(complex_pdb, "wb") as f:
            f.write(complex_pdb_content)

        # ----------------------------------------------------------------------
        # 1. Split Complex into Protein (Chain A) & Nucleic Acid (Chain B)
        # ----------------------------------------------------------------------
        prot_aa_pdb = wd / "protein_aa.pdb"
        na_aa_pdb = wd / "na_aa.pdb"

        prot_lines, na_raw_lines = [], []
        with open(complex_pdb, "r") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    chain = line[21].strip()
                    if chain == "A":
                        prot_lines.append(line)
                    elif chain == "B":
                        na_raw_lines.append(line)

        # Fallback if chains aren't explicitly tagged
        if not prot_lines or not na_raw_lines:
            all_na_res = {
                "DA", "DT", "DG", "DC", "A", "T", "G", "C", "U", "URA", "ADE", "CYT", "GUA",
                "CFZ", "UMS", "UFT", "ADS"
            }
            prot_lines, na_raw_lines = [], []
            with open(complex_pdb, "r") as f:
                for line in f:
                    if line.startswith(("ATOM", "HETATM")):
                        r_name = line[17:20].strip()
                        if r_name in all_na_res:
                            na_raw_lines.append(line[:21] + "B" + line[22:])
                        else:
                            prot_lines.append(line[:21] + "A" + line[22:])

        with open(prot_aa_pdb, "w") as f:
            f.writelines(prot_lines)
            f.write("END\n")

        # ----------------------------------------------------------------------
        # 2. Nucleic Acid Type Detection (DNA vs RNA)
        # ----------------------------------------------------------------------
        # Check for RNA identifiers: Uracil variants, ribose O2', or modified bases (CFZ, UMS, UFT, ADS)
        is_rna = False
        for line in na_raw_lines:
            r_name = line[17:20].strip()
            a_name = line[12:16].strip()
            if r_name in {"U", "URA", "RU", "CFZ", "UMS", "UFT", "ADS"} or a_name in {"O2'", "HO2'", "F2'", "SE2'"}:
                is_rna = True
                break

        # ----------------------------------------------------------------------
        # 3. Sanitize Nucleic Acid for the Appropriate Tool
        # ----------------------------------------------------------------------
        # Alias atom names so centers of mass are computed accurately
        atom_alias = {
            "S2P": "OP1", "S3P": "OP2", "O1P": "OP1", "O2P": "OP2"
        }

        # Residue mappings
        if is_rna:
            # RNA targets canonical A, C, G, U
            res_map = {
                "CFZ": "C", "CYT": "C", "DC": "C", "C": "C",
                "UFT": "U", "UMS": "U", "URA": "U", "DT": "U", "T": "U", "U": "U",
                "ADS": "A", "ADE": "A", "DA": "A", "A": "A",
                "GUA": "G", "DG": "G", "G": "G"
            }
        else:
            # DNA targets canonical DA, DC, DG, DT
            res_map = {
                "CFZ": "DC", "CYT": "DC", "C": "DC", "DC": "DC",
                "UFT": "DT", "UMS": "DT", "URA": "DT", "U": "DT", "DT": "DT", "T": "DT",
                "ADS": "DA", "ADE": "DA", "A": "DA", "DA": "DA",
                "GUA": "DG", "G": "DG", "DG": "DG"
            }

        sanitized_na_lines = []
        for line in na_raw_lines:
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()

            clean_atom = atom_alias.get(atom_name, atom_name)
            clean_res = res_map.get(res_name, "U" if is_rna else "DT")

            fmt_atom = f"{clean_atom:<4}" if len(clean_atom) == 4 or clean_atom[0].isdigit() else f" {clean_atom:<3}"
            fmt_res = f"{clean_res:>3}"
            sanitized_na_lines.append(f"ATOM  {line[6:12]}{fmt_atom} {fmt_res} B{line[22:]}")

        with open(na_aa_pdb, "w") as f:
            f.writelines(sanitized_na_lines)
            f.write("END\n")

        # ----------------------------------------------------------------------
        # 4. Coarse-Grain Nucleic Acid (DNA vs RNA with Safe Fallback)
        # ----------------------------------------------------------------------
        subprocess.run(["gmx", "editconf", "-f", "na_aa.pdb", "-o", "na_aa.gro"], cwd=str(wd), check=True, capture_output=True)

        if is_rna:
            # Check for RNA scripts
            rna_candidates = [
                cg_tools_dir / "martini-nucleotide.py",
                cg_tools_dir / "martinize-nucleotide.py"
            ]
            rna_script = next((s for s in rna_candidates if s.exists()), None)

            if rna_script:
                rna_type = "ss" if dnatype == "ss" else "ds-stiff"
                na_cmd = [
                    "python2", str(rna_script),
                    "-type", rna_type,
                    "-f", "na_aa.gro",
                    "-o", "dna.top",
                    "-x", "dna_CG.pdb"
                ]
            else:
                # Fallback to martinize-dna.py / martini-dna.py if nucleotide script is missing
                dna_candidates = [
                    cg_tools_dir / "martini-dna.py",
                    cg_tools_dir / "martinize-dna.py"
                ]
                fallback_script = next((s for s in dna_candidates if s.exists()), None)
                if not fallback_script:
                    raise FileNotFoundError("No coarse-graining script found in /app/tools/CGMD")

                na_cmd = [
                    "python2", str(fallback_script),
                    "-dnatype", dnatype,
                    "-f", "na_aa.gro",
                    "-o", "dna.top",
                    "-x", "dna_CG.pdb"
                ]
        else:
            # DNA branch
            dna_candidates = [
                cg_tools_dir / "martini-dna.py",
                cg_tools_dir / "martinize-dna.py"
            ]
            dna_script = next((s for s in dna_candidates if s.exists()), None)
            if not dna_script:
                raise FileNotFoundError("Neither martini-dna.py nor martinize-dna.py was found in /app/tools/CGMD")

            na_cmd = [
                "python2", str(dna_script),
                "-dnatype", dnatype,
                "-f", "na_aa.gro",
                "-o", "dna.top",
                "-x", "dna_CG.pdb"
            ]

        subprocess.run(na_cmd, cwd=str(wd), check=True, capture_output=True, text=True)
        subprocess.run(["gmx", "editconf", "-f", "dna_CG.pdb", "-o", "aptamer_CG.gro"], cwd=str(wd), check=True, capture_output=True)

        # ----------------------------------------------------------------------
        # 5. Coarse-Grain Protein (martinize.py)
        # ----------------------------------------------------------------------
        prot_cmd = [
            "python2", str(cg_tools_dir / "martinize.py"),
            "-f", "protein_aa.pdb",
            "-o", "system_prot_temp.top",
            "-x", "protein_CG.pdb",
            "-p", "backbone",
            "-ff", "martini22"
        ]
        if elastic:
            prot_cmd.append("-elastic")
        subprocess.run(prot_cmd, cwd=str(wd), check=True, capture_output=True, text=True)
        subprocess.run(["gmx", "editconf", "-f", "protein_CG.pdb", "-o", "protein_CG.gro"], cwd=str(wd), check=True, capture_output=True)

        # ----------------------------------------------------------------------
        # 6. Generate Merged complex_CG.gro
        # ----------------------------------------------------------------------
        with open(wd / "protein_CG.gro", "r") as f:
            p_gro_lines = f.readlines()
        with open(wd / "aptamer_CG.gro", "r") as f:
            d_gro_lines = f.readlines()

        p_atoms = p_gro_lines[2:-1]
        d_atoms = d_gro_lines[2:-1]
        total_atoms = len(p_atoms) + len(d_atoms)
        box_vector = p_gro_lines[-1] if len(p_gro_lines) > 2 else "  10.00000  10.00000  10.00000\n"

        mol_title = "RNA Aptamer" if is_rna else "DNA Aptamer"
        merged_gro = [f"Complex Coarse-Grained Structure ({mol_title})\n", f"{total_atoms:>5}\n"]
        atom_count = 1
        for line in p_atoms:
            merged_gro.append(f"{line[:15]}{atom_count:>5}{line[20:]}")
            atom_count += 1
        for line in d_atoms:
            merged_gro.append(f"{line[:15]}{atom_count:>5}{line[20:]}")
            atom_count += 1
        merged_gro.append(box_vector)

        with open(wd / "complex_CG.gro", "w") as f:
            f.writelines(merged_gro)

        # ----------------------------------------------------------------------
        # 7. Box Solvation & Ionization (insane.py)
        # ----------------------------------------------------------------------
        if (cg_tools_dir / "insane.py").exists():
            shutil.copy(cg_tools_dir / "insane.py", wd / "insane.py")

        insane_cmd = [
            "python2", "insane.py",
            "-f", "complex_CG.gro",
            "-o", "system.gro",
            "-pbc", "cubic",
            "-box", "10,10,10",
            "-salt", "0.15",
            "-charge", "auto",
            "-sol", "W"
        ]
        insane_proc = subprocess.run(insane_cmd, cwd=str(wd), check=True, capture_output=True, text=True)

        solvent_molecules = []
        for line in insane_proc.stderr.splitlines():
            line_str = line.strip()
            if line_str and not line_str.startswith(";"):
                parts = line_str.split()
                if len(parts) >= 2 and parts[1].isdigit():
                    mol_name_clean = parts[0].replace("NA+", "NA").replace("CL-", "CL")
                    solvent_molecules.append(f"{mol_name_clean:<10} {parts[1]:>7}\n")

        if not solvent_molecules and (wd / "system.gro").exists():
            w_cnt, na_cnt, cl_cnt = 0, 0, 0
            with open(wd / "system.gro", "r") as f:
                for line in f.readlines()[2:-1]:
                    r = line[5:10].strip()
                    if r == "W": w_cnt += 1
                    elif "NA" in r: na_cnt += 1
                    elif "CL" in r: cl_cnt += 1
            if w_cnt > 0: solvent_molecules.append(f"{'W':<10} {w_cnt:>7}\n")
            if na_cnt > 0: solvent_molecules.append(f"{'NA':<10} {na_cnt:>7}\n")
            if cl_cnt > 0: solvent_molecules.append(f"{'CL':<10} {cl_cnt:>7}\n")

        # ----------------------------------------------------------------------
        # 8. Unified system.top File
        # ----------------------------------------------------------------------
        itp_to_copy = ["martini_v2.0_ions.itp", "martini_v2.2.itp"]
        if is_rna:
            itp_to_copy.append("martini_v2.1-dna.itp")
            na_itp_include = '#include "martini_v2.1-dna.itp"\n'
        else:
            itp_to_copy.append("martini_v2.1P-dna.itp")
            na_itp_include = '#include "martini_v2.1P-dna.itp"\n'

        for itp in itp_to_copy:
            if (cg_tools_dir / itp).exists():
                shutil.copy(cg_tools_dir / itp, wd / itp)

        system_top_lines = [
            na_itp_include,
            '#include "martini_v2.0_ions.itp"\n\n',
            '#define RUBBER_BANDS\n\n',
            '#include "Protein_A.itp"\n',
            '#include "Nucleic_A.itp"\n\n',
            '[ system ]\n',
            '; name\n',
            f'Complex System ({mol_title})\n\n',
            '[ molecules ]\n',
            '; name        number\n',
            'Protein_A 	 1\n',
            'Nucleic_A        1\n'
        ]
        system_top_lines.extend(solvent_molecules)

        with open(wd / "system.top", "w") as f:
            f.writelines(system_top_lines)

        # ----------------------------------------------------------------------
        # 9. Backward Mapping Definition (.map)
        # ----------------------------------------------------------------------
        if is_rna:
            map_content = """[ molecule ]
RNA_Aptamer

[ martini ]
BB1 BB2 BB3 SC1 SC2 SC3 SC4

[ mapping ]
charmm36

[ residues ]
G      -> BB1 BB2 BB3 SC1 SC2 SC3 SC4
A      -> BB1 BB2 BB3 SC1 SC2 SC3 SC4
C      -> BB1 BB2 BB3 SC1 SC2 SC3
U      -> BB1 BB2 BB3 SC1 SC2 SC3
DG     -> BB1 BB2 BB3 SC1 SC2 SC3 SC4
DA     -> BB1 BB2 BB3 SC1 SC2 SC3 SC4
CFZ    -> BB1 BB2 BB3 SC1 SC2 SC3
UMS    -> BB1 BB2 BB3 SC1 SC2 SC3
UFT    -> BB1 BB2 BB3 SC1 SC2 SC3
ADS    -> BB1 BB2 BB3 SC1 SC2 SC3 SC4
"""
        else:
            map_content = """[ molecule ]
DNA_Aptamer

[ martini ]
BB1 BB2 BB3 SC1 SC2 SC3 SC4

[ mapping ]
charmm36

[ residues ]
DG     -> BB1 BB2 BB3 SC1 SC2 SC3 SC4
DA     -> BB1 BB2 BB3 SC1 SC2 SC3 SC4
DT     -> BB1 BB2 BB3 SC1 SC2 SC3
DC     -> BB1 BB2 BB3 SC1 SC2 SC3
"""
        with open(wd / "aptamer.map", "w") as f:
            f.write(map_content)

        # ----------------------------------------------------------------------
        # 10. Package Output Archive
        # ----------------------------------------------------------------------
        zip_path = wd / "complex_cg_prep.zip"
        files_to_zip = [
            "system.gro", "system.top", "complex_CG.gro", "protein_CG.gro",
            "aptamer_CG.gro", "protein_CG.pdb", "dna_CG.pdb", "dna.top",
            "Protein_A.itp", "Nucleic_A.itp", "martini_v2.0_ions.itp", "aptamer.map"
        ]
        if is_rna:
            files_to_zip.append("martini_v2.1-dna.itp")
        else:
            files_to_zip.append("martini_v2.1P-dna.itp")

        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as z:
            for fname in files_to_zip:
                if (wd / fname).exists():
                    z.write(wd / fname, arcname=fname)
            for itp in wd.glob("*.itp"):
                if itp.name not in files_to_zip:
                    z.write(itp, arcname=itp.name)

        gro_preview = ""
        if (wd / "system.gro").exists():
            with open(wd / "system.gro", "r") as f:
                gro_preview = "".join(f.readlines()[:40])

        return {
            "status": "success",
            "result_file": str(zip_path),
            "sys_top_preview": "".join(system_top_lines),
            "gro_preview": gro_preview
        }
    except subprocess.CalledProcessError as e:
        err = e.stderr if e.stderr else e.stdout
        return {"status": "failed", "message": f"CG Preparation Error: {err}"}
    except Exception as e:
        return {"status": "failed", "message": str(e)}

# ==============================================================================
# TAB 6: BACKMAPPING WORKFLOW (PROTEIN - NUCLEIC ACID)
# ==============================================================================

INITRAM_SCRIPT = r"""#!/bin/bash

PROGRAM=initram.sh
VERSION=0.7-purepython
AUTHOR="Tsjerk A. Wassenaar, PhD (Updated for GMX 2019+)"

DEPENDENCIES=(backward.py gmx python3)
SDIR=$( [[ $0 != ${0%/*} ]] && cd ${0%/*}; pwd )

INP=
TOP=
OUT=backmapped.gro
OTP=backmapped.top
NDX=backmapped.ndx
RAW=projected.gro
BW=0-backward.gro
CG=martini
AA=charmm36
NP=0
KICK=0.5
KEEP=false
TRJ=false
POSRE=true

# Read steps from environment variables (Python memory feature)
EMSTEPS=${EM_STEPS:-3000}
NBSTEPS=${NB_STEPS:-1000}
MDSTEPS=${MD_STEPS:-500}
DT=0.0002,0.0005,0.001,0.002
MDP=()

while [ -n "$1" ]; do
  case $1 in
        -f)    INP=$2   ; shift 2; continue ;;
        -p)    TOP=$2   ; shift 2; continue ;;
       -po)    OTP=$2   ; shift 2; continue ;;
        -o)    OUT=$2   ; shift 2; continue ;;
     -from)    CG=$2    ; shift 2; continue ;;
       -to)    AA=$2    ; shift 2; continue ;;
     -kick)    KICK=$2  ; shift 2; continue ;;
     -keep)    KEEP=true; shift  ; continue ;;
      -nopr)   POSRE=false; shift; continue ;;
          *)   shift ;;
  esac
done

[[ -z $INP ]] && echo "FATAL ERROR: MISSING INPUT STRUCTURE FILE" && exit 1
[[ -z $TOP ]] && echo "FATAL ERROR: MISSING INPUT TOPOLOGY FILE" && exit 1

$POSRE && MDPDEF=-DPOSRES || MDPDEF=
GARBAGE=()
trash() { for item in $@; do GARBAGE[${#GARBAGE[@]}]=$item; done; }

echo "=========================================================="
echo " Running Initram v0.7 with Pure-Python Spatial De-Clasher"
echo "=========================================================="

GRO=$BW
B="$SDIR/backward.py -f $INP -raw $RAW -o $GRO -kick $KICK -sol -p $TOP -po $OTP -n $NDX -from $CG -to $AA"
echo $B; $B || exit 1
trash $RAW $GRO

# --- ZERO-DEPENDENCY PURE PYTHON SPATIAL DE-CLASHER ---
python3 - << 'EOF'
import math
import random

gro_file = '0-backward.gro'
with open(gro_file, 'r') as f:
    lines = f.readlines()

header = lines[:2]
atom_lines = lines[2:-1]
box_line = lines[-1]

coords = []
for line in atom_lines:
    x = float(line[20:28])
    y = float(line[28:36])
    z = float(line[36:44])
    coords.append([x, y, z])

min_dist = 0.105
min_dist_sq = min_dist * min_dist

for iteration in range(120):
    grid = {}
    clashes_found = 0
    for idx, (x, y, z) in enumerate(coords):
        cell = (int(x / min_dist), int(y / min_dist), int(z / min_dist))
        if cell not in grid:
            grid[cell] = []
        grid[cell].append(idx)

    for cell, indices in grid.items():
        cx, cy, cz = cell
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    neighbor = (cx + dx, cy + dy, cz + dz)
                    if neighbor in grid:
                        for i in indices:
                            for j in grid[neighbor]:
                                if i < j:
                                    x1, y1, z1 = coords[i]
                                    x2, y2, z2 = coords[j]
                                    rx, ry, rz = x2 - x1, y2 - y1, z2 - z1
                                    d2 = rx*rx + ry*ry + rz*rz
                                    if d2 < min_dist_sq:
                                        clashes_found += 1
                                        dist = math.sqrt(d2)
                                        if dist < 1e-5:
                                            rx = random.uniform(-1, 1)
                                            ry = random.uniform(-1, 1)
                                            rz = random.uniform(-1, 1)
                                            dist = math.sqrt(rx*rx + ry*ry + rz*rz)
                                        ux, uy, uz = rx/dist, ry/dist, rz/dist
                                        overlap = (min_dist - dist) / 2.0
                                        coords[i][0] -= ux * overlap
                                        coords[i][1] -= uy * overlap
                                        coords[i][2] -= uz * overlap
                                        coords[j][0] += ux * overlap
                                        coords[j][1] += uy * overlap
                                        coords[j][2] += uz * overlap

    if clashes_found == 0:
        print(f"-> All spatial clashes resolved in {iteration} iterations!")
        break

new_atom_lines = []
for idx, line in enumerate(atom_lines):
    x, y, z = coords[idx]
    new_line = line[:20] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + line[44:]
    new_atom_lines.append(new_line)

with open(gro_file, 'w') as f:
    f.writelines(header + new_atom_lines + [box_line])
EOF

i=1

###############################################
## STEP 2: ENERGY MINIMIZATION PASS 1        ##
###############################################
if [ "$EMSTEPS" -gt 0 ]; then
    BASE=$i-EM
    mdp=$BASE.mdp

cat << __MDP__ > $mdp
define                    = -DFLEXIBLE
integrator                = steep
nsteps                    = $EMSTEPS
emstep                    = 0.001
emtol                     = 100.0
pbc                       = xyz
cutoff-scheme             = Verlet
nstlist                   = 10
rlist                     = 1.0
rcoulomb                  = 1.0
rvdw                      = 1.0
coulombtype               = Cut-off
vdwtype                   = Cut-off
constraints               = none
nstxout-compressed        = 0
__MDP__

    G="gmx grompp -f $mdp -c $GRO -n $NDX -p $OTP -o $BASE -maxwarn 2"
    echo $G; $G || exit 1
    M="gmx mdrun -deffnm $BASE -v -nt $NP"
    echo $M; $M || { echo "FATAL: EM1 failed."; exit 1; }

    trash $BASE.*
    GRO=$((i++))-EM.gro
fi

############################################
## STEP 3: ENERGY MINIMIZATION PASS 2     ##
############################################
if [ "$NBSTEPS" -gt 0 ]; then
    BASE=$i-EM
cat << __MDP__ > $BASE.mdp
define                    = -DFLEXIBLE
integrator                = steep
nsteps                    = $NBSTEPS
emstep                    = 0.005
emtol                     = 10.0
pbc                       = xyz
cutoff-scheme             = Verlet
nstlist                   = 10
rlist                     = 1.0
rcoulomb                  = 1.0
rvdw                      = 1.0
coulombtype               = Cut-off
vdwtype                   = Cut-off
constraints               = none
nstxout-compressed        = 0
__MDP__
    mdp=$BASE.mdp

    G="gmx grompp -f $mdp -c $GRO -n $NDX -p $OTP -o $BASE -maxwarn 2"
    echo $G; $G || exit 1
    M="gmx mdrun -deffnm $BASE -v -nt $NP"
    echo $M; $M || { echo "FATAL: EM2 failed."; exit 1; }

    trash $BASE.*
    GRO=$i-EM.gro
fi

##########################################
## STEP 4: MD EQUILIBRATION CYCLES      ##
##########################################
if [ "$MDSTEPS" -gt 0 ]; then
    ifs=$IFS
    IFS=,
    DT=($DT)
    IFS=$ifs

    for DELTA_T in ${DT[@]}; do
        BASE=$((++i))-mdpr-$DELTA_T
        mdp=$BASE.mdp

cat << __MDP__ > $mdp
define                    = $MDPDEF
integrator                = md
nsteps                    = $MDSTEPS
dt                        = $DELTA_T
pbc                       = xyz
cutoff-scheme             = Verlet
nstlist                   = 20
rcoulomb                  = 1.0
rlist                     = 1.0
rvdw                      = 1.0
tcoupl                    = v-rescale
ref_t                     = 300
tau_t                     = 0.1
tc_grps                   = System
gen_vel                   = yes
gen_temp                  = 300
constraints               = h-bonds
nstxout-compressed        = 0
__MDP__

      G="gmx grompp -f $mdp -c $GRO -r $BW -p $OTP -o $BASE -maxwarn 2"
      echo $G; $G || exit 1
      M="gmx mdrun -deffnm $BASE -v -nt $NP"
      echo $M; $M || { echo "FATAL: Step $BASE (MD) failed."; exit 1; }
      trash $BASE.*
      GRO=$BASE.gro
    done
fi

cp $GRO $OUT
rm -f ${GARBAGE[@]}
echo "Backmapping successfully completed!"
"""

@celery_app.task
def run_backmapping_workflow(
    work_dir: str, protein_pdb: str, aptamer_pdb: str, xtc_file: str,
    tpr_file: str, mdp_file: str, ndx_file: str, system_gro: str,
    aptamer_map: str = None, topol_top: str = None,
    loop_from: int = 0, loop_until: int = 10, solute_group: int = 1
):
    wd = Path(work_dir)
    try:
        # 1. Copy Backmapping scripts and Mapping library
        shutil.copy("/app/tools/Backmapping/backward.py", wd / "backward.py")
        shutil.copytree("/app/tools/Backmapping/Mapping", wd / "Mapping", dirs_exist_ok=True)

        if aptamer_map and (wd / aptamer_map).exists():
            shutil.copy2(wd / aptamer_map, wd / "Mapping" / "aptamer.map")

        # 2. Locate and copy the modified charmm36-jul2022.ff into work directory
        ff_source = Path("/app/tools/Backmapping/charmm36-jul2022.ff")
        if not ff_source.exists():
            ff_source = Path("/app/tools/Backmapping/charmm36.ff")

        if not ff_source.exists():
            raise FileNotFoundError("Custom CHARMM36 force field not found in /app/tools/Backmapping/")

        target_ff_dir = wd / ff_source.name
        if not target_ff_dir.exists():
            shutil.copytree(ff_source, target_ff_dir)

        ff_name = ff_source.stem

        # 3. Inject local residuetypes.dat into the task working directory
        res_types_content = """CFZ     RNA
UFT     RNA
UMS     RNA
ADS     RNA
DA      RNA
DG      RNA
DC      RNA
DT      RNA
A       RNA
C       RNA
G       RNA
U       RNA
"""
        with open(wd / "residuetypes.dat", "w") as f:
            f.write(res_types_content)

        # 4. Write executable initram-v5.sh
        initram_path = wd / "initram-v5.sh"
        with open(initram_path, "w") as f:
            f.write(INITRAM_SCRIPT)
        subprocess.run(["chmod", "+x", "initram-v5.sh", "backward.py"], cwd=wd, check=True)

        # 5. Generate All-Atom topol.top and aa.gro natively with modified charmm36
        topol_path = wd / "topol.top"
        if not topol_path.exists():
            clean_prot_lines = []
            with open(wd / protein_pdb, "r") as f:
                for line in f:
                    if line.startswith(("ATOM", "HETATM")):
                        atom_name = line[12:16].strip()

                        # Skip secondary C-terminal oxygens so pdb2gmx builds COO- cleanly
                        if atom_name in {"OXT", "OT2"}:
                            continue

                        # If the terminal oxygen is named OT1, standardize it to O
                        if atom_name == "OT1":
                            line = f"{line[:12]} O   {line[17:]}"

                        clean_prot_lines.append(f"{line[:21]}A{line[22:]}")

            clean_na_lines = []
            with open(wd / aptamer_pdb, "r") as f:
                for line in f:
                    if line.startswith(("ATOM", "HETATM")):
                        clean_na_lines.append(f"{line[:21]}B{line[22:]}")

            combined_pdb = wd / "combined.pdb"
            with open(combined_pdb, "w") as f:
                f.writelines(clean_prot_lines)
                f.write("TER\n")
                f.writelines(clean_na_lines)
                f.write("TER\nEND\n")

            env = os.environ.copy()
            env["GMXLIB"] = str(wd)

            pdb2gmx_cmd = [
                "gmx", "pdb2gmx",
                "-f", "combined.pdb",
                "-o", "aa.gro",
                "-p", "topol.top",
                "-ff", ff_name,
                "-water", "tip3p",
                "-ignh"
            ]
            subprocess.run(
                pdb2gmx_cmd,
                cwd=wd,
                env=env,
                input="1\n1\n1\n1\n1\n",
                capture_output=True,
                text=True,
                timeout=180,
                check=True
            )

        # 6. Output directories for trajectory reconstruction
        out_xtc_dir = wd / "all_xtc"
        out_gro_dir = wd / "gro"
        out_xtc_dir.mkdir(exist_ok=True)
        out_gro_dir.mkdir(exist_ok=True)

        # 7. Execute Backmapping loop
        for i in range(loop_from, loop_until + 1):
            dump_time = f"{i}000"
            trj_gro = f"trj_{i}ns.gro"
            aa_gro = f"aa_{i}ns.gro"
            aa_xtc = f"aa_{i}ns.xtc"

            # Extract CG frame from dynamic_fitted.xtc
            subprocess.run([
                "gmx", "trjconv", "-f", xtc_file, "-s", system_gro, "-n", ndx_file,
                "-o", trj_gro, "-dump", dump_time
            ], cwd=wd, input=f"{solute_group}\n", capture_output=True, text=True, check=True)

            loop_env = os.environ.copy()
            loop_env["GMXLIB"] = str(wd)
            if i == loop_from:
                loop_env["EM_STEPS"] = "3000"
                loop_env["NB_STEPS"] = "1000"
                loop_env["MD_STEPS"] = "500"
            else:
                loop_env["EM_STEPS"] = "2000"
                loop_env["NB_STEPS"] = "0"
                loop_env["MD_STEPS"] = "0"

            # Execute backmapping projection using charmm36
            subprocess.run([
                "./initram-v5.sh", "-f", trj_gro, "-o", aa_gro,
                "-to", "charmm36", "-p", "topol.top"
            ], cwd=wd, env=loop_env, capture_output=True, text=True, check=True)

            # Update Frame Header cleanly
            with open(wd / trj_gro, "r") as f:
                header_line = f.readline().strip()
            with open(wd / aa_gro, "r") as f:
                aa_lines = f.readlines()
            aa_lines[0] = f"Protein-Nucleic Acid Backmapped {i}ns {header_line}\n"
            with open(wd / aa_gro, "w") as f:
                f.writelines(aa_lines)

            # Convert Frame GRO to XTC
            subprocess.run([
                "gmx", "trjconv", "-f", aa_gro, "-o", aa_xtc
            ], cwd=wd, input="0\n", capture_output=True, text=True, check=True)

            shutil.copy(wd / aa_xtc, out_xtc_dir / aa_xtc)
            shutil.copy(wd / aa_gro, out_gro_dir / aa_gro)

        # 8. Concatenate all generated XTC frames
        xtc_files = sorted([f.name for f in out_xtc_dir.glob("*.xtc")])
        if not xtc_files:
            return {"status": "failed", "message": "No all-atom XTC frames were produced."}

        trjcat_cmd = ["gmx", "trjcat", "-f"] + xtc_files + ["-o", "combine.xtc"]
        subprocess.run(trjcat_cmd, cwd=out_xtc_dir, capture_output=True, text=True, check=True)

        # 9. Assemble all-atomic dynamic.tpr
        earliest_gro = f"aa_{loop_from}ns.gro"
        source_gro = out_gro_dir / earliest_gro
        target_gro = wd / earliest_gro
        if source_gro.exists():
            shutil.copy2(source_gro, target_gro)
        else:
            first_gro = sorted(list(out_gro_dir.glob("*.gro")))[0]
            shutil.copy2(first_gro, target_gro)
            earliest_gro = first_gro.name

        mdp_path = wd / mdp_file
        if mdp_path.exists():
            with open(mdp_path, "r") as f:
                lines = f.readlines()
            new_lines = []
            for l in lines:
                if l.strip().startswith("tc-grps"): new_lines.append("tc-grps = System\n")
                elif l.strip().startswith("tau-t"): new_lines.append("tau-t = 0.1\n")
                elif l.strip().startswith("ref-t"): new_lines.append("ref-t = 300\n")
                else: new_lines.append(l)
            with open(mdp_path, "w") as f:
                f.writelines(new_lines)

        grompp_cmd = [
            "gmx", "grompp",
            "-f", mdp_file,
            "-c", earliest_gro,
            "-p", "topol.top",
            "-o", "dynamic.tpr",
            "-maxwarn", "5"
        ]
        subprocess.run(grompp_cmd, cwd=wd, capture_output=True, text=True, check=True)
        shutil.copy2(wd / "dynamic.tpr", out_xtc_dir / "dynamic.tpr")

        # 10. Center-fit PBC and generate index.ndx + concatenated.pdb
        if (out_xtc_dir / "combine.xtc").exists() and (out_xtc_dir / "dynamic.tpr").exists():
            trjconv_pbc_cmd = [
                "gmx", "trjconv",
                "-s", "dynamic.tpr",
                "-f", "combine.xtc",
                "-o", "all.xtc",
                "-pbc", "mol",
                "-center"
            ]
            subprocess.run(trjconv_pbc_cmd, cwd=out_xtc_dir, input="1\n0\n", capture_output=True, text=True, check=True)
            (out_xtc_dir / "combine.xtc").unlink()

        subprocess.run(["gmx", "make_ndx", "-f", "dynamic.tpr", "-o", "index.ndx"], cwd=out_xtc_dir, input="q\n", capture_output=True, text=True, check=True)
        subprocess.run(["gmx", "trjconv", "-s", "dynamic.tpr", "-n", "index.ndx", "-f", "all.xtc", "-dt", "1000", "-o", "concatenated.pdb"], cwd=out_xtc_dir, input="0\n", capture_output=True, text=True, check=True)

        # 11. Clean up intermediate frame files from the root directory
        for f in wd.glob("aa_*ns.gro"): f.unlink()
        for f in wd.glob("aa_*ns.xtc"): f.unlink()
        for f in wd.glob("trj_*ns.gro"): f.unlink()

        # 12. Package final output archive
        zip_path = wd / "backmapped_results.zip"
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root, dirs, files in os.walk(out_xtc_dir):
                for file in files:
                    file_path = os.path.join(root, file)
                    zipf.write(file_path, arcname=f"all_xtc/{file}")
            for root, dirs, files in os.walk(out_gro_dir):
                for file in files:
                    file_path = os.path.join(root, file)
                    zipf.write(file_path, arcname=f"gro/{file}")

        return {
            "status": "success",
            "message": "Backmapping completed successfully.",
            "result_file": str(zip_path)
        }
    except subprocess.CalledProcessError as e:
        err = e.stderr if e.stderr else e.stdout
        return {"status": "failed", "message": f"Command Failed: {err}"}
    except Exception as e:
        return {"status": "failed", "message": str(e)}

run_backmapping_task = run_backmapping_workflow

# ==============================================================================
# TAB 7: TRAJECTORY ANALYSIS (RMSD, RMSF, Rg)
# ==============================================================================
@celery_app.task
def run_trajectory_analysis_task(job_id: str, rmsd_group: int = 4, rmsf_group: int = 1, rg_group: int = 1):
    import matplotlib.pyplot as plt
    os.environ["MPLBACKEND"] = "Agg"

    job_dir = Path(f"/app/shared_data/{job_id}")
    xtc_dir = job_dir / "Backmapping" / "all_xtc" if (job_dir / "Backmapping" / "all_xtc").exists() else job_dir / "all_xtc"
    analysis_dir = job_dir / "analysis" / "trajectory"
    analysis_dir.mkdir(parents=True, exist_ok=True)

    tpr = xtc_dir / "dynamic.tpr"
    xtc = xtc_dir / "all.xtc"
    ndx = xtc_dir / "index.ndx"

    if not (tpr.exists() and xtc.exists()):
        return {"status": "error", "message": f"Missing dynamic.tpr or all.xtc in {xtc_dir}"}

    try:
        subprocess.run(["gmx", "rms", "-s", str(tpr), "-f", str(xtc), "-n", str(ndx), "-o", str(analysis_dir / "rmsd.xvg"), "-tu", "ns"], input=f"{rmsd_group}\n{rmsd_group}\n", text=True, capture_output=True, cwd=analysis_dir, check=True)
        subprocess.run(["gmx", "rmsf", "-s", str(tpr), "-f", str(xtc), "-n", str(ndx), "-o", str(analysis_dir / "rmsf.xvg"), "-res"], input=f"{rmsf_group}\n", text=True, capture_output=True, cwd=analysis_dir, check=True)
        subprocess.run(["gmx", "gyrate", "-s", str(tpr), "-f", str(xtc), "-n", str(ndx), "-o", str(analysis_dir / "gyration.xvg")], input=f"{rg_group}\n", text=True, capture_output=True, cwd=analysis_dir, check=True)

        def plot_xvg(xvg, title, xl, yl, out_png, scale=1.0):
            if not Path(xvg).exists(): return
            x, y = [], []
            with open(xvg, 'r') as f:
                for line in f:
                    if not line.startswith(('@', '#')):
                        p = line.split()
                        if len(p) >= 2: x.append(float(p[0]) * scale); y.append(float(p[1]))
            plt.figure(figsize=(8, 5))
            plt.plot(x, y, color='royalblue', linewidth=1.5)
            plt.title(title, fontsize=13, fontweight='bold')
            plt.xlabel(xl); plt.ylabel(yl)
            plt.grid(True, linestyle='--', alpha=0.6)
            plt.tight_layout()
            plt.savefig(out_png, dpi=300)
            plt.close()

        plot_xvg(analysis_dir / "rmsd.xvg", "Backbone RMSD", "Time (ns)", "RMSD (nm)", analysis_dir / "RMSD.png")
        plot_xvg(analysis_dir / "rmsf.xvg", "Per-Residue RMSF", "Residue Number", "RMSF (nm)", analysis_dir / "RMSF.png")
        plot_xvg(analysis_dir / "gyration.xvg", "Radius of Gyration", "Time (ns)", "Rg (nm)", analysis_dir / "Gyration.png", scale=0.001)

        zip_path = analysis_dir / "trajectory_results.zip"
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as z:
            for f in analysis_dir.glob("*"):
                if f.name != "trajectory_results.zip": z.write(f, arcname=f.name)

        return {"status": "success", "result_file": str(zip_path)}
    except Exception as e:
        return {"status": "error", "message": str(e)}

# ==============================================================================
# TAB 8: PLIP INTERACTION PROFILING
# ==============================================================================
@celery_app.task
def run_plip_analysis_task(job_id: str):
    os.environ["MPLBACKEND"] = "Agg"
    env = os.environ.copy()
    env["MPLBACKEND"] = "Agg"

    job_dir = Path(f"/app/shared_data/{job_id}")
    work_dir = job_dir / "analysis" / "plip"
    work_dir.mkdir(parents=True, exist_ok=True)
    xtc_dir = job_dir / "Backmapping" / "all_xtc" if (job_dir / "Backmapping" / "all_xtc").exists() else job_dir / "all_xtc"

    tpr = xtc_dir / "dynamic.tpr"
    xtc = xtc_dir / "all.xtc"
    ndx = xtc_dir / "index.ndx"

    if not (tpr.exists() and xtc.exists()):
        return {"status": "error", "message": f"Missing trajectories in {xtc_dir}"}

    try:
        # Step 1: Convert trajectory frames to concatenated PDB
        subprocess.run(
            ["gmx", "trjconv", "-s", str(tpr), "-n", str(ndx), "-f", str(xtc), "-dt", "1000", "-o", "concatenated.pdb"],
            input="0\n", cwd=work_dir, check=True, text=True, capture_output=True, env=env
        )

        analysis_dir = work_dir / "Analysis"
        if analysis_dir.exists(): shutil.rmtree(analysis_dir)
        shutil.copytree("/app/tools/Analysis", analysis_dir)
        shutil.move(work_dir / "concatenated.pdb", analysis_dir / "concatenated.pdb")

        split_dir = analysis_dir / "split_pdbs"
        split_dir.mkdir(parents=True, exist_ok=True)

        # Step 2: Dynamic Path Patching (Replace hardcoded Google Drive paths)
        colab_prefixes = [
            "/content/drive/MyDrive/ColabMD-Edu_Protein-NA/Analysis/split_pdbs",
            "/content/drive/MyDrive/ColabMD-Edu_Protein-Protein/NA/split_pdbs",
            "/content/drive/MyDrive/ColabMD-Edu_Protein-NA/Analysis",
            "/content/drive/MyDrive/ColabMD-Edu_Protein-NA",
            "/content/drive/MyDrive/ColabMD-Edu_Protein-Protein/NA",
            "/content/drive/MyDrive/ColabMD-Edu_Protein-Ligand/Analysis",
            "/content/drive/MyDrive/ColabMD-Edu_Protein-Ligand",
        ]

        for script_file in list(analysis_dir.glob("*.py")) + list(analysis_dir.glob("*.sh")):
            with open(script_file, "r") as f:
                content = f.read()

            for old_p in colab_prefixes:
                if "split_pdbs" in old_p:
                    content = content.replace(old_p, str(split_dir))
                else:
                    content = content.replace(old_p, str(analysis_dir))

            # Remove PyMOL visual rendering flags so PLIP runs headlessly
            if script_file.name == "plip.sh":
                content = content.replace("--pics", "").replace("--pymol", "").replace(" -p", "").replace(" -y", "")

            if script_file.name == "timeline_interaction_color.py":
                content = content.replace("plt.tight_layout()", "plt.xticks(rotation=90, fontsize=7)\n    plt.tight_layout()")

            with open(script_file, "w") as f:
                f.write(content)

        # Step 3: Split the concatenated.pdb into split_pdbs/
        subprocess.run([sys.executable, "split_pdb.py"], cwd=analysis_dir, check=True, capture_output=True, text=True, env=env)

        # Step 4: Copy scripts into split_pdbs (Matches Google Colab Notebook Cell 147)
        for sf in list(analysis_dir.glob("*.py")) + list(analysis_dir.glob("*.sh")):
            shutil.copy2(sf, split_dir / sf.name)

        # Step 5: Run PLIP inside split_pdbs where the split frame PDBs live
        if (split_dir / "plip.sh").exists():
            subprocess.run(["bash", "plip.sh"], cwd=split_dir, check=True, capture_output=True, text=True, env=env)
        elif (analysis_dir / "plip.sh").exists():
            subprocess.run(["bash", "plip.sh"], cwd=analysis_dir, check=True, capture_output=True, text=True, env=env)

        # Step 6: Convert XML to JSON in split_pdbs
        if (split_dir / "convert_xml_to_json.py").exists():
            subprocess.run([sys.executable, "convert_xml_to_json.py"], cwd=split_dir, check=True, capture_output=True, text=True, env=env)
        else:
            subprocess.run([sys.executable, "convert_xml_to_json.py"], cwd=analysis_dir, check=True, capture_output=True, text=True, env=env)

        # Sync extracted_data.json between split_pdbs and analysis_dir
        if (split_dir / "extracted_data.json").exists():
            shutil.copy(split_dir / "extracted_data.json", analysis_dir / "extracted_data.json")
        elif (analysis_dir / "extracted_data.json").exists():
            shutil.copy(analysis_dir / "extracted_data.json", split_dir / "extracted_data.json")
        else:
            raise Exception("PLIP completed, but no interaction XML files were produced. Verify PLIP binary is in PATH.")

        # Step 7: Process JSON into CSV folders
        subprocess.run([sys.executable, "process_json.py"], cwd=analysis_dir, check=True, capture_output=True, text=True, env=env)

        # Sync bond_csv_files across both directories
        bond_csv = analysis_dir / "bond_csv_files"
        split_bond_csv = split_dir / "bond_csv_files"
        if bond_csv.exists():
            if split_bond_csv.exists(): shutil.rmtree(split_bond_csv)
            shutil.copytree(bond_csv, split_bond_csv)
        elif split_bond_csv.exists():
            if bond_csv.exists(): shutil.rmtree(bond_csv)
            shutil.copytree(split_bond_csv, bond_csv)

        # Step 8: Generate interaction figures and summary plots
        plot_scripts = [
            "percentage_interactions.py",
            "type_interactions_color.py",
            "residue_interactions_color_csv_a.py",
            "timeline_interaction_color.py"
        ]
        for script in plot_scripts:
            if (split_dir / script).exists():
                subprocess.run([sys.executable, script], cwd=split_dir, check=True, capture_output=True, text=True, env=env)
            elif (analysis_dir / script).exists():
                subprocess.run([sys.executable, script], cwd=analysis_dir, check=True, capture_output=True, text=True, env=env)

        # Move generated plot images to analysis_dir root to prevent duplicates in split_pdbs
        for ext in ["*.png", "*.jpg"]:
            for f in split_dir.glob(ext):
                dest = analysis_dir / f.name
                if dest.exists():
                    dest.unlink()
                shutil.move(str(f), str(dest))

        # Copy CSV files to root while leaving originals in split_pdbs
        for f in split_dir.glob("*.csv"):
            shutil.copy2(f, analysis_dir / f.name)

        # Step 9: Package final output archive
        zip_path = work_dir / "plip_results.zip"
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as z:
            for root, _, files in os.walk(analysis_dir):
                for file in files:
                    fp = os.path.join(root, file)
                    z.write(fp, arcname=os.path.relpath(fp, analysis_dir))

        return {"status": "success", "result_file": str(zip_path)}

    except subprocess.CalledProcessError as e:
        err = e.stderr if e.stderr else e.stdout
        return {"status": "error", "message": f"PLIP Command Failed: {err}"}
    except Exception as e:
        return {"status": "error", "message": str(e)}

# ==============================================================================
# TAB 9: MM-GBSA BINDING ENERGY (PROTEIN - NUCLEIC ACID)
# ==============================================================================
@celery_app.task
def run_mmpbsa_task(job_id: str, ligand_name: str = "DNA", protein_idx: int = 1, ligand_idx: int = 12, start_frame: int = 1, end_frame: int = 0):
    try:
        os.chdir("/app")
    except Exception:
        pass

    # 1. Enforce distinct indices (Failsafe against 12 12 duplicate errors)
    try:
        p_idx, l_idx = int(protein_idx), int(ligand_idx)
    except Exception:
        p_idx, l_idx = 1, 12

    if p_idx == l_idx:
        p_idx, l_idx = 1, 12

    # 2. Setup PATH for AmberTools
    extra_paths = ["/root/micromamba/bin", "/usr/local/bin", os.path.expanduser("~/.local/bin"), "/opt/conda/bin"]
    os.environ["PATH"] = ":".join(extra_paths) + ":" + os.environ.get("PATH", "")

    mmpbsa_exec = shutil.which("gmx_MMPBSA")
    if not mmpbsa_exec:
        for candidate in ["/root/micromamba/bin/gmx_MMPBSA", "/usr/local/bin/gmx_MMPBSA", "/opt/conda/bin/gmx_MMPBSA"]:
            if os.path.exists(candidate):
                mmpbsa_exec = candidate
                break

    if not mmpbsa_exec:
        return {"status": "error", "message": "gmx_MMPBSA executable not found in Colab PATH."}

    job_dir = Path(f"/app/shared_data/{job_id}")
    analysis_dir = job_dir / "analysis"
    mmgbsa_dir = analysis_dir / "MMGBSA"
    mmgbsa_dir.mkdir(parents=True, exist_ok=True)

    try:
        # 3. Locate inputs
        bm_base = job_dir / "Backmapping" if (job_dir / "Backmapping").exists() else job_dir
        xtc_dir = bm_base / "all_xtc" if (bm_base / "all_xtc").exists() else job_dir / "all_xtc"

        tpr_src = xtc_dir / "dynamic.tpr"
        xtc_src = xtc_dir / "all.xtc"
        ndx_src = xtc_dir / "index.ndx"
        gro_src = bm_base / "aa.gro" if (bm_base / "aa.gro").exists() else job_dir / "aa.gro"
        top_src = bm_base / "topol.top" if (bm_base / "topol.top").exists() else job_dir / "topol.top"

        if not (tpr_src.exists() and xtc_src.exists() and top_src.exists()):
            return {"status": "error", "message": f"Missing input files for MM-GBSA in {bm_base}."}

        shutil.copy2(xtc_src, mmgbsa_dir / "all.xtc")
        shutil.copy2(tpr_src, mmgbsa_dir / "dynamic.tpr")
        if ndx_src.exists(): shutil.copy2(ndx_src, mmgbsa_dir / "index.ndx")
        if gro_src.exists(): shutil.copy2(gro_src, mmgbsa_dir / "aa.gro")
        shutil.copy2(top_src, mmgbsa_dir / "topol.top")

        for search_dir in [job_dir, bm_base]:
            for item in search_dir.glob("*.itp"): shutil.copy2(item, mmgbsa_dir / item.name)
            for item in search_dir.glob("*.prm"): shutil.copy2(item, mmgbsa_dir / item.name)
            if (search_dir / "Mapping").exists() and not (mmgbsa_dir / "Mapping").exists():
                shutil.copytree(search_dir / "Mapping", mmgbsa_dir / "Mapping")

        # 4. Generate mmpbsa.in
        xtc_files = list(xtc_dir.glob("aa_*ns.xtc"))
        final_end_frame = end_frame if end_frame > 0 else (len(xtc_files) if xtc_files else 1)

        mmpbsa_template = f"""# General input for MM-GBSA evaluation
&general
sys_name="Protein-DNA_Complex",
startframe={start_frame},
endframe={final_end_frame},
interval=1,
verbose=2,
/
&gb
igb=5,
saltcon=0.150,
/
&decomp
idecomp=1,
dec_verbose=1,
print_res="within 6.0",
/
"""
        with open(mmgbsa_dir / "mmpbsa.in", "w") as f:
            f.write(mmpbsa_template)

        # 5. INJECT CUSTOM DECOMPOSITION GRAPHING SCRIPT
        plot_script = """import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def parse_decomp():
    dat_file = "FINAL_DECOMP_MMPBSA.dat"
    if not os.path.exists(dat_file):
        return
    
    data = []
    with open(dat_file, 'r') as f:
        text = f.read()
        
    deltas = text.split("DELTAS:")
    if len(deltas) < 2:
        return
        
    lines = deltas[1].splitlines()
    for line in lines:
        if line.startswith("Residue") or line.startswith(",") or not line.strip() or line.startswith("Total"): continue
        parts = line.split(',')
        if len(parts) >= 19:
            res_info = parts[0]
            tokens = res_info.split(':')
            if len(tokens) >= 4:
                loc = tokens[0]        # 'R' or 'L'
                resname = tokens[2]    # 'THR' or 'DA'
                resid = tokens[3]      # '57' or '5'
                try:
                    # Index 16 from left, or -3 from right is the Total Avg
                    total_val = float(parts[-3])
                    data.append({'Residue': f"{resname}{resid}", 'Loc': loc, 'Total': total_val})
                except Exception:
                    pass
                        
    if not data: return
    df = pd.DataFrame(data)
    
    # Filter for favorable negative energies and plot top 40 contributors
    df = df[df['Total'] < -0.1].sort_values('Total')
    if df.empty: return
    df_top = df.head(40)
    
    plt.figure(figsize=(14, 7))
    colors = ['#1f77b4' if loc == 'R' else '#ff7f0e' for loc in df_top['Loc']]
    plt.bar(df_top['Residue'], df_top['Total'], color=colors, edgecolor='black')
    
    plt.axhline(0, color='black', linewidth=1)
    plt.xticks(rotation=60, ha='right', fontsize=10)
    plt.ylabel(r"$\Delta G_{bind}$ Contribution (kcal/mol)", fontsize=12)
    plt.title("Top Favorable Per-Residue Decomposition (Protein & DNA Aptamer)", fontsize=14, fontweight='bold')
    
    from matplotlib.patches import Patch
    plt.legend(handles=[
        Patch(facecolor='#1f77b4', label='Protein (Receptor)'),
        Patch(facecolor='#ff7f0e', label='DNA Aptamer (Ligand)')
    ], fontsize=12)
    
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig("decomp_favorable_interactions.png", dpi=300)

if __name__ == "__main__":
    parse_decomp()
"""
        with open(mmgbsa_dir / "plot_decomp.py", "w") as f:
            f.write(plot_script)

        # 6. Execute gmx_MMPBSA
        mmpbsa_cmd = [
            mmpbsa_exec, "-O", "-i", "mmpbsa.in", "-cs", "dynamic.tpr",
            "-ct", "all.xtc", "-ci", "index.ndx", "-cg", str(p_idx), str(l_idx),
            "-cp", "topol.top", "-o", "FINAL_RESULTS_MMPBSA.dat",
            "-eo", "FINAL_RESULTS_MMPBSA.csv", "-do", "FINAL_DECOMP_MMPBSA.dat",
            "-deo", "FINAL_DECOMP_MMPBSA.csv", "-nogui"
        ]

        try:
            subprocess.run(mmpbsa_cmd, cwd=mmgbsa_dir, check=True, capture_output=True, text=True, env=os.environ)
        except subprocess.CalledProcessError as e:
            err_out = e.stderr if e.stderr else e.stdout
            if "********* " in err_out or "ValueError" in err_out:
                for bad_file in mmgbsa_dir.glob("_GMXMMPBSA_*"):
                    if bad_file.is_file():
                        try:
                            with open(bad_file, "r") as f: text = f.read()
                            if "*********" in text:
                                with open(bad_file, "w") as f: f.write(text.replace("*********", " 9999.999"))
                        except UnicodeDecodeError: pass
                subprocess.run(mmpbsa_cmd + ["--rewrite-output"], cwd=mmgbsa_dir, check=True, capture_output=True, text=True, env=os.environ)
            else:
                raise e

        # 7. Generate Decomposition Graph
        subprocess.run(["python", "plot_decomp.py"], cwd=mmgbsa_dir, check=True, capture_output=True, text=True, env=os.environ)

        # 8. Extract Summary and Package Results
        summary_data = ""
        if (mmgbsa_dir / "FINAL_RESULTS_MMPBSA.dat").exists():
            with open(mmgbsa_dir / "FINAL_RESULTS_MMPBSA.dat", "r") as f:
                summary_data = "".join(f.readlines()[-30:])

        zip_path = analysis_dir / "mmgbsa_results.zip"
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root, _, files in os.walk(mmgbsa_dir):
                for file in files:
                    fp = os.path.join(root, file)
                    zipf.write(fp, arcname=os.path.relpath(fp, mmgbsa_dir))

        return {
            "status": "success",
            "message": "MM-GBSA analysis complete.",
            "result_file": str(zip_path),
            "summary_data": summary_data
        }

    except subprocess.CalledProcessError as e:
        err = e.stderr if e.stderr else e.stdout
        return {"status": "error", "message": f"gmx_MMPBSA Error: {err}"}
    except Exception as e:
        return {"status": "error", "message": str(e)}

# ==============================================================================
# TAB 10: PCA & FEL (NO PYMOL)
# ==============================================================================
@celery_app.task
def run_pca_fel_task(job_id: str, backbone_group: int = 4):
    os.environ["MPLBACKEND"] = "Agg"
    env = os.environ.copy()
    env["MPLBACKEND"] = "Agg"

    job_dir = Path(f"/app/shared_data/{job_id}")
    pca_dir = job_dir / "analysis" / "PCA-FEL"
    pca_dir.mkdir(parents=True, exist_ok=True)
    xtc_dir = job_dir / "Backmapping" / "all_xtc" if (job_dir / "Backmapping" / "all_xtc").exists() else job_dir / "all_xtc"

    tpr = xtc_dir / "dynamic.tpr"
    xtc = xtc_dir / "all.xtc"
    ndx = xtc_dir / "index.ndx"

    if not (tpr.exists() and xtc.exists()):
        return {"status": "error", "message": "Missing trajectory files."}

    try:
        grp = f"{backbone_group}\n{backbone_group}\n"
        subprocess.run(["gmx", "covar", "-f", str(xtc), "-s", str(tpr), "-n", str(ndx), "-o", "eigenval.xvg", "-v", "eigenvec.trr", "-av", "average.pdb", "-l", "covar.log", "-b", "0", "-tu", "ns"], input=grp, text=True, capture_output=True, cwd=pca_dir, check=True)

        for fst, lst in [("1", "2"), ("1", "3"), ("2", "3")]:
            subprocess.run(["gmx", "anaeig", "-v", "eigenvec.trr", "-f", str(xtc), "-s", str(tpr), "-n", str(ndx), "-comp", f"eigcomp{fst}-{lst}.xvg", "-rmsf", f"eigrmsf{fst}-{lst}.xvg", "-2d", f"2dproj{fst}-{lst}.xvg", "-b", "0", "-tu", "ns", "-first", fst, "-last", lst], input=grp, text=True, capture_output=True, cwd=pca_dir, check=True)

        subprocess.run(["gmx", "sham", "-f", "2dproj1-2.xvg", "-notime", "-ls", "gibb1-2.xpm"], cwd=pca_dir, check=True, capture_output=True)

        for s in ["xpm2txt.py", "PCA.py", "FEL.py"]:
            if (Path(f"/app/tools/Analysis/{s}")).exists(): shutil.copy(f"/app/tools/Analysis/{s}", pca_dir / s)

        subprocess.run(["python2", "xpm2txt.py", "-f", "gibb1-2.xpm", "-o", "FEL.dat"], cwd=pca_dir, check=True, capture_output=True)
        subprocess.run([sys.executable, "PCA.py"], cwd=pca_dir, check=True, capture_output=True, env=env)
        subprocess.run([sys.executable, "FEL.py"], cwd=pca_dir, check=True, capture_output=True, env=env)

        zip_path = pca_dir / "pca_fel_results.zip"
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as z:
            for root, _, files in os.walk(pca_dir):
                for file in files:
                    if file != "pca_fel_results.zip":
                        z.write(os.path.join(root, file), arcname=file)

        return {"status": "success", "result_file": str(zip_path)}
    except Exception as e:
        return {"status": "error", "message": f"PCA/FEL Error: {str(e)}"}
