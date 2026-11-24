import streamlit as st
import requests
import time
import zipfile
import io
import os
import re
import pandas as pd
import py3Dmol
import streamlit.components.v1 as components
from Bio.PDB import PDBParser

# ==========================================
# 1. PAGE SETUP (Sets browser tab metadata)
# ==========================================
st.set_page_config(page_title="ColabMD - Protein-Nucleic Acid", page_icon="🧬", layout="wide")

# ==========================================
# 2. 20-USER AUTHENTICATION GATE
# ==========================================
USER_CREDENTIALS = {
    "admin": "colabmd2024",
    "user1": "pass1",   "user2": "pass2",   "user3": "pass3",   "user4": "pass4",
    "user5": "pass5",   "user6": "pass6",   "user7": "pass7",   "user8": "pass8",
    "user9": "pass9",   "user10": "pass10", "user11": "pass11", "user12": "pass12",
    "user13": "pass13", "user14": "pass14", "user15": "pass15", "user16": "pass16",
    "user17": "pass17", "user18": "pass18", "user19": "pass19", "user20": "pass20",
}

if "authenticated" not in st.session_state:
    st.session_state.authenticated = False
if "current_user" not in st.session_state:
    st.session_state.current_user = None

if not st.session_state.authenticated:
    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        st.title("🧬 ColabMD Login")
        st.markdown("Please log in to access the Protein–Nucleic Acid MD Web Service.")
        
        username = st.text_input("Username").strip()
        password = st.text_input("Password", type="password").strip()
        
        if st.button("Login", use_container_width=True):
            if username in USER_CREDENTIALS and USER_CREDENTIALS[username] == password:
                st.session_state.authenticated = True
                st.session_state.current_user = username
                st.rerun()
            else:
                st.error("❌ Invalid username or password.")
    st.stop()

# ==========================================
# 3. SIDEBAR CONFIGURATION
# ==========================================
with st.sidebar:
    # Display the ColabMD-Edu Mascot Image
    st.image("https://raw.githubusercontent.com/tasyriqomar/ColabMD-Edu/main/Mascot_ColabMD-Edu.png", use_container_width=True)
    
    # Authentication Status
    st.success(f"👤 Logged in as: **{st.session_state.current_user}**")
    st.markdown("---")
    
    # Colab Worker Link
    st.header("🚀 Worker Environment")
    st.markdown("""
    - SIGN UP and get ngrok_auth_token with [ngrok](https://ngrok.com)
    - Launch the background worker before running tasks:
    """)
    st.link_button("🔗 Open Google Colab Worker", "https://colab.research.google.com/github/tasyriqomar/ColabMD-Edu_Protein-NA/blob/main/Google_Colab_Worker.ipynb", use_container_width=True)
    st.markdown("---")
    
    # Backend URL Configuration
    st.header("⚙️ Backend Connection")
    API_URL = st.text_input("Colab Ngrok URL", value=st.session_state.get("api_url", "https://YOUR-NGROK-URL.ngrok-free.app/api/v1"))
    st.session_state.api_url = API_URL
    
    if not API_URL.endswith("/api/v1"):
        st.warning("⚠️ URL should end with /api/v1")
        
    st.markdown("---")
    
    # Logout
    if st.button("Logout", use_container_width=True):
        st.session_state.authenticated = False
        st.session_state.current_user = None
        st.rerun()
    
    st.divider()

    if st.button("🔄 Restart (Clear Tabs 1, 3 & 4)", use_container_width=True):
        # 1. Collect any active/completed task IDs so the backend can revoke/clean them
        task_ids = [
            st.session_state.get("tab1_task_id"),
            st.session_state.get("tab3_task_id"),
            st.session_state.get("tab4_task_id"),
        ]
        task_ids = [tid for tid in task_ids if tid]

        # 2. Tell the FastAPI backend to delete Tab 1, 3, 4 files in shared_data & revoke tasks
        api_url = st.session_state.get("api_url")
        username = st.session_state.get("username")
        if api_url and username:
            try:
                requests.post(
                    f"{api_url}/cleanup-workspace",
                    json={"username": username, "tabs": [1, 3, 4], "task_ids": task_ids},
                    timeout=10,
                )
            except Exception as e:
                st.warning(f"Backend cleanup warning: {e}")

        # 3. Preserve login & connection keys; wipe everything else (Tabs 1, 3, 4 data/flags)
        preserve_keys = {"logged_in", "authenticated", "username", "token", "api_url"}
        for key in list(st.session_state.keys()):
            if key not in preserve_keys:
                del st.session_state[key]

        st.success("Tabs 1, 3, and 4 data cleared.")
        st.rerun()

    if st.button("🚪 Logout", use_container_width=True):
        st.session_state.clear()
        st.rerun()

# ==========================================
# 4. STATE PERSISTENCE
# ==========================================
state_keys = [
    "job_id", "aptamer_bytes", "docking_zip_bytes", "docking_scores_df",
    "docking_interactions_df", "rank2_model_name",
    "cg_prep_zip_bytes", "cg_gro_preview", "dna_top_preview", "sys_top_preview",
    "bm_task_id", "traj_zip_bytes", "plip_zip_bytes", "mmgbsa_zip_bytes",
    "pca_fel_zip_bytes", "mmgbsa_summary"
]
for k in state_keys:
    if k not in st.session_state:
        st.session_state[k] = "" if "summary" in k or "job_id" in k or "preview" in k else None

# ==========================================
# 5. VISIBLE FRONT PAGE HEADER
# ==========================================
st.title("🧬 ColabMD - Protein–Nucleic Acid Dynamics")
st.markdown("Automated macromolecular docking, coarse-grained dynamics, backmapping, interaction profiling (PLIP), MM-GBSA free energy calculations, and PCA-FEL analysis.")
st.markdown("---")

# ==========================================
# 6. PIPELINE TABS (1 - 10)
# ==========================================
tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8, tab9, tab10 = st.tabs([
    "1. Aptamer Prep",
    "2. All-atom MD (DNA)",
    "3. Docking Processing",
    "4. Complex CG Prep",
    "5. CGMD Simulation",
    "6. Backmapping",
    "7. Trajectory Analysis",
    "8. PLIP",
    "9. MM-GBSA",
    "10. PCA and FEL"
])

# ------------------------------------------
# TAB 1: APTAMER PREP
# ------------------------------------------
with tab1:
    st.header("Step 1: Aptamer Processing")
    st.markdown("Process DNA/RNA structure and rename residues for GROMACS parameterization.")
    
    st.info("💡 **Don't have a 3D structure yet?**")
    st.markdown("""
    - Predict secondary structure with [RNAfold](http://rna.tbi.univie.ac.at/cgi-bin/RNAWebSuite/RNAfold.cgi)
    - Submit sequence to [RNAComposer](https://rnacomposer.cs.put.poznan.pl/) and download `Aptamer.pdb`
    """)
    st.markdown("---")
        
    aptamer_file = st.file_uploader("Upload Aptamer Structure (PDB)", type=["pdb"], key="apt_upload")
    convert_rna = st.checkbox("Convert RNA (Uracil) to DNA (Thymine)", value=True)
    
    if st.button("Process Aptamer PDB") and aptamer_file:
        st.session_state.aptamer_bytes = None
        with st.spinner("Processing Aptamer on Colab..."):
            files = {"aptamer_pdb": (aptamer_file.name, aptamer_file.getvalue(), "chemical/x-pdb")}
            data = {"convert_rna": convert_rna}
            try:
                res = requests.post(f"{API_URL}/run-aptamer", files=files, data=data)
                if res.status_code == 200:
                    task_id = res.json()["task_id"]
                    status_ph = st.empty()
                    while True:
                        s_res = requests.get(f"{API_URL}/status/{task_id}").json()
                        st_val = s_res["status"]
                        if st_val == "SUCCESS":
                            t_res = s_res.get("result", {})
                            if isinstance(t_res, dict) and t_res.get("status") == "success":
                                dl = requests.get(f"{API_URL}/download-prep-result/{task_id}")
                                if dl.status_code == 200:
                                    st.session_state.aptamer_bytes = dl.content
                                    status_ph.success("✅ Aptamer processed successfully!")
                                    time.sleep(1)
                                    st.rerun()
                            break
                        elif st_val == "FAILED":
                            status_ph.error("❌ Celery Task Failed on Colab.")
                            break
                        time.sleep(2)
            except Exception as e:
                st.error(f"Connection error: {e}")

    if st.session_state.get("aptamer_bytes"):
        st.markdown("---")
        st.download_button(
            label="⬇️ Download aptamer_fixed.pdb",
            data=st.session_state.aptamer_bytes,
            file_name="aptamer_fixed.pdb",
            mime="chemical/x-pdb",
            key="dl_apt_btn"
        )

# ------------------------------------------
# TAB 2: ALL-ATOM MD (DNA)
# ------------------------------------------
with tab2:
    st.header("Step 2: All-atom MD Simulations of DNA Aptamer")
    st.info("Execute DNA/RNA AAMD simulations on Google Colab GPU runtimes.")
    
    st.markdown("""
    ### Pipeline Operations:
    1. **Setup:** Install GROMACS and configure custom DNA/RNA forcefield parameters.
    2. **Preparation:** Generate topology, define the box, solvate, and add ions.
    3. **Equilibration:** Run energy minimization, NVT, and NPT simulations.
    4. **Production:** Execute the main All-Atomic Molecular Dynamics (AAMD) simulation.
    5. **Processing:** Center, fit trajectories, and output final PDB files.
    """)
    
    colab_url_aamd = "https://colab.research.google.com/github/tasyriqomar/ColabMD-Edu_Protein-NA/blob/main/ColabMD_Edu_Protein_NA_Run_AAMD.ipynb"
    st.markdown(f'<a href="{colab_url_aamd}" target="_blank"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab" style="height:38px;"/></a>', unsafe_allow_html=True)
    st.caption("Once production finishes, download `aptamer_relaxed.pdb` and upload into Tab 3")

# ------------------------------------------
# TAB 3: DOCKING PROCESSING (HDOCK + EGNN AI SCORING)
# ------------------------------------------
with tab3:
    st.header("Step 3: Protein–Nucleic Acid Docking (HDOCK + EGNN AI Scoring)")
    st.markdown("Perform ab-initio docking with **HDOCKlite** (Chain A: Protein, Chain B: Aptamer) and re-rank poses using the pretrained **EGNN** scorer.")

    st.info("💡 **Don't have a protein structure?**")
    st.markdown("- Generate one using [AlphaFold2.ipynb](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)")
    st.info("💡 **Renumber residue numbers in PDB file?**")
    st.markdown("- Revise PDB using [www.canoz.com](https://www.canoz.com/sdh/renumberpdbchain.pl)")
    st.markdown("---")

    c1, c2 = st.columns(2)
    with c1:
        prot_dock_file = st.file_uploader("Upload Protein PDB (Receptor → Chain A)", type=["pdb"], key="prot_dock_up")
        nmax_models = st.number_input("Top HDOCK Models to Generate & Score (-nmax)", min_value=10, max_value=100, value=100, step=10)
    with c2:
        dna_dock_file = st.file_uploader("Upload Aptamer PDB (Ligand → Chain B)", type=["pdb"], key="dna_dock_up")
        interaction_cutoff = st.number_input("Interaction Cutoff (Å)", min_value=2.5, max_value=8.0, value=4.5, step=0.5)

    if st.button("Run HDOCK Docking & AI Scoring") and prot_dock_file and dna_dock_file:
        st.session_state.docking_zip_bytes = None
        st.session_state.docking_scores_df = None
        st.session_state.docking_interactions_df = None
        st.session_state.rank2_model_name = ""

        with st.spinner("Running HDOCKlite docking and EGNN AI scoring on Colab..."):
            files = {
                "protein_pdb": (prot_dock_file.name, prot_dock_file.getvalue(), "chemical/x-pdb"),
                "dna_pdb": (dna_dock_file.name, dna_dock_file.getvalue(), "chemical/x-pdb"),
            }
            data = {"nmax": int(nmax_models), "interaction_cutoff": float(interaction_cutoff)}
            try:
                res = requests.post(f"{API_URL}/run-hdock-docking", files=files, data=data)
                if res.status_code == 200:
                    task_id = res.json()["task_id"]
                    st.session_state["tab3_task_id"] = task_id
                    status_ph = st.empty()
                    while True:
                        s_res = requests.get(f"{API_URL}/status/{task_id}").json()
                        st_val = s_res["status"]
                        if st_val == "SUCCESS":
                            t_res = s_res.get("result", {})
                            if isinstance(t_res, dict) and t_res.get("status") == "success":
                                dl = requests.get(f"{API_URL}/download-prep-result/{task_id}")
                                if dl.status_code == 200:
                                    st.session_state.docking_zip_bytes = dl.content
                                    if "scores" in t_res:
                                        st.session_state.docking_scores_df = pd.DataFrame(t_res["scores"])
                                    if "rank2_interactions" in t_res:
                                        st.session_state.docking_interactions_df = pd.DataFrame(t_res["rank2_interactions"])
                                    st.session_state.rank2_model_name = t_res.get("rank2_model", "Rank 2")
                                    status_ph.success("✅ HDOCK docking, EGNN scoring, and Rank 2 interaction profiling complete!")
                                    time.sleep(1)
                                    st.rerun()
                            else:
                                err_msg = t_res.get("message", "Unknown error") if isinstance(t_res, dict) else str(t_res)
                                status_ph.error(f"❌ HDOCK Error: {err_msg}")
                            break
                        elif st_val == "FAILED":
                            status_ph.error("❌ HDOCK Celery task failed on Colab.")
                            break
                        time.sleep(3)
                else:
                    st.error(f"API Error {res.status_code}: {res.text}")
            except Exception as e:
                st.error(f"Connection error: {e}")

    if st.session_state.get("docking_scores_df") is not None:
        st.subheader("🏆 Final EGNN AI Rankings (HDOCK Models)")
        st.dataframe(st.session_state.docking_scores_df, use_container_width=True)

    if st.session_state.get("docking_zip_bytes"):
        st.download_button(
            label="⬇️ Download All HDOCK Results (.zip)",
            data=st.session_state.docking_zip_bytes,
            file_name="hdock_docking_results.zip",
            mime="application/zip",
            key="dl_dock_btn",
        )
        st.info("🔗 **`docked_complex_best.pdb` (Chain A: Protein, Chain B: Aptamer) is automatically loaded into Tab 4.**")

        # --- Visualize Rank 2 Complex & Interactions (py3Dmol) ---
        st.markdown("---")
        rank2_label = st.session_state.get("rank2_model_name", "Rank 2")
        st.subheader(f"🧬 Visualize Rank 2 Complex & Interactions ({rank2_label})")
        try:
            from Bio.PDB import NeighborSearch

            with zipfile.ZipFile(io.BytesIO(st.session_state.docking_zip_bytes)) as z:
                target_pdb = "docked_complex_rank2.pdb" if "docked_complex_rank2.pdb" in z.namelist() else "docked_complex_best.pdb"
                pdb_string = z.read(target_pdb).decode("utf-8")

            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("rank2", io.StringIO(pdb_string))
            model = structure[0]

            prot_atoms = [a for a in model["A"].get_atoms() if a.element != "H"] if "A" in model else []
            na_atoms = [a for a in model["B"].get_atoms() if a.element != "H"] if "B" in model else []

            interacting_prot, interacting_na = {}, {}
            if prot_atoms and na_atoms:
                ns_prot = NeighborSearch(prot_atoms)
                for atom in na_atoms:
                    close_prot = ns_prot.search(atom.get_coord(), float(interaction_cutoff), level="R")
                    if close_prot:
                        na_res = atom.get_parent()
                        na_key = ("B", na_res.id[1])
                        if na_key not in interacting_na:
                            anchor = na_res["P"].get_coord() if "P" in na_res else atom.get_coord()
                            interacting_na[na_key] = (na_res.resname.strip(), anchor)
                        for p_res in close_prot:
                            p_key = ("A", p_res.id[1])
                            if p_key not in interacting_prot:
                                anchor = p_res["CA"].get_coord() if "CA" in p_res else next(p_res.get_atoms()).get_coord()
                                interacting_prot[p_key] = (p_res.resname.strip(), anchor)

            st.caption(f"Found **{len(interacting_prot)}** protein (Chain A) & **{len(interacting_na)}** nucleic acid (Chain B) interacting residues within {interaction_cutoff} Å.")

            view = py3Dmol.view(width=900, height=580)
            view.addModel(pdb_string, "pdb")
            view.setBackgroundColor("white")
            view.setStyle({"chain": "A"}, {"cartoon": {"color": "royalblue", "opacity": 0.85}})
            view.setStyle({"chain": "B"}, {"cartoon": {"color": "darkorange", "opacity": 0.85}})

            for (ch, resi), (resname, coord) in interacting_prot.items():
                view.addStyle({"chain": ch, "resi": resi}, {"stick": {"colorscheme": "cyanCarbon", "radius": 0.2}})
                view.addLabel(
                    f"{resname}{resi}:{ch}",
                    {"position": {"x": float(coord[0]), "y": float(coord[1]), "z": float(coord[2])},
                     "backgroundColor": "#005f73", "fontColor": "white", "fontSize": 11,
                     "backgroundOpacity": 0.85, "borderThickness": 0}
                )

            for (ch, resi), (resname, coord) in interacting_na.items():
                view.addStyle({"chain": ch, "resi": resi}, {"stick": {"colorscheme": "orangeCarbon", "radius": 0.2}})
                view.addLabel(
                    f"{resname}{resi}:{ch}",
                    {"position": {"x": float(coord[0]), "y": float(coord[1]), "z": float(coord[2])},
                     "backgroundColor": "#ae2012", "fontColor": "white", "fontSize": 11,
                     "backgroundOpacity": 0.85, "borderThickness": 0}
                )

            view.zoomTo()
            v_html = view._make_html()
            btn_html = """
            <div style="text-align: center; margin-top: 6px;">
                <button id="download_rank2_btn" style="background-color: #1a73e8; color: white; padding: 8px 16px; font-size: 14px; font-weight: bold; border: none; border-radius: 5px; cursor: pointer;">
                    📸 Download High-Res PNG
                </button>
            </div>
            <script>
            document.getElementById('download_rank2_btn').onclick = function() {
                var c = document.getElementsByTagName('canvas');
                if (c.length > 0) {
                    var link = document.createElement('a');
                    link.download = 'rank2_complex_interactions.png';
                    link.href = c[c.length - 1].toDataURL("image/png").replace("image/png", "image/octet-stream");
                    document.body.appendChild(link);
                    link.click();
                    document.body.removeChild(link);
                }
            };
            </script>
            """
            components.html(v_html + btn_html, height=640)
        except Exception as e:
            st.warning(f"Could not render Rank 2 3D viewer: {e}")

        # --- Protein–Nucleic Acid Interaction Table (CSV) ---
        st.markdown("---")
        st.subheader("📊 Protein–Nucleic Acid Interaction Table (Rank 2)")
        df_int = st.session_state.get("docking_interactions_df")
        if df_int is not None and not df_int.empty:
            st.dataframe(df_int, use_container_width=True)
            csv_bytes = df_int.to_csv(index=False).encode("utf-8")
            st.download_button(
                label="⬇️ Download Interaction Table (CSV)",
                data=csv_bytes,
                file_name="rank2_interactions.csv",
                mime="text/csv",
                key="dl_rank2_csv_btn",
            )
        else:
            st.info(f"No protein–nucleic acid interactions found within {interaction_cutoff} Å.")

# ------------------------------------------
# TAB 4: COMPLEX COARSE-GRAINED PREPARATION
# ------------------------------------------
with tab4:
    st.header("Step 4: Complex Coarse-Grained Preparation")
    st.markdown("Parameterize the best docked complex into coarse-grained representation, solvate, ionize, and construct.")
    
    st.subheader("1. All-Atomic Complex Input")
    has_docking_pose = st.session_state.get("docking_zip_bytes") is not None
    
    input_choice = st.radio(
        "Choose structure source:",
        ["Use best docked pose from Tab 3 (docked_complex_best.pdb)", "Upload custom docked complex (.pdb)"],
        index=0 if has_docking_pose else 1
    )
    
    selected_complex_bytes = None
    if input_choice.startswith("Use best"):
        if has_docking_pose:
            with zipfile.ZipFile(io.BytesIO(st.session_state.docking_zip_bytes)) as z:
                selected_complex_bytes = z.read("docked_complex_best.pdb")
            st.success("✅ `docked_complex_best.pdb` loaded automatically from Tab 3 memory!")
        else:
            st.warning("⚠️ No docked pose found in memory. Run docking in Tab 3 or upload a PDB file.")
    else:
        custom_pdb = st.file_uploader("Upload Complex PDB (Chain A: Protein, Chain B: DNA/RNA)", type=["pdb"], key="custom_complex_cg_up")
        if custom_pdb:
            selected_complex_bytes = custom_pdb.getvalue()

    c1, c2 = st.columns(2)
    with c1:
        dnatype = st.selectbox("Nucleic Acid Strand Type", ["ss", "ds"], index=0, help="ss: Single-stranded, ds: Double-stranded")
    with c2:
        elastic_net = st.checkbox("Enable Elastic Network for Protein (-elastic)", value=True)

    if st.button("Generate Coarse-Grained Complex & Solvated System") and selected_complex_bytes:
        st.session_state.cg_prep_zip_bytes = None
        st.session_state.cg_gro_preview = ""
        st.session_state.sys_top_preview = ""
        
        with st.spinner("Generating CG models, running insane.py (10x10x10 nm box, 0.15 M NaCl), and writing system.top..."):
            files = {"complex_pdb": ("docked_complex_best.pdb", selected_complex_bytes, "chemical/x-pdb")}
            data = {"dnatype": dnatype, "elastic": elastic_net}
            try:
                res = requests.post(f"{API_URL}/run-complex-cg-prep", files=files, data=data)
                if res.status_code == 200:
                    task_id = res.json()["task_id"]
                    status_ph = st.empty()
                    while True:
                        s_res = requests.get(f"{API_URL}/status/{task_id}").json()
                        st_val = s_res["status"]
                        if st_val == "SUCCESS":
                            t_res = s_res.get("result", {})
                            if isinstance(t_res, dict) and t_res.get("status") == "success":
                                dl = requests.get(f"{API_URL}/download-prep-result/{task_id}")
                                if dl.status_code == 200:
                                    st.session_state.cg_prep_zip_bytes = dl.content
                                    st.session_state.cg_gro_preview = t_res.get("gro_preview", "")
                                    st.session_state.sys_top_preview = t_res.get("sys_top_preview", "")
                                    status_ph.success("✅ Coarse-graining, solvation (system.gro), and system.top generation complete!")
                                    time.sleep(1)
                                    st.rerun()
                            else:
                                err_m = t_res.get("message", "Task failed.") if isinstance(t_res, dict) else str(t_res)
                                status_ph.error(f"❌ Error: {err_m}")
                            break
                        elif st_val == "FAILED":
                            status_ph.error("❌ Coarse-graining task failed on Colab.")
                            break
                        time.sleep(3)
            except Exception as e:
                st.error(f"Connection error: {e}")

    # 4.2 Topology Mapping & Previews
    if st.session_state.get("cg_prep_zip_bytes"):
        st.markdown("---")
        st.subheader("2. Generated `system.top` Topology")
        st.code(st.session_state.get("sys_top_preview", ""), language="ini")

        st.markdown("---")
        st.subheader("3. Solvated & Ionized Coordinate File (`system.gro`)")
        st.text_area("system.gro Preview (Header & First Atoms)", value=st.session_state.get("cg_gro_preview", ""), height=220)

        st.markdown("---")
        st.download_button(
            label="⬇️ Download All CG Preparation Files (.zip)",
            data=st.session_state.cg_prep_zip_bytes,
            file_name="complex_cg_prep.zip",
            mime="application/zip",
            key="dl_cg_prep_btn"
        )

# ------------------------------------------
# TAB 5: COARSE-GRAINED MD (CGMD)
# ------------------------------------------
with tab5:
    st.header("Step 5: Coarse-Grained MD (CGMD)")
    st.info("Execute long-timescale CGMD simulations on Google Colab GPU runtimes.")
    
    st.markdown("""
    ### Pipeline Operations:
    1. **Setup:** Install GPU-enabled GROMACS and clone the CGMD workspace.
    2. **Minimization:** Compile topology and run energy minimization.
    3. **Indexing:** Group Protein with DNA/RNA, and Water with Ions.
    4. **Simulation:** Execute equilibration, followed by the 100ns production MD.
    5. **Processing:** Center, fit trajectories to the backbone, and export as a PDB file.
    """)
    
    colab_url_cgmd = "https://colab.research.google.com/github/tasyriqomar/ColabMD-Edu_Protein-NA/blob/main/ColabMD_Edu_Protein_NA_Run_CGMD.ipynb"
    st.markdown(f'<a href="{colab_url_cgmd}" target="_blank"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab" style="height:38px;"/></a>', unsafe_allow_html=True)
    st.caption("After completing CGMD, download `dynamic_fitted.xtc`, `dynamic.tpr`, `dynamic.mdp`, `index.ndx`, and `dynamic.pdb` and upload into Tab 6.")

# ------------------------------------------
# TAB 6: BACKMAPPING (CG TO ALL-ATOM RESOLUTION)
# ==========================================
with tab6:
    st.header("Step 6: Backmapping (CG to All-Atom Resolution)")
    st.markdown("Upload your CGMD simulation files and topologies to convert back to all-atom resolution.")

    col1, col2 = st.columns(2)
    with col1:
        protein_pdb_bm = st.file_uploader("1. Upload protein.pdb (AA)", type=["pdb"], key="prot_bm")
        aptamer_pdb_bm = st.file_uploader("2. Upload aptamer_relaxed.pdb (AA)", type=["pdb"], key="apt_bm")
        xtc_file = st.file_uploader("3. Upload dynamic_fitted.xtc", type=["xtc"], key="xtc_bm")
        tpr_file = st.file_uploader("4. Upload dynamic.tpr (CG)", type=["tpr"], key="tpr_bm")
        mdp_file = st.file_uploader("5. Upload dynamic.mdp", type=["mdp"], key="mdp_bm")
    with col2:
        ndx_file = st.file_uploader("6. Upload index.ndx (CG)", type=["ndx", "txt"], key="ndx_bm")
        system_gro = st.file_uploader("7. Upload system.gro (CGMD)", type=["gro"], key="sys_gro_bm")
        aptamer_map = st.file_uploader("8. Upload aptamer.map (from Tab 4, Optional)", type=["map"], key="map_bm")
        topol_top = st.file_uploader("9. Upload topol.top (AA, Optional)", type=["top"], key="top_bm")

        st.markdown("#### Frame Loop Settings")
        c_l1, c_l2 = st.columns(2)
        with c_l1:
            loop_from = st.number_input("Loop from (ns):", min_value=0, max_value=1000, value=0, step=1, key="bm_loop_from")
        with c_l2:
            loop_until = st.number_input("Loop until (ns):", min_value=1, max_value=1000, value=10, step=1, key="bm_loop_until")
            
        solute_group = st.number_input("CG Solute Group Index (Default: 1)", value=1, step=1, key="bm_solute_group")

    all_files = [protein_pdb_bm, aptamer_pdb_bm, xtc_file, tpr_file, mdp_file, ndx_file, system_gro]

    if st.button("Run Backmapping"):
        st.session_state.bm_task_id = None

        if not all(all_files):
            st.warning("⚠️ Please upload all required files (1–7).")
        elif loop_from > loop_until:
            st.error("❌ 'Loop from' must be less than or equal to 'Loop until'.")
        else:
            with st.spinner(f"Running Backmapping loop from {loop_from} ns to {loop_until} ns on Colab..."):
                files = {
                    "protein_pdb": (protein_pdb_bm.name, protein_pdb_bm.getvalue(), "chemical/x-pdb"),
                    "aptamer_pdb": (aptamer_pdb_bm.name, aptamer_pdb_bm.getvalue(), "chemical/x-pdb"),
                    "xtc_file": (xtc_file.name, xtc_file.getvalue(), "application/octet-stream"),
                    "tpr_file": (tpr_file.name, tpr_file.getvalue(), "application/octet-stream"),
                    "mdp_file": (mdp_file.name, mdp_file.getvalue(), "text/plain"),
                    "ndx_file": (ndx_file.name, ndx_file.getvalue(), "text/plain"),
                    "system_gro": (system_gro.name, system_gro.getvalue(), "text/plain")
                }
                if aptamer_map:
                    files["aptamer_map"] = (aptamer_map.name, aptamer_map.getvalue(), "text/plain")
                if topol_top:
                    files["topol_top"] = (topol_top.name, topol_top.getvalue(), "text/plain")

                data = {
                    "loop_from": loop_from,
                    "loop_until": loop_until,
                    "solute_group": solute_group
                }

                try:
                    response = requests.post(f"{API_URL}/run-backmapping", files=files, data=data)

                    if response.status_code == 200:
                        task_id = response.json()["task_id"]
                        status_placeholder = st.empty()

                        while True:
                            status_res = requests.get(f"{API_URL}/status/{task_id}").json()
                            status = status_res["status"]

                            if status == "SUCCESS":
                                task_result = status_res.get("result", {})
                                if task_result.get("status") == "success":
                                    res_file_path = task_result.get("result_file", "")
                                    if res_file_path:
                                        extracted_job_id = res_file_path.replace('\\', '/').split('/')[-2]

                                        # Auto-propagate Job ID to master state and downstream analysis tabs
                                        st.session_state["job_id"] = extracted_job_id
                                        st.session_state["traj_job_id"] = extracted_job_id
                                        st.session_state["plip_job_id"] = extracted_job_id
                                        st.session_state["mmgbsa_job_id"] = extracted_job_id
                                        st.session_state["pca_job_id"] = extracted_job_id

                                    # Save task ID to keep the download section visible across tab navigation
                                    st.session_state["bm_task_id"] = task_id
                                    time.sleep(1)
                                    st.rerun()
                                else:
                                    status_placeholder.error(f"❌ Failed: {task_result.get('message')}")
                                break
                            elif status == "FAILED":
                                status_placeholder.error("❌ Celery Task Failed on Colab.")
                                break
                            else:
                                status_placeholder.info(f"⏳ Status: {status} (Reconstructing all-atom frames {loop_from} to {loop_until} ns...)")
                                time.sleep(3)
                    else:
                        st.error(f"Error {response.status_code}: {response.text}")
                except Exception as e:
                    st.error(f"Connection error: {e}")

    # --- PERSISTENT DOWNLOAD & LINKING DISPLAY (Rendered outside button click) ---
    if st.session_state.get("bm_task_id") and st.session_state.get("job_id"):
        st.markdown("---")
        st.success(f"✅ Backmapping Complete! **Job ID:** `{st.session_state.job_id}`")

        dl_url = f"{API_URL}/download-prep-result/{st.session_state.bm_task_id}"
        st.markdown(f"### [⬇️ Click Here to Download Backmapped Trajectories (.zip)]({dl_url})")
        st.info("🔗 **Job ID auto-linked to Tabs 7–10!** You can now open Tabs 7 through 10 directly without manually re-entering your Job ID.")

# ------------------------------------------
# TAB 7: TRAJECTORY ANALYSIS
# ------------------------------------------
with tab7:
    st.header("Step 7: Trajectory Analysis (RMSD, RMSF, Rg)")
    
    job_id_traj = st.text_input("Job ID", value=st.session_state.get("job_id", ""), key="traj_job_id")
    
    if not job_id_traj:
        st.info("⏳ Waiting for Backmapping (Tab 6) or paste your Job ID above.")
    else:
        st.session_state.job_id = job_id_traj
        col1, col2, col3 = st.columns(3)
        with col1:
            rmsd_g = st.number_input("RMSD Group (4=Backbone)", value=4, step=1)
        with col2:
            rmsf_g = st.number_input("RMSF Group (1=Protein)", value=1, step=1)
        with col3:
            rg_g = st.number_input("Rg Group (1=Protein)", value=1, step=1)
            
        if st.button("Calculate Trajectory Dynamics"):
            st.session_state.traj_zip_bytes = None
            with st.spinner("Computing RMSD, RMSF, and Radius of Gyration..."):
                payload = {"job_id": job_id_traj, "rmsd_group": rmsd_g, "rmsf_group": rmsf_g, "rg_group": rg_g}
                res = requests.post(f"{API_URL}/analysis/trajectory", json=payload)
                if res.status_code == 200:
                    task_id = res.json()["task_id"]
                    while True:
                        s_res = requests.get(f"{API_URL}/status/{task_id}").json()
                        if s_res["status"] == "SUCCESS":
                            dl = requests.get(f"{API_URL}/download-prep-result/{task_id}")
                            st.session_state.traj_zip_bytes = dl.content
                            time.sleep(1)
                            st.rerun()
                            break
                        time.sleep(3)

    if st.session_state.get("traj_zip_bytes"):
        st.markdown("---")
        st.download_button(
            label="⬇️ Download GROMACS Analysis Results (.zip)",
            data=st.session_state.traj_zip_bytes,
            file_name=f"{job_id_traj}_gromacs_analysis.zip",
            mime="application/zip",
            key="traj_download"
        )
        try:
            with zipfile.ZipFile(io.BytesIO(st.session_state.traj_zip_bytes)) as z:
                pngs = [f for f in z.namelist() if f.endswith('.png')]
                if pngs:
                    cols = st.columns(3)
                    for i, png in enumerate(sorted(pngs)):
                        with cols[i % 3]:
                            st.image(z.read(png), caption=os.path.basename(png), use_container_width=True)
        except Exception as e:
            st.warning(f"Error previewing plots: {e}")

# ------------------------------------------
# TAB 8: PLIP INTERACTION PROFILING
# ------------------------------------------
with tab8:
    st.header("Step 8: PLIP (Protein–Nucleic Acid Interaction Profiler)")
    
    job_id_plip = st.text_input("Job ID", value=st.session_state.get("job_id", ""), key="plip_job_id")
    
    if not job_id_plip:
        st.info("⏳ Waiting for Backmapping (Tab 6) or paste your Job ID above.")
    else:
        st.session_state.job_id = job_id_plip
        if st.button("Run PLIP Analysis Pipeline"):
            st.session_state.plip_zip_bytes = None
            with st.spinner("Processing trajectory frames and extracting chemical interactions..."):
                res = requests.post(f"{API_URL}/analysis/plip", json={"job_id": job_id_plip})
                if res.status_code == 200:
                    task_id = res.json()["task_id"]
                    status_ph = st.empty()
                    while True:
                        s_res = requests.get(f"{API_URL}/status/{task_id}").json()
                        st_val = s_res["status"]
                        
                        if st_val == "SUCCESS":
                            t_res = s_res.get("result", {})
                            
                            # Safely check if the backend Python function succeeded before downloading
                            if isinstance(t_res, dict) and t_res.get("status") in ["success", "SUCCESS"]:
                                dl = requests.get(f"{API_URL}/download-prep-result/{task_id}")
                                
                                # Verify it's not a JSON error response
                                if dl.status_code == 200 and not dl.content.startswith(b'{"error"'):
                                    st.session_state.plip_zip_bytes = dl.content
                                    status_ph.success("✅ PLIP Analysis Complete!")
                                    time.sleep(1)
                                    st.rerun()
                                else:
                                    status_ph.error("❌ Download failed or file was missing on server.")
                            else:
                                err_m = t_res.get("message", "Task failed silently.") if isinstance(t_res, dict) else str(t_res)
                                status_ph.error(f"❌ PLIP Error: {err_m}")
                            break
                        elif st_val == "FAILED":
                            status_ph.error("❌ Celery Task Failed on Colab.")
                            break
                        else:
                            status_ph.info(f"⏳ Status: {st_val} (Running PLIP Analysis...)")
                            time.sleep(3)
                else:
                    st.error(f"API Error {res.status_code}")

    if st.session_state.get("plip_zip_bytes"):
        st.markdown("---")
        st.download_button(
            label="⬇️ Download All PLIP Results (.zip)",
            data=st.session_state.plip_zip_bytes,
            file_name=f"{job_id_plip}_plip_results.zip",
            mime="application/zip",
            key="plip_download"
        )
        try:
            with zipfile.ZipFile(io.BytesIO(st.session_state.plip_zip_bytes)) as z:
                # Deduplicate by basename so each plot is only rendered once
                seen_filenames = set()
                unique_png_paths = []
                
                for f in z.namelist():
                    if f.endswith(('.png', '.jpg')):
                        base_name = os.path.basename(f)
                        if base_name and base_name not in seen_filenames:
                            seen_filenames.add(base_name)
                            unique_png_paths.append(f)
                
                if unique_png_paths:
                    cols = st.columns(2)
                    for i, p in enumerate(sorted(unique_png_paths)):
                        with cols[i % 2]:
                            st.image(z.read(p), caption=os.path.basename(p), use_container_width=True)
        except Exception as e:
            st.error(f"Error previewing interaction plots: {e}")

# ------------------------------------------
# TAB 9: MM-GBSA BINDING ENERGY
# ------------------------------------------
with tab9:
    st.header("Step 9: MM-GBSA & Per-Residue Decomposition")
    
    job_id_mm = st.text_input("Job ID", value=st.session_state.get("job_id", ""), key="mmgbsa_job_id")
    
    if not job_id_mm:
        st.info("⏳ Waiting for Backmapping (Tab 6) or paste your Job ID above.")
    else:
        st.session_state.job_id = job_id_mm
        col1, col2 = st.columns(2)
        with col1:
            prot_idx = st.number_input("Receptor: Protein Group Index (from index.ndx)", min_value=0, max_value=50, value=1, step=1, key="receptor_group_in")
        with col2:
            na_idx = st.number_input("Ligand: DNA Aptamer Group Index (from index.ndx)", min_value=0, max_value=50, value=12, step=1, key="ligand_group_in")

        if st.button("Execute MM-GBSA Calculation"):
            if int(prot_idx) == int(na_idx):
                st.error("❌ Receptor and Ligand groups must be different! Protein is Group 1 and DNA is Group 12.")
            else:
                st.session_state.mmgbsa_zip_bytes = None
                st.session_state.mmgbsa_summary = ""
                with st.spinner("Computing binding free energies via gmx_MMPBSA on Colab..."):
                    payload = {
                        "job_id": job_id_mm,
                        "ligand_name": "DNA",
                        "protein_idx": int(prot_idx),
                        "ligand_idx": int(na_idx)
                    }
                    res = requests.post(f"{API_URL}/analysis/mmpbsa", json=payload)
                    if res.status_code == 200:
                        task_id = res.json()["task_id"]
                        status_ph = st.empty()
                        while True:
                            s_res = requests.get(f"{API_URL}/status/{task_id}").json()
                            st_val = s_res["status"]
                            if st_val == "SUCCESS":
                                t_res = s_res.get("result", {})
                                if isinstance(t_res, dict) and t_res.get("status") == "success":
                                    dl = requests.get(f"{API_URL}/download-prep-result/{task_id}")
                                    if dl.status_code == 200 and dl.content.startswith(b"PK"):
                                        st.session_state.mmgbsa_zip_bytes = dl.content
                                        st.session_state.mmgbsa_summary = t_res.get("summary_data", "")
                                        status_ph.success("✅ MM-GBSA Calculation Complete!")
                                        time.sleep(1)
                                        st.rerun()
                                    else:
                                        status_ph.error("❌ Received invalid file from server.")
                                else:
                                    err_msg = t_res.get("message", "MM-GBSA calculation failed.") if isinstance(t_res, dict) else str(t_res)
                                    status_ph.error(f"❌ MM-GBSA Error: {err_msg}")
                                break
                            elif st_val == "FAILED":
                                status_ph.error("❌ Celery Task Failed on Colab.")
                                break
                            else:
                                status_ph.info(f"⏳ Status: {st_val} (Evaluating enthalpy & solvation free energies...)")
                                time.sleep(3)
                    else:
                        st.error(f"API Error {res.status_code}: {res.text}")

    if st.session_state.get("mmgbsa_summary"):
        st.subheader("📊 Binding Free Energy Summary ($\Delta G_{\\text{bind}}$)")
        st.code(st.session_state.mmgbsa_summary, language="text")

    if st.session_state.get("mmgbsa_zip_bytes"):
        st.markdown("---")
        st.download_button(
            label="⬇️ Download MM-GBSA Results (.zip)",
            data=st.session_state.mmgbsa_zip_bytes,
            file_name=f"{job_id_mm}_mmgbsa_results.zip",
            mime="application/zip",
            key="mmgbsa_download"
        )
        
        # Display the custom Protein-DNA decomposition graph
        try:
            with zipfile.ZipFile(io.BytesIO(st.session_state.mmgbsa_zip_bytes)) as z:
                if "decomp_favorable_interactions.png" in z.namelist():
                    st.subheader("🔬 Protein-DNA Favorable Interactions (Decomposition)")
                    st.image(z.read("decomp_favorable_interactions.png"), use_container_width=True)
        except Exception as e:
            pass

# ------------------------------------------
# TAB 10: PCA AND FEL (NO PYMOL)
# ------------------------------------------
with tab10:
    st.header("Step 10: Principal Component Analysis (PCA) & Free Energy Landscape (FEL)")
    
    job_id_pca = st.text_input("Job ID", value=st.session_state.get("job_id", ""), key="pca_job_id")
    
    if not job_id_pca:
        st.info("⏳ Waiting for Backmapping (Tab 6) or paste your Job ID above.")
    else:
        st.session_state.job_id = job_id_pca
        bb_group = st.number_input("Backbone Index (Default: 4)", value=4, step=1, key="bb_grp_in")
        
        if st.button("Run PCA & FEL Calculation"):
            st.session_state.pca_fel_zip_bytes = None
            with st.spinner("Calculating covariance matrix, 2D projections, and free energy surface..."):
                payload = {"job_id": job_id_pca, "backbone_group": bb_group}
                res = requests.post(f"{API_URL}/analysis/pca-fel", json=payload)
                if res.status_code == 200:
                    task_id = res.json()["task_id"]
                    while True:
                        s_res = requests.get(f"{API_URL}/status/{task_id}").json()
                        if s_res["status"] == "SUCCESS":
                            dl = requests.get(f"{API_URL}/download-prep-result/{task_id}")
                            st.session_state.pca_fel_zip_bytes = dl.content
                            time.sleep(1)
                            st.rerun()
                            break
                        time.sleep(3)

    if st.session_state.get("pca_fel_zip_bytes"):
        st.markdown("---")
        st.download_button(
            label="⬇️ Download PCA & FEL Figures (.zip)",
            data=st.session_state.pca_fel_zip_bytes,
            file_name=f"{job_id_pca}_pca_fel_results.zip",
            mime="application/zip",
            key="pca_fel_download"
        )
        try:
            with zipfile.ZipFile(io.BytesIO(st.session_state.pca_fel_zip_bytes)) as z:
                pngs = [f for f in z.namelist() if f.endswith('.png')]
                if pngs:
                    cols = st.columns(2)
                    for i, p in enumerate(sorted(pngs)):
                        with cols[i % 2]:
                            st.image(z.read(p), caption=os.path.basename(p), use_container_width=True)
        except Exception as e:
            st.warning(f"Error previewing PCA plots: {e}")
