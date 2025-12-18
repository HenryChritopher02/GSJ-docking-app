import streamlit as st
import subprocess
import os
import zipfile
import shutil
from pathlib import Path
import pandas as pd
import re
import plotly.express as px
import py3Dmol
from stmol import showmol

# Giữ lại các import từ file utils cục bộ để tận dụng cấu trúc hiện có
# Lưu ý: Vì logic đơn giản hóa, ta sẽ không dùng hết tất cả biến, nhưng giữ lại import để tránh lỗi
from utils.paths import (
    APP_VERSION, BASE_GITHUB_URL_FOR_DATA, 
    APP_ROOT, VINA_EXECUTABLE_NAME, VINA_PATH_LOCAL,
    RECEPTOR_DIR_LOCAL, CONFIG_DIR_LOCAL,
    LIGAND_PREP_DIR_LOCAL, LIGAND_UPLOAD_TEMP_DIR, ZIP_EXTRACT_DIR_LOCAL,
    DOCKING_OUTPUT_DIR_LOCAL, WORKSPACE_PARENT_DIR
)
from utils.app_utils import (
    initialize_directories, download_file_from_github, 
    check_vina_binary, convert_df_to_csv
)

# --- CẤU HÌNH CÁC MỤC TIÊU TIỂU ĐƯỜNG ---
# Giả định các file này nằm trong thư mục 'receptors' và 'configs' trên GitHub
# Bạn cần đảm bảo tên file trên GitHub khớp với định nghĩa ở đây.
DIABETES_TARGETS = {
    "DPP-4": {
        "pdbqt": "dpp4.pdbqt",
        "config": "dpp4.txt"
    },
    "GLP1-R": {
        "pdbqt": "glp1r.pdbqt",
        "config": "glp1r.txt"
    },
    "PPAR-γ": {
        "pdbqt": "pparg.pdbqt",
        "config": "pparg.txt"
    },
    "SGLT2": {
        "pdbqt": "sglt2.pdbqt",
        "config": "sglt2.txt"
    },
    "SUR1": {
        "pdbqt": "sur1.pdbqt",
        "config": "sur1.txt"
    }
}

def parse_vina_score_from_file(file_path):
    """
    Hàm đọc file output PDBQT và lấy điểm năng lượng liên kết thấp nhất (best affinity).
    """
    best_affinity = None
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('REMARK VINA RESULT'):
                    parts = line.split()
                    # Định dạng thường là: REMARK VINA RESULT: -9.5 0.000 0.000
                    if len(parts) >= 4:
                        best_affinity = float(parts[3])
                    break
    except Exception:
        pass
    return best_affinity

def run_single_docking(vina_path, receptor_path, ligand_path, config_path, output_path):
    """
    Hàm chạy Vina cho 1 cặp Receptor - Ligand.
    """
    cmd = [
        str(vina_path),
        "--receptor", str(receptor_path),
        "--ligand", str(ligand_path),
        "--config", str(config_path),
        "--out", str(output_path),
        "--cpu", "2" # Sử dụng 2 CPU cho mỗi tác vụ để cân bằng
    ]
    
    # Chạy lệnh
    proc = subprocess.run(cmd, capture_output=True, text=True)
    return proc.returncode, proc.stdout, proc.stderr
    
def view_complex(protein_path, ligand_path):
    """
    Generates a 3D visualization of the Protein-Ligand complex
    """
    try:
        with open(protein_path, 'r') as f:
            protein_data = f.read()
        with open(ligand_path, 'r') as f:
            ligand_data = f.read()

        view = py3Dmol.view(width=800, height=500)
        view.addModelsAsFrames(protein_data)
        view.setStyle({'model': -1}, {"cartoon": {'color': 'spectrum'}})
        view.addModelsAsFrames(ligand_data)
        view.setStyle({'model': -1}, {"stick": {'colorscheme': 'greenCarbon'}})
        view.zoomTo()
        showmol(view, height=500, width=800)
    except FileNotFoundError:
        st.error("Could not find PDBQT files for visualization.")
        
def display_diabetes_docking_procedure():
    st.header(f"Diabetes Targets Molecular Docking (v{APP_VERSION})")
    st.image("https://raw.githubusercontent.com/HenryChritopher02/GSJ/main/docking-app.png", use_column_width=True)
    
    # Initialize session state
    if 'docking_results' not in st.session_state:
        st.session_state.docking_results = []
    if 'prepared_ligand_paths' not in st.session_state:
        st.session_state.prepared_ligand_paths = []

    # --- SIDEBAR SETTINGS (Keep this logic) ---
    with st.sidebar:
        vina_ready = check_vina_binary(show_success=False)
        st.subheader("Select Targets")
        st.caption("Choose targets associated with Diabetes Type 2:")
        
        selected_targets_keys = st.multiselect(
            "Select Target(s):",
            options=list(DIABETES_TARGETS.keys()),
            default=[list(DIABETES_TARGETS.keys())[0]]
        )
        
        if st.button("Fetch Selected Targets Data", key="fetch_targets_btn"):
            if not selected_targets_keys:
                st.warning("Please select at least one target.")
            else:
                with st.spinner("Downloading receptor and config files..."):
                    download_count = 0
                    for key in selected_targets_keys:
                        info = DIABETES_TARGETS[key]
                        download_file_from_github(BASE_GITHUB_URL_FOR_DATA, f"targets/{info['pdbqt']}", info['pdbqt'], RECEPTOR_DIR_LOCAL)
                        download_file_from_github(BASE_GITHUB_URL_FOR_DATA, f"configs/{info['config']}", info['config'], CONFIG_DIR_LOCAL)
                        download_count += 1
                    st.success(f"Successfully checked/downloaded data for {download_count} targets.")

    # --- NEW TABS LAYOUT ---
    tab1, tab2, tab3 = st.tabs(["📂 1. Ligand Input", "🚀 2. Run Docking", "📊 3. Analysis & 3D"])

    # --- TAB 1: INPUT ---
    with tab1:
        st.info("Upload ligands in .pdbqt or .zip format.")
        input_method = st.radio("Upload Method:", ("Upload PDBQT File(s)", "Upload ZIP Archive"), horizontal=True)
        new_ligands = []
        
        if input_method == "Upload PDBQT File(s)":
            uploaded_files = st.file_uploader("Select .pdbqt files:", type="pdbqt", accept_multiple_files=True)
            if uploaded_files and st.button("Process PDBQT Files"):
                for up_file in uploaded_files:
                    dest_path = LIGAND_PREP_DIR_LOCAL / up_file.name
                    with open(dest_path, "wb") as f: f.write(up_file.getbuffer())
                    new_ligands.append(str(dest_path))
                st.success(f"Added {len(new_ligands)} ligands.")

        elif input_method == "Upload ZIP Archive":
            uploaded_zip = st.file_uploader("Select .zip file:", type="zip")
            if uploaded_zip and st.button("Process ZIP File"):
                if ZIP_EXTRACT_DIR_LOCAL.exists(): shutil.rmtree(ZIP_EXTRACT_DIR_LOCAL)
                ZIP_EXTRACT_DIR_LOCAL.mkdir(parents=True, exist_ok=True)
                temp_zip_path = LIGAND_UPLOAD_TEMP_DIR / uploaded_zip.name
                with open(temp_zip_path, "wb") as f: f.write(uploaded_zip.getbuffer())
                try:
                    with zipfile.ZipFile(temp_zip_path, 'r') as zip_ref:
                        zip_ref.extractall(ZIP_EXTRACT_DIR_LOCAL)
                    for item in ZIP_EXTRACT_DIR_LOCAL.rglob("*.pdbqt"):
                        dest_path = LIGAND_PREP_DIR_LOCAL / item.name
                        shutil.copy(item, dest_path)
                        new_ligands.append(str(dest_path))
                    st.success(f"Extracted and added {len(new_ligands)} ligands.")
                except Exception as e: st.error(f"Error processing ZIP: {e}")
                finally: 
                    if temp_zip_path.exists(): temp_zip_path.unlink()

        if new_ligands:
            current_paths = set(st.session_state.prepared_ligand_paths)
            for p in new_ligands: current_paths.add(p)
            st.session_state.prepared_ligand_paths = list(current_paths)

        if st.session_state.prepared_ligand_paths:
            st.write(f"**Current Ligands ({len(st.session_state.prepared_ligand_paths)}):**")
            with st.expander("View List"):
                for p in st.session_state.prepared_ligand_paths: st.text(Path(p).name)
            if st.button("Clear Ligand List"):
                st.session_state.prepared_ligand_paths = []
                st.experimental_rerun()

    # --- TAB 2: EXECUTION ---
    with tab2:
        st.write("### Simulation Controls")
        if st.button("Start Screening", type="primary"):
            if not vina_ready: st.error("Vina executable is missing.")
            elif not selected_targets_keys: st.error("No targets selected.")
            elif not st.session_state.prepared_ligand_paths: st.error("No ligands loaded.")
            else:
                targets_ready = []
                for t_key in selected_targets_keys:
                    t_info = DIABETES_TARGETS[t_key]
                    r_path = RECEPTOR_DIR_LOCAL / t_info['pdbqt']
                    c_path = CONFIG_DIR_LOCAL / t_info['config']
                    if r_path.exists() and c_path.exists(): targets_ready.append((t_key, r_path, c_path))
                    else: st.error(f"Files missing for {t_key}.")
                
                if len(targets_ready) == len(selected_targets_keys):
                    st.info(f"Docking {len(st.session_state.prepared_ligand_paths)} ligands vs {len(targets_ready)} targets.")
                    progress_bar = st.progress(0)
                    status_text = st.empty()
                    total_tasks = len(st.session_state.prepared_ligand_paths) * len(targets_ready)
                    completed_tasks = 0
                    results_data = []
                    DOCKING_OUTPUT_DIR_LOCAL.mkdir(parents=True, exist_ok=True)

                    for lig_path_str in st.session_state.prepared_ligand_paths:
                        lig_path = Path(lig_path_str)
                        lig_name = lig_path.stem
                        row_data = {"Ligand": lig_name}
                        
                        for t_name, r_path, c_path in targets_ready:
                            status_text.text(f"Docking {lig_name} against {t_name}...")
                            out_filename = f"{lig_name}_{DIABETES_TARGETS[t_name]['pdbqt'].replace('.pdbqt', '')}_out.pdbqt"
                            out_path = DOCKING_OUTPUT_DIR_LOCAL / out_filename
                            
                            ret_code, stdout, stderr = run_single_docking(VINA_PATH_LOCAL, r_path, lig_path, c_path, out_path)
                            
                            if ret_code == 0 and out_path.exists():
                                score = parse_vina_score_from_file(out_path)
                                row_data[t_name] = score if score is not None else "N/A"
                            else: row_data[t_name] = "Error"
                            completed_tasks += 1
                            progress_bar.progress(completed_tasks / total_tasks)
                        results_data.append(row_data)

                    st.session_state.docking_results = results_data
                    status_text.text("Docking completed!")
                    st.success("Run Finished.")
                    st.balloons()

    # --- TAB 3: ANALYSIS ---
    with tab3:
        if st.session_state.docking_results:
            df_results = pd.DataFrame(st.session_state.docking_results)
            score_cols = [col for col in df_results.columns if col != 'Ligand']
            for col in score_cols: df_results[col] = pd.to_numeric(df_results[col], errors='coerce')

            # 1. HEATMAP TABLE
            st.subheader("🔥 Affinity Heatmap")
            st.dataframe(
                df_results.style.background_gradient(
                    cmap='RdYlGn_r', subset=score_cols, vmin=-12, vmax=-4
                ).format(precision=2, na_rep="N/A"),
                use_container_width=True
            )
            
            col_dl, col_chart = st.columns([1, 2])
            with col_dl:
                csv = convert_df_to_csv(df_results)
                st.download_button("Download CSV", csv, "docking_results.csv", "text/csv")

            # 2. DISTRIBUTION CHART
            st.markdown("---")
            st.subheader("📈 Score Distribution")
            try:
                df_melted = df_results.melt(id_vars=['Ligand'], var_name='Target', value_name='Score')
                df_melted = df_melted.dropna()
                fig = px.box(df_melted, x='Target', y='Score', points="all", color='Target', title="Binding Energy Distribution")
                st.plotly_chart(fig, use_container_width=True)
            except Exception as e:
                st.warning("Not enough data for chart.")

            # 3. 3D VISUALIZATION
            st.markdown("---")
            st.subheader("🧬 3D Complex Visualization")
            
            c1, c2 = st.columns(2)
            with c1:
                selected_ligand = st.selectbox("Select Ligand:", df_results['Ligand'].unique())
            with c2:
                selected_target = st.selectbox("Select Target:", score_cols)

            if st.button("Render 3D Structure"):
                target_info = DIABETES_TARGETS[selected_target]
                receptor_file = RECEPTOR_DIR_LOCAL / target_info['pdbqt']
                
                # Reconstruct output filename
                out_filename = f"{selected_ligand}_{target_info['pdbqt'].replace('.pdbqt', '')}_out.pdbqt"
                docked_ligand_file = DOCKING_OUTPUT_DIR_LOCAL / out_filename

                if receptor_file.exists() and docked_ligand_file.exists():
                    st.write(f"Visualizing: **{selected_ligand}** bound to **{selected_target}**")
                    view_complex(str(receptor_file), str(docked_ligand_file))
                else:
                    st.error(f"Output file not found: {out_filename}. Did the docking finish successfully?")
        else:
            st.info("No docking results to analyze yet. Please run docking in Tab 2.")

def display_about_page():
    st.header("About Diabetes Docking App")
    st.markdown(f"**Diabetes Molecular Docking Suite v{APP_VERSION}**")
    st.markdown("""
    This application is specialized for screening compounds against key therapeutic targets for Type 2 Diabetes.
    
    **Features:**
    - **Focused Targets:** Pre-configured screening against 5 major diabetes-related proteins:
        1. **PTP1B (1X70):** Negative regulator of insulin signaling.
        2. **DPP-4 (4A5S):** Enzyme that degrades incretins.
        3. **Alpha-Glucosidase (3A4A):** Enzyme involved in carbohydrate digestion.
        4. **Alpha-Amylase (1B2Y):** Enzyme involved in starch breakdown.
        5. **PPAR-gamma (2PRG):** Nuclear receptor regulating fatty acid storage and glucose metabolism.
    - **Simplified Input:** Direct upload of `.pdbqt` files or `.zip` archives.
    - **Automated Vina:** Runs AutoDock Vina automatically for all combinations.
    """)

def main():
    st.set_page_config(layout="wide", page_title=f"Diabetes Docking v{APP_VERSION}")
    
    initialize_directories()

    #st.sidebar.image("https://raw.githubusercontent.com/HenryChritopher02/GSJ/main/docking-app.png", width=300)
    st.sidebar.title("Navigation")

    app_mode = st.sidebar.radio(
        "Go to:",
        ("Diabetes Docking", "About"),
    )
    st.sidebar.markdown("---")

    if app_mode == "Diabetes Docking":
        display_diabetes_docking_procedure()
    elif app_mode == "About":
        display_about_page()

if __name__ == "__main__":
    main()
