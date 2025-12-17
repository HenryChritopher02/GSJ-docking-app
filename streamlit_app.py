import streamlit as st
import subprocess
import os
import zipfile
import shutil
from pathlib import Path
import pandas as pd
import re

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
    "PTP1B (1X70)": {
        "pdbqt": "3duy.pdbqt",
        "config": "3duy.txt",
        "desc": "Protein Tyrosine Phosphatase 1B - Insulin signaling regulation"
    },
    "DPP-4 (4A5S)": {
        "pdbqt": "3ine.pdbqt",
        "config": "3ine.txt",
        "desc": "Dipeptidyl peptidase-4 - Glucose homeostasis"
    },
    "Alpha-Glucosidase (3A4A)": {
        "pdbqt": "3inf.pdbqt",
        "config": "3inf.txt",
        "desc": "Alpha-Glucosidase - Carbohydrate digestion"
    },
    "Alpha-Amylase (1B2Y)": {
        "pdbqt": "3inh.pdbqt",
        "config": "3inh.txt",
        "desc": "Alpha-Amylase - Starch hydrolysis"
    },
    "PPAR-gamma (2PRG)": {
        "pdbqt": "3ooz.pdbqt",
        "config": "3ooz.txt",
        "desc": "Peroxisome proliferator-activated receptor gamma - Insulin sensitization"
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

def display_diabetes_docking_procedure():
    st.header(f"Diabetes Targets Molecular Docking (v{APP_VERSION})")
    st.image("https://raw.githubusercontent.com/HenryChritopher02/GSJ/main/docking-app.png", use_column_width=True)
    st.markdown("---")
    
    # Khởi tạo session state
    if 'docking_results' not in st.session_state:
        st.session_state.docking_results = []
    if 'prepared_ligand_paths' not in st.session_state:
        st.session_state.prepared_ligand_paths = []

    # --- SIDEBAR: Cấu hình ---
    with st.sidebar:
        # st.header("⚙️ Settings")
        # st.subheader("1. System Check")
        vina_ready = check_vina_binary(show_success=False)
        
        st.subheader("Select Targets")
        st.caption("Choose targets associated with Diabetes Type 2:")
        
        selected_targets_keys = st.multiselect(
            "Select Target(s):",
            options=list(DIABETES_TARGETS.keys()),
            default=[list(DIABETES_TARGETS.keys())[0]]
        )
        
        # Nút tải dữ liệu Receptor
        if st.button("Fetch Selected Targets Data", key="fetch_targets_btn"):
            if not selected_targets_keys:
                st.warning("Please select at least one target.")
            else:
                with st.spinner("Downloading receptor and config files..."):
                    download_count = 0
                    for key in selected_targets_keys:
                        info = DIABETES_TARGETS[key]
                        # Tải Receptor (.pdbqt)
                        download_file_from_github(
                            BASE_GITHUB_URL_FOR_DATA, 
                            f"targets/{info['pdbqt']}", # Giả định cấu trúc thư mục trên GitHub
                            info['pdbqt'], 
                            RECEPTOR_DIR_LOCAL
                        )
                        # Tải Config (.txt)
                        download_file_from_github(
                            BASE_GITHUB_URL_FOR_DATA, 
                            f"configs/{info['config']}", 
                            info['config'], 
                            CONFIG_DIR_LOCAL
                        )
                        download_count += 1
                    st.success(f"Successfully checked/downloaded data for {download_count} targets.")

    # --- MAIN UI ---
    
    # 1. INPUT LIGAND SECTION
    st.subheader("📂 1. Ligand Input")
    st.info("Please upload your prepared ligands. Accepted formats: `.pdbqt` (single/multiple) or `.zip` (containing .pdbqt files).")
    
    input_method = st.radio("Upload Method:", ("Upload PDBQT File(s)", "Upload ZIP Archive"), horizontal=True)
    
    new_ligands = []
    
    if input_method == "Upload PDBQT File(s)":
        uploaded_files = st.file_uploader("Select .pdbqt files:", type="pdbqt", accept_multiple_files=True)
        if uploaded_files:
            if st.button("Process PDBQT Files"):
                for up_file in uploaded_files:
                    dest_path = LIGAND_PREP_DIR_LOCAL / up_file.name
                    with open(dest_path, "wb") as f:
                        f.write(up_file.getbuffer())
                    new_ligands.append(str(dest_path))
                st.success(f"Added {len(new_ligands)} ligands.")

    elif input_method == "Upload ZIP Archive":
        uploaded_zip = st.file_uploader("Select .zip file:", type="zip")
        if uploaded_zip:
            if st.button("Process ZIP File"):
                # Xóa thư mục giải nén cũ nếu có
                if ZIP_EXTRACT_DIR_LOCAL.exists():
                    shutil.rmtree(ZIP_EXTRACT_DIR_LOCAL)
                ZIP_EXTRACT_DIR_LOCAL.mkdir(parents=True, exist_ok=True)
                
                # Lưu file zip tạm
                temp_zip_path = LIGAND_UPLOAD_TEMP_DIR / uploaded_zip.name
                with open(temp_zip_path, "wb") as f:
                    f.write(uploaded_zip.getbuffer())
                
                # Giải nén
                try:
                    with zipfile.ZipFile(temp_zip_path, 'r') as zip_ref:
                        zip_ref.extractall(ZIP_EXTRACT_DIR_LOCAL)
                    
                    # Tìm tất cả file pdbqt trong thư mục giải nén
                    for item in ZIP_EXTRACT_DIR_LOCAL.rglob("*.pdbqt"):
                        dest_path = LIGAND_PREP_DIR_LOCAL / item.name
                        shutil.copy(item, dest_path)
                        new_ligands.append(str(dest_path))
                    
                    st.success(f"Extracted and added {len(new_ligands)} ligands from ZIP.")
                except Exception as e:
                    st.error(f"Error processing ZIP: {e}")
                finally:
                    if temp_zip_path.exists(): temp_zip_path.unlink()

    # Cập nhật danh sách ligand vào session state
    if new_ligands:
        current_paths = set(st.session_state.prepared_ligand_paths)
        for p in new_ligands:
            current_paths.add(p)
        st.session_state.prepared_ligand_paths = list(current_paths)

    # Hiển thị trạng thái Ligand hiện tại
    if st.session_state.prepared_ligand_paths:
        with st.expander(f"✅ {len(st.session_state.prepared_ligand_paths)} Ligands Ready for Docking", expanded=False):
            for p in st.session_state.prepared_ligand_paths:
                st.text(Path(p).name)
        
        if st.button("Clear Ligand List", key="clear_lig_list"):
            st.session_state.prepared_ligand_paths = []
            st.experimental_rerun()
    else:
        st.warning("No ligands loaded yet.")

    st.markdown("---")

    # 2. DOCKING EXECUTION SECTION
    st.subheader("🚀 2. Run Docking")
    
    if st.button("Start Screening", type="primary"):
        # Kiểm tra điều kiện
        if not vina_ready:
            st.error("Vina executable is missing.")
        elif not selected_targets_keys:
            st.error("No targets selected.")
        elif not st.session_state.prepared_ligand_paths:
            st.error("No ligands loaded.")
        else:
            # Kiểm tra file Target có tồn tại không
            targets_ready = []
            for t_key in selected_targets_keys:
                t_info = DIABETES_TARGETS[t_key]
                r_path = RECEPTOR_DIR_LOCAL / t_info['pdbqt']
                c_path = CONFIG_DIR_LOCAL / t_info['config']
                if r_path.exists() and c_path.exists():
                    targets_ready.append((t_key, r_path, c_path))
                else:
                    st.error(f"Files missing for {t_key}. Please click 'Fetch Selected Targets Data' in sidebar.")
            
            if len(targets_ready) == len(selected_targets_keys):
                st.info(f"Starting docking: {len(st.session_state.prepared_ligand_paths)} ligands vs {len(targets_ready)} targets.")
                
                # Tạo thanh tiến trình
                progress_bar = st.progress(0)
                status_text = st.empty()
                
                total_tasks = len(st.session_state.prepared_ligand_paths) * len(targets_ready)
                completed_tasks = 0
                results_data = []

                # Tạo thư mục output
                DOCKING_OUTPUT_DIR_LOCAL.mkdir(parents=True, exist_ok=True)

                for lig_path_str in st.session_state.prepared_ligand_paths:
                    lig_path = Path(lig_path_str)
                    lig_name = lig_path.stem
                    
                    row_data = {"Ligand": lig_name}
                    
                    for t_name, r_path, c_path in targets_ready:
                        status_text.text(f"Docking {lig_name} against {t_name}...")
                        
                        # Định nghĩa tên file output
                        out_filename = f"{lig_name}_{DIABETES_TARGETS[t_name]['pdbqt'].replace('.pdbqt', '')}_out.pdbqt"
                        out_path = DOCKING_OUTPUT_DIR_LOCAL / out_filename
                        
                        # Chạy Docking
                        ret_code, stdout, stderr = run_single_docking(
                            VINA_PATH_LOCAL, r_path, lig_path, c_path, out_path
                        )
                        
                        if ret_code == 0 and out_path.exists():
                            score = parse_vina_score_from_file(out_path)
                            row_data[t_name] = score if score is not None else "N/A"
                        else:
                            row_data[t_name] = "Error"
                            print(f"Error docking {lig_name} vs {t_name}: {stderr}")

                        completed_tasks += 1
                        progress_bar.progress(completed_tasks / total_tasks)
                    
                    results_data.append(row_data)

                st.session_state.docking_results = results_data
                status_text.text("Docking completed!")
                st.success("Docking Run Finished.")
                st.balloons()

    # 3. RESULTS SECTION
    if st.session_state.docking_results:
        st.markdown("---")
        st.subheader("📊 3. Docking Results Summary (kcal/mol)")
        
        df_results = pd.DataFrame(st.session_state.docking_results)
        
        # --- BẮT ĐẦU SỬA LỖI ---
        # 1. Xác định các cột chứa điểm số (tất cả trừ cột 'Ligand')
        score_cols = [col for col in df_results.columns if col != 'Ligand']
        
        # 2. Chuyển đổi dữ liệu cột điểm số sang dạng số thực (float)
        # Các giá trị "N/A" hoặc "Error" sẽ bị biến thành NaN để không gây lỗi khi so sánh
        for col in score_cols:
            df_results[col] = pd.to_numeric(df_results[col], errors='coerce')

        # 3. Hiển thị DataFrame với Style đã sửa
        # subset=score_cols: Chỉ tô màu các cột điểm số
        # na_rep="N/A": Hiển thị NaN (lỗi) thành chữ "N/A" cho đẹp
        st.dataframe(
            df_results.style.highlight_min(
                axis=1, 
                color='lightgreen', 
                subset=score_cols 
            ).format(precision=2, na_rep="N/A")
        )
        # --- KẾT THÚC SỬA LỖI ---
        
        # Download button
        csv = convert_df_to_csv(df_results)
        st.download_button(
            label="Download Results as CSV",
            data=csv,
            file_name="diabetes_docking_results.csv",
            mime="text/csv",
        )
        
        st.caption("Lower scores indicate better binding affinity.")

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
