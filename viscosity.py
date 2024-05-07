##conda install -c bioconda anarci -y
##pip install biopython
import streamlit as st
import subprocess
from Bio import SeqIO

# 定义辅助函数
def extract_sequences_from_fasta(fasta_file):
    """
    从FASTA文件中提取氨基酸序列。
    """
    with open(fasta_file) as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    return [str(record.seq) for record in records]

def calculate_hydrophobicity_index(sequence, amino_acid_props):
    """
    计算给定序列的疏水指数。
    """
    amino_acid_counts = {aa: sequence.count(aa) for aa in amino_acid_props}
    sum_hydrophobic = sum(amino_acid_counts[aa] * props for aa, props in amino_acid_props.items() if aa in 'ACFILPVWY')
    sum_hydrophilic = sum(amino_acid_counts[aa] * props for aa, props in amino_acid_props.items() if aa in 'DEGHKMNPQRST')
    return - (sum_hydrophobic / sum_hydrophilic) if sum_hydrophilic != 0 else 0

def calculate_net_charge(sequence, amino_acid_pKa, pH):
    """
    计算序列在特定pH下的净电荷。
    """
    net_charge = sum((1 / (1 + 10 ** (pKa - pH))) * (
        1 if aa in 'DERK' else
        (1 - 1 if aa == 'H' else 1)  # 修正了组氨酸的处理
    ) for aa, pKa in ((aa, amino_acid_pKa[aa]) for aa in sequence if aa in amino_acid_pKa))
    return net_charge

# 定义氨基酸属性
amino_acid_hydrophobicity = {
    'A': 0.62, 'C': 0.29, 'F': 1.19, 'I': 1.38, 'L': 1.06, 'P': 0.74, 'V': 1.08, 'W': 0.81, 'Y': 0.26,
    'D': -0.90, 'E': -0.74, 'G': -0.72, 'H': -0.40, 'K': -1.50, 'M': 0.64, 'N': -0.78, 'Q': -0.85, 'R': -2.53, 'S': -0.18, 'T': -0.05
}
amino_acid_pKa = {
    'D': 3.9, 'E': 4.3, 'H': 6.0, 'C': 8.3, 'Y': 10.1, 'K': 10.5, 'R': 12.0
}

#定义上传的fasta文件
def load_user_fasta():
    uploaded_file = st.file_uploader("选择您的FASTA文件", type=["fasta"])
    if uploaded_file is not None:
        fasta_content = uploaded_file.getvalue().decode("utf-8")
        temp_file = "temp.fasta"
        with open(temp_file, "w") as f:
            f.write(fasta_content)
        return temp_file
    else:
        return None

# 设置Streamlit应用
st.set_page_config(layout="wide")
st.title('Antibody Fv Sequences Analysis')

# ANARCI命令
anarci_input_file = 'test.fasta'
anarci_output_file = 'output.fasta'
anarci_command = ["ANARCI", "-i", anarci_input_file, "-o", anarci_output_file]

# 运行ANARCI命令
try:
    subprocess.run(anarci_command, check=True)
    st.success("ANARCI command executed successfully.")
except subprocess.CalledProcessError as e:
    st.error(f"ANARCI command failed with error: {e}")
    st.stop()
except FileNotFoundError:
    st.error("ANARCI executable not found. Please ensure that ANARCI is installed and available in your PATH.")
    st.stop()

# 加载和处理ANARCI序列
try:
    light_chain_fv, heavy_chain_fv = extract_sequences_from_fasta(anarci_output_file)
except FileNotFoundError:
    st.error("ANARCI output file not found.")
    st.stop()

# 显示提取的序列
st.subheader('Extracted Sequences')
with st.expander('View Sequences'):
    st.code(f"Light chain:\n{light_chain_fv}\n\nHeavy chain:\n{heavy_chain_fv}")

# 计算和显示参数
pH = st.slider("Select pH value", min_value=0.0, max_value=14.0, value=7.4, step=0.1)
phi = calculate_hydrophobicity_index(light_chain_fv + heavy_chain_fv, amino_acid_hydrophobicity)
net_charge = calculate_net_charge(light_chain_fv + heavy_chain_fv, amino_acid_pKa, pH)
vh_net_charge = calculate_net_charge(heavy_chain_fv, amino_acid_pKa, pH)
vl_net_charge = calculate_net_charge(light_chain_fv, amino_acid_pKa, pH)
fv_csp = vh_net_charge * vl_net_charge

# 显示计算参数
st.subheader('Calculated Parameters')
parameters_table = {
    "Hydrophobicity Index (phi)": phi,
    f"Net Charge at pH {pH}": net_charge,
    "Net Charge VH Domain": vh_net_charge,
    "Net Charge VL Domain": vl_net_charge,
    "FvCSP": fv_csp
}
st.table(parameters_table)

# 计算和显示最终结果
final_result = 10**(0.15 + 1.26 * phi - 0.043 * net_charge - 0.020 * fv_csp)

st.subheader('Final Result')
st.markdown(f"**η,cP(180 mg/mL, 25°C):** {final_result}")