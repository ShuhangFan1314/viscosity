import streamlit as st
import subprocess
from Bio import SeqIO
import pandas as pd
import os

class AminoAcidProperties:
    def __init__(self):
        self.hydrophobicity = {
            'A': 0.62, 'C': 0.29, 'F': 1.19, 'I': 1.38, 'L': 1.06, 'P': 0.74, 'V': 1.08, 'W': 0.81, 'Y': 0.26,
            'D': -0.90, 'E': -0.74, 'G': -0.72, 'H': -0.40, 'K': -1.50, 'M': 0.64, 'N': -0.78, 'Q': -0.85, 'R': -2.53, 'S': -0.18, 'T': -0.05
        }
        self.pKa_values = {
            'D': 3.9, 'E': 4.3, 'H': 6.0, 'C': 8.3, 'Y': 10.1, 'K': 10.5, 'R': 12.0
        }

    def calculate_hydrophobicity_index(self, sequence):
        return self._calculate_index(sequence, self.hydrophobicity)

    def calculate_net_charge(self, sequence, pH):
        return self._calculate_net_charge(sequence, self.pKa_values, pH)

    def _calculate_index(self, sequence, property_dict):
        amino_acid_counts = {aa: sequence.count(aa) for aa in property_dict}
        sum_hydrophobic = sum(amino_acid_counts[aa] * prop for aa, prop in property_dict.items() if aa in 'ACFILPVWY')
        sum_hydrophilic = sum(amino_acid_counts[aa] * prop for aa, prop in property_dict.items() if aa in 'DEGHKMNPQRST')
        return - (sum_hydrophobic / sum_hydrophilic) if sum_hydrophilic else 0

    def _calculate_net_charge(self, sequence, pKa_dict, pH):
        net_charge = sum((1 / (1 + 10 ** (pKa - pH))) * (
            1 if aa in 'DERK' else
            (1 - 1 if aa == 'H' else 1)
        ) for aa, pKa in ((aa, pKa_dict[aa]) for aa in sequence if aa in pKa_dict))
        return net_charge

def extract_sequences_from_fasta(fasta_file):
    with open(fasta_file) as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    return [str(record.seq).replace('-', '') for record in records]

def run_anarci(input_file, output_file):
    anarci_command = ["ANARCI", "-i", input_file, "-o", output_file]
    try:
        subprocess.run(anarci_command, check=True)
        return True
    except subprocess.CalledProcessError:
        st.error("ANARCI命令执行失败，请检查是否已安装ANARCI并配置好环境变量。")
        return False
    except FileNotFoundError:
        st.error("ANARCI可执行文件未找到，请确保ANARCI已安装并在系统PATH中。")
        return False

def analyze_antibody_fv_sequences(fasta_sequences, pH):
    aac = AminoAcidProperties()
    fv_sequences = ''.join(fasta_sequences)
    
    phi = aac.calculate_hydrophobicity_index(fv_sequences)
    net_charge = aac.calculate_net_charge(fv_sequences, pH)
    
    # 假设我们有从fasta序列中分离出的VH和VL序列
    vh_sequence, vl_sequence = fasta_sequences
    vh_net_charge = aac.calculate_net_charge(vh_sequence, pH)
    vl_net_charge = aac.calculate_net_charge(vl_sequence, pH)
    fv_csp = vh_net_charge * vl_net_charge

    return phi, net_charge, vh_net_charge, vl_net_charge, fv_csp

def display_analysis_results(results, pH):
    phi, net_charge, vh_net_charge, vl_net_charge, fv_csp = results

    st.subheader("计算参数")
    params_data = [
        {"参数": "疏水性指数 (phi)", "值": phi},
        {"参数": f"在pH {pH} 下的净电荷", "值": net_charge},
        {"参数": "VH区域的净电荷", "值": vh_net_charge},
        {"参数": "VL区域的净电荷", "值": vl_net_charge},
        {"参数": "FvCSP", "值": fv_csp}
    ]
    st.table(pd.DataFrame(params_data))

def load_user_fasta():
    uploaded_file = st.file_uploader("请选择您的FASTA文件", type=["fasta"])
    if uploaded_file is not None:
        fasta_content = uploaded_file.getvalue().decode("utf-8")
        temp_file = "temp.fasta"
        with open(temp_file, "w") as f:
            f.write(fasta_content)
        return temp_file
    else:
        return None

if __name__ == "__main__":
    st.set_page_config(layout="wide", page_title="Antibody Fv Sequences Analysis Tool")

    st.markdown("""
    # 抗体 Fv 序列分析工具
   粘度是影响抗体药物开发的重要因素，临床上抗体往往需要静脉内或皮下给药，需要高浓度的抗体溶液（>100mg/mL）才能以小剂量注射获得与治疗相关的剂量，但是高浓度的抗体往往表现出高粘度，这对抗体药物的开发，制造和给药提出了挑   战。研究发现，抗体序列是决定抗体粘度的关键因素，文献报道抗体粘度与 Fv 区域的电荷、VH 和 VL 区域电荷的不对称性 FvCSP 和 Fv 区域的疏水指数HI存在相关性，基于抗体序列预测抗体粘度是一个有效方法。

    此工具参考了以下研究：
    
    Data obtained from the Sharma Vikas K，[In silico selection of therapeutic antibodies for development: viscosity, clearance, and chemical stability](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4284567/). ***PNAS*** 111,52 (2014): 18601-6.

    **应用目的：**
    本工具旨在帮助制剂科学家快速评估抗体的粘度，从而更好地进行早期开发。
    """)

    uploaded_fasta_file = load_user_fasta()

    if uploaded_fasta_file:
        if run_anarci(uploaded_fasta_file, 'output.fasta'):
            fasta_sequences = extract_sequences_from_fasta('output.fasta')
            if fasta_sequences:
                pH = st.sidebar.slider("选择pH值", min_value=0.0, max_value=14.0, value=7.4, step=0.1)
                results = analyze_antibody_fv_sequences(fasta_sequences, pH)
                display_analysis_results(results, pH)
            else:
                st.warning("未能成功从ANARCI输出中提取序列。")
        else:
            st.warning("未能成功运行ANARCI，请确认上传的文件是有效的FASTA格式。")
    else:
        st.info("请上传您的FASTA文件以进行分析。")

    # 清理临时文件
    if uploaded_fasta_file and os.path.exists(uploaded_fasta_file):
        os.remove(uploaded_fasta_file)

    st.write("\n**注意：** 本工具所作的计算基于经验尺度和简化的假设，如需更精确预测，请使用专门的生物信息学工具，并查阅相关原始研究论文。")