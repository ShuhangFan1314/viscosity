import os
import subprocess
import tempfile
from Bio import SeqIO
import streamlit as st
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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
        sum_hydrophobic = sum(amino_acid_counts[aa] * property_dict[aa] for aa in 'ACFILPVWY')
        sum_hydrophilic = sum(amino_acid_counts[aa] * property_dict[aa] for aa in 'DEGHKMNPQRST')
        hydrophobicity_index = - (sum_hydrophobic / sum_hydrophilic)
        return hydrophobicity_index

    def _calculate_net_charge(self, sequence, pKa_dict, pH):
        net_charge = 0
        for aa in sequence:
            if aa in pKa_dict:
                pKa = pKa_dict[aa]
                charge = 1 / (1 + 10**(pKa - pH))
                if aa in 'DERK':
                    net_charge += charge
                elif aa in 'HY':
                    net_charge += charge - 1  # Histidine can donate and accept two protons
        return net_charge

def validate_sequence(sequence):
    # Check if the sequence contains only valid amino acids
    valid_amino_acids = set(AminoAcidProperties().hydrophobicity.keys())
    return all(aa in valid_amino_acids for aa in sequence)

def remove_invalid_characters(sequence):
    valid_amino_acids = set(AminoAcidProperties().hydrophobicity.keys())
    return ''.join(aa for aa in sequence if aa in valid_amino_acids)

def run_anarci_and_extract_fv(input_file="test.fasta", output_file="output.fasta"):
    # 使用ANARCI对输入文件进行处理
    anarci_command = ["ANARCI", "-i", input_file, "-o", output_file]
    subprocess.run(anarci_command, check=True)
    # 在实际脚本中执行此命令（此处仅为演示，不实际执行）
    # os.system(anarci_command)

    # 从ANARCI输出的FASTA文件中提取轻链和重链FV序列
    def extract_fv_from_record(record):
        fv_regions = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3']
        fv_sequence = ''
        for feature in record.features:
            if feature.type == 'chain' and (feature.qualifiers['region'][0] in fv_regions):
                fv_sequence += str(record.seq[feature.location.start:feature.location.end])
        # 移除任何可能出现的连字符（-）
        fv_sequence = fv_sequence.replace('-', '')
        return fv_sequence

    light_chain_fv = None
    heavy_chain_fv = None

    with open(output_file, 'r') as f:
        records = list(SeqIO.parse(f, "fasta"))

    for record in records:
        if record.id.startswith('IGHV'):
            heavy_chain_fv = extract_fv_from_record(record)
        elif record.id.startswith('IGLV'):
            light_chain_fv = extract_fv_from_record(record)

    if light_chain_fv is not None and heavy_chain_fv is not None:
        return light_chain_fv, heavy_chain_fv
    else:
        raise ValueError("无法从ANARCI输出中成功提取轻链和重链的可变区序列")


def extract_sequences_from_input():
    light_chain_str = st.text_input("请输入轻链序列（不含分隔符）")
    heavy_chain_str = st.text_input("请输入重链序列（不含分隔符）")

    if light_chain_str and heavy_chain_str:
        light_chain_str = remove_invalid_characters(light_chain_str)
        heavy_chain_str = remove_invalid_characters(heavy_chain_str)

        if validate_sequence(light_chain_str) and validate_sequence(heavy_chain_str):
            temp_light_chain_file = tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".fasta")
            temp_heavy_chain_file = tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".fasta")

            temp_light_chain_file.write(f">{light_chain_str}\n{light_chain_str}\n")
            temp_heavy_chain_file.write(f">{heavy_chain_str}\n{heavy_chain_str}\n")

            temp_light_chain_file.close()
            temp_heavy_chain_file.close()

            light_chain_fv, heavy_chain_fv = run_anarci_and_extract_fv(temp_light_chain_file.name, temp_heavy_chain_file.name)

            os.remove(temp_light_chain_file.name)
            os.remove(temp_heavy_chain_file.name)

            if light_chain_fv and heavy_chain_fv:
                return light_chain_fv, heavy_chain_fv
            else:
                st.error("无法从输入序列中提取Fv区域。")
                return None, None
        else:
            st.error("请确保轻链和重链序列只包含有效的氨基酸字符。")
            return None, None
    else:
        st.error("请确保轻链和重链序列均已被填写。")
        return None, None

def validate_and_extract_sequences(light_chain_str, heavy_chain_str):
    if light_chain_str and heavy_chain_str:
        light_chain_str = remove_invalid_characters(light_chain_str)
        heavy_chain_str = remove_invalid_characters(heavy_chain_str)
    
        if validate_sequence(light_chain_str) and validate_sequence(heavy_chain_str):
            return SeqRecord(Seq(light_chain_str), id="light_chain"), SeqRecord(Seq(heavy_chain_str), id="heavy_chain")
        else:
            st.error("输入的序列包含无效的氨基酸字符。")
            return None, None
    else:
        st.error("请确保轻链和重链序列均已被填写。")
        return None, None

def analyze_antibody_sequences(light_chain_fv, heavy_chain_fv):
    aac = AminoAcidProperties()
    fv_sequences = str(light_chain_fv) + str(heavy_chain_fv)
    
    phi = aac.calculate_hydrophobicity_index(fv_sequences)
    net_charge = aac.calculate_net_charge(fv_sequences, 7.4)
    
    # Calculate the net charge for VH and VL
    vh_net_charge = aac.calculate_net_charge(str(heavy_chain_fv), 7.4)
    vl_net_charge = aac.calculate_net_charge(str(light_chain_fv), 7.4)
    
    # Calculate FvCSP
    fv_csp = vh_net_charge * vl_net_charge
    
    # Calculate viscosity using the formula from the reference
    viscosity = 10**(0.15 + 1.26 * phi - 0.043 * net_charge - 0.020 * fv_csp)

    return viscosity

def display_final_result(viscosity):
    st.subheader('Viscosity')
    st.markdown(f"**η,cP(180 mg/mL, 25°C):** {viscosity:.2f}")

def main():
    st.set_page_config(layout="wide", page_title="抗体Fv序列分析工具")

    st.markdown("""
     # 抗体 Fv 序列分析工具
    粘度是影响抗体药物开发的重要因素，临床上抗体往往需要静脉内或皮下给药，需要高浓度的抗体溶液（>100mg/mL）才能以小剂量注射获得与治疗相关的剂量，但是高浓度的抗体往往表现出高粘度，这对抗体药物的开发，制造和给药提出了挑战。研究发现，抗体粘度与 Fv 区域的电荷、VH 和 VL 区域电荷的不对称性 FvCSP 和 Fv 区域的疏水指数HI存在相关性，基于抗体序列预测抗体粘度是一个有效方法。

    此工具参考了以下研究：
    
    Data obtained from the Sharma Vikas K，[In silico selection of therapeutic antibodies for development: viscosity, clearance, and chemical stability](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4284567/). ***PNAS*** 111,52 (2014): 18601-6.

    **应用目的：**
    本工具旨在帮助制剂科学家快速评估抗体的粘度，只需直接输入轻链和重链序列即可进行分析。
    """)

    light_chain_fv, heavy_chain_fv = extract_sequences_from_input()

    if light_chain_fv and heavy_chain_fv:
        calculate_button = st.button("开始计算")
        if calculate_button:
            viscosity = analyze_antibody_sequences(light_chain_fv, heavy_chain_fv)
            display_final_result(viscosity)
    else:
        st.info("请分别输入轻链和重链序列以进行分析。")
   
    st.write("\n**注意：** 本工具所作的计算基于经验尺度和简化的假设，如需更精确预测，请使用专门的生物信息学工具，并查阅相关原始研究论文。")

if __name__ == "__main__":
    main()