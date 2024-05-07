import streamlit as st
import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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
        net_charge = sum((1 / (1 + 10 ** (pKa - 7.4))) * (
            1 if aa in 'DERK' else
            (1 - 1 if aa == 'H' else 1)
        ) for aa, pKa in ((aa, pKa_dict[aa]) for aa in sequence if aa in pKa_dict))
        return net_charge

def extract_sequences_from_input():
    light_chain_seq = st.text_input("请输入轻链序列（不含分隔符）")
    heavy_chain_seq = st.text_input("请输入重链序列（不含分隔符）")

    if light_chain_seq and heavy_chain_seq:
        return SeqRecord(Seq(light_chain_seq), id="light_chain"), \
               SeqRecord(Seq(heavy_chain_seq), id="heavy_chain")
    else:
        st.error("请确保轻链和重链序列均已被填写。")
        return None, None

def analyze_antibody_sequences(light_chain, heavy_chain):
    aac = AminoAcidProperties()
    fv_sequences = ''.join([str(light_chain.seq), str(heavy_chain.seq)])
    
    phi = aac.calculate_hydrophobicity_index(fv_sequences)
    net_charge = aac.calculate_net_charge(fv_sequences, 7.4)
    
    vh_sequence = str(heavy_chain.seq)
    vl_sequence = str(light_chain.seq)
    vh_net_charge = aac.calculate_net_charge(vh_sequence, 7.4)
    vl_net_charge = aac.calculate_net_charge(vl_sequence, 7.4)
    fv_csp = vh_net_charge * vl_net_charge
    viscosity = 10**(0.15 + 1.26 * phi - 0.043 * net_charge - 0.020 * fv_csp)

    return viscosity

def display_final_result(viscosity):
    st.subheader('Viscosity')
    st.markdown(f"**η,cP(180 mg/mL, 25°C):** {viscosity:.2f}")

def main():
    st.set_page_config(layout="wide", page_title="抗体Fv序列分析工具")

    st.markdown("""
     # 抗体 Fv 序列分析工具
    粘度是影响抗体药物开发的重要因素，临床上抗体往往需要静脉内或皮下给药，需要高浓度的抗体溶液（>100mg/mL）才能以小剂量注射获得与治疗相关的剂量，但是高浓度的抗体往往表现出高粘度，这对抗体药物的开发，制造和给药提出了挑   战。研究发现，抗体序列是决定抗体粘度的关键因素，文献报道抗体粘度与 Fv 区域的电荷、VH 和 VL 区域电荷的不对称性 FvCSP 和 Fv 区域的疏水指数HI存在相关性，基于抗体序列预测抗体粘度是一个有效方法。

    此工具参考了以下研究：
    
    Data obtained from the Sharma Vikas K，[In silico selection of therapeutic antibodies for development: viscosity, clearance, and chemical stability](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4284567/). ***PNAS*** 111,52 (2014): 18601-6.

    **应用目的：**
    本工具旨在帮助制剂科学家快速评估抗体的粘度，只需直接输入轻链和重链序列即可进行分析。
    """)

    light_chain, heavy_chain = extract_sequences_from_input()

    if light_chain and heavy_chain:
        calculate_button = st.button("开始计算")
        if calculate_button:
            viscosity = analyze_antibody_sequences(light_chain, heavy_chain)
            if viscosity is not None:
                display_final_result(viscosity)
            else:
                st.warning("未能成功分析序列。")
    else:
        st.info("请分别输入轻链和重链序列以进行分析。")
   
    st.write("\n**注意：** 本工具所作的计算基于经验尺度和简化的假设，如需更精确预测，请使用专门的生物信息学工具，并查阅相关原始研究论文。")

if __name__ == "__main__":
    main()