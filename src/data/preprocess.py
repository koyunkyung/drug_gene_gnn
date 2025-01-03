import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from cmapPy.pandasGEXpress.parse import parse

compound_info = pd.read_csv("~/projects/drug_gene_gnn/data/raw/compound_info.txt", sep="\t")
gene_info = pd.read_csv("~/projects/drug_gene_gnn/data/raw/gene_info.txt", sep="\t")
gctx_data = parse("~/projects/drug_gene_gnn/data/raw/level5_beta_trt_cp.gctx")

## 'compound_info '(약물 정보) 전처리 ##
filtered_compounds = compound_info[['pert_id', 'canonical_smiles', 'moa']].dropna()
# SMILES 벡터화 함수
def smiles_to fingerprint(smiles, n_bits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(n_bits)
    return np.array(AllChem.GetMorganFingerpringAsBitVect(mol, 2, nBits=n_bits))
# 각 약물('pert_id')에 대한 2048 차원의 분자 특징 벡터 생성
filtered_compounds['features'] = filtered_compounds['canonical_smiles'].apply(smiles_to_fingerprint)


## 'gene_info'(유전자 정보) 전처리 ##
filtered_genes = gene_info[['gene_id', 'gene_symbol', 'gene_type']].dropna()


## 'trt_cp.gctx'(약물-유전자 관계 정보) 전처리 ##
expression_data = gctx_data.data_df
edges = []
for col in expression_data.columns:
    pert_id = col.split(":")[0]  # 약물 ID 추출
    top_genes = expression_data[col].nlargest(10).index  # 상위 10개 유전자
    for gene_id in top_genes:
        edges.append((pert_id, gene_id))


## 데이터프레임 생성 후 데이터 저장 ##
# 엣지 데이터프레임
edges_df = pd.DataFrame(edges, columns=["drug", "gene"])
edges_df.to_csv("~/projects/drug_gene_gnn/data/processed/edges.csv", index=False)

# 약물 노트 특징 데이터프레임
drug_feature = pd.DataFrame(
    list(filtered_compounds['features']),
    index=filtered_compounds['pert_id'],
    columns=[f"feature_{i}" for i in range(filtered_compounds['features'][0].shape[0])]
)
drug_feaures.reset_index(inplace=True)
drug_features.rename(columns={'index': 'node_id'}, inplace=True)
drug_features.to_csv("~/projects/drug_gene_gnn/data/processed/node_features.csv", index=False)

# 유전자 노드 특징 데이터프레임
gene_features = pd.get_dummies(filtered_genes['gene_type'])
gene_features.insert(0, 'node_id', filtered_genes['gene_id'])
gene_features.to_csv("~/projects/drug_gene_gnn/data/processed/node_features.csv", mode='a', header=False, index=False)

# 노드 메타데이터
metadata = pd.concat([
    filtered_compounds[['pert_id', 'moa']].rename(columns={'pert_id': 'node_id'}),
    filtered_genes[['gene_id', 'gene_symbol']].rename(columns={'gene_id': 'node_id'})
], ignore_index=True)
metadata.to_csv("~/projects/drug_gene_gnn/data/processed/metadata.csv", index=False)