import glob
import pyreadr
import pandas as pd
from Bio.Seq import Seq
from utils.functions import *
import pickle

#input
files=glob.glob("data/*.rds")

df_list=[]
for f in files:
    result = pyreadr.read_r(f)[None]
    df_list.append(result)
df=pd.concat(df_list, ignore_index=True)

#load ASC map
asc_map=pd.read_csv("assets/ASC_mapping_table.csv")
asc_map["ASC_G"]=asc_map[~((asc_map["ASC"].str.contains("_"))|(asc_map["Allele"].str.contains("_")))]["ASC"].str.split("*").str.get(0).str.split("-").str.get(1)
asc_map=asc_map.dropna()
asc_map=asc_map[["ASC_G","Allele"]].set_index("Allele").to_dict()["ASC_G"]

#compute heavy chain mutation positions
df["gaps"]=df["heavy_sequence"].apply(lambda x: locate_gaps(x))
print("gaps computed")
df["heavy_germline_w_gaps"]=df.apply(lambda row: transfer_gaps(row['heavy_germline'],row['gaps']), axis=1)
print("gaps transfered to germline")
df["heavy_parent_eq_germline"]=df["heavy_parent"]==df["heavy_germline_w_gaps"]

df["heavy_sequence_mut_pos"]=df.apply(lambda row: detect_mut_per_position(row['heavy_sequence'],row['heavy_germline']), axis=1)
print("detected mutation NUC positions")

df["heavy_germline_w_gaps_AA"]=df["heavy_germline_w_gaps"].apply(lambda x: tokenize_and_translate(x))
print("translated heavy germline")

df["heavy_sequence_w_gaps_AA"]=df["heavy_sequence"].apply(lambda x: tokenize_and_translate(x))
print("translated heavy sequence")

df["heavy_sequence_mut_pos_AA"]=df.apply(lambda row: detect_mut_per_position(row['heavy_sequence_w_gaps_AA'],row['heavy_germline_w_gaps_AA']), axis=1)
print("detected mutation AA positions")

#output mutation positions and translated sequences
df.to_csv("ALL_data_translated_and_mutation_indexed.csv",sep=",")

#map to ASC
df["ASC"]=df["v_call_heavy"].map(asc_map)
df=df[~df["ASC"].isna()]

#create ASC index dictionary
ASCindex_dict={}
for a in set(df["ASC"].values.tolist()):
    ASCindex=list(df["ASC"]==a)
    ASCindex_dict[a]=ASCindex

#retrieve mutations
mut_AA_list=get_mutation_AA(df["heavy_sequence_mut_pos_AA"].values.tolist())

#create mutation matrix
mut_AA_matrix=pd.DataFrame(create_mut_matrix(mut_AA_list, germline=False))

#create PSSMs
PSSM_dict={}
for a in ASCindex_dict.keys():
    ASCindex=ASCindex_dict[a]
    pssm_df=create_pssm(mut_AA_matrix, ASCindex).T
    PSSM_dict[a]=pssm_df

#output PSSMs
with open("HeavyChain_ASC_PSSMs.pkl", "wb") as f:
    pickle.dump(PSSM_dict, f)