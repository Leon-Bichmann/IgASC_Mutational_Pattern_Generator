all_aas=["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]

#helper functions
def locate_gaps(x):
    locs=[i for i,c in enumerate(x) if c=="."]
    return locs

def transfer_gaps(template,gaps):
    it=template
    for i in gaps:
        it=it[:i]+"."+it[i:]
    return it

def detect_mut_per_position(template,germline):
    muts=[germline[i]+">"+template[i] if template[i]!=germline[i] else "." for i in range(0,len(germline))]
    return muts

def tokenize_and_translate(template):
    template=list(template)
    grouped_words = [str(Seq(''.join(template[i: i + 3])).translate()) if "." not in template[i:i+3] else "..." for i in range(0, len(template), 3)]
    return grouped_words

def get_mutation_AA(heavy_sequence_mut_pos_AA):
    mut_AA_list=[]
    for muts in heavy_sequence_mut_pos_AA:
        mut_AA_list.append(muts.strip("[").strip("]").split(","))
    return mut_AA_list

def create_mut_matrix(mut_AA_list,germline=False):
    i = 0 if germline else 1
    mut_AA_matrix=[]
    for m in mut_AA_list:
        L=[]
        for n in m:
            if ">" in n:
                L.append(n.strip(" ").strip("'").split(">")[i])
            else:
                L.append(".")
        mut_AA_matrix.append(L)
    return mut_AA_matrix

def create_pssm(mut_AA_matrix, ASCindex):
    pssm={}
    max_position=105
    for i in range(0,max_position):
        pssm[i]=[]
        for aa in all_aas:
            freqs=mut_AA_matrix[ASCindex][i].value_counts()
            if aa in freqs.index.tolist():
                count_aa=mut_AA_matrix[ASCindex][i].value_counts()[aa]
                count_all=mut_AA_matrix[ASCindex][i].value_counts()["."]
                perc=count_aa/count_all
            else:
                perc=0
            pssm[i].append(perc)
    pssm_df=pd.DataFrame(pssm).T
    pssm_df.columns=all_aas
    return pssm_df