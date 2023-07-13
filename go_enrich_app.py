import os
import streamlit as st
import pandas as pd
import numpy as np
from scipy import sparse, stats
from statsmodels.stats import multitest


data_path = os.path.join(os.path.dirname(__file__), 'data')

@st.cache_data
def get_file_with_cache(filename):
    df = pd.read_csv(os.path.join(data_path, filename))
    return df


bp_df = get_file_with_cache('bp_go_annotation_df.csv.gz')
cc_df = get_file_with_cache('cc_go_annotation_df.csv.gz')
mf_df = get_file_with_cache('mf_go_annotation_df.csv.gz')
go_dictionary_df = get_file_with_cache('go_dictionary_df.csv')
default_list = get_file_with_cache('default_ranking.txt')


def calculate_auroc(sum_of_ranks: int, n_total: int, n_pos: int, n_neg: int) -> float:
    if n_pos == 0:
      return 0
    else:
      average_rank = sum_of_ranks / n_pos
      min_average_rank = (n_pos + 1)/2
      shifted_average_rank = average_rank - min_average_rank
      auroc = shifted_average_rank / n_neg
      return auroc

def calculate_pval(n_pos: int, n_neg: int, auroc: float) -> float:   
    U = auroc * n_pos * n_neg
    Z = (np.abs(U -(n_pos * n_neg / 2))) / np.sqrt(n_pos * n_neg *(n_pos + n_neg + 1) / 12)
    p = stats.norm.sf(Z)
    return p

def gene_table(gene_list, go_annotation_df, go_dictionary_df):
    columns_only_in_ranked_list_df = set(gene_list).intersection(go_annotation_df.columns)
    #print(len(columns_only_in_ranked_list_df))
    #filter the ranked list for the intersection, keep order intact
    ranked_list_df = [x for x in gene_list if x in columns_only_in_ranked_list_df]
    ranked_list_df = pd.DataFrame(ranked_list_df)
    ranked_list_df = ranked_list_df.rename({0: 'gene'}, axis=1)
    ranked_list_df['index'] = list(reversed((ranked_list_df.index + 1).tolist()))
    ranked_list_df = ranked_list_df.set_index('gene')
    ranked_list_df = ranked_list_df.T

    #this operates locally so the full annotation matrix will be untouched for the next call
    go_annotation_df = go_annotation_df[list(columns_only_in_ranked_list_df)]

    go_annotation_df = go_annotation_df.reindex(columns=ranked_list_df.columns)
    go_annotation_df = go_annotation_df.T
    analysis_df = ranked_list_df @ go_annotation_df

    analysis_df = analysis_df.rename(index={'index': 'sum_of_ranks'})

    analysis_df = analysis_df.T
    
    columns = ['n_total', 'n_pos', 'n_neg', 'auroc', 'pval', 'name']
    for column in columns: analysis_df[column] = np.nan

    n_total = go_annotation_df.shape[0]
    for index, row in analysis_df.iterrows():
      n_pos = sum(go_annotation_df[index])
      n_neg = n_total - n_pos
      analysis_df.loc[index, 'n_total'], analysis_df.loc[index, 'n_pos'], analysis_df.loc[index, 'n_neg'] = n_total, n_pos, n_neg
      analysis_df.loc[index, 'auroc'] = calculate_auroc(analysis_df.loc[index, 'sum_of_ranks'], n_total, n_pos, n_neg)
      analysis_df.loc[index, 'pval'] = calculate_pval(n_pos, n_neg, analysis_df.loc[index, 'auroc'])
      analysis_df.loc[index, 'name'] = go_dictionary_df.loc[index, 'name']

    sorted_analysis_df = analysis_df
    sorted_analysis_df = sorted_analysis_df.drop(['sum_of_ranks', 'n_total', 'n_neg'], axis=1)


    sorted_analysis_df = sorted_analysis_df.reset_index()
    columns = sorted_analysis_df.columns.tolist()
    columns = [columns[-1]] + columns[:-1]
    sorted_analysis_df = sorted_analysis_df[columns]

    #sorted_analysis_df.reset_index(inplace=True, drop=True)
    #sorted_analysis_df.index = [""] * len(sorted_analysis_df)

#
    sorted_analysis_df['signed_log_p'] = np.sign(sorted_analysis_df['auroc']-0.5) * -1 * np.log(sorted_analysis_df['pval'])
    sorted_analysis_df['pvalue_fdr'] = multitest.multipletests(sorted_analysis_df['pval'].tolist(), method="fdr_bh")[1]
    sorted_analysis_df = sorted_analysis_df.sort_values('pval', ascending=True)
    sorted_analysis_df = sorted_analysis_df.rename(columns = {"n_pos" : "genes"})

    # Display the dataframe as a table
    st.write(sorted_analysis_df.style.format({'auroc' : "{:.2f}", 'genes':"{:.0f}", "pval": "{:.2g}", "pvalue_fdr": "{:.2g}"}))

    st.download_button(
        "Download table",
        sorted_analysis_df.to_csv(index=False).encode('utf-8'),
        "AUROC_table_genome_wide_predictions.csv",
        "text/csv",
        key='table-download-csv'
    )



st.title('Simple GO enrichment tester')


gene_list = st.sidebar.text_area("Ranked genes (one gene symbol per line)", '\n'.join(default_list.gene_symbol), height=100)
pathway_database = st.sidebar.selectbox("Pathway Database",["All", "Biological Process", "Cellular Component", "Molecular Function"])

if pathway_database == "Biological Process":
   go_annotation_df = bp_df
elif pathway_database == "Cellular Component":
   go_annotation_df = cc_df
elif pathway_database == "Molecular Function":
   go_annotation_df = mf_df
else:
    go_annotation_df = pd.concat([bp_df, cc_df, mf_df], axis=0)
    
go_annotation_df = go_annotation_df.set_index('go_id').copy()
go_dictionary_df = go_dictionary_df.set_index('go_id').copy()



# Split the input into a list of genes
gene_list = gene_list.split('\n')

# Remove any empty entries from the list
gene_list = [gene.strip() for gene in gene_list if gene.strip()]

# Generate and display the gene table
gene_table(gene_list, go_annotation_df, go_dictionary_df)


