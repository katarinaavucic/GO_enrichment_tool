import os
import streamlit as st
import pandas as pd
import numpy as np
from scipy import sparse, stats

data_path = os.path.join(os.path.dirname(__file__), 'data')

@st.cache
def get_file_with_cache(filename):
    df = pd.read_csv(os.path.join(data_path, filename))
    return df


raw_go_annotation_df = get_file_with_cache('filtered_go_annotation_df.csv.gz')
go_dictionary_df = get_file_with_cache('go_dictionary_df.csv')
default_list = get_file_with_cache('default_ranking.txt')

go_annotation_df = raw_go_annotation_df.set_index('go_id').copy()
go_dictionary_df = go_dictionary_df.set_index('go_id').copy()

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

    go_annotation_df = go_annotation_df[columns_only_in_ranked_list_df]
    # Create a sample dataframe with gene names
    
    go_annotation_df = go_annotation_df.reindex(columns=ranked_list_df.columns)
    go_annotation_df = go_annotation_df.T
    analysis_df = ranked_list_df @ go_annotation_df
    analysis_df = analysis_df.rename(index={0: 'sum_of_ranks'})

    df = pd.DataFrame(analysis_df.head())

    # Display the dataframe as a table
    st.dataframe(df)

st.title('Simple GO enrichment tester')


gene_list = st.sidebar.text_area("Ranked genes input (gene symbol per line)", '\n'.join(default_list.gene_symbol), height=100)


# Create an input field for the gene list
#gene_list = st.text_area('Enter gene names (one per line)')

# Split the input into a list of genes
gene_list = gene_list.split('\n')

# Remove any empty entries from the list
gene_list = [gene.strip() for gene in gene_list if gene.strip()]

# Generate and display the gene table
gene_table(gene_list, go_annotation_df, go_dictionary_df)

