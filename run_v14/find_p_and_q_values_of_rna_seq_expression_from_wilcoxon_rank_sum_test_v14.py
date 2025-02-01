import pandas as pd
import numpy as np
from qtl import annotation
from scipy.stats import ranksums
from statsmodels.stats import multitest


# The annotation object that we specifically want inputted here is annotation_object = annotation.Annotation("gencode.v36.annotation.gtf") for CTSP_DLBCL project
# The annotation object that we specifically want inputted here is annotation_object = annotation.Annotation("gencode.v19.annotation.gtf") for any PCAWG project (e.g. PCAWG_lymphoma)
# We specifically want to use gencode v36 since that was what was used for the batch-corrected RNA-seq data
def produce_p_values_for_genes_enriched_and_depleted(annotation_object, SV_list_path, which_data):
    annotated_SV_list_df = pd.read_csv(SV_list_path, sep='\t')
    
    name_s = pd.Series([g.name for g in annotation_object.genes if g.type == 'protein_coding']).value_counts()
    exclude_names = set(name_s[name_s == 2].index)
    name2id_dict = {g.name:g.id for g in annotation_object.genes if g.name not in exclude_names}
    id2name_dict = {g.id:g.name for g in annotation_object.genes if g.name not in exclude_names}
    # transcript_id2name_dict = {}
    # for g in annotation_object.genes:
    #     if not g.name in exclude_names:
    #         for t in g.transcripts:
    #             transcript_id2name_dict[t.id] = g.name

    if which_data == 'CTSP_DLBCL': # In this case the SV_list_path would be "All_Pairs.dRanger_etc.filtered_SV.109_pairs.with_tumor_submitter_id.16_low_purity_samples_removed.tsv" (for 3/14/2023)
    
        # We have no batch-corrected RNA-seq for this patient so I am taking them out unfortunately:
        annotated_SV_list_df = annotated_SV_list_df[annotated_SV_list_df['individual'] != '50e2b8b0-306f-4f9e-9888-9c16e15e2caa_b6a46a7d-48c2-450d-b239-6a2c39a974a1']

        rna_seq_data_df = pd.read_csv('Corrected_all_RNA_Seq_space_separated.tsv',' ',index_col=0)
        pairs_to_submitter_ids_df = pd.read_csv('pairs_to_submitter_ids.tsv',sep='\t')

        new_column_names = []
        for col_name in list(rna_seq_data_df.columns):
            if "CTSP" in col_name and col_name.replace(".", "-") in list(pairs_to_submitter_ids_df['tumor_submitter_id']):
                new_column_names.append(pairs_to_submitter_ids_df[pairs_to_submitter_ids_df['tumor_submitter_id'] == col_name.replace(".", "-")].iloc[0]['entity:pair_id'])
            else:
                new_column_names.append(col_name)

        rna_seq_data_df.columns = new_column_names
    
    elif which_data == 'PCAWG_lymphoma': # In this case the SV_list_path would be "PCAWG_data/merged_1.6.1.PCAWG_SV_list.subset_lymphoma.lifted_over.annotated.tsv" (for 3/14/2023)
        
        lymphoma_PCAWG_SV_df = pd.read_csv("PCAWG_data/merged_1.6.1.PCAWG_SV_list.subset_lymphoma.csv")
        # pcawg_rnaseq_exp_fpkm_df = pd.read_csv("PCAWG_data/pcawg.rnaseq.transcript.expr.fpkm.tsv",sep='\t')
        # rnaseq_metadata_aliquot_df = pd.read_csv("PCAWG_data/rnaseq.extended.metadata.aliquot_id.V4.tsv",sep='\t')
        pcawg_sample_sheet_df = pd.read_csv("PCAWG_data/pcawg_sample_sheet.tsv",sep='\t')
        pcawg_histology_df = pd.read_csv("PCAWG_data/pcawg_specimen_histology_August2016_v6.tsv",sep='\t')
        
        pcawg_histology_df_subset_lymphoma = pcawg_histology_df[(pcawg_histology_df['histology_abbreviation'].str.startswith('Lymph')) & (pcawg_histology_df['histology_abbreviation'] != "Lymph-CLL") & (pcawg_histology_df['donor_wgs_included_excluded'] == 'Included')]
        donor_unique_ids_107 = list(set(pcawg_histology_df_subset_lymphoma['donor_unique_id'])) # There should be 107 unique donors for lymphoma cases in PCAWG
        rna_sample_unique_ids_107 = list(map(lambda x: pcawg_sample_sheet_df[(pcawg_sample_sheet_df['donor_unique_id'] == x) & (pcawg_sample_sheet_df['library_strategy'] == 'RNA-Seq')]['aliquot_id'].iloc[0] \
                                                   if pcawg_sample_sheet_df[(pcawg_sample_sheet_df['donor_unique_id'] == x) & (pcawg_sample_sheet_df['library_strategy'] == 'RNA-Seq')].shape[0] == 1 
                                                   else "N/A", donor_unique_ids_107))
        rna_sample_unique_ids_105_to_donor_unique_ids_105 = {} # There should be 105 unique donors with RNA-seq samples for lymphoma cases in PCAWG
        for i in range(len(rna_sample_unique_ids_107)):
            if rna_sample_unique_ids_107[i] != "N/A":
                rna_sample_unique_ids_105_to_donor_unique_ids_105[rna_sample_unique_ids_107[i]] = donor_unique_ids_107[i]
        # rna_sample_unique_ids_105 = list(filter(lambda x: x != "N/A", rna_sample_unique_ids_107))
        # rna_sample_unique_ids_105_to_sids_105 = {}
        # rna_sample_unique_ids_105_to_donor_unique_ids_105 = {}
        # for rid in rna_sample_unique_ids_105:
        #     rna_sample_unique_ids_105_to_sids_105[rid] = rnaseq_metadata_aliquot_df[rnaseq_metadata_aliquot_df['rna_seq_aliquot_id'] == rid]['wgs_aliquot_id'].iloc[0]
        #     print("rid:", rid)
        #     print("rna_sample_unique_ids_105_to_sids_105[rid]:", rna_sample_unique_ids_105_to_sids_105[rid])
        #     rna_sample_unique_ids_105_to_donor_unique_ids_105[rid] = lymphoma_PCAWG_SV_df[lymphoma_PCAWG_SV_df["sid"] == rna_sample_unique_ids_105_to_sids_105[rid]]['donor_unique_id'].iloc[0]
        # donor_unique_ids_with_rna_seq_samples = list(lymphoma_PCAWG_SV_df[lymphoma_PCAWG_SV_df['sid'].isin(rna_sample_unique_ids_105_non_NA)]['donor_unique_id'])
        annotated_SV_list_df = annotated_SV_list_df[annotated_SV_list_df["individual"].isin(list(rna_sample_unique_ids_105_to_donor_unique_ids_105.values()))]
        
        rna_seq_data_df = pd.read_csv("PCAWG_data/tophat_star_fpkm.v2_aliquot_gl.tsv",sep='\t')
        rna_seq_data_df = rna_seq_data_df[["feature"] + list(rna_sample_unique_ids_105_to_donor_unique_ids_105.keys())]
        rna_seq_data_df.columns = ["Gene_ID"] + list(map(lambda x: rna_sample_unique_ids_105_to_donor_unique_ids_105[x], list(rna_seq_data_df.columns)[1:]))
        
        # rna_seq_data_df['ID'] = list(map(lambda x: transcript_id2name_dict[x] if x in transcript_id2name_dict else x, list(rna_seq_data_df['ID']))).copy()
        # rna_seq_data_df.set_index('ID',inplace=True)
        
    else:
        raise("Input for which_data has not been accounted for in this method")
    
    rna_seq_data_df['Gene_ID'] = list(map(lambda x: id2name_dict[x] if x in id2name_dict else x, list(rna_seq_data_df['Gene_ID'])))
    # rna_seq_data_df.rename(columns={'Gene_ID':'ID'},inplace=True)
    rna_seq_data_df.set_index('Gene_ID',inplace=True)
    
    unique_genes_list = list(set(list(annotated_SV_list_df['gene1']) + list(annotated_SV_list_df['gene2'])))
    dict_of_genes_to_corresponding_SV_samples = {}
    for unique_gene in unique_genes_list:
        dict_of_genes_to_corresponding_SV_samples[unique_gene] = []
    for index, row in annotated_SV_list_df.iterrows():
        dict_of_genes_to_corresponding_SV_samples[row['gene1']].append(row['individual'])
        dict_of_genes_to_corresponding_SV_samples[row['gene2']].append(row['individual'])
        dict_of_genes_to_corresponding_SV_samples[row['gene1']] = list(set(dict_of_genes_to_corresponding_SV_samples[row['gene1']]))
        dict_of_genes_to_corresponding_SV_samples[row['gene2']] = list(set(dict_of_genes_to_corresponding_SV_samples[row['gene2']]))

    rna_seq_data_df_columns_list_for_relevant_samples = list(filter(lambda x: x in set(list(annotated_SV_list_df['individual'])), list(rna_seq_data_df.columns)))
    genes_not_in_rna_seq_data = []
    genes_and_p_values = []
    genes_and_p_values_enriched = []
    genes_and_p_values_depleted = []
    for gene in dict_of_genes_to_corresponding_SV_samples:
        if gene in list(rna_seq_data_df.index):
            list_of_rna_seq_values_for_samples_with_gene_SV = list(rna_seq_data_df.loc[gene][dict_of_genes_to_corresponding_SV_samples[gene]])
            list_of_rna_seq_values_for_samples_without_gene_SV = list(rna_seq_data_df.loc[gene][list(set(rna_seq_data_df_columns_list_for_relevant_samples) - set(dict_of_genes_to_corresponding_SV_samples[gene]))])
            statistic_enriched, p_value_enriched = ranksums(list_of_rna_seq_values_for_samples_with_gene_SV, list_of_rna_seq_values_for_samples_without_gene_SV, alternative='greater')
            statistic_depleted, p_value_depleted = ranksums(list_of_rna_seq_values_for_samples_with_gene_SV, list_of_rna_seq_values_for_samples_without_gene_SV, alternative='less')
            genes_and_p_values_enriched.append([gene, float(p_value_enriched)])
            genes_and_p_values_depleted.append([gene, float(p_value_depleted)])
        else:
            genes_not_in_rna_seq_data.append(gene)

    data_genes_and_p_values_enriched = {'genes': np.array(genes_and_p_values_enriched)[:,0], 'p-values': np.array(genes_and_p_values_enriched)[:,1]}
    genes_and_p_values_enriched_df = pd.DataFrame(data=data_genes_and_p_values_enriched)
    genes_and_p_values_enriched_df['p-values'] = pd.to_numeric(genes_and_p_values_enriched_df['p-values'])
    genes_and_p_values_enriched_df.sort_values(by=['p-values'],inplace=True,ascending=True)
    genes_and_p_values_enriched_df['q-values'] = list(multitest.multipletests(genes_and_p_values_enriched_df['p-values'] , method = "fdr_bh")[1])

    data_genes_and_p_values_depleted = {'genes': np.array(genes_and_p_values_depleted)[:,0], 'p-values': np.array(genes_and_p_values_depleted)[:,1]}
    genes_and_p_values_depleted_df = pd.DataFrame(data=data_genes_and_p_values_depleted)
    genes_and_p_values_depleted_df['p-values'] = pd.to_numeric(genes_and_p_values_depleted_df['p-values'])
    genes_and_p_values_depleted_df.sort_values(by=['p-values'],inplace=True,ascending=True)
    genes_and_p_values_depleted_df['q-values'] = list(multitest.multipletests(genes_and_p_values_depleted_df['p-values'] , method = "fdr_bh")[1])

    # Also including the number of SVs and number of breakpoints in the dataframes we return:
    orig_annotated_SV_list_df = pd.read_csv(SV_list_path, sep='\t')
    genes_and_p_values_enriched_df['Number of patients'] = list(map(lambda x: len(set(orig_annotated_SV_list_df[(orig_annotated_SV_list_df['gene1'] == x) | (orig_annotated_SV_list_df['gene2'] == x)]['individual'])), list(genes_and_p_values_enriched_df['genes'])))
    genes_and_p_values_enriched_df['Number of SVs'] = list(map(lambda x: orig_annotated_SV_list_df[(orig_annotated_SV_list_df['gene1'] == x) | (orig_annotated_SV_list_df['gene2'] == x)].shape[0], list(genes_and_p_values_enriched_df['genes'])))
    genes_and_p_values_enriched_df['Number of breakpoints'] = list(map(lambda x: (list(orig_annotated_SV_list_df['gene1']) + list(orig_annotated_SV_list_df['gene2'])).count(x), list(genes_and_p_values_enriched_df['genes'])))
    genes_and_p_values_depleted_df['Number of patients'] = list(map(lambda x: len(set(orig_annotated_SV_list_df[(orig_annotated_SV_list_df['gene1'] == x) | (orig_annotated_SV_list_df['gene2'] == x)]['individual'])), list(genes_and_p_values_depleted_df['genes'])))
    genes_and_p_values_depleted_df['Number of SVs'] = list(map(lambda x: orig_annotated_SV_list_df[(orig_annotated_SV_list_df['gene1'] == x) | (orig_annotated_SV_list_df['gene2'] == x)].shape[0], list(genes_and_p_values_depleted_df['genes'])))
    genes_and_p_values_depleted_df['Number of breakpoints'] = list(map(lambda x: (list(orig_annotated_SV_list_df['gene1']) + list(orig_annotated_SV_list_df['gene2'])).count(x), list(genes_and_p_values_depleted_df['genes'])))

    return genes_and_p_values_enriched_df, genes_and_p_values_depleted_df

