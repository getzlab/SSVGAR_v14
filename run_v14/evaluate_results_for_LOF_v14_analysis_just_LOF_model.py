import pandas as pd
import numpy as np
import random
import argparse
from statsmodels.stats import multitest
# from find_p_and_q_values_of_rna_seq_expression_from_wilcoxon_rank_sum_test_v14 import produce_p_values_for_genes_enriched_and_depleted
from SV_analysis_loss_of_function_events_v14 import load_and_prepare_R_mat_file_into_df
from SV_analysis_loss_of_function_events_v14 import convert_chr_xavi


# # I don't think we need to anonymize gene names anymore:
# def anonymize_genes(gene_list):
#     the_gene_list = sorted(gene_list)
#     random.seed(100) # Chose random seed number arbitrarily/randomly; it doesn't matter what it is as long as it remains the same for our analysis
#     random.shuffle(the_gene_list)
#     anonymized_ids = list(map(lambda x: 'anonymized_gene_' + str(x), list(range(0,len(the_gene_list)))))
#     genes_to_anonymized_ids = {}
#     for i in range(len(the_gene_list)):
#         genes_to_anonymized_ids[the_gene_list[i]] = anonymized_ids[i]
#     return genes_to_anonymized_ids

def produce_final_hit_list_table(cohort, distance_threshold=1000000, min_num_patients_threshold_for_FDR=5):
    
    # Load in model results based off cohort and distance threshold parameter specified (default distance threshold parameter is 1000000bp):
    LOF_v14_analysis_df = pd.read_csv("../model_results/" + cohort + "/LOF_model_results_for_all_genes_that_have_SVs_across_at_least_two_patients_threshold_of_proximity_" + str(distance_threshold) + "bp_v14.tsv",sep='\t')
    
    R_mat_file_path='../reference_files/annotation_file/hg38_R_with_hugo_symbols_with_DUX4L1_HMGN2P46_MALAT1.mat'
    R_mat_loaded_and_prepared_df = load_and_prepare_R_mat_file_into_df(R_mat_file_path)

    genes_to_cytobands_df = pd.read_csv("../reference_files/genes_to_cytobands/genes_to_cytobands.tsv",sep='\t')
    genes_to_cytobands_dict = {}
    for i,r in genes_to_cytobands_df.iterrows():
        genes_to_cytobands_dict[r['gene']] = r['cytoband']

    # # I don't think we need to anonymize gene names anymore:
    # gene_list = sorted(list(set(R_mat_loaded_and_prepared_df['gene'])))
    # genes_to_anonymized_ids = anonymize_genes(gene_list)
    
    CGC_df = pd.read_csv("../reference_files/CGC_gene_list/Census_allThu_Jun_20_12_28_25_2024_hg38.tsv",sep='\t')
    list_of_known_TSG_genes_from_CGC = list(CGC_df[~CGC_df['Role in Cancer'].isna()][CGC_df[~CGC_df['Role in Cancer'].isna()]['Role in Cancer'].str.contains('TSG')]['Gene Symbol'])
    list_of_known_fusion_genes_from_CGC = list(CGC_df[~CGC_df['Role in Cancer'].isna()][CGC_df[~CGC_df['Role in Cancer'].isna()]['Role in Cancer'].str.contains('fusion')]['Gene Symbol'])
    
    filtered_LOF_v14_analysis_df = LOF_v14_analysis_df[LOF_v14_analysis_df['num_total_events_used_for_test'] > 0].copy()
    
    filtered_LOF_v14_analysis_df['gene_cytoband'] = list(map(lambda g: genes_to_cytobands_dict[g], list(filtered_LOF_v14_analysis_df['gene'])))
    filtered_LOF_v14_analysis_df['is_TSG_gene_in_CGC'] = list(map(lambda g: "yes" if g in list_of_known_TSG_genes_from_CGC else "no", list(filtered_LOF_v14_analysis_df['gene'])))
    filtered_LOF_v14_analysis_df['is_fusion_gene_in_CGC'] = list(map(lambda g: "yes" if g in list_of_known_fusion_genes_from_CGC else "no", list(filtered_LOF_v14_analysis_df['gene'])))
    filtered_LOF_v14_analysis_df = filtered_LOF_v14_analysis_df[filtered_LOF_v14_analysis_df['num_patients_used_for_test'] >= min_num_patients_threshold_for_FDR].copy()
    filtered_LOF_v14_analysis_df['q_value_LOF'] = multitest.multipletests(list(filtered_LOF_v14_analysis_df['p_value_LOF']), method = "fdr_bh")[1]
    filtered_LOF_v14_analysis_df = filtered_LOF_v14_analysis_df[['gene','is_TSG_gene_in_CGC','is_fusion_gene_in_CGC','num_patients_used_for_test','num_patients_with_at_least_one_LOF_event',
                                                                'num_total_events_used_for_test','num_intragenic_events','num_intergenic_events','num_LOF_events','prob_of_LOF','p_value_LOF','q_value_LOF','gene_cytoband']]
    
    filtered_LOF_v14_analysis_df.sort_values(by=['p_value_LOF'],inplace=True)
    FDR_corrected_filtered_LOF_v14_analysis_df_path = "../model_results/" + cohort + "/LOF_model_results_for_all_genes_that_have_SVs_across_at_least_two_patients_FDR_correction_applied_min_num_patients_" + str(min_num_patients_threshold_for_FDR) + "_threshold_of_proximity_" + str(distance_threshold) + "bp_v14.tsv"
    filtered_LOF_v14_analysis_df.to_csv(FDR_corrected_filtered_LOF_v14_analysis_df_path,sep='\t',index=False)
    final_hist_list_df = filtered_LOF_v14_analysis_df[filtered_LOF_v14_analysis_df['q_value_LOF'] <= 0.25].sort_values(by=['q_value_LOF'])
    print("Final FDR-corrected hit list (windowing parameter = " + str(distance_threshold) + "; minimum number of patients with corresponding breakpoint necessary to qualify for FDR = " + str(min_num_patients_threshold_for_FDR) + ") with q-value ≤ 0.25:")
    print(final_hist_list_df.to_string(index=False))
    FDR_corrected_filtered_LOF_v14_analysis_final_hit_list_df_path = "../model_results/" + cohort + "/LOF_model_results_filtered_via_FDR_correction_applied_min_num_patients_" + str(min_num_patients_threshold_for_FDR) + "_threshold_of_proximity_" + str(distance_threshold) + "bp_v14.tsv"
    final_hist_list_df.to_csv(FDR_corrected_filtered_LOF_v14_analysis_final_hit_list_df_path,sep='\t',index=False)
    print('\n' + 'FDR-corrected gene list table saved to: "' + FDR_corrected_filtered_LOF_v14_analysis_df_path + '"')
    print('FDR-corrected gene list table filtered by q-value ≤ 0.25 saved to: "' + FDR_corrected_filtered_LOF_v14_analysis_final_hit_list_df_path + '"')
    

if __name__ == '__main__':
    ### NEED TO SPECIFY PARAMETERS FOR WHICH YOU ARE RUNNING VIA THE COMMAND LINE
    parser = argparse.ArgumentParser()
    parser.add_argument('--cohort', type=str, required=True, help="Project/cohort name")
    parser.add_argument('--distance_threshold', type=int, default=1000000, help="Windowing parameter in bp from gene edges used for SSVGAR (default: 1000000 for 1Mbp)")
    parser.add_argument('--min_num_patients_with_SV_for_gene_cutoff_for_FDR', type=int, default=5, help="Minimum number of patients which must have an SV corresponding to a gene for that gene to be considered for the FDR correction")
    args = parser.parse_args()
    produce_final_hit_list_table(args.cohort, args.distance_threshold, args.min_num_patients_with_SV_for_gene_cutoff_for_FDR)
