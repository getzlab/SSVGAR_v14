import pandas as pd
from SV_analysis_loss_of_function_events_v14 import produce_annotations_from_r_mat_file
from SV_analysis_loss_of_function_events_v14 import find_surrounding_annotations_for_gene_with_parameters_for_list_of_input_genes_post
from SV_analysis_loss_of_function_events_v14 import load_and_prepare_R_mat_file_into_df


def producing_surrounding_annotations_for_relevant_genes_v14(project):
    all_genes_with_at_least_2_patients_with_corresponding_SVs = []

    file_path_of_candidate_genes = "../lists_of_genes_annotated_to_SVs_recurring_across_at_least_2_patients_for_cohorts/v14/" + project + "/list_of_genes_annotated_to_SVs_recurring_across_at_least_2_patients." + project + ".txt"
    with open(file_path_of_candidate_genes) as list_file:
        for line in list_file:
            all_genes_with_at_least_2_patients_with_corresponding_SVs.append(line.strip())

    # Can run manually run these in parallel in separate tmux sessions (with indexing going up to however many genes there are within all_genes_with_at_least_2_patients_with_corresponding_SVs:
    # find_surrounding_annotations_for_gene_with_parameters_for_list_of_input_genes_pre(all_genes_with_at_least_2_patients_with_corresponding_SVs[0:50],10)
    produce_annotations_from_r_mat_file(selected_genes = all_genes_with_at_least_2_patients_with_corresponding_SVs)

    # And then need to run this after:
    R_mat_file_path = "../reference_files/annotation_file/hg38_R_with_hugo_symbols_with_DUX4L1_HMGN2P46_MALAT1.mat"
    R_mat_loaded_and_prepared_df = load_and_prepare_R_mat_file_into_df(R_mat_file_path)
    dict_of_genes_and_gene_orientations = {}
    dict_of_genes_and_gene_chrs = {}
    dict_of_genes_and_gene_starts = {}
    dict_of_genes_and_gene_ends = {}
    for g in all_genes_with_at_least_2_patients_with_corresponding_SVs:
        # Every gene has only one orientation/strand
        dict_of_genes_and_gene_orientations[g] = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == g]['strand'].iloc[0]
        dict_of_genes_and_gene_chrs[g] = sorted(list(set(R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == g]['chr'])))
        dict_of_genes_and_gene_starts[g] = list(map(lambda c: R_mat_loaded_and_prepared_df[(R_mat_loaded_and_prepared_df['gene'] == g) & (R_mat_loaded_and_prepared_df['chr'] == c)]['gene_start'].min(), dict_of_genes_and_gene_chrs[g]))
        dict_of_genes_and_gene_ends[g] = list(map(lambda c: R_mat_loaded_and_prepared_df[(R_mat_loaded_and_prepared_df['gene'] == g) & (R_mat_loaded_and_prepared_df['chr'] == c)]['gene_end'].max(), dict_of_genes_and_gene_chrs[g]))
    surrounding_window_lengths_one_way_of_gene_for_analysis = [1000000]
    step_size = 1
    for d in surrounding_window_lengths_one_way_of_gene_for_analysis:
        print('d:',d)
        find_surrounding_annotations_for_gene_with_parameters_for_list_of_input_genes_post(all_genes_with_at_least_2_patients_with_corresponding_SVs,step_size,d,dict_of_genes_and_gene_orientations,dict_of_genes_and_gene_chrs,dict_of_genes_and_gene_starts,dict_of_genes_and_gene_ends)


