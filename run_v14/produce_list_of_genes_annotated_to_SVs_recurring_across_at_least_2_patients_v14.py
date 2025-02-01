import pandas as pd
import os

def produce_list_of_genes_annotated_to_SVs_recurring_across_at_least_2_patients_v14(project, path_to_SV_list_file):

    directory_to_save_lists_of_genes_annotated_to_SVs_recurrent_across_at_least_2_patients = "../lists_of_genes_annotated_to_SVs_recurring_across_at_least_2_patients_for_cohorts/v14/" + project
    if not os.path.isdir(directory_to_save_lists_of_genes_annotated_to_SVs_recurrent_across_at_least_2_patients):
        os.makedirs(directory_to_save_lists_of_genes_annotated_to_SVs_recurrent_across_at_least_2_patients)

    filename = directory_to_save_lists_of_genes_annotated_to_SVs_recurrent_across_at_least_2_patients + "/list_of_genes_annotated_to_SVs_recurring_across_at_least_2_patients." + project + ".txt"
    SV_df = pd.read_csv(path_to_SV_list_file,sep='\t')
    
    unique_genes = list(set(list(SV_df['gene1']) + list(SV_df['gene2'])))
    counts_of_patients_corresponding_to_unique_genes = list(map(lambda x: len(set(SV_df[(SV_df['gene1'] == x) | (SV_df['gene2'] == x)]['individual'])), unique_genes))
    genes_to_counts_of_patients = {}
    for i in range(len(unique_genes)):
        genes_to_counts_of_patients[unique_genes[i]] = counts_of_patients_corresponding_to_unique_genes[i]
    ordered_genes_that_recur_across_at_least_2_patients = list(map(lambda z: z[0], list(filter(lambda y: y[1] > 1, sorted(genes_to_counts_of_patients.items(),key = lambda x: x[1],reverse=True)))))

    # Produce .txt file:
    file_to_create = open(filename, "w")
    for gene in ordered_genes_that_recur_across_at_least_2_patients:
        file_to_create.write(gene + "\n")
    file_to_create.close()

