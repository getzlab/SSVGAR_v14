import pandas as pd
import os
from qtl import annotation
from find_p_and_q_values_of_rna_seq_expression_from_wilcoxon_rank_sum_test_v7 import produce_p_values_for_genes_enriched_and_depleted


# I set the low parameter to what I want it to be:
which_dataset = "PCAWG_lymphoma" # possible datasets right now: only CTSP_DLBCL, because that is the only dataset where we have RNA-seq data available

if which_dataset == "CTSP_DLBCL":
    # SV_list_file = 'All_Pairs.dRanger_etc.filtered_SV.109_pairs.with_tumor_submitter_id.16_low_purity_samples_removed.tsv'
    SV_list_file = '____set_this_to_appropriate_SV_list_path____'
    annot_file = annotation.Annotation("gencode.v36.annotation.gtf")
elif which_dataset == "PCAWG_lymphoma":
    # SV_list_file = 'PCAWG_data/merged_1.6.1.PCAWG_SV_list.subset_lymphoma.lifted_over.annotated.tsv'
    SV_list_file = '____set_this_to_appropriate_SV_list_path____'
    annot_file = annotation.Annotation("gencode.v19.annotation.gtf")
# elif which_dataset == "PCAWG_lymphoid":
#     SV_list_file = 'PCAWG_data/merged_1.6.1.PCAWG_SV_list.subset_lymphoid.annotated.tsv'
else:
    raise("Unknown dataset")


SV_df = pd.read_csv(SV_list_file,sep='\t')
the_enriched_rna_df, the_depleted_rna_df = produce_p_values_for_genes_enriched_and_depleted(annot_file, SV_list_file, which_dataset)

directory_to_save_file = "../rna_wilcoxon_tsvs/" + which_dataset
the_enriched_rna_df.to_csv(directory_to_save_file + "/enriched_rna_wilcoxon_tsv_v14.tsv",sep='\t',index=False)
the_depleted_rna_df.to_csv(directory_to_save_file + "/depleted_rna_wilcoxon_tsv_v14.tsv",sep='\t',index=False)
