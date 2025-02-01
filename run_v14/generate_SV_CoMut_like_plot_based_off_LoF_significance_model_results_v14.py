import numpy as np
import pandas as pd
import os
import glob
import argparse
from tqdm import tqdm
from statsmodels.stats import multitest
import matplotlib.pyplot as plt
from matplotlib import patches



### THIS PART OF THE CODE IS WHERE WE PRODUCE THE DATAFRAME FOR THE SV COMUT-LIKE PLOT:
def produce_comut_like_plot_dataframe(cohort, surrounding_window, min_patient_thresh, path_to_save_dataframe):
    # Get list of cases/individuals in the cohort:
    cohort_filename_list = glob.glob("../model_results/" + cohort + "/input_SV_list_copy/*")
    if len(cohort_filename_list) != 1:
        raise("There is not only one list in the ../model_results/" + cohort + "/input_SV_list_copy directory for this cohort as there should be")
    else:
        cohort_filename = glob.glob("../model_results/" + cohort + "/input_SV_list_copy/*")[0]
    individuals_in_cohort = list(set(pd.read_csv(cohort_filename,sep='\t')['individual']))
    # Obtain results of SV LoF significance model after multiple test correction:
    results_df = pd.read_csv("../model_results/" + cohort + "/LOF_model_results_filtered_via_FDR_correction_applied_min_num_patients_" + str(min_patient_thresh) + "_threshold_of_proximity_" + str(surrounding_window) + "bp_v14.tsv", sep='\t')
    # Set up dataframe for CoMut plot with empty values but appropriate row and column names:
    dataframe_for_plot = pd.concat([pd.Series(map(lambda x: "", list(range(len(results_df['gene'])))))] * len(individuals_in_cohort),axis=1)
    dataframe_for_plot.set_index(pd.Index(list(results_df['gene'])),inplace=True)
    dataframe_for_plot.columns = individuals_in_cohort
    # Fill in dataframe for CoMut plot as appropriate based off whether there is at least one SV for that patient and gene as well as whether or not at least one is LoF:
    for g in list(results_df['gene']):
        SV_sublist_for_gene_df = pd.read_csv("../SV_list_subsets_of_regions_within_distance_of_gene/v14/" + cohort + "/annotations_relative_to_given_gene_and_boundary_threshold_with_indication_of_LOF_or_not/SV_list_filtered_" + g + "_within_" + str(surrounding_window) + "bp.reannotated_with_LOF_indication.tsv",sep='\t')
        individuals_with_at_least_one_LoF_event_for_gene = list(set(SV_sublist_for_gene_df[SV_sublist_for_gene_df['LOF_or_not'] == 1]['individual']))
        individuals_with_no_LoF_event_for_gene = list(set(SV_sublist_for_gene_df['individual']) - set(individuals_with_at_least_one_LoF_event_for_gene))
        dataframe_for_plot.loc[g, individuals_with_at_least_one_LoF_event_for_gene] = "At least one SV present; at least one LoF"
        dataframe_for_plot.loc[g, individuals_with_no_LoF_event_for_gene] = "At least one SV present; no LoF"
    # Order the patients (x-axis) as appropriate:
    print("Setting up dataframe for CoMut-like plot (setting up axis for each gene)...")
    reordered_individuals = []
    not_yet_reordered_individuals = set(dataframe_for_plot.columns)
    for gene_index in tqdm(range(len(results_df))):
        for rev_index in range(len(results_df),gene_index,-1):
            individuals_to_append = list(filter(lambda i: (dataframe_for_plot.iloc[gene_index:rev_index][i] != "").all(), not_yet_reordered_individuals))
            if not individuals_to_append:
                continue
            else:
                individuals_to_append_reordered = []
                for x in range(gene_index,rev_index):
                    for indiv in individuals_to_append:
                        if (dataframe_for_plot.iloc[x][indiv] == "At least one SV present; no LoF") and (not indiv in individuals_to_append_reordered):
                            individuals_to_append_reordered = [indiv] + individuals_to_append_reordered
                individuals_to_append_reordered = list(set(individuals_to_append) - set(individuals_to_append_reordered)) + individuals_to_append_reordered
                reordered_individuals += individuals_to_append_reordered
                not_yet_reordered_individuals -= set(individuals_to_append_reordered)
    dataframe_for_plot = dataframe_for_plot[reordered_individuals + list(not_yet_reordered_individuals)]
    dataframe_for_plot = dataframe_for_plot.loc[list(dataframe_for_plot.index)[::-1]]
    dataframe_for_plot.to_csv(path_to_save_dataframe,sep='\t')
    print('Created CoMut-like plot dataframe saved to: "' + path_to_save_dataframe + '"')
    return dataframe_for_plot


### THIS PART OF THE CODE IS WHERE WE PRODUCE THE ACTUAL PLOT:
def produce_comut_like_plot_plotting(cohort, dataframe_for_plot,path_to_save_plot):
    labels_to_colors = {"At least one SV present; at least one LoF" : (0.749, 0.047, 0.047),
                        "At least one SV present; no LoF" : (0.031, 0.549, 0.388)}
    figsize = (20,6)
    x_padding = 0.08
    y_padding = 0.08
    width = 1-2*x_padding
    height = 1-2*y_padding
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot()
    print("Producing CoMut-like plot (plotting for each gene)...")
    for g in tqdm(range(len(dataframe_for_plot.index))):
        gene = dataframe_for_plot.index[g]
        for i in range(len(dataframe_for_plot.columns)):
            x_base = i + x_padding
            y_base = g + y_padding
            indiv = dataframe_for_plot.columns[i]
            if dataframe_for_plot.loc[gene,indiv]:
                patch_options = {'facecolor': labels_to_colors[dataframe_for_plot.loc[gene,indiv]]}
            else:
                patch_options = {'facecolor': (1, 1, 1)}
            rect = patches.Rectangle((x_base,y_base),width,height,**patch_options)
            ax.add_patch(rect)
    ax.set_ylim([0, len(dataframe_for_plot.index) + y_padding])
    ax.set_xlim([0, len(dataframe_for_plot.columns) + x_padding])
    ax.set_yticks(np.arange(0.5, len(dataframe_for_plot.index) + 0.5))
    ax.set_yticklabels(dataframe_for_plot.index, style='italic',size=5)
    ax.tick_params(axis='both',which='both',bottom=False,top=False,length=0)
    for loc in ['top', 'right', 'bottom', 'left']:
        ax.spines[loc].set_visible(False)
    # Create legend:
    legend_patches = [patches.Patch(color=color, label=label) for label, color in labels_to_colors.items()]
    ax.legend(handles=legend_patches, loc="upper right", bbox_to_anchor=(0.8, 1), fontsize=8)
    # Label axes:
    ax.set_xlabel("Patients", fontsize=12)
    ax.set_ylabel("Genes", fontsize=12)
    # Give plot a title:
    ax.set_title("Co-occurrences of driver LoF SVs across " + cohort + " cohort", fontsize=18)
    plt.savefig(path_to_save_plot)
    print('Created CoMut-like plot saved to: "' + path_to_save_plot + '"')
    


if __name__ == '__main__':
    ### NEED TO SPECIFY PARAMETERS FOR WHICH YOU ARE RUNNING VIA THE COMMAND LINE
    parser = argparse.ArgumentParser()
    parser.add_argument('--cohort', type=str, required=True, help="Project/cohort name")
    parser.add_argument('--distance_threshold', type=int, default=1000000, help="Windowing parameter in bp from gene edges used for SSVGAR (default: 1000000 for 1Mbp)")
    parser.add_argument('--min_num_patients_with_SV_for_gene_cutoff_for_FDR', type=int, default=5, help="Minimum number of patients which must have an SV corresponding to a gene for that gene to be considered for the FDR correction")
    args = parser.parse_args()
    
    if not os.path.exists('../model_results/' + args.cohort + '/CoMut_like_plots'):
        os.mkdir('../model_results/' + args.cohort + '/CoMut_like_plots')
    
    path_to_save_comut_like_plot = "../model_results/" + args.cohort + "/CoMut_like_plots/CoMut_like_plot_LOF_model_" + args.cohort + "_thresh_of_proximity_min_num_patients_" + str(args.min_num_patients_with_SV_for_gene_cutoff_for_FDR) + "_threshold_of_proximity_" + str(args.distance_threshold) + 'bp_v14.png'
    path_to_save_comut_like_plot_dataframe = "../model_results/" + args.cohort + "/CoMut_like_plots/CoMut_like_plot_LOF_model_" + args.cohort + "_thresh_of_proximity_min_num_patients_" + str(args.min_num_patients_with_SV_for_gene_cutoff_for_FDR) + "_threshold_of_proximity_" + str(args.distance_threshold) + 'bp_v14_dataframe_table.tsv'
    comut_like_dataframe = produce_comut_like_plot_dataframe(args.cohort, args.distance_threshold, args.min_num_patients_with_SV_for_gene_cutoff_for_FDR, path_to_save_comut_like_plot_dataframe)
    produce_comut_like_plot_plotting(args.cohort, comut_like_dataframe,path_to_save_comut_like_plot)
    