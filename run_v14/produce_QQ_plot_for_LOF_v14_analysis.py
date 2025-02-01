import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
from mpl_toolkits.mplot3d import Axes3D
from statsmodels.stats import multitest
from adjustText import adjust_text
import scipy.stats
from statsmodels.stats import multitest

# This is me modifying Julian's version of his QQ plots function for my analysis:
def custom_QQ(pvalues, pmidvalues = None, labels = None, sig_thresh = 0.1, near_sig_thresh = 0.25, fnum = None, neg_log10_cutoff = None): # neg_log10_cutoff is a new parameter I introduced
    si = np.argsort(-np.log10(pvalues))
    logp = -np.log10(pvalues)[si]
    
    if pmidvalues:
        si_pmid = np.argsort(-np.log10(pmidvalues))
        logp_pmid = -np.log10(pmidvalues)[si_pmid]
        assert len(logp) == len(logp_pmid)
    
    n_p = len(logp)

    # expected quantiles
    x = -np.log10(np.r_[n_p:0:-1]/(n_p + 1))

    # FDR
    _, q, _, _ = multitest.multipletests(pvalues[si], method = "fdr_bh")
    sig_idx = q < sig_thresh
    n_sig_idx = q < near_sig_thresh

    fdr_color_idx = np.c_[sig_idx, n_sig_idx]@np.r_[2, 1]

    #                      nonsig     near sig   <padding>     sig
    fdr_colors = np.array([[0, 0, 1], [0, 1, 1], [-1, -1, -1], [1, 0, 0]])

    #
    # plot
    f = plt.figure(fnum); plt.clf()
    if not pmidvalues:
        plt.scatter(x, logp, c = fdr_colors[fdr_color_idx])
    else:
        plt.scatter(x, logp_pmid, c = fdr_colors[fdr_color_idx])

    # 1:1 line
    plt.plot(plt.xlim(), plt.xlim(), color = 'k', linestyle = ':')
    
    # Modification I made -Xavi:
    if neg_log10_cutoff:
        orig_bottom,orig_top = plt.ylim()
        plt.ylim(orig_bottom,neg_log10_cutoff)

    #
    # labels (if given)
    if labels is not None:
        if len(labels) != len(x):
            raise ValueError("Length of labels must match length of p-value array")
        labels = labels[si]
        label_plt = [plt.text(x, y, l) for x, y, l in zip(x[n_sig_idx], logp[n_sig_idx], labels[n_sig_idx])]
        adjust_text(label_plt, arrowprops = { 'color' : 'k', "arrowstyle" : "-" })

    plt.xlabel("Expected quantile (-log10)")
    plt.ylabel("Observed quantile (-log10)")

    return f

def produce_final_hit_list_table(cohort, distance_threshold=1000000, min_num_patients_for_multiple_test_correction=5, use_p_mid=False, neg_log10_cutoff_for_plot=None):
    LOF_v14_analysis_df = pd.read_csv("../model_results/" + cohort + "/LOF_model_results_for_all_genes_that_have_SVs_across_at_least_two_patients_threshold_of_proximity_" + str(distance_threshold) + "bp_v14.tsv", sep='\t')
    np.random.seed(min_num_patients_for_multiple_test_correction * 10**8 + distance_threshold)
    filtered_LOF_v14_analysis_df = LOF_v14_analysis_df[LOF_v14_analysis_df['num_total_events_used_for_test'] > 0].copy()
    filtered_LOF_v14_analysis_df = filtered_LOF_v14_analysis_df[filtered_LOF_v14_analysis_df['num_patients_used_for_test'] >= min_num_patients_for_multiple_test_correction].copy()
    if not use_p_mid:
        fig_to_save = custom_QQ(np.array(list(filtered_LOF_v14_analysis_df['p_value_LOF'])),neg_log10_cutoff=neg_log10_cutoff_for_plot)
        additional_clause_indicating_p_mid = ""
    else:
        list_of_p_mids = []
        for i,row in filtered_LOF_v14_analysis_df.iterrows():
            x = row['num_LOF_events']
            n = row['num_total_events_used_for_test']
            p = row['prob_of_LOF']
            if x == n:
                list_of_p_mids.append(row['p_value_LOF'] / 2)
            else:
                pval = scipy.stats.binom_test(x, n=n, p=p, alternative='greater')
                low_pval = scipy.stats.binom_test(x+1, n=n, p=p, alternative='greater')
                list_of_p_mids.append(np.random.uniform(low_pval, pval))
        filtered_LOF_v14_analysis_df['p_value_mid_LOF'] = list_of_p_mids
        fig_to_save = custom_QQ(np.array(list(filtered_LOF_v14_analysis_df['p_value_mid_LOF'])),neg_log10_cutoff=neg_log10_cutoff_for_plot)
        additional_clause_indicating_p_mid = ".pmid"
    filtered_LOF_v14_analysis_df['q_value_LOF'] = multitest.multipletests(list(filtered_LOF_v14_analysis_df['p_value_LOF']), method = "fdr_bh")[1]
    print("25 most significant hits in QQ plot (windowing parameter = " + str(distance_threshold) + "; minimum number of patients with corresponding breakpoint necessary to qualify for FDR = " + str(min_num_patients_for_multiple_test_correction) + "):")
    print("(May be less than 25 if less than 25 genes qualified to be FDR-corrected overall)")
    pd.set_option('display.max_rows', 25)
    print(filtered_LOF_v14_analysis_df.sort_values(by=['p_value_LOF']).head(25).to_string(index=False))
    if not os.path.exists('../model_results/' + cohort + '/QQ_plots'):
        os.mkdir('../model_results/' + cohort + '/QQ_plots')
    if neg_log10_cutoff_for_plot:
        additional_clause_indicating_neg_log10_cutoff_for_plot = '.cutoff_observed_quantile_at_1e-' + str(neg_log10_cutoff_for_plot)
    else:
        additional_clause_indicating_neg_log10_cutoff_for_plot = ''
    path_of_saved_QQ_plot = '../model_results/' + cohort + '/QQ_plots/QQ_plot_LOF_model_' + cohort + '_thresh_of_proximity_min_num_patients_' + str(min_num_patients_for_multiple_test_correction) + '_threshold_of_proximity_' + str(distance_threshold) + 'bp_' + 'v14' + additional_clause_indicating_neg_log10_cutoff_for_plot + additional_clause_indicating_p_mid + '.png'
    print('\n' + 'Created QQ plot saved to: "' + path_of_saved_QQ_plot + '"')
    fig_to_save.savefig(path_of_saved_QQ_plot)
    
if __name__ == '__main__':
    ### NEED TO SPECIFY PARAMETERS FOR WHICH YOU ARE RUNNING VIA THE COMMAND LINE
    parser = argparse.ArgumentParser()
    parser.add_argument('--cohort', type=str, required=True, help="Project/cohort name")
    parser.add_argument('--distance_threshold', type=int, default=1000000, help="Windowing parameter in bp from gene edges used for SSVGAR (default: 1000000 for 1Mbp)")
    parser.add_argument('--min_num_patients_with_SV_for_gene_cutoff_for_FDR', type=int, default=5, help="Minimum number of patients which must have an SV corresponding to a gene for that gene to be considered for the FDR correction")
    parser.add_argument('--use_p_mid', action='store_true', help="Use p-mid correction for producing QQ plots")
    parser.add_argument('--neg_log10_cutoff_for_plot', type=float, default=None, help="Specify a specific negative-log10 cutoff (such as 1, 2, 3, 4, 2.5, 3.6, etc.) the QQ plot's y-axis if you so desire,")
    args = parser.parse_args()
    produce_final_hit_list_table(args.cohort, args.distance_threshold, args.min_num_patients_with_SV_for_gene_cutoff_for_FDR, args.use_p_mid, args.neg_log10_cutoff_for_plot)

