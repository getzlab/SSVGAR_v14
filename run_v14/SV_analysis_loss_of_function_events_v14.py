import pandas as pd
import numpy as np
import scipy.io as sio
import os
import mat4py
from scipy import stats
import time
from tqdm import tqdm
import subprocess
from copy import deepcopy

# This is the path to the R_mat file we will always use:
R_MAT_FILE_PATH = "../reference_files/annotation_file/hg38_R_with_hugo_symbols_with_DUX4L1_HMGN2P46_MALAT1.mat"

# # Helper function to flatten a list used in completely_flatten_list
# def flatten_list(l):
#     return [item for sublist in l for item in sublist]

# Helper function to completely flatten a list (i.e. remove all nested lists within a list) used 
# in find_overall_exonic_and_intronic_portions_of_genome_setup to access R.mat file data on exon 
# starts and exon ends
def completely_flatten_list(l):
    if isinstance(l, list):
        if len(l) == 0:
            return []
        elif len(l) == 1:
            return completely_flatten_list(l[0])
        else:
            return completely_flatten_list(l[0]) + completely_flatten_list(l[1:])
    else:
        return [l]

def convert_chr_xavi(chr_val):
    # This is my own chr conversion function that takes in a chr value and converts it to an int appropriately -Xavi
    # I also made this function assuming that a human genome is being used. May have to account for other species later.
    def can_convert_to_int(some_value):
        try:
            int(some_value)
            return True
        except ValueError:
            return False
    if chr_val == 'chrX':
        return 23
    elif chr_val == 'chrY':
        return 24
    elif chr_val == 'chrM':
        return 0
    elif isinstance(chr_val, str):
        if chr_val.startswith('chr'):
            return int(chr_val[3:])
        else:
            if can_convert_to_int(chr_val):
                return int(chr_val)
            else:
                raise("Invalid input value for chr")
    elif isinstance(chr_val, float):
        return int(chr_val)
    elif isinstance(chr_val, int):
        return chr_val

def load_and_prepare_R_mat_file_into_df(R_mat_file_path):
    # Below is the code for loading in the R.mat file
    # The assumption I am making is that there are 14 different columns in the following order:
    # id, transcript, chr, strand, tx_start, tx_end, code_start, code_end, n_exons, exon_starts, exon_ends, gene, exon_frames, version
    R_mat_contents = sio.loadmat(R_mat_file_path)
    data_for_R_mat_dataframe = {}
    data_for_R_mat_dataframe['id'] = np.squeeze(R_mat_contents['R'][0][0][0])
    data_for_R_mat_dataframe['transcript'] = np.array(list(map(lambda x: x[0], np.squeeze(R_mat_contents['R'][0][0][1]))))
    data_for_R_mat_dataframe['chr'] = np.array(list(map(lambda x: x[0], np.squeeze(R_mat_contents['R'][0][0][2]))))
    data_for_R_mat_dataframe['strand'] = np.array(list(map(lambda x: x[0], np.squeeze(R_mat_contents['R'][0][0][3]))))
    data_for_R_mat_dataframe['tx_start'] = np.squeeze(R_mat_contents['R'][0][0][4])
    data_for_R_mat_dataframe['tx_end'] = np.squeeze(R_mat_contents['R'][0][0][5])
    data_for_R_mat_dataframe['code_start'] = np.squeeze(R_mat_contents['R'][0][0][6])
    data_for_R_mat_dataframe['code_end'] = np.squeeze(R_mat_contents['R'][0][0][7])
    data_for_R_mat_dataframe['n_exons'] = np.squeeze(R_mat_contents['R'][0][0][8])
    data_for_R_mat_dataframe['exon_starts'] = np.array(list(map(lambda x: np.squeeze(x,axis=1), np.squeeze(R_mat_contents['R'][0][0][9]))),dtype=object)
    data_for_R_mat_dataframe['exon_ends'] = np.array(list(map(lambda x: np.squeeze(x,axis=1), np.squeeze(R_mat_contents['R'][0][0][10]))),dtype=object)
    data_for_R_mat_dataframe['gene'] = np.array(list(map(lambda x: x[0], np.squeeze(R_mat_contents['R'][0][0][11]))))
    data_for_R_mat_dataframe['exon_frames'] = np.array(list(map(lambda x: np.squeeze(x,axis=1), np.squeeze(R_mat_contents['R'][0][0][12]))),dtype=object)
    data_for_R_mat_dataframe['version'] = np.array(list(map(lambda x: x[0], np.squeeze(R_mat_contents['R'][0][0][13]))))
    # The contents below are part of load_refseq.m
    # I have not (and maybe might not do so anyway) implemented the conditionals for the build for load_refseq
    # The way I am implementing things right now is assuming that the build input is the path to an actual R.mat file, and
    # not a value like 'hg18', 'hg19', 'hg38', etc. But what's below is necessary regardless:
    # %% collapse small nucleolar RNA subtypes
    data_for_R_mat_dataframe['gene'] = np.array(list(map(lambda x: x.split('-')[0] if (x.startswith('SNORD') and '-' in x) else x, data_for_R_mat_dataframe['gene'])))
    # %% add txlen, cdlen,. protlen
    tx_len_values = np.zeros(len(R_mat_contents['R'][0][0][0]))
    code_len_values = np.zeros(len(R_mat_contents['R'][0][0][0]))
    for i in range(0,len(R_mat_contents['R'][0][0][0])):
        for e in range(0,data_for_R_mat_dataframe['n_exons'][i]):
            st = data_for_R_mat_dataframe['exon_starts'][i][e]
            en = data_for_R_mat_dataframe['exon_ends'][i][e]
            tx_len_values[i] += (en-st+1)
            if (en < data_for_R_mat_dataframe['code_start'][i]) | (st > data_for_R_mat_dataframe['code_end'][i]):
                continue
            if st < data_for_R_mat_dataframe['code_start'][i]:
                st = data_for_R_mat_dataframe['code_start'][i]
            if en > data_for_R_mat_dataframe['code_end'][i]:
                en = data_for_R_mat_dataframe['code_end'][i]
            code_len_values[i] += (en-st+1)
    n_codons_values = code_len_values / 3
    data_for_R_mat_dataframe['tx_len'] = tx_len_values
    data_for_R_mat_dataframe['code_len'] = code_len_values
    data_for_R_mat_dataframe['n_codons'] = n_codons_values

    data_for_R_mat_dataframe['chr'] = list(map(convert_chr_xavi, data_for_R_mat_dataframe['chr']))

    data_for_R_mat_dataframe['gene_start'] = data_for_R_mat_dataframe['tx_start'].copy()
    data_for_R_mat_dataframe['gene_end'] = data_for_R_mat_dataframe['tx_end'].copy()
    # if P_parameters_dictionary['impute_promoters']:
    if True:
        imputed_promoter_size = 3000
        idx = np.squeeze(np.argwhere(data_for_R_mat_dataframe['strand']=='+'))
        # data_for_R_mat_dataframe['gene_start'][idx] -= P_parameters_dictionary['imputed_promoter_size']
        data_for_R_mat_dataframe['gene_start'][idx] -= imputed_promoter_size
        idx = np.squeeze(np.argwhere(data_for_R_mat_dataframe['strand']=='-'))
        # data_for_R_mat_dataframe['gene_end'][idx] += P_parameters_dictionary['imputed_promoter_size']
        data_for_R_mat_dataframe['gene_end'][idx] += imputed_promoter_size
        
    R_mat_loaded_and_prepared_df = pd.DataFrame(data_for_R_mat_dataframe)
    return R_mat_loaded_and_prepared_df


# Making this function to be able to account for pseudoautosomal genes
def is_within_specified_range_of_gene(R_mat_loaded_and_prepared_df,gene_name,distance_threshold_within_gene,bp_pos,bp_chr):
    time_start = time.time()
    R_mat_df_subset_by_chr = R_mat_loaded_and_prepared_df[(R_mat_loaded_and_prepared_df['gene'] == gene_name) & (R_mat_loaded_and_prepared_df['chr'] == bp_chr)]
    # print("TIME POINT 1:", time.time()-time_start)
    if R_mat_df_subset_by_chr.empty:
        return False
    # to_return_temp = (R_mat_df_subset_by_chr['gene_start'].min()-distance_threshold_within_gene <= bp_pos <= R_mat_df_subset_by_chr['gene_end'].max()+distance_threshold_within_gene)
    # print("TIME POINT 2:", time.time()-time_start)
    # return to_return_temp
    return (R_mat_df_subset_by_chr['gene_start'].min()-distance_threshold_within_gene <= bp_pos <= R_mat_df_subset_by_chr['gene_end'].max()+distance_threshold_within_gene)


# Does the setup (creates a file called exons_from_updated_R_mat_file_start_and_end_coords.tsv) that acts as the setup file for 
# find_overall_exonic_and_intronic_portions_of_genome_results. In the end, this is all just for calculating the proportion of the 
# genome that is exonic and the proportion of the genome that is intronic.
def find_overall_exonic_and_intronic_portions_of_genome_setup(R_mat_file_path=R_MAT_FILE_PATH):
    if not os.path.exists("../intermediate_files/v14"):
        os.makedirs("../intermediate_files/v14")
    data = mat4py.loadmat(R_mat_file_path)
    chr_vals = np.array(data['R']['chr'])
    start_coords = np.array(data['R']['code_start']).astype(float)
    end_coords = np.array(data['R']['code_end']).astype(float)
    exon_start_coords = list(map(lambda x: completely_flatten_list(x), (data['R']['exon_starts'])))
    exon_end_coords = list(map(lambda x: completely_flatten_list(x), (data['R']['exon_ends'])))
    genes = list(np.array(data['R']['gene']).flatten())
    gene_orientations = list(np.array(data['R']['strand']).flatten())
    # As it turns out the same gene can occur multiple times in R.mat, but I am going with the assumption that the orientation for a given gene (and its transcripts) is always the same
    dict_of_genes_and_orientations = {}
    for g_i in range(0, len(genes)):
        dict_of_genes_and_orientations[genes[g_i]] = gene_orientations[g_i]

    # This is for producing the intervals of genes:
    chr_and_start_and_end_coords = np.concatenate((chr_vals, start_coords, end_coords), axis=1)
    d = {'chr': chr_and_start_and_end_coords[:,0], 'start': chr_and_start_and_end_coords[:,1], 'end': chr_and_start_and_end_coords[:,2]}
    chr_and_start_and_end_coords_df = pd.DataFrame(data=d)
    # Need to drop any columns with nan start/end coordinates:
    chr_and_start_and_end_coords_df = chr_and_start_and_end_coords_df[(chr_and_start_and_end_coords_df['start'] != 'nan') & (chr_and_start_and_end_coords_df['end'] != 'nan')]
    chr_and_start_and_end_coords_df['start'] = np.array(chr_and_start_and_end_coords_df.copy()['start']).astype(float).astype(int)
    chr_and_start_and_end_coords_df['end'] = np.array(chr_and_start_and_end_coords_df.copy()['end']).astype(float).astype(int)
    chr_and_start_and_end_coords_df.to_csv('../intermediate_files/v14/genes_from_updated_R_mat_file_start_and_end_coords.tsv',header=False,sep='\t',index=False)

    # This is for producing the intervals of exons:
    exon_chr_vals = []
    exon_start_coord_vals = []
    exon_end_coord_vals = []
    for i in range(0,len(chr_vals)):
        for j in range(0,len(exon_start_coords[i])):
            exon_chr_vals.append(chr_vals[i][0])
            exon_start_coord_vals.append(exon_start_coords[i][j])
            exon_end_coord_vals.append(exon_end_coords[i][j])
    d_exon = {'chr': exon_chr_vals, 'start': exon_start_coord_vals, 'end': exon_end_coord_vals}
    exon_chr_and_start_and_end_coords_df = pd.DataFrame(data=d_exon)
    exon_chr_and_start_and_end_coords_df.to_csv('../intermediate_files/v14/exons_from_updated_R_mat_file_start_and_end_coords.tsv',header=False,sep='\t',index=False)
    # After this function (find_overall_exonic_and_intronic_portions_of_genome_setup) is run, need to run the commands from the following two files:
    # - commands_i_ran_to_get_merged_intervals_for_genes_from_R_mat_file.txt
    # - commands_i_ran_to_get_merged_intervals_for_exons_from_R_mat_file.txt
    # This allows me to work with the output of these commands for running the next find_overall_exonic_and_intronic_portions_of_genome_results


# Finds the proportion of the genome that is exonic and the proportion that is intronic. This is used downstream to calculate 
# the probability of LOF for the intergenic case.
def find_overall_exonic_and_intronic_portions_of_genome_results(path_to_coding_intervals_tsv='../intermediate_files/v14/genes_from_updated_R_mat_file_start_and_end_coords_sorted_merged.tsv', path_to_exon_intervals_tsv='../intermediate_files/v14/exons_from_updated_R_mat_file_start_and_end_coords_sorted_merged.tsv'):
    merged_gene_intervals_df = pd.read_csv(path_to_coding_intervals_tsv, sep='\t',names=["chr", "start", "end"])
    merged_gene_intervals_df['interval width'] = merged_gene_intervals_df['end'].subtract(merged_gene_intervals_df['start'])
    span_that_gene_coding_regions_take_up = merged_gene_intervals_df['interval width'].sum()
    merged_exon_intervals_df = pd.read_csv(path_to_exon_intervals_tsv, sep='\t',names=["chr", "start", "end"])
    merged_exon_intervals_df['interval width'] = merged_exon_intervals_df['end'].subtract(merged_exon_intervals_df['start'])
    span_that_exonic_regions_take_up = merged_exon_intervals_df['interval width'].sum()
    # This is the total length of the genome for males (meaning Y is included) according to hg38:
    length_of_genome_male = 3088269832
    # This is the total length of the genome for females (meaning Y is excluded) according to hg38:
    length_of_genome_female = 3031042417
    # This is the average length of the genome (if approximately 50% of the population is male vs. female)
    length_of_genome = int(np.round(length_of_genome_male * 0.5 + length_of_genome_female * 0.5))
    proportion_of_genome_thats_intronic_or_exonic_interval = span_that_gene_coding_regions_take_up / length_of_genome
    proportion_of_genome_thats_exonic_interval = span_that_exonic_regions_take_up / length_of_genome
    proportion_of_genome_thats_intronic_interval = proportion_of_genome_thats_intronic_or_exonic_interval - proportion_of_genome_thats_exonic_interval
    return proportion_of_genome_thats_exonic_interval, proportion_of_genome_thats_intronic_interval


# Finds the amount of the genome that is covered by an annotated gene in the R_mat_file for hg38. This is used to find the average distance between genes 
# for the method find_average_distance_between_genes.
def find_amount_of_genome_covered_by_genes_in_bp(R_mat_file_path=R_MAT_FILE_PATH):
    R_mat_loaded_and_prepared_df = load_and_prepare_R_mat_file_into_df(R_mat_file_path)
    total_bp_covered_by_genome = 0
    for chr in range(1,25):
        R_mat_subset_chr_df = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['chr'] == chr]
        gene_starts_and_ends = []
        for i,r in R_mat_subset_chr_df.iterrows():
            gene_starts_and_ends.append([r['gene_start'],r['gene_end']])
        gene_starts_and_ends.sort()
        merged_intervals = [gene_starts_and_ends[0]]
        for interv in gene_starts_and_ends[1:]:
            if merged_intervals[-1][0] <= interv[0] <= merged_intervals[-1][-1]:
                merged_intervals[-1][-1] = max(merged_intervals[-1][-1], interv[-1])
            else:
                merged_intervals.append(interv)
        total_bp_covered_by_genome += np.sum(list(map(lambda x: np.abs(x[1]-x[0]), merged_intervals)))
    return total_bp_covered_by_genome 


# Find the average distance between genes in the genome based off the R_mat_file. We can find a windowing parameter for the LoF enrichment model
# based off this.
def find_average_distance_between_genes(R_mat_file_path=R_MAT_FILE_PATH):
    total_bp_of_genome_covered_by_genes = find_amount_of_genome_covered_by_genes_in_bp(R_mat_file_path=R_mat_file_path)
    # This is the total length of the genome for males (meaning Y is included) according to hg38:
    length_of_genome_male = 3088269832
    # This is the total length of the genome for females (meaning Y is excluded) according to hg38:
    length_of_genome_female = 3031042417
    # This is the average length of the genome (if approximately 50% of the population is male vs. female)
    length_of_genome = int(np.round(length_of_genome_male * 0.5 + length_of_genome_female * 0.5))
    total_bp_of_genome_not_covered_by_genes = length_of_genome - total_bp_of_genome_covered_by_genes
    R_mat_loaded_and_prepared_df = load_and_prepare_R_mat_file_into_df(R_mat_file_path)
    average_distance_between_genes = total_bp_of_genome_not_covered_by_genes / (len(set(R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['chr'].isin(list(range(1,25)))]['gene'])) + 24)
    return int(np.round(average_distance_between_genes))


### NOTE: The method find_probability_of_intergenic_vs_intragenic_gene_for_cohort didn't actually account for chromosome matching for eligible transcripts before,
###       and also didn't properly account for pseudoautosomal genes. It is now defunct!
# # Finds the overall probability of an SV being intergenic vs. intragenic based on the cohort. This is just based off the counts of intergenic events and intragenic events that is returned by this function.
# def find_probability_of_intergenic_vs_intragenic_gene_for_cohort(dict_of_genes_to_reannotated_filtered_SV_lists, threshold_of_proximity, R_mat_file_path=R_MAT_FILE_PATH):
#     print("Finding probability of an event being intergenic vs. intragenic for cohort...")
#     R_mat_loaded_and_prepared_df = load_and_prepare_R_mat_file_into_df(R_mat_file_path)
#     counts_of_each_type_of_event = {'intragenic': 0, 'intergenic': 0}
#     for gene in dict_of_genes_to_reannotated_filtered_SV_lists:
#         SV_sublist_reannotated_df = pd.read_csv(dict_of_genes_to_reannotated_filtered_SV_lists[gene],sep='\t')
#         gene_start = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]['gene_start'].min()
#         gene_end = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]['gene_end'].max()
#         gene_chr = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]['chr'].iloc[0]
#         num_intragenic_this_gene = 0
#         num_intergenic_this_gene = 0
#         for index, row in SV_sublist_reannotated_df.iterrows():
#             num_events_not_within_threshold_of_gene = 0
#             if not ((gene_start-threshold_of_proximity <= row['pos1'] <= gene_end+threshold_of_proximity) and (row['chr1'] == gene_chr)):
#                 num_events_not_within_threshold_of_gene += 1
#             if not ((gene_start-threshold_of_proximity <= row['pos2'] <= gene_end+threshold_of_proximity) and (row['chr2'] == gene_chr)):
#                 num_events_not_within_threshold_of_gene += 1
#             if num_events_not_within_threshold_of_gene == 0:
#                 counts_of_each_type_of_event['intragenic'] += 1
#                 num_intragenic_this_gene += 1
#             elif num_events_not_within_threshold_of_gene == 1:
#                 counts_of_each_type_of_event['intergenic'] += 1
#                 num_intergenic_this_gene += 1
#             else:
#                 raise("This shouldn't be possible, or else something is wrong")
#         print("gene:", gene)
#         print("number of intragenic events:", num_intragenic_this_gene)
#         print("number of total events:", num_intragenic_this_gene+num_intergenic_this_gene)
#     return counts_of_each_type_of_event


def filter_dataframe_to_have_only_events_with_bps_within_threshold_of_gene(gene_name,threshold_within_gene_edges_bp,SV_list_df,R_mat_loaded_and_prepared_df,which_project):
    # print("gene:", gene_name)
    if not os.path.exists("../SV_list_subsets_of_regions_within_distance_of_gene/v14/" + which_project + "/original_annotations"):
        os.makedirs("../SV_list_subsets_of_regions_within_distance_of_gene/v14/" + which_project + "/original_annotations")
    filtered_df_file_path = "../SV_list_subsets_of_regions_within_distance_of_gene/v14/" + which_project + "/original_annotations/"+"SV_list_filtered_"+gene_name+"_within_"+str(threshold_within_gene_edges_bp)+"bp.orig.tsv"
    if not os.path.isfile(filtered_df_file_path):
        
        # time_start = time.time()
        
        # print("* filter_dataframe_to_have_only_events_with_bps_within_threshold_of_gene Time point 1:", time.time()-time_start)
        
        ###### I think we can try filtering based off the chromosomes corresponding to the gene,
        ###### and then keep everything else the same
        # R_mat_loaded_and_prepared_df_subset = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene_name]
        R_mat_loaded_and_prepared_df_subset = R_mat_loaded_and_prepared_df.iloc[np.where(R_mat_loaded_and_prepared_df['gene'].to_numpy() == gene_name)[0]]
        # print("* filter_dataframe_to_have_only_events_with_bps_within_threshold_of_gene Time point 2:", time.time()-time_start)
        gene_chrs = set(R_mat_loaded_and_prepared_df_subset['chr'])
        # print("* filter_dataframe_to_have_only_events_with_bps_within_threshold_of_gene Time point 3:", time.time()-time_start)
        # print("* filter_dataframe_to_have_only_events_with_bps_within_threshold_of_gene Time point 4:", time.time()-time_start)

        # gene_chr = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene_name]['chr'].iloc[0]
        # gene_start = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene_name]['gene_start'].min()
        # gene_end = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene_name]['gene_end'].max()
        # SV_filtered_df = SV_df[((SV_df['chr1'] == gene_chr) & (SV_df['pos1'].between(gene_start-threshold_within_gene_edges_bp,gene_end+threshold_within_gene_edges_bp))) | \
        #                        ((SV_df['chr2'] == gene_chr) & (SV_df['pos2'].between(gene_start-threshold_within_gene_edges_bp,gene_end+threshold_within_gene_edges_bp)))].copy()
        
        # print("pd.Series(list(map(lambda coords: (is_within_specified_range_of_gene(R_mat_loaded_and_prepared_df,gene_name,threshold_within_gene_edges_bp,coords[0],coords[1]) or \
        #                                                 is_within_specified_range_of_gene(R_mat_loaded_and_prepared_df,gene_name,threshold_within_gene_edges_bp,coords[2],coords[3])), \
        #                                 zip(list(SV_df['pos1']),list(SV_df['chr1']),list(SV_df['pos2']),list(SV_df['chr2']))))).astype('bool'):")
        # print(pd.Series(list(map(lambda coords: (is_within_specified_range_of_gene(R_mat_loaded_and_prepared_df,gene_name,threshold_within_gene_edges_bp,coords[0],coords[1]) or \
        #                                                 is_within_specified_range_of_gene(R_mat_loaded_and_prepared_df,gene_name,threshold_within_gene_edges_bp,coords[2],coords[3])), 
        #                                 zip(list(SV_df['pos1']),list(SV_df['chr1']),list(SV_df['pos2']),list(SV_df['chr2']))))).astype('bool'))
        # SV_filtered_df = SV_df[pd.Series(list(map(lambda coords: (is_within_specified_range_of_gene(R_mat_loaded_and_prepared_df,gene_name,threshold_within_gene_edges_bp,coords[0],coords[1]) or \
        #                                                 is_within_specified_range_of_gene(R_mat_loaded_and_prepared_df,gene_name,threshold_within_gene_edges_bp,coords[2],coords[3])), \
        #                                 zip(list(SV_df['pos1']),list(SV_df['chr1']),list(SV_df['pos2']),list(SV_df['chr2']))))).astype('bool')].copy()
        
        # This step was done like this to reduce runtime:
        if len(gene_chrs) <= 1:
            gene_chr = R_mat_loaded_and_prepared_df_subset['chr'].iat[0]
            gene_start = np.amin(R_mat_loaded_and_prepared_df_subset['gene_start'].to_numpy())
            gene_end = np.amax(R_mat_loaded_and_prepared_df_subset['gene_end'].to_numpy())
            SV_list_df = SV_list_df.iloc[np.where((SV_list_df['chr1'].to_numpy() == gene_chr) | (SV_list_df['chr2'].to_numpy() == gene_chr))]
            R_mat_loaded_and_prepared_df_subset = R_mat_loaded_and_prepared_df.iloc[np.where(R_mat_loaded_and_prepared_df['gene'].to_numpy() == gene_name)[0]]
            # SV_filtered_df = SV_list_df[((SV_list_df['chr1'] == gene_chr) & (SV_list_df['pos1'].between(gene_start-threshold_within_gene_edges_bp,gene_end+threshold_within_gene_edges_bp))) | \
            #                             ((SV_list_df['chr2'] == gene_chr) & (SV_list_df['pos2'].between(gene_start-threshold_within_gene_edges_bp,gene_end+threshold_within_gene_edges_bp)))].copy()
            SV_filtered_df = SV_list_df.iloc[np.where(((SV_list_df['chr1'].to_numpy() == gene_chr) & (SV_list_df['pos1'].to_numpy() >= gene_start-threshold_within_gene_edges_bp) & (SV_list_df['pos1'].to_numpy() <= gene_end+threshold_within_gene_edges_bp)) | \
                                                      ((SV_list_df['chr2'].to_numpy() == gene_chr) & (SV_list_df['pos2'].to_numpy() >= gene_start-threshold_within_gene_edges_bp) & (SV_list_df['pos2'].to_numpy() <= gene_end+threshold_within_gene_edges_bp)))]
        else:
            SV_list_df = SV_list_df[SV_list_df['chr1'].isin(gene_chrs) | SV_list_df['chr2'].isin(gene_chrs)].copy()
            SV_filtered_df = SV_list_df[list(map(lambda coords: bool(is_within_specified_range_of_gene(R_mat_loaded_and_prepared_df_subset,gene_name,threshold_within_gene_edges_bp,coords[0],coords[1]) or \
                                                                     is_within_specified_range_of_gene(R_mat_loaded_and_prepared_df_subset,gene_name,threshold_within_gene_edges_bp,coords[2],coords[3])), \
                                                    zip(list(SV_list_df['pos1']),list(SV_list_df['chr1']),list(SV_list_df['pos2']),list(SV_list_df['chr2']))))].copy()

        
        # # This is how things were originally (most recently):
        # SV_filtered_df = SV_list_df[list(map(lambda coords: bool(is_within_specified_range_of_gene(R_mat_loaded_and_prepared_df_subset,gene_name,threshold_within_gene_edges_bp,coords[0],coords[1]) or \
        #                                                 is_within_specified_range_of_gene(R_mat_loaded_and_prepared_df_subset,gene_name,threshold_within_gene_edges_bp,coords[2],coords[3])), \
        #                                 zip(list(SV_list_df['pos1']),list(SV_list_df['chr1']),list(SV_list_df['pos2']),list(SV_list_df['chr2']))))].copy()
        # print("* filter_dataframe_to_have_only_events_with_bps_within_threshold_of_gene Time point 5:", time.time()-time_start)
        SV_filtered_df.to_csv(filtered_df_file_path, sep='\t',index=False)
        # print("* filter_dataframe_to_have_only_events_with_bps_within_threshold_of_gene Time point 6:", time.time()-time_start)
    else:
        pass
    return filtered_df_file_path

def find_intron_in_merged_transcript(gene_name, merged_transcript, pos, multi_chr_gene, chrs_to_R_mat_subset_dfs_dict_for_multi_chr_gene):
    intronic_regions = [merged_transcript[i] for i in range(len(merged_transcript)) if merged_transcript[i][2][:6] == "intron"]
    
    
    #for pseudoautosomal - if PAR1, then merged transcript is right. if PAR2 and on the X (23) chromosome (around pos 155,000,000, for Y it's around 57,000,000) then we have to translate the coordinate to the saved merged transcript on chrY. The boundary is not very specific so I just picked the round number of 100,000,000.
    if multi_chr_gene:
        chr23_subset = chrs_to_R_mat_subset_dfs_dict_for_multi_chr_gene[23].reset_index()
        chr24_subset = chrs_to_R_mat_subset_dfs_dict_for_multi_chr_gene[24].reset_index()
        if pos > 100000000 and int(intronic_regions[0][0]) < 100000000: # need to change
            pos -= chr23_subset.loc[0, "tx_start"] - chr24_subset.loc[0, "tx_start"]
        elif pos < 100000000 and int(intronic_regions[0][0]) > 100000000: # need to change
            pos += chr23_subset.loc[0, "tx_start"] - chr24_subset.loc[0, "tx_start"]
        # else: it's all ok
        # print("New position:", pos)

    for a_region in intronic_regions:
        if pos >= int(a_region[0]) and pos <= int(a_region[1]):
            return int(a_region[2][6:])
    raise("intron not found")


def modified_dRanger_annotate_sites_driver_CGC_genes_prioritized(X_df_path, gene, distance_threshold_within_gene, R_mat_loaded_and_prepared_df, chrs_to_R_mat_subset_dfs_dict, which_project): # Remember X_df should be a subset of the original SV_list dataframe for just the sites corresponding to a given gene
        
    def bp2str(b):
        possible_strings = []
        possible_strings.append(str(int(b)) + 'bp')
        if b >= 1000:
            possible_strings.append(str(int(np.round(b / 1000))) + 'Kb')
        if b >= 1000000:
            possible_strings.append(str(int(np.round(b / 1000000))) + 'Mb')
        ind = np.argmin(list(map(len, possible_strings)))
        return possible_strings[ind]
    
    time_start = time.time()
    
    # print("Gene:", gene)
    
    R_mat_loaded_and_prepared_df_subset = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]
    chrs_for_gene = set(R_mat_loaded_and_prepared_df_subset['chr'])
    multi_chr_gene = bool(len(chrs_for_gene) > 1)
    
    chrs_to_R_mat_subset_dfs_dict_for_multi_chr_gene = {}
    if not multi_chr_gene:
        R_mat_subset_for_gene_and_chr_assuming_one_chr_for_gene = chrs_to_R_mat_subset_dfs_dict[list(chrs_for_gene)[0]][chrs_to_R_mat_subset_dfs_dict[list(chrs_for_gene)[0]]['gene'] == gene]
        # Additional variables to track added by me (Xavi) to make things a bit easier:
        # gene_chr = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]['chr'].iloc[0]
        # gene_start = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]['gene_start'].min()
        # gene_end = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]['gene_end'].max()
        gene_chr = R_mat_loaded_and_prepared_df_subset['chr'].iat[0]
        gene_start = np.amin(R_mat_loaded_and_prepared_df_subset['gene_start'].to_numpy())
        gene_end = np.amax(R_mat_loaded_and_prepared_df_subset['gene_end'].to_numpy())
    else:
        for c in chrs_for_gene:
            chrs_to_R_mat_subset_dfs_dict_for_multi_chr_gene[c] = chrs_to_R_mat_subset_dfs_dict[c][chrs_to_R_mat_subset_dfs_dict[c]['gene'] == gene]
    
    # print(chrs_to_R_mat_subset_dfs_dict_for_multi_chr_gene)

    # Modification by me (Xavi): Going to keep this list (list_of_cgc_driver_transloc_genes) empty. That's fine, because when we reannotate it doesn't matter, because the reannotations will only be done relative to a specific gene's transcripts.
    # So we don't have to prioritize certain driver genes' transcripts over other genes' transcripts
    list_of_cgc_driver_transloc_genes = []
    inferred_threshold_within_gene_edges_bp = int(X_df_path.split('bp.orig.tsv')[0].split('_')[-1])
    if not os.path.isdir("../SV_list_subsets_of_regions_within_distance_of_gene/v14/" + which_project + "/annotations_relative_to_given_gene_and_boundary_threshold"):
        os.makedirs("../SV_list_subsets_of_regions_within_distance_of_gene/v14/" + which_project + "/annotations_relative_to_given_gene_and_boundary_threshold")
    reannotated_file_path = "../SV_list_subsets_of_regions_within_distance_of_gene/v14/" + which_project + "/annotations_relative_to_given_gene_and_boundary_threshold/"+"SV_list_filtered_"+gene+"_within_"+str(inferred_threshold_within_gene_edges_bp)+"bp.reannotated.tsv"
    # print("* modified_dRanger_annotate_sites_driver_CGC_genes_prioritized Time point 1:", time.time()-time_start)
    # Modification I (Xavi) made: Only runs/produces a newly annotated file if the output file as specified by reannotated_file_path does not already exist
    if not os.path.isfile(reannotated_file_path):
        X_df = pd.read_csv(X_df_path,sep='\t')
        # print("* modified_dRanger_annotate_sites_driver_CGC_genes_prioritized Time point 2:", time.time()-time_start)
        ## Now the actual start of the logic of the dRanger annotator:
        X_df_gene1_list = []
        X_df_site1_list = []
        X_df_gene2_list = []
        X_df_site2_list = []
        X_df_fusion_list = []
        # % process each rearrangement
        nx = X_df.shape[0]
        # print("* modified_dRanger_annotate_sites_driver_CGC_genes_prioritized Time point 3:", time.time()-time_start)
        for x in range(1,nx+1):
            time_start_annotator = time.time()
            if x % 1000 == 0:
                print(str(x) + '/' + str(nx))
            # % first: identify each end
            intronnum = np.nan
            intronframe = np.nan
            # print("** annotator Time point 1:", time.time()-time_start_annotator)
            X_df_row = X_df.iloc[x-1]
            for ee in [1,2]:

                # print("*** annotator nested Time point 1:", time.time()-time_start_annotator)
                if ee == 1:
                    chr_ = X_df_row['chr1']
                    pos_ = X_df_row['pos1']
                    str_ = X_df_row['str1']
                else:
                    chr_ = X_df_row['chr2']
                    pos_ = X_df_row['pos2']
                    str_ = X_df_row['str2']
                # # original:
                # overlaps_df = R_mat_loaded_and_prepared_df[(R_mat_loaded_and_prepared_df['chr'] == chr_) & (R_mat_loaded_and_prepared_df['gene_start'] <= pos_) & (R_mat_loaded_and_prepared_df['gene_end'] >= pos_)]
                # modified for our purposes:
                # print("*** annotator nested Time point 2:", time.time()-time_start_annotator)
                name = '' # having the name variable start off with an empty value
                a = '' # having the a variable start off with an empty value
                zone = -1 # having the zone variable start off with a -1 value
                strand = '' # having the a variable start off with an empty value
                # print("*** annotator nested Time point 3:", time.time()-time_start_annotator)
                boolean_if_overlaps_df_is_def_empty = False
                if chr_ in chrs_for_gene:
                    if not multi_chr_gene:
                        # overlaps_df = R_mat_subset_for_gene_and_chr_assuming_one_chr_for_gene[(R_mat_subset_for_gene_and_chr_assuming_one_chr_for_gene['gene_start'] <= pos_) & (R_mat_subset_for_gene_and_chr_assuming_one_chr_for_gene['gene_end'] >= pos_)]
                        overlaps_df_indices = list(set(np.nonzero(np.array(R_mat_subset_for_gene_and_chr_assuming_one_chr_for_gene['gene_start'])<=pos_)[0]).intersection(set(np.nonzero(np.array(R_mat_subset_for_gene_and_chr_assuming_one_chr_for_gene['gene_end'])>=pos_)[0])))
                        # print("overlaps_df_indices:", overlaps_df_indices)
                        if not overlaps_df_indices:
                            # print("OKOKOKOK")
                            boolean_if_overlaps_df_is_def_empty = True
                        overlaps_df = R_mat_subset_for_gene_and_chr_assuming_one_chr_for_gene.iloc[overlaps_df_indices]
                    else:
                        overlaps_df = chrs_to_R_mat_subset_dfs_dict[chr_][(chrs_to_R_mat_subset_dfs_dict[chr_]['gene'] == gene) & (chrs_to_R_mat_subset_dfs_dict[chr_]['gene_start'] <= pos_) & (chrs_to_R_mat_subset_dfs_dict[chr_]['gene_end'] >= pos_)]
                        if len(overlaps_df) == 0:
                            boolean_if_overlaps_df_is_def_empty = True
                else:
                    boolean_if_overlaps_df_is_def_empty = True
                # print("*** annotator nested Time point 4:", time.time()-time_start_annotator)
                # if boolean_if_overlaps_df_is_def_empty or overlaps_df.empty: 
                if boolean_if_overlaps_df_is_def_empty: 
                    # it's intergenic
                    # Making this into an if-else clause to make my code more efficient runtime-wise
                    # print("SUBTEST TIME POINT 1:", time.time()-time_start_annotator)
                    if not multi_chr_gene:
                        # print("SUBTEST TIME POINT 2:", time.time()-time_start_annotator)
                        if (chr_ == gene_chr) and (gene_start-distance_threshold_within_gene <= pos_ <= gene_end+distance_threshold_within_gene):
                            # print("SUBTEST TRUE TIME POINT 3:", time.time()-time_start_annotator)
                            dist_before = R_mat_subset_for_gene_and_chr_assuming_one_chr_for_gene['tx_start'].subtract(pos_).abs()
                            # print("SUBTEST TRUE TIME POINT 4:", time.time()-time_start_annotator)
                            dist_after = (-R_mat_subset_for_gene_and_chr_assuming_one_chr_for_gene['tx_end']).add(pos_).abs()
                            # print("SUBTEST TRUE TIME POINT 5:", time.time()-time_start_annotator)
                        else:
                            # print("SUBTEST FALSE TIME POINT 3:", time.time()-time_start_annotator)
                            dist_before = chrs_to_R_mat_subset_dfs_dict[chr_]['tx_start'].subtract(pos_).abs()
                            # print("SUBTEST FALSE TIME POINT 4:", time.time()-time_start_annotator)
                            dist_after = (-chrs_to_R_mat_subset_dfs_dict[chr_]['tx_end']).add(pos_).abs()
                            # print("SUBTEST FALSE TIME POINT 5:", time.time()-time_start_annotator)
                    else:
                        if is_within_specified_range_of_gene(R_mat_loaded_and_prepared_df,gene,distance_threshold_within_gene,pos_,chr_):
                            dist_before = chrs_to_R_mat_subset_dfs_dict_for_multi_chr_gene[chr_]['tx_start'].subtract(pos_).abs()
                            dist_after = (-chrs_to_R_mat_subset_dfs_dict_for_multi_chr_gene[chr_]['tx_end']).add(pos_).abs()
                        else:
                            dist_before = chrs_to_R_mat_subset_dfs_dict[chr_]['tx_start'].subtract(pos_).abs()
                            dist_after = (-chrs_to_R_mat_subset_dfs_dict[chr_]['tx_end']).add(pos_).abs()
                    # print("SUBTEST TIME POINT 1:", time.time()-time_start_annotator)
                    # if is_within_specified_range_of_gene(R_mat_loaded_and_prepared_df,gene,distance_threshold_within_gene,pos_,chr_):
                    # # if (chr_ == gene_chr) and (gene_start-distance_threshold_within_gene <= pos_ <= gene_end+distance_threshold_within_gene):
                    #     print("SUBTEST TRUE TIME POINT 2:", time.time()-time_start_annotator)
                    #     dist_before = chrs_to_R_mat_subset_dfs_dict[chr_][chrs_to_R_mat_subset_dfs_dict[chr_]['gene'] == gene]['tx_start'].subtract(pos_).abs()
                    #     print("SUBTEST TRUE TIME POINT 3:", time.time()-time_start_annotator)
                    #     dist_after = (-chrs_to_R_mat_subset_dfs_dict[chr_][chrs_to_R_mat_subset_dfs_dict[chr_]['gene'] == gene]['tx_end']).add(pos_).abs()
                    #     print("SUBTEST TRUE TIME POINT 4:", time.time()-time_start_annotator)
                    # else:
                    #     print("SUBTEST FALSE TIME POINT 2:", time.time()-time_start_annotator)
                    #     dist_before = chrs_to_R_mat_subset_dfs_dict[chr_]['tx_start'].subtract(pos_).abs()
                    #     print("SUBTEST FALSE TIME POINT 3:", time.time()-time_start_annotator)
                    #     dist_after = (-chrs_to_R_mat_subset_dfs_dict[chr_]['tx_end']).add(pos_).abs()
                    #     print("SUBTEST FALSE TIME POINT 4:", time.time()-time_start_annotator)
                    # print("*** annotator nested EMPTY CASE Time point 5:", time.time()-time_start_annotator)
                    for y in list(np.arange(1,3,0.3)):
                        # print("****** ITERATIVE DISTANCE Time point 1:", time.time()-time_start_annotator)
                        if np.any(dist_before.to_numpy() <= 1000**y):
                            # print("****** ITERATIVE DISTANCE dist_before Time point 1:", time.time()-time_start_annotator)
                            # idx = list(dist_before[dist_before <= 1000**y].index)
                            # print("****** ITERATIVE DISTANCE dist_before Time point 2:", time.time()-time_start_annotator)
                            # i = dist_before[idx].idxmin()
                            i = dist_before.idxmin()
                            # tmp = dist_before.loc[i]
                            # print("****** ITERATIVE DISTANCE dist_before Time point 3:", time.time()-time_start_annotator)
                            # print('gene:',gene)
                            # print('chr_:',chr_)
                            # print('pos_:',pos_)
                            # print('dist_before:\n',dist_before)
                            name = chrs_to_R_mat_subset_dfs_dict[chr_].at[i,'gene']
                            # print("****** ITERATIVE DISTANCE dist_before Time point 4:", time.time()-time_start_annotator)
                            strand = chrs_to_R_mat_subset_dfs_dict[chr_].at[i,'strand']
                            # print("****** ITERATIVE DISTANCE dist_before Time point 5:", time.time()-time_start_annotator)
                            zone = 0
                            # print("****** ITERATIVE DISTANCE dist_before Time point 6:", time.time()-time_start_annotator)
                            a = 'IGR: ' + bp2str(dist_before.at[i]) + ' before ' + name +  '(' + strand + ')'
                            # print("****** ITERATIVE DISTANCE dist_before Time point 7:", time.time()-time_start_annotator)
                            break
                        # print("****** ITERATIVE DISTANCE Time point 2:", time.time()-time_start_annotator)
                        if np.any(dist_after.to_numpy() <= 1000**y):
                            # print("****** ITERATIVE DISTANCE dist_after Time point 1:", time.time()-time_start_annotator)
                            # idx = list(dist_after[dist_after <= 1000**y].index)
                            # print("****** ITERATIVE DISTANCE dist_after Time point 2:", time.time()-time_start_annotator)
                            # i = dist_after[idx].idxmin()
                            i = dist_after.idxmin()
                            # tmp = dist_after.loc[i]
                            # print("****** ITERATIVE DISTANCE dist_after Time point 3:", time.time()-time_start_annotator)
                            name = chrs_to_R_mat_subset_dfs_dict[chr_].at[i,'gene']
                            # print("****** ITERATIVE DISTANCE dist_after Time point 4:", time.time()-time_start_annotator)
                            strand = chrs_to_R_mat_subset_dfs_dict[chr_].at[i,'strand']
                            # print("****** ITERATIVE DISTANCE dist_after Time point 5:", time.time()-time_start_annotator)
                            zone = 0
                            # print("****** ITERATIVE DISTANCE dist_after Time point 6:", time.time()-time_start_annotator)
                            a = 'IGR: ' + bp2str(dist_after.at[i]) + ' after ' + name +  '(' + strand + ')'
                            # print("****** ITERATIVE DISTANCE dist_after Time point 7:", time.time()-time_start_annotator)
                            break
                        # print("****** ITERATIVE DISTANCE Time point 3:", time.time()-time_start_annotator)
                    # print("*** annotator nested EMPTY CASE Time point 6:", time.time()-time_start_annotator)
                else:
                    # it's within a transcript: determine promoter/UTR/intron/exon
                    no = overlaps_df.shape[0]
                    # Note: TBTGFC means "Transcript belongs to gene from COSMIC"
                    c = np.full(no, np.nan) # zone: 1=TBTGFC+exon, 2=TBTGFC+intron, 3=TBTGFC+3'-UTR, 4=TBTGFC+5'-UTR, 5=TBTGFC+promoter, 6=exon, 7=intron, 8=3'-UTR, 9=5'-UTR, 10=promoter
                    d = np.full(no, np.nan) # for exons: which one, and how far
                    e = np.full(no, np.nan) # for exons: which one, and how far
                    d1 = np.full(no, np.nan) # for introns: between which exons and how far?
                    d2 = np.full(no, np.nan) # for introns: between which exons and how far?
                    e1 = np.full(no, np.nan) # for introns: between which exons and how far?
                    e2 = np.full(no, np.nan) # for introns: between which exons and how far?
                    f = np.full(no, np.nan) # for introns: how many bases in the partially completed codon?
                    # print("*** annotator nested NON-EMPTY CASE Time point 5:", time.time()-time_start_annotator)
                    for j in range(0,no):
                        i = overlaps_df.index[j]
                        relevant_transcript_row = R_mat_loaded_and_prepared_df.loc[i]
                        gene_name = relevant_transcript_row['gene']
                        # is contained in the list of COSMIC driver genes with differing transcripts?
                        addend_for_if_not_TBTGFC = 0
                        if not gene_name in list_of_cgc_driver_transloc_genes:
                            addend_for_if_not_TBTGFC = 5
                        # in promoter?
                        if pos_ < relevant_transcript_row['tx_start']:
                            c[j] = 5 + addend_for_if_not_TBTGFC
                            d[j] = relevant_transcript_row['tx_start'] - pos_
                            continue
                        if pos_ > relevant_transcript_row['tx_end']:
                            c[j] = 5 + addend_for_if_not_TBTGFC
                            d[j] = pos_ - relevant_transcript_row['tx_end']
                            continue
                        # in UTR?
                        if relevant_transcript_row['strand'] == '+':
                            if relevant_transcript_row['code_start'] > pos_:
                                c[j] = 4 + addend_for_if_not_TBTGFC
                                d[j] = relevant_transcript_row['code_start'] - pos_
                                continue
                            if relevant_transcript_row['code_end'] < pos_:
                                c[j] = 3 + addend_for_if_not_TBTGFC
                                d[j] = pos_-relevant_transcript_row['code_end']
                                continue
                        else: # (-)
                            if relevant_transcript_row['code_start'] > pos_:
                                c[j] = 3 + addend_for_if_not_TBTGFC
                                d[j] = relevant_transcript_row['code_start'] - pos_
                                continue
                            if relevant_transcript_row['code_end'] < pos_:
                                c[j] = 4 + addend_for_if_not_TBTGFC
                                d[j] = pos_-relevant_transcript_row['code_end']
                                continue
                        # in exon(s)?
                        in_e = list(np.array(range(1,len(relevant_transcript_row['exon_starts'])+1))[(relevant_transcript_row['exon_starts'] <= pos_) & (relevant_transcript_row['exon_ends'] >= pos_)])
                        if in_e:
                            c[j] = 1 + addend_for_if_not_TBTGFC
                            e[j] = np.array(in_e)
                            continue
                        # otherwise: in intron
                        c[j] = 2 + addend_for_if_not_TBTGFC
                        for k in range(0, relevant_transcript_row['n_exons']-1):
                            if (relevant_transcript_row['exon_ends'][k] < pos_) and (relevant_transcript_row['exon_starts'][k+1] > pos_):
                                if relevant_transcript_row['strand'] == '+':
                                    f[j] = relevant_transcript_row['exon_frames'][k+1]
                                else:
                                    f[j] = relevant_transcript_row['exon_frames'][k]
                                e1[j] = k+1
                                e2[j] = k+2
                                d1[j] = pos_ - relevant_transcript_row['exon_ends'][k]
                                d2[j] = relevant_transcript_row['exon_starts'][k+1] - pos_
                                d[j] = min(d1[j], d2[j])
                                break
                    # print("*** annotator nested NON-EMPTY CASE Time point 6:", time.time()-time_start_annotator)
                    # % find the transcript in the highest-priority class
                    zone = -1
                    for cidx in list(range(1,11)):
                        # print("**** annotator nested NON-EMPTY CASE Time point 6.1:", time.time()-time_start_annotator)
                        idx = [z for z in range(len(c)) if c[z] == cidx]
                        # print("**** annotator nested NON-EMPTY CASE Time point 6.2:", time.time()-time_start_annotator)
                        if not idx:
                            continue
                        if cidx == 1 or cidx == 6:
                            j = idx[0]
                        else:
                            k = np.argmin(d[np.array(idx)])
                            j = idx[k]
                        # print("**** annotator nested NON-EMPTY CASE Time point 6.3:", time.time()-time_start_annotator)
                        i = overlaps_df.index[j]
                        relevant_transcript_row = R_mat_loaded_and_prepared_df.loc[i]
                        name = relevant_transcript_row['gene']
                        strand = relevant_transcript_row['strand']
                        # print("**** annotator nested NON-EMPTY CASE Time point 6.4:", time.time()-time_start_annotator)
                        if strand == '-':
                            e[j] = relevant_transcript_row['n_exons'] - e[j] + 1
                            e1[j] = relevant_transcript_row['n_exons'] - e1[j] + 1
                            e2[j] = relevant_transcript_row['n_exons'] - e2[j] + 1
                        zone = cidx
                        # print("**** annotator nested NON-EMPTY CASE Time point 6.5:", time.time()-time_start_annotator)
                        break
                    # print("*** annotator nested NON-EMPTY CASE Time point 7:", time.time()-time_start_annotator)
                    if zone == 1: # exon + TBTGFC
                        a = "Exon %d of %s(%s)"%(e[j],name,strand)
                    elif zone == 2: # intron + TBTGFC
                        if strand == '+':
                            if d1[j] < d2[j]:
                                a = "Intron of %s(%s): %s after exon %d" % (name,strand,bp2str(d1[j]),e1[j])
                            else:
                                a = "Intron of %s(%s): %s before exon %d" % (name,strand,bp2str(d2[j]),e2[j])
                        else: # (-)
                            if d1[j] < d2[j]:
                                a = "Intron of %s(%s): %s before exon %d" % (name,strand,bp2str(d1[j]),e1[j])
                            else:
                                a = "Intron of %s(%s): %s after exon %d" % (name,strand,bp2str(d2[j]),e2[j])
                        intronnum = e1[j]
                        intronframe = f[j]
                    elif zone == 3: # 3'-UTR + TBTGFC
                        a = "3'-UTR of %s(%s): %s after coding stop" % (name,strand,bp2str(d[j]))
                    elif zone == 4: # 5'-UTR + TBTGFC
                        a = "5'-UTR of %s(%s): %s before coding start" % (name,strand,bp2str(d[j]))
                    elif zone == 5: # promoter + TBTGFC
                        a = "Promoter of %s(%s): %s from tx start" % (name,strand,bp2str(d[j]))
                    elif zone == 6: # exon
                        a = "Exon %d of %s(%s)"%(e[j],name,strand)
                    elif zone == 7: # intron
                        if strand == '+':
                            if d1[j] < d2[j]:
                                a = "Intron of %s(%s): %s after exon %d" % (name,strand,bp2str(d1[j]),e1[j])
                            else:
                                a = "Intron of %s(%s): %s before exon %d" % (name,strand,bp2str(d2[j]),e2[j])
                        else: # (-)
                            if d1[j] < d2[j]:
                                a = "Intron of %s(%s): %s before exon %d" % (name,strand,bp2str(d1[j]),e1[j])
                            else:
                                a = "Intron of %s(%s): %s after exon %d" % (name,strand,bp2str(d2[j]),e2[j])
                        intronnum = e1[j]
                        intronframe = f[j]
                    elif zone == 8: # 3'-UTR
                        a = "3'-UTR of %s(%s): %s after coding stop" % (name,strand,bp2str(d[j]))
                    elif zone == 9: # 5'-UTR
                        a = "5'-UTR of %s(%s): %s before coding start" % (name,strand,bp2str(d[j]))
                    elif zone == 10: # promoter
                        a = "Promoter of %s(%s): %s from tx start" % (name,strand,bp2str(d[j]))
                    elif zone == -1:
                        a = 'unexpected error: zone = -1'
                    # print("*** annotator nested NON-EMPTY CASE Time point 8:", time.time()-time_start_annotator)
            
                if ee == 1:
                    gene1 = name
                    X_df_gene1_list.append(gene1)
                    X_df_site1_list.append(a)
                    zone1 = zone
                    strandgene1 = strand
                    str1 = str_
                    intronnum1 = intronnum #not used anymore -J
                    # print("GETTING REGION1")
                    if zone1 in [2, 7]:
                        with open("../surrounding_annotations_for_each_gene/v14/merged_transcript/gene_" + gene1 + ".step_size_1bp.merged_transcript.tsv") as f:
                            merged_transcript = f.readline().strip("\n").split()
                        merged_transcript = [[merged_transcript[i], merged_transcript[i+1], merged_transcript[i+2]] for i in range(0, len(merged_transcript), 3)]
                        # print(merged_transcript)
                        # print(merged_transcript, pos_, intronnum1, gene1)
                        intronic_region1 = find_intron_in_merged_transcript(gene1, merged_transcript, pos_, multi_chr_gene, chrs_to_R_mat_subset_dfs_dict_for_multi_chr_gene)
                    intronframe1 = intronframe
                else:
                    gene2 = name
                    X_df_gene2_list.append(gene2)
                    X_df_site2_list.append(a)
                    zone2 = zone
                    strandgene2 = strand
                    str2 = str_
                    intronnum2 = intronnum #not used anymore -J
                    # print("GETTING REGION2")
                    if zone2 in [2, 7]:
                        with open("../surrounding_annotations_for_each_gene/v14/merged_transcript/gene_" + gene2 + ".step_size_1bp.merged_transcript.tsv") as f:
                            merged_transcript = f.readline().strip("\n").split()
                        merged_transcript = [[merged_transcript[i], merged_transcript[i+1], merged_transcript[i+2]] for i in range(0, len(merged_transcript), 3)]
                        # print(merged_transcript)
                        # print(merged_transcript, pos_, intronnum2, gene2)
                        intronic_region2 = find_intron_in_merged_transcript(gene2, merged_transcript, pos_, multi_chr_gene, chrs_to_R_mat_subset_dfs_dict_for_multi_chr_gene)
                    intronframe2 = intronframe
                
                # print("*** annotator nested Time point 9:", time.time()-time_start_annotator)

            # % next end of this rearrangement

            # % does this rearrangement generate an interesting fusion product?
            # %    zone: 0=IGR 1=exon+TBTGFC, 2=intron+TBTGFC, 3=3'-UTR+TBTGFC, 4=5'-UTR,+TBTGFC 5=promoter+TBTGFC, 6=exon, 7=intron, 8=3'-UTR, 9=5'-UTR, 10=promoter

            strandmatch1 = (strandgene1 == '+' and str1 == 0) or (strandgene1 == '-' and str1 == 1)
            strandmatch2 = (strandgene2 == '+' and str2 == 0) or (strandgene2 == '-' and str2 == 1)
            txactive1 = zone1 > 0 and zone1 != 3 and zone1 != 8
            txactive2 = zone2 > 0 and zone2 != 3 and zone2 != 8
            # print("** annotator Time point 2:", time.time()-time_start_annotator)
            txt = '-'
            if (txactive1 and strandmatch1) and (txactive2 and strandmatch2):
                txt = 'Antisense fusion'
            elif (txactive1 and strandmatch1) or (txactive2 and strandmatch2):
                if txactive1 and txactive2:
                    if gene1 == gene2: # within-gene event
                        if str1 == 0 and str2 == 1:
                            class_ = 'Deletion'
                        elif str1==1 and str2==0:
                            class_ = 'Duplication'
                        else:
                            class_ = 'Inversion'
                        if (zone1==2 or zone1==7) and (zone2==2 or zone2==7):
                            # if intronnum1==intronnum2:
                            if intronic_region1 == intronic_region2:
                                txt = class_ + ' within intron'
                            else:
                                # if intronnum1 == intronnum2: print("intronic region and num disagree!!!")
                                # numexons = int(abs(intronnum2-intronnum1))
                                numexons = int(abs(intronic_region2 - intronic_region1))
                                txt = class_ + ' of ' + str(numexons) + ' exon'
                                if numexons > 1:
                                    txt = txt + 's'
                                if not class_.lower() == 'Inversion'.lower():
                                    if intronframe1 == intronframe2:
                                        txt = txt + ': in frame'
                                    else:
                                        txt = txt + ': out of frame'
                        else:
                            txt = class_ + ' within transcript'
                            if zone1==1 or zone2==1 or zone1==6 or zone2==6:
                                txt = txt + ': mid-exon'
                    else: # inter-gene event
                        if strandmatch1:
                            fusname = gene1 + '-' + gene2
                        else:
                            fusname = gene2 + '-' + gene1
                        if (zone1==2 or zone1==7) and (zone2==2 or zone2==7):
                            if intronframe1==intronframe2:
                                txt = 'Protein fusion: in frame'
                            else:
                                txt = 'Protein fusion: out of frame'
                        elif zone1 == 1 or zone2 == 1 or zone1==6 or zone2==6:
                            txt = 'Protein fusion: mid-exon'
                        else:
                            txt = 'Transcript fusion'
                        txt = txt + ' (' + fusname + ')'
                        
            X_df_fusion_list.append(txt)
            # print("** annotator Time point 3:", time.time()-time_start_annotator)
            
        # print("* modified_dRanger_annotate_sites_driver_CGC_genes_prioritized Time point 4:", time.time()-time_start)
        
        X_df['gene1'] = X_df_gene1_list
        X_df['site1'] = X_df_site1_list
        X_df['gene2'] = X_df_gene2_list
        X_df['site2'] = X_df_site2_list
        X_df['fusion'] = X_df_fusion_list

        X_df['span'] = pd.array(list(X_df['span']), dtype="Int64") # Need to always do this to ensure span has int values
        # When X_df gets loaded in the span column values are converted to floats to account for the NaNs in the column
        
        # print("* modified_dRanger_annotate_sites_driver_CGC_genes_prioritized Time point 5:", time.time()-time_start)
        
        X_df.to_csv(reannotated_file_path, sep='\t',index=False)
        
        # print("* modified_dRanger_annotate_sites_driver_CGC_genes_prioritized Time point 6:", time.time()-time_start)
        
        # print('\nDone\n')
    
    else:
        pass
    
    return reannotated_file_path

def find_proportions_of_1_bp_2_bp_events_and_deletions_tandemdups_inversions_for_intragenic_all_genes(dict_of_genes_to_reannotated_filtered_SV_lists, R_mat_loaded_and_prepared_df, threshold_of_proximity):
    
    time_start = time.time()
    
    num_total_events = 0
    num_total_intragenic_events = 0
    num_deletions = 0
    num_tandemdups = 0
    num_inversions = 0
    
    # print("find_proportions_of_deletions_tandemdups_inversions_for_intragenic_all_genes Time point 1:", time.time()-time_start)

    genes = list(dict_of_genes_to_reannotated_filtered_SV_lists.keys())
    
    # for gene in dict_of_genes_to_reannotated_filtered_SV_lists:
    for i in tqdm(range(len(genes))):
        gene = genes[i]
        # print("find_proportions_of_deletions_tandemdups_inversions_for_intragenic_all_genes for-loop Time point 1:", time.time()-time_start)
        reannotated_df = pd.read_csv(dict_of_genes_to_reannotated_filtered_SV_lists[gene],sep='\t')
        
        # print("gene:",gene)
        # print("reannotated_df:")
        # print(reannotated_df)
        # print("len(reannotated_df):")
        # print(len(reannotated_df))
        # reannotated_df_intragenic = reannotated_df[(reannotated_df['gene1'] == gene) & (reannotated_df['gene2'] == gene)]

        # gene_chr_in_R_mat = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]['chr'].iloc[0]
        # gene_start_in_R_mat = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]['gene_start'].min()
        # gene_end_in_R_mat = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]['gene_end'].max()

        # reannotated_df_intragenic = reannotated_df[((reannotated_df['chr1'] == gene_chr_in_R_mat) & (reannotated_df['pos1'] <= gene_end_in_R_mat+threshold_of_proximity) & (reannotated_df['pos1'] >= gene_start_in_R_mat-threshold_of_proximity)) | \
        #                                            ((reannotated_df['chr2'] == gene_chr_in_R_mat) & (reannotated_df['pos2'] <= gene_end_in_R_mat+threshold_of_proximity) & (reannotated_df['pos2'] >= gene_start_in_R_mat-threshold_of_proximity))].copy()
        
        # reannotated_df_intragenic = reannotated_df[pd.Series(list(map(lambda coords: (is_within_specified_range_of_gene(R_mat_loaded_and_prepared_df,gene,threshold_of_proximity,coords[0],coords[1]) or \
        #                                                                     is_within_specified_range_of_gene(R_mat_loaded_and_prepared_df,gene,threshold_of_proximity,coords[2],coords[3])), 
        #                                                     zip(list(reannotated_df['pos1']),list(reannotated_df['chr1']),list(reannotated_df['pos2']),list(reannotated_df['chr2']))))).astype('bool')].copy()
        
        # reannotated_df_intragenic = reannotated_df[pd.Series(list(map(lambda coords: (is_within_specified_range_of_gene(R_mat_loaded_and_prepared_df,gene,threshold_of_proximity,coords[0],coords[1]) and \
        #                                                                     is_within_specified_range_of_gene(R_mat_loaded_and_prepared_df,gene,threshold_of_proximity,coords[2],coords[3])), 
        #                                                     zip(list(reannotated_df['pos1']),list(reannotated_df['chr1']),list(reannotated_df['pos2']),list(reannotated_df['chr2']))))).astype('bool')].copy()
        
        # print("find_proportions_of_deletions_tandemdups_inversions_for_intragenic_all_genes for-loop Time point 2:", time.time()-time_start)
        
        gene_chrs = set(R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]['chr'])
        if len(gene_chrs) <= 1:
            gene_chr = list(gene_chrs)[0] # It should be impossible for this indexing to fail. Otherwise the gene could not be used for annotation in the first place.
            gene_start_in_R_mat = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]['gene_start'].min()
            gene_end_in_R_mat = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]['gene_end'].max()
            # Bug fix here (Xavi): Previously I used | instead of & for combining these two booleans
            reannotated_df_intragenic = reannotated_df[((reannotated_df['chr1'] == gene_chr) & (reannotated_df['pos1'] <= gene_end_in_R_mat+threshold_of_proximity) & (reannotated_df['pos1'] >= gene_start_in_R_mat-threshold_of_proximity)) & \
                                                       ((reannotated_df['chr2'] == gene_chr) & (reannotated_df['pos2'] <= gene_end_in_R_mat+threshold_of_proximity) & (reannotated_df['pos2'] >= gene_start_in_R_mat-threshold_of_proximity))].copy()
            ### Adding to try to replicate results of previous version of code (v7):
            # num_total_intragenic_events_technically_correct_just_to_replicate_v7_results += reannotated_df[((reannotated_df['chr1'] == gene_chr) & (reannotated_df['pos1'] <= gene_end_in_R_mat+threshold_of_proximity) & (reannotated_df['pos1'] >= gene_start_in_R_mat-threshold_of_proximity)) & \
            #                                            ((reannotated_df['chr2'] == gene_chr) & (reannotated_df['pos2'] <= gene_end_in_R_mat+threshold_of_proximity) & (reannotated_df['pos2'] >= gene_start_in_R_mat-threshold_of_proximity))].shape[0]
            # num_total_intragenic_events_technically_correct_just_to_replicate_v7_results += reannotated_df[((reannotated_df['pos1'] <= gene_end_in_R_mat+threshold_of_proximity) & (reannotated_df['pos1'] >= gene_start_in_R_mat-threshold_of_proximity)) & \
            #                                            ((reannotated_df['pos2'] <= gene_end_in_R_mat+threshold_of_proximity) & (reannotated_df['pos2'] >= gene_start_in_R_mat-threshold_of_proximity))].shape[0]

        else:
            # gene_start_in_R_mat = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]['gene_start'].min()
            # gene_end_in_R_mat = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]['gene_end'].max()
            # Bug fix here (Xavi): Previously I used or instead of and for combining these two booleans 
            reannotated_df_intragenic = reannotated_df[pd.Series(list(map(lambda coords: (is_within_specified_range_of_gene(R_mat_loaded_and_prepared_df,gene,threshold_of_proximity,coords[0],coords[1]) and \
                                                                                          is_within_specified_range_of_gene(R_mat_loaded_and_prepared_df,gene,threshold_of_proximity,coords[2],coords[3])), 
                                                                zip(list(reannotated_df['pos1']),list(reannotated_df['chr1']),list(reannotated_df['pos2']),list(reannotated_df['chr2']))))).astype('bool')].copy()
            ### Adding to try to replicate results of previous version of code (v7):
            # num_total_intragenic_events_technically_correct_just_to_replicate_v7_results += reannotated_df[((reannotated_df['pos1'] <= gene_end_in_R_mat+threshold_of_proximity) & (reannotated_df['pos1'] >= gene_start_in_R_mat-threshold_of_proximity)) & \
            #                                            ((reannotated_df['pos2'] <= gene_end_in_R_mat+threshold_of_proximity) & (reannotated_df['pos2'] >= gene_start_in_R_mat-threshold_of_proximity))].shape[0]
            # num_total_intragenic_events_technically_correct_just_to_replicate_v7_results += reannotated_df[((reannotated_df['chr1'] == gene_chr) & (reannotated_df['pos1'] <= gene_end_in_R_mat+threshold_of_proximity) & (reannotated_df['pos1'] >= gene_start_in_R_mat-threshold_of_proximity)) & \
            #                                            ((reannotated_df['chr2'] == gene_chr) & (reannotated_df['pos2'] <= gene_end_in_R_mat+threshold_of_proximity) & (reannotated_df['pos2'] >= gene_start_in_R_mat-threshold_of_proximity))].shape[0]
            # num_total_intragenic_events_technically_correct_just_to_replicate_v7_results += reannotated_df[pd.Series(list(map(lambda coords: (is_within_specified_range_of_gene(R_mat_loaded_and_prepared_df,gene,threshold_of_proximity,coords[0],coords[1]) or \
            #                                                                                                                                   is_within_specified_range_of_gene(R_mat_loaded_and_prepared_df,gene,threshold_of_proximity,coords[2],coords[3])), 
            #                                                     zip(list(reannotated_df['pos1']),list(reannotated_df['chr1']),list(reannotated_df['pos2']),list(reannotated_df['chr2']))))).astype('bool')].shape[0]
        
        # Another bug that is fixed with this code but wasn't explicitly fixed before:
        # Before (SVGAR_sf v7) I used to not determine which events were 2-breakpoint and 1-breakpoint for a gene appropriately because I wouldn't take into account the chr of the breakpoint.
        # However, I've now fixed it so that I do.
        
        # print("find_proportions_of_deletions_tandemdups_inversions_for_intragenic_all_genes for-loop Time point 3:", time.time()-time_start)
        
        # if len(reannotated_df_intragenic) != len(reannotated_df_intragenic_2):
        #     print("gene:",gene)
        #     print("reannotated_df_intragenic:")
        #     print(reannotated_df_intragenic)
        #     print("reannotated_df_intragenic_2:")
        #     print(reannotated_df_intragenic_2)
        
        num_total_events += reannotated_df.shape[0]
        num_total_intragenic_events += reannotated_df_intragenic.shape[0]
        
        # print("find_proportions_of_deletions_tandemdups_inversions_for_intragenic_all_genes for-loop Time point 4:", time.time()-time_start)
    
        # print("reannotated_df_intragenic:")
        # print(reannotated_df_intragenic)
        
        num_deletions += reannotated_df_intragenic[(reannotated_df_intragenic['str1'] == 0) & (reannotated_df_intragenic['str2'] == 1)].shape[0]
        num_tandemdups += reannotated_df_intragenic[(reannotated_df_intragenic['str1'] == 1) & (reannotated_df_intragenic['str2'] == 0)].shape[0]
        num_inversions += reannotated_df_intragenic[((reannotated_df_intragenic['str1'] == 0) & (reannotated_df_intragenic['str2'] == 0)) | ((reannotated_df_intragenic['str1'] == 1) & (reannotated_df_intragenic['str2'] == 1))].shape[0]
        
        # print("gene:", gene)
        # print("number of intragenic events:", reannotated_df_intragenic.shape[0])
        # print("number of total events:", reannotated_df.shape[0])
        
        # print("find_proportions_of_deletions_tandemdups_inversions_for_intragenic_all_genes for-loop Time point 5:", time.time()-time_start)
    
    # print("find_proportions_of_deletions_tandemdups_inversions_for_intragenic_all_genes Time point 2:", time.time()-time_start)
    
    # print("num_total_events:", num_total_events)
    # print("num_total_intragenic_events:", num_total_intragenic_events)
    
    proportion_intergenic_events = (num_total_events-num_total_intragenic_events)/num_total_events
    proportion_intragenic_events = num_total_intragenic_events/num_total_events
    # proportion_intergenic_events = 0.5429756220155818
    # proportion_intragenic_events = 0.4570243779844182
    
    proportion_deletions = num_deletions/num_total_intragenic_events
    proportion_tandemdups = num_tandemdups/num_total_intragenic_events
    proportion_inversions = num_inversions/num_total_intragenic_events
    
    print("~PARAMETERS USED FOR MODEL~")
    
    print("Proportion of 1-breakpoint events:", proportion_intergenic_events)
    print("Proportion of 2-breakpoint events:", proportion_intragenic_events)
    
    print("Proportion of deletions among two-breakpoint events:", str(proportion_deletions))
    print("Proportion of tandem duplications among two-breakpoint events:", str(proportion_tandemdups))
    print("Proportion of inversions among two-breakpoint events:", str(proportion_inversions))
    
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    
    # print("find_proportions_of_deletions_tandemdups_inversions_for_intragenic_all_genes Time point 3:", time.time()-time_start)
    
    return [proportion_intergenic_events, proportion_intragenic_events, proportion_deletions, proportion_tandemdups, proportion_inversions]
        


def insert_introns_between_exons(data, addend_for_if_not_TBTGFC, gene_name):
    # print("\nData to insert:", data)
    new_data = []
    for i in range(len(data)-1):
        #add exon
        new_data.append([data[i][0], data[i][1], 1 + addend_for_if_not_TBTGFC])

        #add intron between exon and next exon
        new_data.append([data[i][1]+1, data[i+1][0]-1, 2 + addend_for_if_not_TBTGFC])
    new_data.append([data[-1][0], data[-1][1], 1 + addend_for_if_not_TBTGFC])

    # print("Result of insert:\n", new_data)
    # if gene_name == "SH3YL1": print("intronned new data", new_data)
    return new_data

# for [[a, b], [c, d], ...] lists
def slicing_start_end_lists(low, high, data, gene_name):
    # print("\nData to slice:", low, high, data)
    low_found = False

    index_cut_off_for_low = 0
    index_cut_off_for_high = 0

    for i in range(len(data)):
        if not low_found:
            if data[i][1] >= low:
                index_cut_off_for_low = i
                # print("got low", index_cut_off_for_low)
                low_found = True
        # keep setting high until you can't
        if data[i][0] <= high:
            index_cut_off_for_high = i
            # print("got high", index_cut_off_for_high)

    # print("low", index_cut_off_for_low)
    # print("high", index_cut_off_for_high)

    new_data = deepcopy(data[index_cut_off_for_low : index_cut_off_for_high + 1])

    # if gene_name == "SH3YL1": print("newdata", new_data)
    new_data[0][0] = low
    new_data[-1][1] = high
    # if gene_name == "SH3YL1": print("newdata", new_data)
    
    # print("Result of slice:\n", new_data)
    return new_data

def converts_rmat_file_into_gene_format(R_mat_loaded_and_prepared_df, max_transcripts):
    #format: will always be inclusive start, exclusive end. so I can have [0, 10], [10, 20], etc

    
    
    # R_mat_loaded_and_prepared_df.to_csv("rmat.tsv", sep = "\t")
    # gene_format_df = pd.DataFrame({"transcript" : [], "data" : [], "name" : [], "chr" : [], "gene_start" : [], "gene_end" : []})
    chr_transcripts_dict = {}


    print("Finding CGC Genes:")
    # get cgc list
    list_of_cgc_driver_transloc_genes = []
    with open("../reference_files/CGC_gene_list/COSMIC_Cancer_Gene_Census_genes_to_matched_Hugo_symbols_Dec_13_2021_only_translocations.txt") as list_of_driver_transloc_file:
        for line in list_of_driver_transloc_file:
            list_of_cgc_driver_transloc_genes.append(line.strip())

    #exon locations are inclusive for both sides! [5] [10] means 1 exon with a length of 6.
    
    transcript_count_for_displaying_progress = 0
    transcripts_with_nan = []

    gene_chr_start_end_strand = {}


    print("Converting rows:")
    for index, row in tqdm(R_mat_loaded_and_prepared_df.iterrows()):
        
        transcript_count_for_displaying_progress += 1

        #quit early
        if transcript_count_for_displaying_progress > max_transcripts: break

        # if transcript_count_for_displaying_progress % 1000 == 0: print(transcript_count_for_displaying_progress)

        gene_name = row["gene"]
        transcript_name = row["transcript"]
        gene_chr = row["chr"]
        gene_strand = row["strand"]
        

        # print("Running transcript", transcript_name)
        exon_data = list(zip(row["exon_starts"], row["exon_ends"]))
        for i in range(len(exon_data)): exon_data[i] = list(exon_data[i])

        # print("\nExon data:", exon_data)
        
                                    # so lowest number = highest priority
        addend_for_if_not_TBTGFC = 0 # zone: 1=TBTGFC+exon, 2=TBTGFC+intron, 3=TBTGFC+3'-UTR, 4=TBTGFC+5'-UTR, 5=TBTGFC+promoter, 6=exon, 7=intron, 8=3'-UTR, 9=5'-UTR, 10=promoter
        if not gene_name in list_of_cgc_driver_transloc_genes:
            addend_for_if_not_TBTGFC = 5 #changed from 0 to 5

        # concated_row_with_new_transcript = [row[transcript], 0, 0, ]
        
        transcript_new_data = []


        # print("type", type(row["code_start"]))
        # print("bool", bool(row["code_start"]))
        # print("it", row["code_start"])
        # print(pd.notnull(row["code_start"]))
        # print("name -", transcript_name)
        #start with just positive case, will copy over for negative after
        if row["strand"] == "+":
            # if not (pd.notnull(row["gene_start"]) and pd.notnull(row["tx_start"]) and pd.notnull(row["code_start"]) and pd.notnull(row["code_end"]) and pd.notnull(row["tx_end"])): 
            #     transcripts_with_nan.append(row["transcript"])
                
            #     print("NaN found")
            #     continue #if nan values

            gene_start = int(row["gene_start"])
            tx_start = int(row["tx_start"])
            
            if pd.notnull(row["code_start"]):
                code_start = int(row["code_start"])
            else: #just applies to ROBO2
                code_start = int(row["tx_start"])
            

            code_end = int(row["code_end"])
            tx_and_gene_end = int(row["tx_end"])



            #ASSUMING that inclusive, higher priority is more inside, so code: 2000-3000, tx: 1000, 4000 means that whatever corresponds to code-tx is 1000-1999, 3001-4000. also matches priority system.
            #promoter
            transcript_new_data = [gene_name, [gene_start, tx_start - 1, 5 + addend_for_if_not_TBTGFC]]
            if tx_start != code_start: # in case there is no 5'UTR
                transcript_new_data.append([tx_start, code_start - 1, 4 + addend_for_if_not_TBTGFC])
            else:
                pass
                # no 5'UTR

            transcript_new_data.extend(insert_introns_between_exons(slicing_start_end_lists(code_start, code_end, exon_data, gene_name), addend_for_if_not_TBTGFC, gene_name))

            if code_end != tx_and_gene_end: # in case there is no 3'UTR
                transcript_new_data.append([code_end + 1, tx_and_gene_end, 3 + addend_for_if_not_TBTGFC])
            else:
                pass
                # no 3'UTR
        # strand is -
                
        else:
            if not (pd.notnull(row["tx_start"]) and pd.notnull(row["code_start"]) and pd.notnull(row["code_end"]) and pd.notnull(row["tx_end"]) and pd.notnull(row["gene_end"])):
                transcripts_with_nan.append(row["transcript"])
                print("NaN found")
                continue #if nan values

            tx_and_gene_start = int(row["tx_start"])
            code_start = int(row["code_start"])

            code_end = int(row["code_end"])
            tx_end = int(row["tx_end"])
            gene_end = int(row["gene_end"])

        

            #ASSUMING that inclusive, higher priority is more inside, so code: 2000-3000, tx: 1000, 4000 means that whatever corresponds to code-tx is 1000-1999, 3001-4000. also matches priority system.
            transcript_new_data = [gene_name]
            if tx_and_gene_start != code_start: # in case there is no 3'UTR
                transcript_new_data.append([tx_and_gene_start, code_start - 1, 3 + addend_for_if_not_TBTGFC])
            else:
                pass
                # no 3'UTR

            transcript_new_data.extend(insert_introns_between_exons(slicing_start_end_lists(code_start, code_end, exon_data, gene_name), addend_for_if_not_TBTGFC, gene_name))

            if code_end != tx_end: # in case there is no 5'UTR
                transcript_new_data.append([code_end + 1, tx_end, 4 + addend_for_if_not_TBTGFC])
            else:
                pass
                # no 5'UTR

            #promoter
            transcript_new_data.append([tx_end + 1, gene_end, 5 + addend_for_if_not_TBTGFC])

        # print(transcript_new_data)
        if gene_chr not in chr_transcripts_dict.keys():
            chr_transcripts_dict[gene_chr] = [transcript_new_data]
        else:
            chr_transcripts_dict[gene_chr].append(transcript_new_data)


        # R_mat_loaded_and_prepared_df_subset = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene_name]
        # chrs_for_gene = set(R_mat_loaded_and_prepared_df_subset['chr'])
        # multi_chr_gene = bool(len(chrs_for_gene) > 1)

        # #noting where each gene is:
        if gene_name not in gene_chr_start_end_strand:
            gene_chr_start_end_strand[gene_name] = [gene_chr, int(row["gene_start"]), int(row["gene_end"]), gene_strand]
        
        #just avoids pseudoautosomal case on diff (X vs Y) chromosomes. We're only saving one tx, but could end up having [gene_start, gene_end] = [55mil, 157mil] which is a problem
        elif abs(gene_chr_start_end_strand[gene_name][1] - int(row["gene_start"])) < 10000000:
                gene_chr_start_end_strand[gene_name] = [gene_chr, min(gene_chr_start_end_strand[gene_name][1], int(row["gene_start"])), max(gene_chr_start_end_strand[gene_name][2], int(row["gene_end"])), gene_strand]
            #same, but w being in Y chr and setting the max
    
                
        

    #                                       first exon   first intron
    #format for chr_transcripts_dict: {1: [[[0, 100, 1], [101, 200, 2]], first transcript
    #                                      [[5, 99, 1], [100, 207, 2]]]} second transcript

    for a_chr in chr_transcripts_dict.keys():

        chr_transcripts_dict[a_chr].sort(key = lambda x : x[1][0])
        # print(a_chr)
    
    print("has nan", transcripts_with_nan)
    # print(gene_chr_start_end_strand["IL9R"])
    return gene_chr_start_end_strand, chr_transcripts_dict


        # lst.sort(key=lambda x: x[2])

    # for gene in gene_different_transcripts_dict.keys():

def get_transcripts_in_range(gene_name, gene_start, gene_end, list_of_transcripts):
    transcripts_in_range = []

    # print(gene_start, gene_end, list_of_transcripts)
    #add every transcript where any part of it is in the range - so either min or max is between gene start and end
    for i in range(len(list_of_transcripts)):
        if gene_name == list_of_transcripts[i][0] and (gene_start <= list_of_transcripts[i][1][0] <= gene_end or gene_start <= list_of_transcripts[i][-1][1] <= gene_end):
            # print(list_of_transcripts[i], i)
            transcripts_in_range.append(list_of_transcripts[i])

    # if gene_name == "SH3YL1": print("\n\nTRANSCRIPTS:\n\n", transcripts_in_range)
        #transcripts are sorted by this index, so if the gene ends before the start of this transcript, none of the rest will either.
        # if gene_end < list_of_transcripts[i][0][0]: break
    return transcripts_in_range

def slen(gene_section):
    return gene_section[1]-gene_section[0]+1

def consolidate_merged(merged_transcript):
    new_merged = [merged_transcript[0]]
    for i in range(1, len(merged_transcript)):
        #push back end
        if slen(merged_transcript[i]) == 0: continue
        
        #if intergenic region: (2 completely non-overlapping transcripts, should be half "before gene" half "after gene", 11 and 12 respectively) 
        if merged_transcript[i][0] > new_merged[-1][1] + 1:
            avg_of_intergenic_region = ((new_merged[-1][1] + 1) + (merged_transcript[i][0] - 1)) // 2
            
            new_merged.append([new_merged[-1][1] + 1, avg_of_intergenic_region, 12]) # because it's after the first tx, then before the second tx. so after_gene precedes before_gene
            new_merged.append([avg_of_intergenic_region + 1, merged_transcript[i][0] - 1, 11])

        if merged_transcript[i][2] == new_merged[-1][2]:
            new_merged[-1][1] = merged_transcript[i][1]
        else:
            new_merged.append(merged_transcript[i])

    #fix promoter/intergenic region situation
    for i in range(1, len(new_merged)):
        if new_merged[i][2] == 10:
            promoter_length = new_merged[i][1] - new_merged[i][0]
            if promoter_length > 3000: #then it's overwriting the intergenic region, half before_gene half after_gene. 11 and 12 respectively, lowest priority
                new_merged[i][1] = new_merged[i][0] + 2999
                new_merged.insert(i+1, [new_merged[i][0] + 3000, new_merged[i][0] + 3000 + (promoter_length-3000)//2, 11])
                new_merged.insert(i+2, [new_merged[i][0] + 3000 + (promoter_length-3000)//2 + 1, new_merged[i][0] + promoter_length, 12])
    return new_merged

def find_merged_transcript(transcripts_in_range, gene_name):

    # print(transcripts_in_range)
    merged_transcript = deepcopy(transcripts_in_range[0][1:]) #have start. Don't need the gene name anymore

    for i in range(1, len(transcripts_in_range)): #max like a couple dozen (t) transcripts
        # print("currently:", merged_transcript)
        # print("Length of transcript", i, len(transcripts_in_range[i]))
        # print("transcript to add", transcripts_in_range[i])
        for j in range(1, len(transcripts_in_range[i])): # n sections, max like 300 but usually around 10
            # print("transcript section to add", transcripts_in_range[i][j])
            
            lower_limit_index_merged = 0
            upper_limit_index_merged = 0
            #go up to index of merged transcript where it actually matters 

            #go until the start of the merged index is above the start of the section, then go one earlier. so until the merged max is above or equal to the current min. so while it's less.
            while lower_limit_index_merged < len(merged_transcript)-1 and merged_transcript[lower_limit_index_merged][1] <= transcripts_in_range[i][j][0]:
                lower_limit_index_merged += 1

            # print("lower limit", lower_limit_index_merged, merged_transcript[lower_limit_index_merged])
            upper_limit_index_merged = lower_limit_index_merged

            # print("math for upper", merged_transcript[upper_limit_index_merged][0], transcripts_in_range[i][j][1])

            #until the next merged loc's low is above the current transcript's high. so while it's at or below
            while upper_limit_index_merged < len(merged_transcript) and merged_transcript[upper_limit_index_merged][0] <= transcripts_in_range[i][j][1]:
                upper_limit_index_merged += 1
                
            #because we stopped when it stopped being the case, the upper limit must be exclusive

            # print("upper limit", upper_limit_index_merged, merged_transcript[upper_limit_index_merged])
            


            priority_of_current_loc_in_current_transcript = transcripts_in_range[i][j][2]
            # print("prio of current loc", priority_of_current_loc_in_current_transcript)


            #now "paste" transcript on merged. go in reverse so i can insert.
            for k in range(upper_limit_index_merged - 1, lower_limit_index_merged - 1, -1):
                # print("index of merged k", k)
                #merged is fully contained - can be overwritten
                priority_of_current_loc_in_merged = merged_transcript[k][2]
                # print("prio of merged loc", priority_of_current_loc_in_merged)
                # print("len of merged loc", merged_transcript[k][0], "to", merged_transcript[k][1], merged_transcript[k][1] - merged_transcript[k][0])
                
                
                if priority_of_current_loc_in_current_transcript < priority_of_current_loc_in_merged:
                    # print("comparing merged tx", merged_transcript[k], "with current tx", transcripts_in_range[i][j])
                    # print("overwrite!")
                    #merged section is fully in transcript section - change to [new]
                    if transcripts_in_range[i][j][0] <= merged_transcript[k][0] and merged_transcript[k][1] <= transcripts_in_range[i][j][1]:
                        # print("merged is fully in current")
                        #just swap
                        merged_transcript[k][2] = priority_of_current_loc_in_current_transcript 
                        # if slen(merged_transcript[k]) <= 0: print("produced section of len", slen(merged_transcript[k]))
                        # print("result", merged_transcript[k])
                    
                    #change to [old, new, old]
                    #transcript section is fully in merged section, such that 3 sections need to be produced (if it's just 2, then it's "to the left/right". [[  5  ]  7  ] is to the left, [ 7  [  5 ] the rest of the 7  ] is fully in)
                    elif merged_transcript[k][0] < transcripts_in_range[i][j][0] and transcripts_in_range[i][j][1] < merged_transcript[k][1]:
                        # print("current is fully in merged")
                        ### can find a better way to do it if it's too slow, to avoid "insert"
                        #unchanged earlier part - bring max lower

                        save_location = merged_transcript[k][1]
                        merged_transcript[k][1] = transcripts_in_range[i][j][0] - 1
                        # print("new k:k+3", merged_transcript[k:k+3])

                        #new middle part
                        merged_transcript.insert(k+1, transcripts_in_range[i][j])
                        # print("new k:k+3", merged_transcript[k:k+3])
                        #new after part

                        merged_transcript.insert(k+2, [transcripts_in_range[i][j][1] + 1, save_location, merged_transcript[k][2]])
                        # print("new k:k+3", merged_transcript[k:k+3])
                        # if slen(merged_transcript[k]) <= 0: print("produced section k of len", slen(merged_transcript[k]))
                        # if slen(merged_transcript[k+1]) <= 0: print("produced section k+1 of len", slen(merged_transcript[k+1]))
                        # if slen(merged_transcript[k+2]) <= 0: print("produced section k+2 of len", slen(merged_transcript[k+2]))
                        # print("result", merged_transcript[k:k+3])
                    
                    #change to new, old. so old has extra to the right
                    #current transcript section is just to the left of the merged section
                    elif transcripts_in_range[i][j][1] < merged_transcript[k][1]:
                        # print("current is partially to the left of merged")
                        save_location = merged_transcript[k][0]
                        merged_transcript[k][0] = transcripts_in_range[i][j][1] + 1

                        #add merged part
                        merged_transcript.insert(k, [save_location, transcripts_in_range[i][j][1], transcripts_in_range[i][j][2]])
                        # print("result", merged_transcript[k:k+2])
                        # if slen(merged_transcript[k]) <= 0: print("produced section k of len", slen(merged_transcript[k]))
                        # if slen(merged_transcript[k+1]) <= 0: print("produced section k+1 of len", slen(merged_transcript[k+1]))

                    #change to old, new. so old has extra to the left.
                    #current transcript section is just to the right of the merged section
                    elif transcripts_in_range[i][j][0] > merged_transcript[k][0]:
                        # print("current is partially to the right of merged")
                        save_location = merged_transcript[k][1]
                        merged_transcript[k][1] = transcripts_in_range[i][j][0] - 1

                        #add merged part
                        merged_transcript.insert(k+1, [transcripts_in_range[i][j][0], save_location, transcripts_in_range[i][j][2]])
                        
                        # if slen(merged_transcript[k]) <= 0: print("produced section k of len", slen(merged_transcript[k]))
                        # if slen(merged_transcript[k+1]) <= 0: print("produced section k+1 of len", slen(merged_transcript[k+1]))
                        # print("result", merged_transcript[k:k+2])
                    else: print("Merging transcript - None of the above? This shouldn't happen")
                # else: print("don't overwrite")
                # if merged_transcript[k][0] > merged_transcript[k][1]: print("Issue in merged! Negative length!", k, merged_transcript[k])
            
            #for if the transcript section starts before the merged transcript does, add that at the start
            if transcripts_in_range[i][j][0] < merged_transcript[0][0]:
                # print('sth before!') 
                # print("inserting", [transcripts_in_range[i][j][0], merged_transcript[0][0]-1, transcripts_in_range[i][j][2]])
                
                # all the way before, for example if there's 2 completely non-overlapping transcripts. then be itself, leave a gap.
                if False and transcripts_in_range[i][j][1] < merged_transcript[0][0]:
                    merged_transcript.insert(0, [transcripts_in_range[i][j][0], transcripts_in_range[i][j][1], transcripts_in_range[i][j][2]])
                else:
                    merged_transcript.insert(0, [transcripts_in_range[i][j][0], merged_transcript[0][0]-1, transcripts_in_range[i][j][2]])
                # print(merged_transcript[:2])
            #for if the transcript section ends after the merged transcript does, add that at the end
            if transcripts_in_range[i][j][1] > merged_transcript[-1][1]:
                # print('sth after!')
                # print("appending", [merged_transcript[-1][1]+1, transcripts_in_range[i][j][1], transcripts_in_range[i][j][2]])
                
                # all the way after, for example if there's 2 completely non-overlapping transcripts. then be itself, leave a gap.
                if False and transcripts_in_range[i][j][0] > merged_transcript[-1][1]:
                    merged_transcript.append([transcripts_in_range[i][j][0], transcripts_in_range[i][j][1], transcripts_in_range[i][j][2]])
                else:
                    merged_transcript.append([merged_transcript[-1][1]+1, transcripts_in_range[i][j][1], transcripts_in_range[i][j][2]])
                # print(merged_transcript[-2:])
    # print("merged", merged_transcript)

    # print(merged_transcript)
    consolidated_merged_transcript = consolidate_merged(merged_transcript)
    # print(consolidated_merged_transcript)
    # if gene_name == "SH3YL1": print("CONSOLIDATED", consolidated_merged_transcript)
    return consolidated_merged_transcript

            #add to merged transcript



    # index_of_each_transcript = [] #start with all 0s
    # for i in range(len(transcripts_in_range)): index_of_each_transcript.append(0)

    # #keep going until at the last index of all transcripts - subtract len(transcripts_in_range) because the index can be max len-1 for each transcript
    
    # while sum(index_of_each_transcript) != sum([len(ts) for ts in transcripts_in_range]) - len(transcripts_in_range):
    #     index_in_index_list_with_lowest_value = 0
    #     lowest_val = 1000000000

    #     for i in range(len(transcripts_in_range)):
    #         if transcripts_in_range[i][index_in_index_list_with_lowest_value[i]] < lowest_val:
        
    #     find_min_with_list_of_index(index_of_each_transcript, transcripts_in_range)


def write_file_using_merged_transcript(pre_annotations_filename, merged_transcript_filename, gene_name, gene_start, gene_end, gene_strand, merged_transcript, R_mat_loaded_and_prepared_df):

    R_mat_loaded_and_prepared_df_subset = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene_name]
    chrs_for_gene = set(R_mat_loaded_and_prepared_df_subset['chr'])
    multi_chr_gene = bool(len(chrs_for_gene) > 1)

    column_names = ["exon","3'UTR","5'UTR","promoter","before_gene","after_gene"]
    #get the before/after gene data right
    #exon, 3'UTR, 5'UTR, promoter, before_gene, after_gene, and then the introns after
    column_data = [0, 0, 0, 0, merged_transcript[0][0] - gene_start, gene_end - merged_transcript[-1][1]]

    # print("start and end", gene_start, gene_end)
    # print("range", merged_transcript[0][0], merged_transcript[-1][1])

    intron_data = []
    exon_data = []

    #adding up sections
    for i in range(len(merged_transcript)):
        if merged_transcript[i][2] in [1, 6]:
            # print(merged_transcript[i][1] - merged_transcript[i][0] + 1)
            column_data[0] += merged_transcript[i][1] - merged_transcript[i][0] + 1
        elif merged_transcript[i][2] in [2, 7]:
            # print(merged_transcript[i][1] - merged_transcript[i][0] + 1)
            intron_data.append(merged_transcript[i][1] - merged_transcript[i][0] + 1)
        elif merged_transcript[i][2] in [3, 8]:
            # print(merged_transcript[i][1] - merged_transcript[i][0] + 1)
            column_data[1] += merged_transcript[i][1] - merged_transcript[i][0] + 1
            if multi_chr_gene: column_data[1] *= 2
        elif merged_transcript[i][2] in [4, 9]:
            # print(merged_transcript[i][1] - merged_transcript[i][0] + 1)
            column_data[2] += merged_transcript[i][1] - merged_transcript[i][0] + 1
            if multi_chr_gene: column_data[2] *= 2
        elif merged_transcript[i][2] in [5, 10]:
            # print(merged_transcript[i][1] - merged_transcript[i][0] + 1)
            column_data[3] += merged_transcript[i][1] - merged_transcript[i][0] + 1
            if multi_chr_gene: column_data[3] *= 2
        elif merged_transcript[i][2] == 11:
            # print("adding11", merged_transcript[i][1] - merged_transcript[i][0] + 1)
            column_data[4] += merged_transcript[i][1] - merged_transcript[i][0] + 1
            if multi_chr_gene: column_data[4] *= 2
        elif merged_transcript[i][2] == 12:
            # print("adding12", merged_transcript[i][1] - merged_transcript[i][0] + 1)
            column_data[5] += merged_transcript[i][1] - merged_transcript[i][0] + 1
            if multi_chr_gene: column_data[5] *= 2
        else:
            raise("invalid section")

    intron_column_names = ["intron" + str(x) for x in range(1, len(intron_data) + 1)]

    count_of_exons = 0
    for a_section in merged_transcript:
        if a_section[2] in [1, 6]:
            count_of_exons += 1
    
    current_exon_number = 1
    count_up_or_down = 1
    if gene_strand == "-":
        current_exon_number = count_of_exons
        count_up_or_down = -1

    flat_merged_tx = []
    for a_section in merged_transcript:
        if a_section[2] in [1, 6]: #exon
            flat_merged_tx.extend([str(a_section[0]), str(a_section[1]), "exon" + str(current_exon_number)])
        elif a_section[2] in [2, 7]: #intron
            flat_merged_tx.extend([str(a_section[0]), str(a_section[1]), "intron" + str(current_exon_number)])
            current_exon_number += count_up_or_down
        elif a_section[2] in [3, 8]: #5'UTR
            flat_merged_tx.extend([str(a_section[0]), str(a_section[1]), "5'UTR"])
        elif a_section[2] in [4, 9]: #3'UTR
            flat_merged_tx.extend([str(a_section[0]), str(a_section[1]), "3'UTR"])
        elif a_section[2] in [5, 10]: #promoter
            flat_merged_tx.extend([str(a_section[0]), str(a_section[1]), "promoter"])
        elif a_section[2] == 11:
            flat_merged_tx.extend([str(a_section[0]), str(a_section[1]), "before_gene"])
        elif a_section[2] == 12:
            flat_merged_tx.extend([str(a_section[0]), str(a_section[1]), "after_gene"])        
        else:
            raise("invalid gene section type")

    with open(merged_transcript_filename, 'w') as f:
        f.write("\t".join(flat_merged_tx))
        f.write("\nnot on COSMIC list" if int(merged_transcript[0][2])>5 else "\non COSMIC list") #check if known driver
    # flat_merged_tx.to_csv(merged_transcript_filename)

    #prep merged tx columns:
    if gene_strand == "-":
        intron_column_names = intron_column_names[::-1]
        
    # elif gene_strand == "-":
    #     #added 8/1/24
    #     intron_column_names = ["intron" + str(x) for x in range(len(intron_data), 0)]

    


    
    
    column_names.extend(intron_column_names)
    column_data.extend(intron_data)

    df = pd.DataFrame(columns = column_names)
    df.loc[len(df.index)] = column_data    

    df.to_csv(pre_annotations_filename, sep = "\t", index = False)
    # print(merged_transcript)
    # print(df)


def produce_annotations_from_r_mat_file(max_transcripts = 10000000, selected_genes = []):
    OVERALL_START = time.time()

    if not os.path.isdir('../surrounding_annotations_for_each_gene/v14/pre_annotations'):
        os.makedirs('../surrounding_annotations_for_each_gene/v14/pre_annotations')
    if not os.path.isdir('../surrounding_annotations_for_each_gene/v14/merged_transcript'):
        os.makedirs('../surrounding_annotations_for_each_gene/v14/merged_transcript')
    R_mat_file_path = R_MAT_FILE_PATH
    # R_mat_loaded_and_prepared_df = load_and_prepare_R_mat_file_into_df(R_mat_file_path)
    print("Loading RMAT File")
    R_mat_loaded_and_prepared_df = load_and_prepare_R_mat_file_into_df(R_mat_file_path)

    gene_chr_start_end_strand, chr_transcripts_dict = converts_rmat_file_into_gene_format(R_mat_loaded_and_prepared_df, max_transcripts)

    pre_annotations_foldername = "../surrounding_annotations_for_each_gene/v14/pre_annotations/"
    merged_transcript_foldername = "../surrounding_annotations_for_each_gene/v14/merged_transcript"
    
    # amount_of_genes_annotated = 0

    for a_gene in tqdm(selected_genes if len(selected_genes) else gene_chr_start_end_strand.keys()):
        pre_annotations_filename = pre_annotations_foldername + "gene_" + a_gene + ".step_size_1bp.pre_annotations.tsv"
        merged_transcript_filename = merged_transcript_foldername + "/gene_" + a_gene + ".step_size_1bp.merged_transcript.tsv"
        # filename = "../surrounding_bugfix_10bp_annotations_for_each_gene/v7/pre_annotations/gene_" + a_gene + ".step_size_1bp.bp.pre_annotations.tsv"

        # print("gene start end", gene_chr_start_end_strand_n_exons[a_gene])
        # print(a_gene)
        # amount_of_genes_annotated += 1
        # if amount_of_genes_annotated % 500 == 0:
        #     print("annotation progress:", amount_of_genes_annotated, "genes")
        gene_chr, gene_start, gene_end, gene_strand = gene_chr_start_end_strand[a_gene]

        transcripts_in_range = get_transcripts_in_range(a_gene, gene_start, gene_end, chr_transcripts_dict[gene_chr])
        # if a_gene == "UBE2D2": print("\ntranscripts for", a_gene); 
        # pprint(transcripts_in_range)
        merged_transcript = find_merged_transcript(transcripts_in_range, a_gene)
        # print("MERGED", merged_transcript)
        
        write_file_using_merged_transcript(pre_annotations_filename, merged_transcript_filename, a_gene, gene_start, gene_end, gene_strand, merged_transcript, R_mat_loaded_and_prepared_df)
        
    print("total time:", time.time() - OVERALL_START)

# This function basically finds how a gene/the region near a gene would be annotated every X (specified by step_size) bp. A resulting file
# gets saved. This file is used to find probability of LOF for a given gene. 
# def find_surrounding_annotations_for_gene_with_parameters_pre(input_gene_name, step_size):
#     if not os.path.isdir('../surrounding_annotations_for_each_gene/v14/pre_annotations'):
#         os.makedirs('../surrounding_annotations_for_each_gene/v14/pre_annotations')
#     file_name_based_off_input_parameters = "../surrounding_annotations_for_each_gene/v14/pre_annotations/gene_" + input_gene_name + ".step_size_" + str(step_size) + "bp.pre_annotations.tsv"
#     # check if file exists:
#     if os.path.isfile(file_name_based_off_input_parameters):
#         pass
#     else:
#         sizes_of_each_chromosome_hg38 = {1: 248956422,
#                                         2: 242193529,
#                                         3: 198295559,
#                                         4: 190214555,
#                                         5: 181538259,
#                                         6: 170805979,
#                                         7: 159345973,
#                                         8: 145138636,
#                                         9: 138394717,
#                                         10: 133797422,
#                                         11: 135086622,
#                                         12: 133275309,
#                                         13: 114364328,
#                                         14: 107043718,
#                                         15: 101991189,
#                                         16: 90338345,
#                                         17: 83257441,
#                                         18: 80373285,
#                                         19: 58617616,
#                                         20: 64444167,
#                                         21: 46709983,
#                                         22: 50818468,
#                                         23: 156040895, # chrX
#                                         24: 57227415} # chrY
#         R_mat_file_path = R_MAT_FILE_PATH
#         R_mat_loaded_and_prepared_df = load_and_prepare_R_mat_file_into_df(R_mat_file_path)
#         if R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == input_gene_name].empty:
#             raise("Error: given gene_name is not found within R.mat file!")
#         else:
#             # gene_chr = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == input_gene_name]['chr'].iloc[0]
#             # gene_start = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == input_gene_name]['gene_start'].min()
#             # gene_end = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == input_gene_name]['gene_end'].max()
#             # Making there be multiple gene_chr's, gene_starts, and gene_ends to account for pseudoautosomal genes:
#             gene_chrs = sorted(list(set(R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == input_gene_name]['chr'])))
#             gene_starts = list(map(lambda c: R_mat_loaded_and_prepared_df[(R_mat_loaded_and_prepared_df['gene'] == input_gene_name) & (R_mat_loaded_and_prepared_df['chr'] == c)]['gene_start'].min(), gene_chrs))
#             gene_ends = list(map(lambda c: R_mat_loaded_and_prepared_df[(R_mat_loaded_and_prepared_df['gene'] == input_gene_name) & (R_mat_loaded_and_prepared_df['chr'] == c)]['gene_end'].max(), gene_chrs))
#         starts_of_annotation_search = []
#         for i in range(0,len(gene_chrs)):
#             if gene_starts[i] < 0:
#                 raise("We have a negative gene_start; this shouldn't be possible")
#             else:
#                 starts_of_annotation_search.append(gene_starts[i])
#         ends_of_annotation_search = []
#         for i in range(0,len(gene_chrs)):
#             if gene_ends[i] > sizes_of_each_chromosome_hg38[gene_chrs[i]]:
#                 raise("We have a gene_end that exceeds the size of the chromosome; this shouldn't be possible")
#             else:
#                 ends_of_annotation_search.append(gene_ends[i])
#         # if gene_start < 0:
#         #     raise("Gene start is negative; this shouldn't be possible")
#         # else:
#         #     start_of_annotation_search = gene_start
#         # if gene_end > sizes_of_each_chromosome_hg38[gene_chr]:
#         #     raise("Gene end exceeds the size of the chromosome; this shouldn't be possible")
#         # else:
#         #     end_of_annotation_search = gene_end
#         # Loading in the file for Cancer Gene Census genes (only translocations) and making a list out of the genes on it
#         list_of_cgc_driver_transloc_genes = []
#         with open("../reference_files/CGC_gene_list/COSMIC_Cancer_Gene_Census_genes_to_matched_Hugo_symbols_Dec_13_2021_only_translocations.txt") as list_of_driver_transloc_file:
#             for line in list_of_driver_transloc_file:
#                 list_of_cgc_driver_transloc_genes.append(line.strip())
#         # R_mat_loaded_and_prepared_df_subset_chr = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['chr'] == gene_chr]
#         R_mat_loaded_and_prepared_df_subset_chr = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['chr'].isin(gene_chrs)]
#         counts_of_annotations_for_each_point_df = pd.DataFrame(np.zeros((1,6),dtype=int))
#         counts_of_annotations_for_each_point_df.index = [input_gene_name]
#         counts_of_annotations_for_each_point_df.columns = ["exon","3'UTR","5'UTR","promoter","before_gene","after_gene"]
#         for gene_chr,start_of_annotation_search,end_of_annotation_search in zip(gene_chrs,starts_of_annotation_search,ends_of_annotation_search):
#             for pos in range(start_of_annotation_search, end_of_annotation_search, step_size):
#                 if (pos - start_of_annotation_search) % 100000 == 0:
#                     print('chr: ' + str(gene_chr) + ', pos: ' + str(pos))
#                 # # This was the original overlaps_df where we're not subsetting based off the name of the gene:
#                 # overlaps_df = R_mat_loaded_and_prepared_df_subset_chr[(R_mat_loaded_and_prepared_df_subset_chr['gene_start'] <= pos) & (R_mat_loaded_and_prepared_df_subset_chr['gene_end'] >= pos)]
#                 # This is the new overlaps_df where we now subset based off the name of the gene. We want to annotate based only off the gene's transcripts:
#                 overlaps_df = R_mat_loaded_and_prepared_df_subset_chr[(R_mat_loaded_and_prepared_df_subset_chr['gene'] == input_gene_name) & (R_mat_loaded_and_prepared_df_subset_chr['gene_start'] <= pos) & (R_mat_loaded_and_prepared_df_subset_chr['gene_end'] >= pos)]
#                 gene_orientation = R_mat_loaded_and_prepared_df_subset_chr[R_mat_loaded_and_prepared_df_subset_chr['gene'] == input_gene_name]['strand'].iat[0]
#                 # it's within a transcript: determine promoter/UTR/intron/exon
#                 no = overlaps_df.shape[0]
#                 if not no > 0:
#                     # This should happen only if there is a region within a gene that does not correspond to transcript/promoter of a transcript
#                     # This would happen if there is a region between the gene_start's and gene_end's for a given gene, and this does happen sometimes
#                     # In this case, I am going to annotate like it were before or after the gene relative to the nearest transcript
#                     dist_before = R_mat_loaded_and_prepared_df_subset_chr['tx_start'].subtract(pos).where(lambda x: x >= 0).dropna()
#                     dist_after = (-R_mat_loaded_and_prepared_df_subset_chr['tx_end']).add(pos).where(lambda x: x >= 0).dropna()
#                     if dist_before.size > 0 and dist_after.size > 0:
#                         if dist_before.min() <= dist_after.min():
#                             if gene_orientation == '+':
#                                 counts_of_annotations_for_each_point_df.at[input_gene_name,"before_gene"] += 1
#                             else:
#                                 counts_of_annotations_for_each_point_df.at[input_gene_name,"after_gene"] += 1
#                         else:
#                             if gene_orientation == '+':
#                                 counts_of_annotations_for_each_point_df.at[input_gene_name,"after_gene"] += 1
#                             else:
#                                 counts_of_annotations_for_each_point_df.at[input_gene_name,"before_gene"] += 1
#                     else:
#                         raise("Something went wrong. There should be valid dist_before and dist_after minimums.")
#                 else:
#                     # Note: TBTGFC means "Transcript belongs to gene from COSMIC"
#                     c = np.full(no, np.nan) # zone: 1=TBTGFC+exon, 2=TBTGFC+intron, 3=TBTGFC+3'-UTR, 4=TBTGFC+5'-UTR, 5=TBTGFC+promoter, 6=exon, 7=intron, 8=3'-UTR, 9=5'-UTR, 10=promoter
#                     d = np.full(no, np.nan) # for exons: which one, and how far
#                     d1 = np.full(no, np.nan) # for introns: between which exons and how far?
#                     d2 = np.full(no, np.nan) # for introns: between which exons and how far?
#                     e1 = np.full(no, np.nan) # for introns: between which exons and how far?
#                     for j in range(0,no):
#                         i = overlaps_df.index[j]
#                         gene_name = R_mat_loaded_and_prepared_df.at[i,'gene']
#                         # is contained in the list of COSMIC driver genes with differing transcripts?
#                         addend_for_if_not_TBTGFC = 0
#                         if not gene_name in list_of_cgc_driver_transloc_genes:
#                             addend_for_if_not_TBTGFC = 5
#                         # in promoter?
#                         if pos < R_mat_loaded_and_prepared_df.at[i,'tx_start']:
#                             c[j] = 5 + addend_for_if_not_TBTGFC
#                             d[j] = R_mat_loaded_and_prepared_df.at[i,'tx_start'] - pos
#                             continue
#                         if pos > R_mat_loaded_and_prepared_df.at[i,'tx_end']:
#                             c[j] = 5 + addend_for_if_not_TBTGFC
#                             d[j] = pos - R_mat_loaded_and_prepared_df.at[i,'tx_end']
#                             continue
#                         # in UTR?
#                         if R_mat_loaded_and_prepared_df.at[i,'strand'] == '+':
#                             if R_mat_loaded_and_prepared_df.at[i,'code_start'] > pos:
#                                 c[j] = 4 + addend_for_if_not_TBTGFC
#                                 d[j] = R_mat_loaded_and_prepared_df.at[i,'code_start'] - pos
#                                 continue
#                             if R_mat_loaded_and_prepared_df.at[i,'code_end'] < pos:
#                                 c[j] = 3 + addend_for_if_not_TBTGFC
#                                 d[j] = pos-R_mat_loaded_and_prepared_df.at[i,'code_end']
#                                 continue
#                         else: # (-)
#                             if R_mat_loaded_and_prepared_df.at[i,'code_start'] > pos:
#                                 c[j] = 3 + addend_for_if_not_TBTGFC
#                                 d[j] = R_mat_loaded_and_prepared_df.at[i,'code_start'] - pos
#                                 continue
#                             if R_mat_loaded_and_prepared_df.at[i,'code_end'] < pos:
#                                 c[j] = 4 + addend_for_if_not_TBTGFC
#                                 d[j] = pos-R_mat_loaded_and_prepared_df.at[i,'code_end']
#                                 continue
#                         # in exon(s)?
#                         in_e = list(np.array(range(1,len(R_mat_loaded_and_prepared_df.at[i,'exon_starts'])+1))[(R_mat_loaded_and_prepared_df.at[i,'exon_starts'] <= pos) & (R_mat_loaded_and_prepared_df.at[i,'exon_ends'] >= pos)])
#                         if in_e:
#                             c[j] = 1 + addend_for_if_not_TBTGFC
#                             continue
#                         # otherwise: in intron
#                         c[j] = 2 + addend_for_if_not_TBTGFC
#                         for k in range(0, R_mat_loaded_and_prepared_df.at[i,'n_exons']-1):
#                             if (R_mat_loaded_and_prepared_df.at[i,'exon_ends'][k] < pos) and (R_mat_loaded_and_prepared_df.at[i,'exon_starts'][k+1] > pos):
#                                 e1[j] = k+1
#                                 d1[j] = pos - R_mat_loaded_and_prepared_df.at[i,'exon_ends'][k]
#                                 d2[j] = R_mat_loaded_and_prepared_df.at[i,'exon_starts'][k+1] - pos
#                                 d[j] = min(d1[j], d2[j])
#                                 break
                    
#                     # % find the transcript in the highest-priority class
#                     zone = -1
#                     for cidx in list(range(1,11)):
#                         idx = [z for z in range(len(c)) if c[z] == cidx]
#                         if not idx:
#                             continue
#                         if cidx == 1 or cidx == 6:
#                             j = idx[0]
#                         else:
#                             k = np.argmin(d[np.array(idx)])
#                             j = idx[k]
#                         i = overlaps_df.index[j]
#                         name = R_mat_loaded_and_prepared_df.at[i,'gene']
#                         strand = R_mat_loaded_and_prepared_df.at[i,'strand']
#                         if strand == '-':
#                             e1[j] = R_mat_loaded_and_prepared_df.at[i,'n_exons'] - e1[j] + 1
#                         zone = cidx
#                         if zone == 1 or zone == 6:
#                             which_column = "exon"
#                         elif zone == 2 or zone == 7:
#                             # This code for intronnum is actually different from the original annotator:
#                             if strand == '+':
#                                 intronnum = e1[j]
#                             elif strand == '-':
#                                 intronnum = e1[j] - 1
#                             else:
#                                 raise("strand variable has an invalid value")
#                             which_column = "intron" + str(int(intronnum))
#                         elif zone == 3 or zone == 8:
#                             which_column = "3'UTR"
#                         elif zone == 4 or zone == 9:
#                             which_column = "5'UTR"
#                         elif zone == 5 or zone == 10:
#                             which_column = "promoter"
#                         else:
#                             raise('zone variable is not what it should be')
#                         if name == input_gene_name:
#                             if not which_column in list(counts_of_annotations_for_each_point_df.columns):
#                                 counts_of_annotations_for_each_point_df[which_column] = 0
#                             counts_of_annotations_for_each_point_df.at[input_gene_name,which_column] += 1
#                             # # FOLLOWING LINE OF CODE SHOULD BE DELETED/COMMENTED OUT WHEN ACTUALLY RUNNING THE ANNOTATIONS:
#                             # print("pos: " + str(pos) + ", type of annotation: " + which_column)
#                         else:
#                             raise("There should never be anything not annotated to our given gene")
#                             #counts_of_annotations_for_each_point_df.loc[input_gene_name,"not_annotated_to_gene"] += 1
                            
#                             # # FOLLOWING LINE OF CODE SHOULD BE DELETED/COMMENTED OUT WHEN ACTUALLY RUNNING THE ANNOTATIONS:
#                             # print("pos: " + str(pos) + ", type of annotation: " + "not_annotated_to_gene")
#                         break
#         counts_of_annotations_for_each_point_df.to_csv(file_name_based_off_input_parameters,sep='\t',index=False)


def find_surrounding_annotations_for_gene_with_parameters_post(input_gene_name, step_size, surrounding_window_length_one_way_of_gene_for_analysis, 
                                                               dict_of_genes_and_gene_orientations, dict_of_genes_and_gene_chrs, dict_of_genes_and_gene_starts, dict_of_genes_and_gene_ends):
    
    sizes_of_each_chromosome_hg38 = {1: 248956422, 2: 242193529, 3: 198295559, 4: 190214555, 5: 181538259, 6: 170805979,
                                7: 159345973, 8: 145138636, 9: 138394717, 10: 133797422, 11: 135086622, 12: 133275309,
                                13: 114364328, 14: 107043718, 15: 101991189, 16: 90338345, 17: 83257441, 18: 80373285,
                                19: 58617616, 20: 64444167, 21: 46709983, 22: 50818468, 23: 156040895, # chrX
                                24: 57227415} # chrY
        
    pre_annotation_counts_of_annotations_for_each_point_df = pd.read_csv("../surrounding_annotations_for_each_gene/v14/pre_annotations/gene_" + input_gene_name + ".step_size_" + str(step_size) + "bp.pre_annotations.tsv",sep='\t')
    
    # Every gene has only one orientation and not multiple:
    gene_orientation = dict_of_genes_and_gene_orientations[input_gene_name]
    # Accounting for how there can be genes belonging to multiple chromosomes:
    gene_chrs = dict_of_genes_and_gene_chrs[input_gene_name]
    gene_starts = dict_of_genes_and_gene_starts[input_gene_name]
    gene_ends = dict_of_genes_and_gene_ends[input_gene_name]
    
    for gene_chr, gene_start, gene_end in zip(gene_chrs, gene_starts, gene_ends):
        if gene_orientation == '+':
            if gene_start - surrounding_window_length_one_way_of_gene_for_analysis < 0:
                pre_annotation_counts_of_annotations_for_each_point_df.at[0,'before_gene'] += int(np.floor(gene_start / step_size))
            else:
                pre_annotation_counts_of_annotations_for_each_point_df.at[0,'before_gene'] += int(np.floor(surrounding_window_length_one_way_of_gene_for_analysis / step_size))
            if gene_end + surrounding_window_length_one_way_of_gene_for_analysis > sizes_of_each_chromosome_hg38[gene_chr]:
                pre_annotation_counts_of_annotations_for_each_point_df.at[0,'after_gene'] += int(np.floor((sizes_of_each_chromosome_hg38[gene_chr] - gene_end) / step_size))
            else:
                pre_annotation_counts_of_annotations_for_each_point_df.at[0,'after_gene'] += int(np.floor(surrounding_window_length_one_way_of_gene_for_analysis / step_size))
        else:
            if gene_start - surrounding_window_length_one_way_of_gene_for_analysis < 0:
                pre_annotation_counts_of_annotations_for_each_point_df.at[0,'after_gene'] += int(np.floor(gene_start / step_size))
            else:
                pre_annotation_counts_of_annotations_for_each_point_df.at[0,'after_gene'] += int(np.floor(surrounding_window_length_one_way_of_gene_for_analysis / step_size))
            if gene_end + surrounding_window_length_one_way_of_gene_for_analysis > sizes_of_each_chromosome_hg38[gene_chr]:
                pre_annotation_counts_of_annotations_for_each_point_df.at[0,'before_gene'] += int(np.floor((sizes_of_each_chromosome_hg38[gene_chr] - gene_end) / step_size))
            else:
                pre_annotation_counts_of_annotations_for_each_point_df.at[0,'before_gene'] += int(np.floor(surrounding_window_length_one_way_of_gene_for_analysis / step_size))
        
    pre_annotation_counts_of_annotations_for_each_point_df.to_csv("../surrounding_annotations_for_each_gene/v14/gene_" + input_gene_name + ".surr_dist_" + str(surrounding_window_length_one_way_of_gene_for_analysis) + \
                                                               "bp" + ".step_size_" + str(step_size) + "bp.tsv",sep='\t',index=False)


# This just applies the function that finds pre-annotations for a gene every X bp (find_surrounding_annotations_for_gene_with_parameters) to a list of genes
# def find_surrounding_annotations_for_gene_with_parameters_for_list_of_input_genes_pre(list_of_input_genes, step_size):
#     for g in list_of_input_genes:
#         print('Gene:', g)
#         find_surrounding_annotations_for_gene_with_parameters_pre(g, step_size)
        

# This just applies the function that finds pre-annotations for a gene every X bp (find_surrounding_annotations_for_gene_with_parameters) to a list of genes
def find_surrounding_annotations_for_gene_with_parameters_for_list_of_input_genes_post(list_of_input_genes, step_size, surrounding_window_length_one_way_of_gene_for_analysis, \
                                                                                       dict_of_genes_and_gene_orientations, dict_of_genes_and_gene_chrs, dict_of_genes_and_gene_starts, dict_of_genes_and_gene_ends):
    for g in list_of_input_genes:
        print('Gene:', g)
        find_surrounding_annotations_for_gene_with_parameters_post(g, step_size, surrounding_window_length_one_way_of_gene_for_analysis, \
                                                                   dict_of_genes_and_gene_orientations, dict_of_genes_and_gene_chrs, dict_of_genes_and_gene_starts, dict_of_genes_and_gene_ends)



# This finds the probability of LOF for a given gene based off the annotations for a gene and the surrounding window of the gene based off an input parameter
def find_probability_of_LOF_for_a_given_gene(gene, \
                                             gene_orientation, \
                                             path_for_annotations_of_gene, \
                                             surrounding_window_length_one_way_of_gene_for_analysis, \
                                             overall_proportion_of_genome_thats_exonic_interval, \
                                             overall_proportion_of_genome_thats_intronic_interval, \
                                             prob_of_intragenic, \
                                             prob_of_intergenic, \
                                             prob_of_deletion, \
                                             prob_of_tandemdup, \
                                             prob_of_inversion, \
                                            #  table_for_intragenic_unlikely_vs_likely_LOF_path="table_of_likely_LOF_vs_unlikely_LOF_for_intragenic_events_v14.tsv"):
                                             table_for_intragenic_unlikely_vs_likely_LOF_df):
    # # prob_of_intragenic and prob_of_intergenic can be calculated via the following 5 lines of code:
    # counts_of_each_type_of_event = find_probability_of_intergenic_vs_intragenic_gene_for_cohort(path_to_final_SV_list_for_analysis, threshold_of_proximity, R_mat_file_path)
    # overall_count_of_intragenic = counts_of_each_type_of_event['intragenic']
    # overall_count_of_intergenic = counts_of_each_type_of_event['intergenic']
    # prob_of_intragenic = overall_count_of_intragenic / (overall_count_of_intragenic + overall_count_of_intergenic)
    # prob_of_intergenic = overall_count_of_intergenic / (overall_count_of_intragenic + overall_count_of_intergenic)

    time_start = time.time()
    
    if gene != ".".join(path_for_annotations_of_gene.split('/')[-1].split('.')[:-3]).split('_')[1]:
    # print("path_for_annotations_of_gene:", path_for_annotations_of_gene)
    # if gene != ".".join(path_for_annotations_of_gene.split('/')[-1].split('.')[:-2]).split('_')[3]:
        raise("Gene doesn't correspond to the path for the surrounding annotations .tsv file!")

    gene_regions_mini_df = pd.read_csv(path_for_annotations_of_gene,sep='\t')
    # print("(**) find_probability_of_LOF_for_a_given_gene Time point 0.5:", time.time()-time_start)
    # gene_regions_series = gene_regions_mini_df.iloc[0]
    gene_regions_dict = gene_regions_mini_df.iloc[0].to_dict()
    
    # print("(**) find_probability_of_LOF_for_a_given_gene Time point 1:", time.time()-time_start)

    intron_columns = list(filter(lambda x: x.startswith("intron"), list(gene_regions_mini_df.columns)))
    total_length_annotated_to_gene_within_threshold = gene_regions_dict["exon"] + gene_regions_dict["3'UTR"] + gene_regions_dict["5'UTR"] + gene_regions_dict["promoter"] + \
                                                      gene_regions_dict['before_gene'] + gene_regions_dict['after_gene']
    # print("(**) find_probability_of_LOF_for_a_given_gene Time point 2:", time.time()-time_start)
    total_num_annotations_for_introns = 0
    total_num_annotations_for_introns_excl_intron_1 = 0
    for intr_col in intron_columns:
        total_length_annotated_to_gene_within_threshold += gene_regions_dict[intr_col]
        total_num_annotations_for_introns += gene_regions_dict[intr_col]
        if intr_col != 'intron1':
            total_num_annotations_for_introns_excl_intron_1 += gene_regions_dict[intr_col]
    # print("(**) find_probability_of_LOF_for_a_given_gene Time point 3:", time.time()-time_start)
    # Find probability of LOF for intragenic events:
    # table_of_intragenic_events_labeled_unlikely_vs_likely_LOF_df = pd.read_csv(table_for_intragenic_unlikely_vs_likely_LOF_path,sep='\t')
    # print("(**) find_probability_of_LOF_for_a_given_gene Time point 3.1:", time.time()-time_start)
    # table_for_intragenic_unlikely_vs_likely_LOF_df.set_index('Placement of breakpoints',inplace=True)
    # print("(**) find_probability_of_LOF_for_a_given_gene Time point 3.2:", time.time()-time_start)
    table_of_intragenic_events_probabilities_df = pd.DataFrame(data=np.zeros((table_for_intragenic_unlikely_vs_likely_LOF_df.shape[0],table_for_intragenic_unlikely_vs_likely_LOF_df.shape[1])),\
                                                                index=table_for_intragenic_unlikely_vs_likely_LOF_df.index,
                                                                columns=table_for_intragenic_unlikely_vs_likely_LOF_df.columns)
    # print("(**) find_probability_of_LOF_for_a_given_gene Time point 3.3:", time.time()-time_start)
    # print("(**) find_probability_of_LOF_for_a_given_gene Time point 4:", time.time()-time_start)
    probability_of_events_happening_with_same_intron_given_both_bp_in_intron = np.sum(list(map(lambda intr_col: (gene_regions_dict[intr_col] / total_num_annotations_for_introns)**2, intron_columns)))
    # print("(**) find_probability_of_LOF_for_a_given_gene Time point 5:", time.time()-time_start)
    for index, row in table_for_intragenic_unlikely_vs_likely_LOF_df.iterrows():
        # print("(**) find_probability_of_LOF_for_a_given_gene for-loop Time point 5.0:", time.time()-time_start)
        # print('index:', index)
        probability_reduction_factor = 1
        if '(in-frame)' in index:
            probability_reduction_factor = 1/3
        elif '(out-of-frame)' in index:
            probability_reduction_factor = 2/3
        # print("(**) find_probability_of_LOF_for_a_given_gene for-loop Time point 5.1:", time.time()-time_start)
        first_region, second_region = index.split(' ')[0].split('-')
        if first_region == 'intron' and second_region == 'intron':
            probability_reduction_factor *= (1-probability_of_events_happening_with_same_intron_given_both_bp_in_intron)
        if first_region.startswith('intron'):
            factor1_for_product_for_probability = total_num_annotations_for_introns / total_length_annotated_to_gene_within_threshold
        else:
            factor1_for_product_for_probability = gene_regions_dict[first_region] / total_length_annotated_to_gene_within_threshold
        if second_region.startswith('intron'):
            factor2_for_product_for_probability = total_num_annotations_for_introns / total_length_annotated_to_gene_within_threshold
        else:
            factor2_for_product_for_probability = gene_regions_dict[second_region] / total_length_annotated_to_gene_within_threshold
        if first_region == second_region:
            probability_of_event_no_matter_bp_orientations = factor1_for_product_for_probability * factor2_for_product_for_probability
        else:
            probability_of_event_no_matter_bp_orientations = 2 * factor1_for_product_for_probability * factor2_for_product_for_probability
        # print("(**) find_probability_of_LOF_for_a_given_gene for-loop Time point 5.2:", time.time()-time_start)
        if not '(within same intron)' in index:
            table_of_intragenic_events_probabilities_df.at[index, 'Deletion cases'] = probability_of_event_no_matter_bp_orientations * prob_of_deletion * probability_reduction_factor
            table_of_intragenic_events_probabilities_df.at[index, 'Tandem duplication cases'] = probability_of_event_no_matter_bp_orientations * prob_of_tandemdup * probability_reduction_factor
            table_of_intragenic_events_probabilities_df.at[index, 'Inversions cases'] = probability_of_event_no_matter_bp_orientations * prob_of_inversion * probability_reduction_factor
        else:
            table_of_intragenic_events_probabilities_df.at[index, 'Deletion cases'] = probability_of_event_no_matter_bp_orientations * prob_of_deletion * probability_of_events_happening_with_same_intron_given_both_bp_in_intron
            table_of_intragenic_events_probabilities_df.at[index, 'Tandem duplication cases'] = probability_of_event_no_matter_bp_orientations * prob_of_tandemdup * probability_of_events_happening_with_same_intron_given_both_bp_in_intron
            table_of_intragenic_events_probabilities_df.at[index, 'Inversions cases'] = probability_of_event_no_matter_bp_orientations * prob_of_inversion * probability_of_events_happening_with_same_intron_given_both_bp_in_intron
        # print("(**) find_probability_of_LOF_for_a_given_gene for-loop Time point 5.3:", time.time()-time_start)
    # print("(**) find_probability_of_LOF_for_a_given_gene Time point 6:", time.time()-time_start)
    p_LOF_intragenic = 0
    for index, row in table_of_intragenic_events_probabilities_df.iterrows():
        for col_name in list(table_of_intragenic_events_probabilities_df.columns):
            if table_for_intragenic_unlikely_vs_likely_LOF_df.at[index,col_name] == 'LL':
                p_LOF_intragenic += row[col_name]
    # print("(**) find_probability_of_LOF_for_a_given_gene Time point 7:", time.time()-time_start)

    # Find probability of LOF for intergenic events:
    # # This was the code for p_LOF_intergenic before:
    # p_LOF_intergenic = (gene_regions_mini_df.iloc[0]["exon"] / total_length_annotated_to_gene_within_threshold) * (1 - (1/6)*(overall_proportion_of_genome_thats_exonic_interval)) + (total_num_annotations_for_introns / total_length_annotated_to_gene_within_threshold) * (1 - (1/6)*(overall_proportion_of_genome_thats_intronic_interval))
    # # This is the code for p_LOF_intergenic discounting intron 1 but without factoring in antisense/out-of-frame fusions that could happen within intron 1:
    # p_LOF_intergenic = (gene_regions_mini_df.iloc[0]["exon"] / total_length_annotated_to_gene_within_threshold) * (1 - (1/6)*(overall_proportion_of_genome_thats_exonic_interval)) + (total_num_annotations_for_introns_excl_intron_1 / total_length_annotated_to_gene_within_threshold) * (1 - (1/6)*(overall_proportion_of_genome_thats_intronic_interval))
    # This is the code for p_LOF_intergenic discounting intron 1 that actually factors in antisense/out-of-frame fusions that could happen within intron 1:
    p_LOF_intergenic = (gene_regions_dict["exon"] / total_length_annotated_to_gene_within_threshold) * (1 - (1/6)*(overall_proportion_of_genome_thats_exonic_interval)) + \
                       (total_num_annotations_for_introns_excl_intron_1 / total_length_annotated_to_gene_within_threshold) * (1 - (1/6)*(overall_proportion_of_genome_thats_intronic_interval)) + \
                       ((total_num_annotations_for_introns - total_num_annotations_for_introns_excl_intron_1) / total_length_annotated_to_gene_within_threshold) * ((5/6)*(overall_proportion_of_genome_thats_intronic_interval) + (overall_proportion_of_genome_thats_exonic_interval))
    # print("(**) find_probability_of_LOF_for_a_given_gene Time point 8:", time.time()-time_start)
    
    # print('gene:', gene)
    # print('p_LOF_intragenic:', p_LOF_intragenic)
    # print('prob_of_intragenic:', prob_of_intragenic)
    # print('p_LOF_intergenic:', p_LOF_intergenic)
    # print('prob_of_intergenic:', prob_of_intergenic)
    
    return p_LOF_intragenic * prob_of_intragenic + p_LOF_intergenic * prob_of_intergenic, table_of_intragenic_events_probabilities_df

# This computes the LOF significance for a gene by going through all the events annotated to a gene, classifying the type of 
# event that each one is, determining if each event corresponds to LOF, and from that getting a count of LOF vs. not LOF events
# and computing a p-value using a binomial test.
# def compute_LOF_significance_for_gene(gene, path_to_final_SV_list_for_analysis, probability_of_LOF_for_gene_for_binomial_test, R_mat_loaded_and_prepared_df, R_mat_loaded_and_prepared_df_subset_for_gene, threshold_of_proximity, which_project, table_for_intragenic_unlikely_vs_likely_LOF_path="table_of_likely_LOF_vs_unlikely_LOF_for_intragenic_events_v14.tsv"):
def compute_LOF_significance_for_gene(gene, path_to_final_SV_list_for_analysis, probability_of_LOF_for_gene_for_binomial_test, R_mat_loaded_and_prepared_df, R_mat_loaded_and_prepared_df_subset_for_gene, threshold_of_proximity, which_project, table_for_intragenic_unlikely_vs_likely_LOF_df):
    time_start = time.time()
    SV_list_df_original = pd.read_csv(path_to_final_SV_list_for_analysis,sep='\t')
    # print("(**) compute_LOF_significance_for_gene Time point 0.1:", time.time()-time_start)
    # R_mat_loaded_and_prepared_df_subset = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]
    # gene_start = R_mat_loaded_and_prepared_df_subset['gene_start'].min()
    # gene_end = R_mat_loaded_and_prepared_df_subset['gene_end'].max()
    # print("(**) compute_LOF_significance_for_gene Time point 0.2:", time.time()-time_start)
    gene_orientation = R_mat_loaded_and_prepared_df_subset_for_gene['strand'].iat[0]
    # print("(**) compute_LOF_significance_for_gene Time point 0.3:", time.time()-time_start)
    is_gene_pseudoautosomal = len(set(R_mat_loaded_and_prepared_df_subset_for_gene['chr'])) > 1
    # print("(**) compute_LOF_significance_for_gene Time point 0.4:", time.time()-time_start)
    # print("(**) compute_LOF_significance_for_gene Time point 1:", time.time()-time_start)
    if not is_gene_pseudoautosomal:
        gene_chr_in_R_mat = R_mat_loaded_and_prepared_df_subset_for_gene['chr'].iat[0]
        gene_start_in_R_mat = R_mat_loaded_and_prepared_df_subset_for_gene['gene_start'].min()
        gene_end_in_R_mat = R_mat_loaded_and_prepared_df_subset_for_gene['gene_end'].max()
    num_patients_events_used_for_test = len(set(SV_list_df_original['individual']))
    # print("(**) compute_LOF_significance_for_gene Time point 2:", time.time()-time_start)
    list_of_patients_with_at_least_one_LOF_event = []
    num_events_total_within_threshold_of_proximity = 0
    num_intragenic_events = 0
    num_intergenic_events = 0
    num_events_LOF = 0
    list_of_indicators_for_LOF_events = []
    # print("(**) compute_LOF_significance_for_gene Time point 3:", time.time()-time_start)
    for index,row in SV_list_df_original.iterrows():
        # print('index:', index)
        pos_that_should_correspond_to_gene = [True,True]
        if not is_gene_pseudoautosomal:
            if (row['chr1'] != gene_chr_in_R_mat) or (row['pos1'] < gene_start_in_R_mat - threshold_of_proximity) or (row['pos1'] > gene_end_in_R_mat + threshold_of_proximity):
                pos_that_should_correspond_to_gene[0] = False
            if (row['chr2'] != gene_chr_in_R_mat) or (row['pos2'] < gene_start_in_R_mat - threshold_of_proximity) or (row['pos2'] > gene_end_in_R_mat + threshold_of_proximity):
                pos_that_should_correspond_to_gene[1] = False
        else:
            if not is_within_specified_range_of_gene(R_mat_loaded_and_prepared_df,gene,threshold_of_proximity,row['pos1'],row['chr1']):
                pos_that_should_correspond_to_gene[0] = False
            if not is_within_specified_range_of_gene(R_mat_loaded_and_prepared_df,gene,threshold_of_proximity,row['pos2'],row['chr2']):
                pos_that_should_correspond_to_gene[1] = False
        # # if (row['chr1'] != gene_chr_in_R_mat) or (row['pos1'] < gene_start_in_R_mat - threshold_of_proximity) or (row['pos1'] > gene_end_in_R_mat + threshold_of_proximity):
        # if not is_within_specified_range_of_gene(R_mat_loaded_and_prepared_df,gene,threshold_of_proximity,row['pos1'],row['chr1']):
        #     pos_that_should_correspond_to_gene[0] = False
        # # if (row['chr2'] != gene_chr_in_R_mat) or (row['pos2'] < gene_start_in_R_mat - threshold_of_proximity) or (row['pos2'] > gene_end_in_R_mat + threshold_of_proximity):
        # if not is_within_specified_range_of_gene(R_mat_loaded_and_prepared_df,gene,threshold_of_proximity,row['pos2'],row['chr2']):
        #     pos_that_should_correspond_to_gene[1] = False
        # print("(***) compute_LOF_significance_for_gene for-loop Time point 3.1:", time.time()-time_start)
        if pos_that_should_correspond_to_gene[0] == False and pos_that_should_correspond_to_gene[1] == False: # should not be counted as an event for this gene
            raise("This should not be possible. All the events we're looking at should have at least one breakpoint within the given range threshold of the gene.")
        elif pos_that_should_correspond_to_gene[0] == False or pos_that_should_correspond_to_gene[1] == False: # intergenic case
            time_start_intergenic = time.time()
            num_events_total_within_threshold_of_proximity += 1
            num_intergenic_events += 1
            if pos_that_should_correspond_to_gene[0]:
                s1 = str(1)
                s2 = str(2)
            elif pos_that_should_correspond_to_gene[1]:
                s1 = str(2)
                s2 = str(1)
            if ((row['site'+s1].lower().startswith('exon')) and (not "in frame" in row['fusion'].lower())):
                num_events_LOF += 1
                list_of_patients_with_at_least_one_LOF_event.append(row['individual'])
                list_of_indicators_for_LOF_events.append(1)
                continue
            elif ((row['site'+s1].lower().startswith('intron')) and (not "in frame" in row['fusion'].lower())):
                if row['site'+s1].lower().endswith('before exon 2') or row['site'+s1].lower().endswith('after exon 1'):
                    # Needed to manually make own condition for antisense in three lines of code below (which assumes an intron-intron SV event). Realized antisense annotations by 
                    # dRanger aren't always correct (or at least there are conditions not considered antisense by the dRanger annotator where the downstream parts of transcripts 
                    # have opposite orientations and this would result in loss of function that we would need to account for)
                    do_gene_orientations_match = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == row['gene1']]['strand'].iat[0] == \
                                                 R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == row['gene2']]['strand'].iat[0]
                    manual_antisense_condition = (do_gene_orientations_match and row['str1'] == row['str2']) or ((not do_gene_orientations_match) and row['str1'] != row['str2'])
                    if (row['site'+s2].lower().startswith('intron')) and (manual_antisense_condition or "out of frame" in row['fusion'].lower()):
                        num_events_LOF += 1
                        list_of_patients_with_at_least_one_LOF_event.append(row['individual'])
                        list_of_indicators_for_LOF_events.append(1)
                        continue
                    elif (row['site'+s2].lower().startswith('exon')):
                        num_events_LOF += 1
                        list_of_patients_with_at_least_one_LOF_event.append(row['individual'])
                        list_of_indicators_for_LOF_events.append(1)
                        continue
                else:
                    num_events_LOF += 1
                    list_of_patients_with_at_least_one_LOF_event.append(row['individual'])
                    list_of_indicators_for_LOF_events.append(1)
                    continue
            # print("(****) compute_LOF_significance_for_gene for-loop intergenic case Time taken:", time.time()-time_start_intergenic)
        else: # intragenic case
            time_start_intragenic = time.time()
            num_events_total_within_threshold_of_proximity += 1
            num_intragenic_events += 1
            annotation_for_each_bp = []
            # print("(****) compute_LOF_significance_for_gene for-loop intragenic case Time point 1 taken:", time.time()-time_start_intragenic)
            for i in [1,2]:
                site_annotation = row["site"+str(i)]
                # print('site_annotation:')
                # print(site_annotation)
                if site_annotation.lower().startswith("exon"):
                    annotation_for_each_bp.append("exon")
                elif site_annotation.lower().startswith("intron"):
                    annotation_for_each_bp.append("intron")
                elif site_annotation.lower().startswith("3'-utr"):
                    annotation_for_each_bp.append("3'UTR")
                elif site_annotation.lower().startswith("5'-utr"):
                    annotation_for_each_bp.append("5'UTR")
                elif site_annotation.lower().startswith("promoter"):
                    annotation_for_each_bp.append("promoter")
                # # This is what I had previously:
                # # I am assuming that if our position here is before any gene start in the R.mat file for that gene that it cannot take place after the gene end for the transcript that is used for the position's actual annotation
                # elif (site_annotation.lower().startswith("igr") and (row['pos'+str(i)] <= R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]['gene_start']).any() and R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]['strand'].iloc[0] == '+') or \
                #      (site_annotation.lower().startswith("igr") and (row['pos'+str(i)] >= R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]['gene_end']).any() and R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]['strand'].iloc[0] == '-'):
                #     print('YEAH, THIS CONDITION 1 WAS CARRIED OUT!!!!')
                #     annotation_for_each_bp.append("before_gene")
                # elif (site_annotation.lower().startswith("igr") and (row['pos'+str(i)] >= R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]['gene_end']).any() and R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]['strand'].iloc[0] == '+') or \
                #      (site_annotation.lower().startswith("igr") and (row['pos'+str(i)] <= R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]['gene_start']).any() and R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]['strand'].iloc[0] == '-'):
                #     print('YEAH, THIS CONDITION 2 WAS CARRIED OUT!!!!')
                #     annotation_for_each_bp.append("after_gene")
                elif site_annotation.lower().startswith("igr"):
                    # print("(*****) compute_LOF_significance_for_gene for-loop intragenic case Time point 1.1 taken:", time.time()-time_start_intragenic)
                    # # dist_before = R_mat_loaded_and_prepared_df_subset['tx_start'].subtract(row["pos"+str(i)]).where(lambda x: x >= 0).dropna()
                    # dist_before_with_negatives = np.array(R_mat_loaded_and_prepared_df_subset['tx_start']) - row["pos"+str(i)]
                    # dist_before = dist_before_with_negatives[dist_before_with_negatives >= 0]
                    # print("(*****) compute_LOF_significance_for_gene for-loop intragenic case Time point 1.2 taken:", time.time()-time_start_intragenic)
                    # # dist_after = (-R_mat_loaded_and_prepared_df_subset['tx_end']).add(row["pos"+str(i)]).where(lambda x: x >= 0).dropna()
                    # dist_after_with_negatives = -np.array(R_mat_loaded_and_prepared_df_subset['tx_end']) + row["pos"+str(i)]
                    # dist_after = dist_after_with_negatives[dist_after_with_negatives >= 0]
                    # print("(*****) compute_LOF_significance_for_gene for-loop intragenic case Time point 1.3 taken:", time.time()-time_start_intragenic)
                    if not is_gene_pseudoautosomal:
                        gene_start = gene_start_in_R_mat
                        gene_end = gene_end_in_R_mat
                    else:
                        gene_start = R_mat_loaded_and_prepared_df_subset_for_gene[R_mat_loaded_and_prepared_df_subset_for_gene["chr"] == row["chr"+str(i)]]['gene_start'].min()
                        gene_end = R_mat_loaded_and_prepared_df_subset_for_gene[R_mat_loaded_and_prepared_df_subset_for_gene["chr"] == row["chr"+str(i)]]['gene_end'].max()
                    if row["pos"+str(i)] <= gene_start:
                        if gene_orientation == '+':
                            annotation_for_each_bp.append("before_gene")
                        else:
                            annotation_for_each_bp.append("after_gene")
                    elif row["pos"+str(i)] >= gene_end:
                        if gene_orientation == '+':
                            annotation_for_each_bp.append("after_gene")
                        else:
                            annotation_for_each_bp.append("before_gene")
                    else: # This is just for the case where there could be an annotation for IGR at a point in between isoforms
                        # dist_before_with_negatives = np.array(R_mat_loaded_and_prepared_df_subset['tx_start']) - row["pos"+str(i)]
                        # dist_before = dist_before_with_negatives[dist_before_with_negatives >= 0]
                        # dist_after_with_negatives = -np.array(R_mat_loaded_and_prepared_df_subset['tx_end']) + row["pos"+str(i)]
                        # dist_after = dist_after_with_negatives[dist_after_with_negatives >= 0]
                        dist_before = R_mat_loaded_and_prepared_df_subset_for_gene['tx_start'].subtract(row["pos"+str(i)]).where(lambda x: x >= 0).dropna()
                        dist_after = (-R_mat_loaded_and_prepared_df_subset_for_gene['tx_end']).add(row["pos"+str(i)]).where(lambda x: x >= 0).dropna()
                        if dist_before.size > 0 and dist_after.size > 0:
                            if dist_before.min() <= dist_after.min():
                                if gene_orientation == '+':
                                    annotation_for_each_bp.append("before_gene")
                                else:
                                    annotation_for_each_bp.append("after_gene")
                            else:
                                if gene_orientation == '+':
                                    annotation_for_each_bp.append("after_gene")
                                else:
                                    annotation_for_each_bp.append("before_gene")
                        else:
                            raise("Something went wrong. There should be valid dist_before and dist_after minimums.")
                    # print("(*****) compute_LOF_significance_for_gene for-loop intragenic case Time point 1.4 taken:", time.time()-time_start_intragenic)
                else:
                    raise("Seemingly we have a site annotation outside of all the possibilities we accounted for... something is wrong here.")
            # print("(****) compute_LOF_significance_for_gene for-loop intragenic case Time point 2 taken:", time.time()-time_start_intragenic)
            # annotation_for_each_bp[0] + "-" + annotation_for_each_bp[1]
            # table_for_intragenic_unlikely_vs_likely_LOF_df = pd.read_csv(table_for_intragenic_unlikely_vs_likely_LOF_path,sep='\t')
            # print("table_for_intragenic_unlikely_vs_likely_LOF_df:")
            # print(table_for_intragenic_unlikely_vs_likely_LOF_df)
            table_for_intragenic_unlikely_vs_likely_LOF_df_binary = table_for_intragenic_unlikely_vs_likely_LOF_df.replace('UL',0).replace('LL',1)
            sub_df = table_for_intragenic_unlikely_vs_likely_LOF_df_binary[table_for_intragenic_unlikely_vs_likely_LOF_df_binary.index.to_series().str.contains(annotation_for_each_bp[0] + "-" + annotation_for_each_bp[1]) | table_for_intragenic_unlikely_vs_likely_LOF_df_binary.index.to_series().str.contains(annotation_for_each_bp[1] + "-" + annotation_for_each_bp[0])]
            if row['str1'] == 1 and row['str2'] == 0:
                column_label_for_type_of_event = 'Tandem duplication cases' # tandem duplication case
            elif row['str1'] == 0 and row['str2'] == 1:
                column_label_for_type_of_event = 'Deletion cases' # deletion case
            else:
                column_label_for_type_of_event = 'Inversions cases' # inversion case
            if sub_df.shape[0] == 0:
                raise("No corresponding events found in table of unlikely vs. likely LOF")
            if sub_df.shape[0] == 1:
                addend = sub_df[column_label_for_type_of_event].iat[0]
                num_events_LOF += addend
                if addend > 0:
                    list_of_patients_with_at_least_one_LOF_event.append(row['individual'])
                    list_of_indicators_for_LOF_events.append(1)
                    continue
            else:
                if annotation_for_each_bp[0] + "-" + annotation_for_each_bp[1] == "exon-intron" or annotation_for_each_bp[1] + "-" + annotation_for_each_bp[0] == "exon-intron":
                    if "in frame" in row['fusion']:
                        addend = table_for_intragenic_unlikely_vs_likely_LOF_df_binary[table_for_intragenic_unlikely_vs_likely_LOF_df_binary.index.to_series() == 'exon-intron (in-frame)'][column_label_for_type_of_event].iat[0]
                        num_events_LOF += addend
                        if addend > 0:
                            list_of_patients_with_at_least_one_LOF_event.append(row['individual'])
                            list_of_indicators_for_LOF_events.append(1)
                            continue
                    else:
                        addend = table_for_intragenic_unlikely_vs_likely_LOF_df_binary[table_for_intragenic_unlikely_vs_likely_LOF_df_binary.index.to_series() == 'exon-intron (out-of-frame)'][column_label_for_type_of_event].iat[0]
                        num_events_LOF += addend
                        if addend > 0:
                            list_of_patients_with_at_least_one_LOF_event.append(row['individual'])
                            list_of_indicators_for_LOF_events.append(1)
                            continue
                if annotation_for_each_bp[0] + "-" + annotation_for_each_bp[1] == "exon-exon":
                    if "in frame" in row['fusion']:
                        addend = table_for_intragenic_unlikely_vs_likely_LOF_df_binary[table_for_intragenic_unlikely_vs_likely_LOF_df_binary.index.to_series() == 'exon-exon (in-frame)'][column_label_for_type_of_event].iat[0]
                        num_events_LOF += addend
                        if addend > 0:
                            list_of_patients_with_at_least_one_LOF_event.append(row['individual'])
                            list_of_indicators_for_LOF_events.append(1)
                            continue
                    else:
                        addend = table_for_intragenic_unlikely_vs_likely_LOF_df_binary[table_for_intragenic_unlikely_vs_likely_LOF_df_binary.index.to_series() == 'exon-exon (out-of-frame)'][column_label_for_type_of_event].iat[0]
                        num_events_LOF += addend
                        if addend > 0:
                            list_of_patients_with_at_least_one_LOF_event.append(row['individual'])
                            list_of_indicators_for_LOF_events.append(1)
                            continue
                if annotation_for_each_bp[0] + "-" + annotation_for_each_bp[1] == "intron-intron":
                    bp1_before_or_after, _, bp1_exon_num = row['site1'].split(" ")[-3:]
                    bp2_before_or_after, _, bp2_exon_num = row['site2'].split(" ")[-3:]
                    condition_for_manually_checking_if_within_same_intron = (bp1_before_or_after == bp2_before_or_after and bp1_exon_num == bp2_exon_num) or \
                                                                            (bp1_before_or_after == 'before' and bp2_before_or_after == 'after' and int(bp1_exon_num) == int(bp2_exon_num)+1) or \
                                                                            (bp1_before_or_after == 'after' and bp2_before_or_after == 'before' and int(bp2_exon_num) == int(bp1_exon_num)+1)
                    if "within intron" in row['fusion'] or condition_for_manually_checking_if_within_same_intron:
                        addend = table_for_intragenic_unlikely_vs_likely_LOF_df_binary[table_for_intragenic_unlikely_vs_likely_LOF_df_binary.index.to_series() == 'intron-intron (within same intron)'][column_label_for_type_of_event].iat[0]
                        num_events_LOF += addend
                        if addend > 0:
                            list_of_patients_with_at_least_one_LOF_event.append(row['individual'])
                            list_of_indicators_for_LOF_events.append(1)
                            continue
                    elif "in frame" in row['fusion']:
                        addend = table_for_intragenic_unlikely_vs_likely_LOF_df_binary[table_for_intragenic_unlikely_vs_likely_LOF_df_binary.index.to_series() == 'intron-intron (in-frame)'][column_label_for_type_of_event].iat[0]
                        num_events_LOF += addend
                        if addend > 0:
                            list_of_patients_with_at_least_one_LOF_event.append(row['individual'])
                            list_of_indicators_for_LOF_events.append(1)
                            continue
                    else:
                        addend = table_for_intragenic_unlikely_vs_likely_LOF_df_binary[table_for_intragenic_unlikely_vs_likely_LOF_df_binary.index.to_series() == 'intron-intron (out-of-frame)'][column_label_for_type_of_event].iat[0]
                        num_events_LOF += addend
                        if addend > 0:
                            list_of_patients_with_at_least_one_LOF_event.append(row['individual'])
                            list_of_indicators_for_LOF_events.append(1)
                            continue
            # print("(****) compute_LOF_significance_for_gene for-loop intragenic case Time point 3 taken:", time.time()-time_start_intragenic)
        list_of_indicators_for_LOF_events.append(0)
        # print("(***) compute_LOF_significance_for_gene for-loop Time point 3.2:", time.time()-time_start)
        
    # print("(**) compute_LOF_significance_for_gene Time point 4:", time.time()-time_start)
    
    SV_list_df_original['LOF_or_not'] = list_of_indicators_for_LOF_events
    file_path_for_SV_list_df_with_LOF_or_not_indicator = "../SV_list_subsets_of_regions_within_distance_of_gene/v14/" + which_project + "/annotations_relative_to_given_gene_and_boundary_threshold_with_indication_of_LOF_or_not/"+"SV_list_filtered_"+gene+"_within_"+str(threshold_of_proximity)+"bp.reannotated_with_LOF_indication.tsv"
    if not os.path.exists("../SV_list_subsets_of_regions_within_distance_of_gene/v14/" + which_project + "/annotations_relative_to_given_gene_and_boundary_threshold_with_indication_of_LOF_or_not"):
        os.makedirs("../SV_list_subsets_of_regions_within_distance_of_gene/v14/" + which_project + "/annotations_relative_to_given_gene_and_boundary_threshold_with_indication_of_LOF_or_not")
    SV_list_df_original.to_csv(file_path_for_SV_list_df_with_LOF_or_not_indicator,sep='\t',index=False)
    
    # print("(**) compute_LOF_significance_for_gene Time point 5:", time.time()-time_start)
    
    p_value_for_LOF_test = stats.binom_test(num_events_LOF, n=num_events_total_within_threshold_of_proximity, p=probability_of_LOF_for_gene_for_binomial_test, alternative='greater')
    
    # print("Gene:", gene)
    # print("Number of intragenic events:", num_intragenic_events)
    # print("Number of intergenic events:", num_intergenic_events)
    # print("Number of total events (with at least one breakpoint within threshold of proximity of gene):", num_events_total_within_threshold_of_proximity)
    # print("Number of LOF events by criteria of model:", num_events_LOF)
    # print("Computed probability of LOF used for binomial test:", probability_of_LOF_for_gene_for_binomial_test)
    # print("p-value for LOF:", p_value_for_LOF_test)
    
    # print("(**) compute_LOF_significance_for_gene Time point 6:", time.time()-time_start)
    
    # # Returning as a list because it is more efficient runtime-wise
    # # I am returning things in the order of the following column names:
    # # gene, num_patients_used_for_test, num_patients_with_at_least_one_LOF_event, num_total_events_used_for_test, num_intragenic_events, num_intergenic_events, num_LOF_events, prob_of_LOF, p_value_LOF
    # return_list = []
    # return_list.append(gene)
    # return_list.append(num_patients_events_used_for_test)
    # return_list.append(len(set(list_of_patients_with_at_least_one_LOF_event)))
    # return_list.append(num_events_total_within_threshold_of_proximity)
    # return_list.append(num_intragenic_events)
    # return_list.append(num_intergenic_events)
    # return_list.append(num_events_LOF)
    # return_list.append(probability_of_LOF_for_gene_for_binomial_test)
    # return_list.append(p_value_for_LOF_test)
    # print("(**) compute_LOF_significance_for_gene Time point 7:", time.time()-time_start)
    # return return_list
    
    return pd.DataFrame(data={'gene': [gene], \
                              'num_patients_used_for_test': [num_patients_events_used_for_test], \
                              'num_patients_with_at_least_one_LOF_event': [len(set(list_of_patients_with_at_least_one_LOF_event))], \
                              'num_total_events_used_for_test': [num_events_total_within_threshold_of_proximity], \
                              'num_intragenic_events': [num_intragenic_events], \
                              'num_intergenic_events': [num_intergenic_events], \
                              'num_LOF_events': [num_events_LOF], \
                              'prob_of_LOF': [probability_of_LOF_for_gene_for_binomial_test], \
                              'p_value_LOF': [p_value_for_LOF_test]})
    

# This just takes a list of genes and relevant parameters and runs the LOF test for each gene, and returns an appropriate dataframe with
# the results at the end of the function.
def perform_LOF_test_for_list_of_genes(list_of_genes, \
                                       path_to_final_SV_list_for_analysis, \
                                       which_project, \
                                       path_to_directory_for_annotations_surrounding_gene, \
                                       threshold_of_proximity, \
                                       step_size_for_annotating_gene_regions=10, \
                                       R_mat_file_path="../reference_files/annotation_file/hg38_R_with_hugo_symbols_with_DUX4L1_HMGN2P46_MALAT1.mat", \
                                       table_for_intragenic_unlikely_vs_likely_LOF_path="table_of_likely_LOF_vs_unlikely_LOF_for_intragenic_events_v14.tsv"):
    
    time_start = time.time()
    
    SV_df = pd.read_csv(path_to_final_SV_list_for_analysis,sep='\t')
    R_mat_loaded_and_prepared_df = load_and_prepare_R_mat_file_into_df(R_mat_file_path)
    table_for_intragenic_unlikely_vs_likely_LOF_df = pd.read_csv(table_for_intragenic_unlikely_vs_likely_LOF_path,sep='\t')
    table_for_intragenic_unlikely_vs_likely_LOF_df.set_index('Placement of breakpoints',inplace=True)
    chrs_to_R_mat_subset_dfs_dict = {}
    for chr in set(R_mat_loaded_and_prepared_df['chr']):
        chrs_to_R_mat_subset_dfs_dict[chr] = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['chr'] == chr]
    print("Time point 1:", time.time()-time_start)
    overall_proportion_of_genome_thats_exonic_interval, overall_proportion_of_genome_thats_intronic_interval = find_overall_exonic_and_intronic_portions_of_genome_results()
    print("Time point 2:", time.time()-time_start)
    print("Producing SV lists for each gene containing events only within a certain distance threshold of the gene (" + str(len(list_of_genes)) + " genes total):")
    filtered_SV_list_gene_subset_paths = []
    for i in tqdm(range(len(list_of_genes))):
        filtered_SV_list_gene_subset_paths.append(filter_dataframe_to_have_only_events_with_bps_within_threshold_of_gene(list_of_genes[i], threshold_of_proximity, SV_df, R_mat_loaded_and_prepared_df, which_project))
    # filtered_SV_list_gene_subset_paths = list(map(lambda g: filter_dataframe_to_have_only_events_with_bps_within_threshold_of_gene(g, threshold_of_proximity, SV_df, R_mat_loaded_and_prepared_df, which_project), list_of_genes))
    print("Time point 3:", time.time()-time_start)
    print("Reannotating SV lists for each gene relative to that gene itself (" + str(len(list_of_genes)) + " genes total):")
    zipped_genes_and_sublist_paths = list(zip(list_of_genes, filtered_SV_list_gene_subset_paths))
    reannotated_filtered_SV_list_gene_subset_paths = []
    for i in tqdm(range(len(zipped_genes_and_sublist_paths))):
        reannotated_filtered_SV_list_gene_subset_paths.append(modified_dRanger_annotate_sites_driver_CGC_genes_prioritized(zipped_genes_and_sublist_paths[i][1], zipped_genes_and_sublist_paths[i][0], threshold_of_proximity, R_mat_loaded_and_prepared_df, chrs_to_R_mat_subset_dfs_dict, which_project))
    # reannotated_filtered_SV_list_gene_subset_paths = list(map(lambda g_and_f: modified_dRanger_annotate_sites_driver_CGC_genes_prioritized(g_and_f[1], g_and_f[0], threshold_of_proximity, R_mat_loaded_and_prepared_df, chrs_to_R_mat_subset_dfs_dict, which_project), zip(list_of_genes, filtered_SV_list_gene_subset_paths)))
    print("Time point 4:", time.time()-time_start)
    dict_of_genes_to_reannotated_filtered_SV_lists = {}
    for i in range((len(list_of_genes))):
        dict_of_genes_to_reannotated_filtered_SV_lists[list_of_genes[i]] = reannotated_filtered_SV_list_gene_subset_paths[i]
    print("Time point 5:", time.time()-time_start)
    print("Finding proportions of 1-breakpoint and 2-breakpoint events, as well as deletions/tandem duplications/inversions for two-breakpoint events across all genes (" + str(len(dict_of_genes_to_reannotated_filtered_SV_lists)) + " total) being tested for cohort:")
    prob_of_intergenic, prob_of_intragenic, proportion_deletions_across_cohort_for_intragenic, proportion_tandemdups_across_cohort_for_intragenic, proportion_inversions_across_cohort_for_intragenic = \
                        find_proportions_of_1_bp_2_bp_events_and_deletions_tandemdups_inversions_for_intragenic_all_genes(dict_of_genes_to_reannotated_filtered_SV_lists, R_mat_loaded_and_prepared_df, threshold_of_proximity)
    print("Time point 6:", time.time()-time_start)
    # # This is how we calculated counts_of_each_type_of_event for previous versions:
    # counts_of_each_type_of_event = find_probability_of_intergenic_vs_intragenic_gene_for_cohort(path_to_final_SV_list_for_analysis, threshold_of_proximity=threshold_of_proximity, R_mat_file_path=R_mat_file_path)
    # This is how we calculate counts_of_each_type_of_event now:
    # counts_of_each_type_of_event = find_probability_of_intergenic_vs_intragenic_gene_for_cohort(dict_of_genes_to_reannotated_filtered_SV_lists, threshold_of_proximity=threshold_of_proximity, R_mat_file_path=R_mat_file_path)
    # overall_count_of_intragenic = counts_of_each_type_of_event['intragenic']
    # overall_count_of_intergenic = counts_of_each_type_of_event['intergenic']
    # print("overall_count_of_intragenic:",overall_count_of_intragenic)
    # print("overall_count_of_intragenic+overall_count_of_intergenic:",overall_count_of_intragenic+overall_count_of_intergenic)
    # prob_of_intragenic = overall_count_of_intragenic / (overall_count_of_intragenic + overall_count_of_intergenic)
    # prob_of_intergenic = overall_count_of_intergenic / (overall_count_of_intragenic + overall_count_of_intergenic)
    # print("This is what we're setting to be the probability of an event being intragenic:",prob_of_intragenic)
    # print("This is what we're setting to be the probability of an event being intergenic:",prob_of_intergenic)
    
    print("Time point 7:", time.time()-time_start)
    
    # lists_of_values_of_rows_per_final_results_df = []
    dfs_to_be_concatted = []
    index_for_reindexing = 0
    print("Preparing (appropriate p parameter) and running binomial test on all " + str(len(list_of_genes)) + " genes being tested for cohort:")
    for i in tqdm(range(len(list_of_genes))):
    # for gene in list_of_genes:
        gene = list_of_genes[i]
        # print("Gene:", gene)
        # time_start_SVGAR_sf_test = time.time()
        path_for_annotations_surrounding_gene = path_to_directory_for_annotations_surrounding_gene + "/gene_" + gene + ".surr_dist_" + str(threshold_of_proximity) + "bp.step_size_" + str(step_size_for_annotating_gene_regions) + "bp.tsv"
        # print("(*) SVGAR_sf statistical test for one gene Time point 1:", time.time()-time_start_SVGAR_sf_test)
        # path_for_annotations_surrounding_gene = dict_of_genes_to_reannotated_filtered_SV_lists[gene]
        R_mat_loaded_and_prepared_df_subset_for_gene = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]
        gene_orientation = np.array(R_mat_loaded_and_prepared_df_subset_for_gene['strand'])[0]
        # print("(*) SVGAR_sf statistical test for one gene Time point 2:", time.time()-time_start_SVGAR_sf_test)
        # # Needed this for slides:
        # prob_LOF_based_off_gene, a = find_probability_of_LOF_for_a_given_gene(gene, \
        prob_LOF_based_off_gene, _ = find_probability_of_LOF_for_a_given_gene(gene, \
                                                 gene_orientation, \
                                                 path_for_annotations_surrounding_gene, \
                                                #  SV_df, \
                                                 threshold_of_proximity, \
                                                 overall_proportion_of_genome_thats_exonic_interval=overall_proportion_of_genome_thats_exonic_interval, \
                                                 overall_proportion_of_genome_thats_intronic_interval=overall_proportion_of_genome_thats_intronic_interval, \
                                                 prob_of_intragenic=prob_of_intragenic, \
                                                 prob_of_intergenic=prob_of_intergenic, \
                                                 prob_of_deletion=proportion_deletions_across_cohort_for_intragenic, \
                                                 prob_of_tandemdup=proportion_tandemdups_across_cohort_for_intragenic, \
                                                 prob_of_inversion=proportion_inversions_across_cohort_for_intragenic, \
                                                 table_for_intragenic_unlikely_vs_likely_LOF_df=table_for_intragenic_unlikely_vs_likely_LOF_df)
        
        # print("(*) SVGAR_sf statistical test for one gene Time point 3:", time.time()-time_start_SVGAR_sf_test)

        # # Needed this for slides:
        # if gene == 'CIITA':
        #     a.to_csv('CIITA_prob_matrix_for_1Mbp_surr_window.tsv',sep='\t')
        
        # print("table_for_intragenic_unlikely_vs_likely_LOF_df:")
        # print(table_for_intragenic_unlikely_vs_likely_LOF_df)
        
        row_for_LOF_test_for_gene_df = compute_LOF_significance_for_gene(gene, \
                                                                         dict_of_genes_to_reannotated_filtered_SV_lists[gene], \
                                                                         prob_LOF_based_off_gene, \
                                                                         R_mat_loaded_and_prepared_df, \
                                                                         R_mat_loaded_and_prepared_df_subset_for_gene, \
                                                                         threshold_of_proximity, \
                                                                         which_project, \
                                                                        #  table_for_intragenic_unlikely_vs_likely_LOF_path="table_of_likely_LOF_vs_unlikely_LOF_for_intragenic_events_v14.tsv")
                                                                         table_for_intragenic_unlikely_vs_likely_LOF_df=table_for_intragenic_unlikely_vs_likely_LOF_df)
        
        # print("(*) SVGAR_sf statistical test for one gene Time point 4:", time.time()-time_start_SVGAR_sf_test)
        
        row_for_LOF_test_for_gene_df.index = [index_for_reindexing]
        dfs_to_be_concatted.append(row_for_LOF_test_for_gene_df)
        index_for_reindexing += 1
        
        # lists_of_values_of_rows_per_final_results_df.append(list_of_values_for_row_for_LOF_test_for_gene_df)
        
        # print("(*) SVGAR_sf statistical test for one gene Time point 5:", time.time()-time_start_SVGAR_sf_test)
        
    
    print("Time point 8:", time.time()-time_start)
    
    concatted_df = pd.concat(dfs_to_be_concatted).sort_values(by=['num_total_events_used_for_test'])
    # final_results_df = pd.DataFrame(lists_of_values_of_rows_per_final_results_df)
    # final_results_df.columns = ["gene", 
    #                             "num_patients_used_for_test", 
    #                             "num_patients_with_at_least_one_LOF_event", 
    #                             "num_total_events_used_for_test", 
    #                             "num_intragenic_events", 
    #                             "num_intergenic_events", 
    #                             "num_LOF_events", 
    #                             "prob_of_LOF", 
    #                             "p_value_LOF"]
    # final_results_df.sort_values(by=['num_total_events_used_for_test'],inplace=True)
    
    print("Time point 9:", time.time()-time_start)
    
    # return 
    return concatted_df


def SV_analysis_loss_of_function_events_v14(version, which_dataset, SV_list_file, distance_threshold):
    
    time_start = time.time()
    
    gene_list_file = "../lists_of_genes_annotated_to_SVs_recurring_across_at_least_2_patients_for_cohorts/v14/" + which_dataset + "/list_of_genes_annotated_to_SVs_recurring_across_at_least_2_patients." + which_dataset + ".txt"

    all_genes_with_at_least_2_patients_with_corresponding_SVs = []
    with open(gene_list_file) as list_file:
        for line in list_file:
            all_genes_with_at_least_2_patients_with_corresponding_SVs.append(line.strip())
            
    # # Remember to comment this out later:
    # all_genes_with_at_least_2_patients_with_corresponding_SVs = ["ASMT","WWOX","SPRY2","P2RY8","MYC"]

    # for d in [0,50000,100000,200000,500000,1000000]:
    # for d in list(range(0,5*10**6+1,100000)):
    # for d in list(range(0,5*10**6+1,1000000)):
    # for d in list(range(0,500000+1,100000)):
    # for d in list(range(500000,1000000+1,100000)):
    # for d in list(range(1000000,1500000+1,100000)):
    # for d in list(range(1500000,2000000+1,100000)):
    # for d in list(range(2000000,2500000+1,100000)):
    # for d in list(range(2500000,3000000+1,100000)):
    # for d in list(range(3000000,3500000+1,100000)):
    # for d in list(range(3500000,4000000+1,100000)):
    # for d in list(range(4000000,4500000+1,100000)):
    # for d in list(range(4500000,5000000+1,100000)):
    # for d in [47825]:
    # for d in [95649]:
    # for d in [944000]:
    
    # Steps to help prepare running the background model for SSVGAR:
    find_overall_exonic_and_intronic_portions_of_genome_setup()
    try:
        subprocess.check_output(["sort -k1,1 -k2,2n ../intermediate_files/v14/genes_from_updated_R_mat_file_start_and_end_coords.tsv > ../intermediate_files/v14/genes_from_updated_R_mat_file_start_and_end_coords_sorted.tsv"],shell=True)
    except subprocess.CalledProcessError as e:
        raise(Exception('subprocess command to sort gene start and end coordinates failed'))
    try:
        subprocess.check_output(["bedtools merge -i ../intermediate_files/v14/genes_from_updated_R_mat_file_start_and_end_coords_sorted.tsv > ../intermediate_files/v14/genes_from_updated_R_mat_file_start_and_end_coords_sorted_merged.tsv"],shell=True)
    except subprocess.CalledProcessError as e:
        raise(Exception('subprocess command to merge overlapping gene coordinates failed... Do you not have bedtools installed properly?'))
    try:
        subprocess.check_output(["sort -k1,1 -k2,2n ../intermediate_files/v14/exons_from_updated_R_mat_file_start_and_end_coords.tsv > ../intermediate_files/v14/exons_from_updated_R_mat_file_start_and_end_coords_sorted.tsv"],shell=True)
    except subprocess.CalledProcessError as e:
        raise(Exception('subprocess command to exon gene start and end coordinates failed'))
    try:
        subprocess.check_output(["bedtools merge -i ../intermediate_files/v14/exons_from_updated_R_mat_file_start_and_end_coords_sorted.tsv > ../intermediate_files/v14/exons_from_updated_R_mat_file_start_and_end_coords_sorted_merged.tsv"],shell=True)
    except subprocess.CalledProcessError as e:
        raise(Exception('subprocess command to merge overlapping exon coordinates failed... Do you not have bedtools installed properly?'))

    print("Time before running perform_LOF_test_for_list_of_genes:", time.time()-time_start)
    
    for d in [distance_threshold]:
        print("distance_thresold d:", d)
        results_df = perform_LOF_test_for_list_of_genes(list_of_genes=all_genes_with_at_least_2_patients_with_corresponding_SVs, \
                                        path_to_final_SV_list_for_analysis=SV_list_file, \
                                        which_project=which_dataset, \
                                        path_to_directory_for_annotations_surrounding_gene='../surrounding_annotations_for_each_gene/' + version, \
                                        threshold_of_proximity=d, \
                                        step_size_for_annotating_gene_regions = 1)
        if not os.path.exists("../model_results/" + which_dataset):
            os.makedirs("../model_results/" + which_dataset)
        results_df.to_csv("../model_results/" + which_dataset + "/LOF_model_results_for_all_genes_that_have_SVs_across_at_least_two_patients_threshold_of_proximity_" + str(d) + "bp" + "_" + version + ".tsv", sep='\t',index=False)

    print("Finished running.")
