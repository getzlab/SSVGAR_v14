import pandas as pd
import numpy as np
import scipy.io as sio

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




cytoband_df = pd.read_csv('cytoBand.txt',sep='\t')
cytoband_df['chr'] = list(map(lambda x: convert_chr_xavi(x), list(cytoband_df['chr'])))

R_mat_file_path="../annotation_file/hg38_R_with_hugo_symbols_with_DUX4L1_HMGN2P46_MALAT1.mat"
R_mat_loaded_and_prepared_df = load_and_prepare_R_mat_file_into_df(R_mat_file_path)

genes_to_cytobands_dict = {}
iteration_counter = 0
for gene in set(R_mat_loaded_and_prepared_df['gene']):
    if iteration_counter % 1000 == 0:
        print("iteration_counter: " + str(iteration_counter) + '/' + str(len(set(R_mat_loaded_and_prepared_df['gene']))))
    iteration_counter += 1
    gene_chr = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]['chr'].iloc[0]
    gene_start = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]['gene_start'].min()
    gene_end = R_mat_loaded_and_prepared_df[R_mat_loaded_and_prepared_df['gene'] == gene]['gene_end'].max()
    cytoband_df_subset = cytoband_df[(cytoband_df['chr'] == gene_chr) & (((cytoband_df['start'] <= gene_start) & (cytoband_df['end'] >= gene_start)) | ((cytoband_df['start'] <= gene_end) & (cytoband_df['end'] >= gene_end)))]
    if len(cytoband_df_subset) == 0:
        genes_to_cytobands_dict[gene] = 'N/A'
    elif len(cytoband_df_subset) == 1:
        genes_to_cytobands_dict[gene] = str(cytoband_df_subset['chr'].iloc[0]) + str(cytoband_df_subset['band'].iloc[0])
    elif len(cytoband_df_subset) == 2:
        if cytoband_df_subset.iloc[0]['end']-gene_start >= gene_end-cytoband_df_subset.iloc[1]['start']:
            genes_to_cytobands_dict[gene] = str(cytoband_df_subset['chr'].iloc[0]) + str(cytoband_df_subset['band'].iloc[0])
        else:
            genes_to_cytobands_dict[gene] = str(cytoband_df_subset['chr'].iloc[1]) + str(cytoband_df_subset['band'].iloc[1])
    else:
        raise("This shouldn't be possible")


genes_to_cytobands_df = pd.DataFrame({"gene":list(genes_to_cytobands_dict.keys()),"cytoband":list(genes_to_cytobands_dict.values())})
genes_to_cytobands_df.to_csv("genes_to_cytobands.tsv",sep='\t',index=False)
