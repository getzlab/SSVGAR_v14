import numpy as np
import pandas as pd

CGC_df = pd.read_csv('Census_allTue_Feb_14_19_31_20_2023.tsv',sep='\t')
boolean_list_based_off_somatic_tumor_type = list(map(lambda x: x == x and ("other cancer" in x or "other tumour types" in x or "NHL" in x or "DLBCL" in x), list(CGC_df['Tumour Types(Somatic)'])))
boolean_list_based_off_role_in_cancer = list(map(lambda x: x == x and ("fusion" in x or "TSG" in x), list(CGC_df['Role in Cancer'])))
filtered_CGC_df = CGC_df[np.array(boolean_list_based_off_somatic_tumor_type) & np.array(boolean_list_based_off_role_in_cancer)]
list_of_genes_from_CGC_with_connection_to_NHL_or_DLBCL = list(filtered_CGC_df['Gene Symbol'])
list_of_genes_from_CGC_with_connection_to_NHL_or_DLBCL.remove('IGH')
list_of_genes_from_CGC_with_connection_to_NHL_or_DLBCL.remove('IGK')

f = open('../genes_from_CGC_with_connection_to_NHL_or_DLBCL.txt',mode='w')
f.write('# To obtain this list, I downloaded the content of the Cancer Gene Census on 2/14/2023, saved it as scripts_and_files_for_processing/Census_allTue_Feb_14_19_31_20_2023.tsv,\n# and ran come_up_with_list_of_genes_from_CGC_with_connection_to_NHL_or_DLBCL.py\n')
for gene in list_of_genes_from_CGC_with_connection_to_NHL_or_DLBCL:
    f.write(gene + '\n')
