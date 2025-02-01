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
import argparse
import subprocess

from produce_list_of_genes_annotated_to_SVs_recurring_across_at_least_2_patients_v14 import produce_list_of_genes_annotated_to_SVs_recurring_across_at_least_2_patients_v14
from producing_surrounding_annotations_for_relevant_genes_v14 import producing_surrounding_annotations_for_relevant_genes_v14
from SV_analysis_loss_of_function_events_v14 import produce_annotations_from_r_mat_file
from SV_analysis_loss_of_function_events_v14 import find_surrounding_annotations_for_gene_with_parameters_for_list_of_input_genes_post
from SV_analysis_loss_of_function_events_v14 import load_and_prepare_R_mat_file_into_df
from SV_analysis_loss_of_function_events_v14 import SV_analysis_loss_of_function_events_v14


### NEED TO SPECIFY PARAMETERS FOR WHICH YOU ARE RUNNING VIA THE COMMAND LINE
version = "v14" # any versions before v12 have to be run manually
parser = argparse.ArgumentParser()
parser.add_argument('--cohort', type=str, required=True, help="Project/cohort name")
parser.add_argument('--SV_list_path', type=str, required=True, help="Path (absolute or relative to the current directory) of the SV list to run SSVGAR on")
parser.add_argument('--distance_threshold', type=int, default=1000000, help="Windowing parameter in bp from gene edges used for SSVGAR (default: 1000000 for 1Mbp)")
args = parser.parse_args()
produce_list_of_genes_annotated_to_SVs_recurring_across_at_least_2_patients_v14(args.cohort, args.SV_list_path)
producing_surrounding_annotations_for_relevant_genes_v14(args.cohort)
SV_analysis_loss_of_function_events_v14(version, args.cohort, args.SV_list_path, args.distance_threshold)
subprocess.run(["mkdir ../model_results/" + args.cohort + "/input_SV_list_copy"], shell=True)
subprocess.run(["cp " + args.SV_list_path + " ../model_results/" + args.cohort + "/input_SV_list_copy"], shell=True)
