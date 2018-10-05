from corr_db import *
from gather_data import *
yamlfn = "./gather_012fm.yaml"
input_dict = readin_stream(yamlfn) 
data_dict, meta_dict = gather_dataset(input_dict) 