from corr_db import *
from gather_data import *

yaml_file = "./gather_012fm.yaml"
default_dict = readin_stream(yaml_file)
data_dict, meta_dict = gather_dataset(default_dict)

for i in range(10):
    for key in data_dict:
        print (data_dict[key])[i]
        print (meta_dict[key])[i]
        break # Only print the first correlator