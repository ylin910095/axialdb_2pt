# axialdb_2pt
Scripts for reading staggered baryon two point correlators from databases. 


## Data
The correlators are generated with __Coulomb gauge fixed__ HISQ ensembles:

l3248f211b580m002426m06730m8447

l4864f211b600m001907m05252m6382

l6496f211b630m0012m0363m432

The correlators are based on the operators found TABLE VIII (or equivalently for isospin doublet nucleon, TABLE VI,) of https://arxiv.org/abs/hep-lat/0611023

In particular, the source/sink operators we constructed are characterized by irreducible representations (16, 8^+, or 8^-) and classes

## Usage

### Prerequisites 
- python 2.7 (tested)
- numpy
- sqlalchemy
- pyyaml

The main interface is the function `gather_dataset` in gather_data.py. It accepts one input argument, `input_dict`, which is a python dictionary that contains various parameters to query the database. There are serval import keys that need to be specificed for the routine to work

- 'db_name': location of the database
- 'data_dir': location of the data cache files
- 'overwrite': True or False. Do you want to always overwrite the data cache files?
- 'ensemble':  one of three ensembles above 
- 'mass': mass of light quarks
- 'op_irrep': Irrep for the source/sink operators
- 'sink_class_list': List of classes of sink operators 
- 'src_class_list': List of classes of source operators
- 'blocking': Number of consecutive blocking of the raw correlators (blocking = 1 is no blocking)
- 'avg_tsrc': True or False. Do you want to average all the time sources for a given gauge configuration?

Usually, these parameters are put into an yaml file and can be read to python dictionary using `readin_stream` function found in corr_db.py. For an example of yaml file, see gather_012fm.yaml

`gather_dataset` will return two python dictionaries. Both dictionaries have keys given by the returned string of `generate_tag_baryon` according to the correlators you query. These keys are called datatags and are used extensively to identity the correlators within the program. 

For a given key, the first dictionary will return a list of raw correlator queried for a given configuration from the database __AFTER__ blocking and time source averaging. The configuration information is given by the second dictionary. The configuration is a string that contains series, trajectory, and time source for the given correlator. 

An example will be 'a00110_t036+a00115_t038'. If blocking or time source averaging, the string can be separated by '+' character. In this case, we are blocking two configurations: series a, trajectory 110, time source 36 and series a, trajectory 115, time source 38. 

### Example
A typical usage will look something like

```python
from corr_db import *
from gather_data import *

yamlfn = "./location/of/yaml/test.yaml" 
input_dict = readin_stream(yamlfn) 
data_dict, meta_dict = gather_dataset(input_dict) 

# Do something else...
```
