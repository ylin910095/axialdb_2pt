import time
import yaml
import os
import pickle
import numpy as np
from corr_db import *
import sys

#ordered direction of corner wall source to be used later by other modules
cw_dir_list = ['0','x','y','xy','z','zx','yz','xyz']

def _generate_correlator_keys_baryon(datatag, input_dict):

    """
    Generate list of correlator keys for all possible time sources and
    trajectories from the database given.
    This works for baryon correlators.
    """
    meta_dict = parse_tag_baryon(datatag)
    src_irrep  = meta_dict['src_irrep']
    sink_irrep = meta_dict['sink_irrep']
    src_class  = meta_dict['src_class']
    sink_class = meta_dict['sink_class']
    mom        = meta_dict['mom']
    mass       = meta_dict['mass']
    ensemble   = meta_dict['ensemble']

    if src_irrep != sink_irrep:
        #I should only give one irrep choice..
        raise ValueError('source irrep must be identical to sink irrep!')

    #this is fine as long as we dont exceed 9 units of mom in any directions.
    mom_int = (int(momi) for momi in mom)
    tot_mom = str(sum(mom_int)).zfill(3)

    #we always want to gather both 16m and 16p together anyway
    if src_irrep == '16m' or src_irrep == '16p' or sink_irrep == '16m' or \
       sink_irrep == '16p':
        raise ValueError('Source/sink irrep cannot be 16m or 16p. Use 16 instead')

    pprefix = "nd_b"

    if src_irrep != '16':
        prefix = pprefix + "_%s_s_%s_%s_s_%s_" %(
                src_irrep, src_class, sink_irrep, sink_class
                )
    else:
        prefix_p = pprefix + "_%s_s_%s_%s_s_%s_" %(
                '16p', src_class, '16p', sink_class
                )
        prefix_m = pprefix + "_%s_s_%s_%s_s_%s_" %(
                '16m', src_class, '16m', sink_class
                )
    suffix = "d_d_d_m%s_m%s_m%s" %(str(mass), str(mass), str(mass))

    #zero mom, old convention
    midfix = "cw0_cw0_cw0_"
    if src_irrep != '16':
        corr_template_list = [prefix + midfix + suffix]
    else:
        corr_template_list = [prefix_p + midfix + suffix,
                              prefix_m + midfix + suffix]
    key_list = []
    for corr_template in corr_template_list:
        key_list.append(corr_template)
    return key_list

##############################################################
#########################Main routine#########################
##############################################################

def main():
    if len(sys.argv) != 2:
        print 'usage:', sys.argv[0], 'input_processed.yaml'
        sys.exit(1)

    # Read in param files
    default_dict = readin_stream(sys.argv[1])
    gather_dataset(default_dict)

def gather_data(datatag, input_dict, out_format="gpl"):
    """
    Gather the set of data given by datatag.
    If data cache is found at directory `output_dir` and `input_dict`
    requires no overwrite, it will read it directly without queuing 
    the database; otherwise, it will dump files to 
    `output_dir`. The out_format accepts either `gpl` or `pickle`
    that determines the output format.

    Output:
        dictionary with raw data with datatags as keys 
    """
    if out_format != "gpl" and out_format != "pickle":
        raise ValueError("out_format needs to be gpl or pickle!")

    #tagging data for internal identification within the fitter
    key_list = _generate_correlator_keys_baryon(datatag, input_dict)
    print "datatag: %s" %(datatag)

    #save name for the dataset
    #TODO: add this data back to database to have a more coherent data
    #management.
    output_dir = input_dict['data_dir']
    save_name = output_dir + '/' + 'raw_' + datatag + '_' + 'baryon' +\
                "_tsrcavg" + str(int(input_dict['avg_tsrc'])) + "_blocking" + str(input_dict['blocking']) +\
                ".%s"%out_format #use datatag as file name
    meta_save_name = output_dir + '/' + 'meta_' + datatag + '_' + 'baryon' +\
                     "_tsrcavg" + str(int(input_dict['avg_tsrc'])) + "_blocking" + str(input_dict['blocking']) +\
                     '.%s'%out_format 
    #overwriting warning
    if (os.path.isfile(save_name) is False or os.path.isfile(meta_save_name) is False) or (input_dict['overwrite'] is True):
        if input_dict['overwrite'] and os.path.isfile(save_name):
            print 'WARNING: Overwriting existing file %s' %save_name

        tot = 0 #count number of configurations found
        dlist = [] #container for all data for given correlator
        if input_dict['op_irrep'] == '16':
            dlist_temp = []
        configid_list = [] #for safety check later
        start_time = time.time()

        # find blocking number
        try:
            blockno = int(input_dict['blocking'])
        except:
            blockno = 1
        if blockno == 1:
            print 'No blocking data!'
        else:
            print 'Block data by %s consecutive trajectories' %blockno
        for corr_name in key_list:
            print corr_name
            meta_info = []
            # Gather entries from database
            npt = Lattice_Corrlator(input_dict['db_name'], corr_name, datatag,
                                            'baryon', verbose=True)
            npt.block(block_no=blockno, avg_tsrc=input_dict['avg_tsrc'])

            configid_list.append(npt.configId)
            if meta_info == []: # Only need to do this once since they are all the same if done correctly
                meta_info = npt.configId
            if input_dict['op_irrep'] == '16':
                # For data with only multiple correlator names
                dlist_temp.append([])
                (dlist_temp[-1]) = npt.output_data
            else:
                # For data with only one correlator name only
                dlist.append(npt.output_data)

            try:
                # cast to numpy all_correlator_keys
                dlist_temp[-1] = np.array(dlist_temp[-1])
            except:
                pass
            tot += npt.nconf
        if input_dict['op_irrep'] == '16':
            # check if all these correlators we average come from the same gauge configurations
            if configid_list.count(configid_list[0]) != len(configid_list):
                raise ValueError('Error in gathering data! Possible errors in generating data!')

            tot = tot/len(dlist_temp)
            dlist = np.sum(np.array(dlist_temp),axis=0)/len(dlist_temp) # sum the raw value for meson and 16 = 16m + 16p
        print 'Total unique configurations (average over tsrc): %d' %tot
        if tot == 0:
            raise ValueError('No configurations found!')
        print "time to query: %.1fs" %((time.time()-start_time))

        data_dict = dict()
        meta_dict = dict()
        data_dict[datatag] = list(dlist)
        meta_dict[datatag] = meta_info
        
        if len(meta_info) != len(list(dlist)):
            raise ValueError("Inconsistent length between meta_info and dlist!")
        fio = open(save_name, 'wb')
        fio_meta = open(meta_save_name, 'wb')
        # Dump to pickle cache file
        if out_format == "pickle":
            pickle.dump(data_dict, fio)
            pickle.dump(meta_dict, fio_meta)
        elif out_format == "gpl":
            dump_gpl(data_dict, meta_dict, fio, fio_meta)
        fio.close()
        fio_meta.close()
        print 'data file saved: %s' %(save_name)
        print 'meta file saved: %s' %(meta_save_name)
    else: #overwrite == False
        fio = open(save_name,'r')
        fio_meta = open(meta_save_name,'r')
        #load the file if found in directory
        print 'Loading existing file: %s ... ' %save_name
        print 'Loading existing meta: %s ... ' %meta_save_name
        if out_format == "pickle":
            data_dict =  pickle.load(fio)
            meta_dict =  pickle.load(fio_meta)
        elif out_format == "gpl":
            data_dict, meta_dict = load_gpl(fio, fio_meta)
        fio.close()
        fio_meta.close()
    return data_dict, meta_dict

def gather_dataset(input_dict):
    """
    Gather a set data according to keywords in input_dict.
    If data cache is found at directory `output_dir` and `input_dict`
    requires no overwrite, it will read it directly without queuing 
    the database; otherwise, it will dump pickle files to 
    `output_dir`

    Output:
        dictionary with raw data with datatags as keys 
    """
    dlist_dict = dict() 
    meta_dict_all = dict()
    # First construct all datatags based on input_dict
    input_dict['datatag_list'] = []
    for src_class in input_dict['src_class_list']:
        for sink_class in input_dict['sink_class_list']:
            datatag = generate_tag_baryon(input_dict['op_irrep'], 
                                          input_dict['op_irrep'],
                                          src_class, sink_class,
                                          '000', mass=input_dict['mass'],
                                          ensemble=input_dict['ensemble'])
            input_dict['datatag_list'].append(datatag)

    # Gather all data
    for datatag in input_dict['datatag_list']:
        data_dict, meta_dict = gather_data(datatag, input_dict)
        dlist_dict[datatag] = data_dict[datatag]
        meta_dict_all[datatag] = meta_dict[datatag]

    return dlist_dict, meta_dict_all

def dump_gpl(data_dict, meta_dict, fio, fio_meta):
    """
    Dump correlators and meta information to `fio` and `fio_meta` text files
    that have white space separating each timeslice and newline separating
    each measurements
    """
    if len(data_dict) != 1:
        raise ValueError("data_dict should only have one key!")
    for datatag in data_dict:
        # Write data
        for iconfig in data_dict[datatag]:
            fio.write(datatag + " ") # write datatag
            strlist = ["%.15e"%x for x in iconfig]
            strout = " ".join(strlist) + "\n"
            fio.write(strout)
        # Write meta information
        for imeta in meta_dict[datatag]:
            fio_meta.write("%s\n"%imeta)

def load_gpl(fio, fio_meta):
    """
    Load the text files created by `dump_txt`. Return a `data_dict` and 
    `meta_dict`
    """
    dat_lines = [(line.rstrip('\n')).split(" ") for line in fio]
    meta_lines = [(line.rstrip('\n')).split(" ") for line in fio_meta]
    if len(dat_lines) != len(meta_lines):
        raise ValueError("Mistmatch in number of lines for data and metadata!")
    datatag = (dat_lines[0])[0] # all files should contain one datatag
    data_dict = {datatag:[]}
    meta_dict = {datatag:[]}
    for datl,metal in zip(dat_lines, meta_lines):
        if datl[0] != datatag:
            raise ValueError("There are more than one datatag for gpl files!")
        data_dict[datatag].append([float(x) for x in datl[1:]])
        meta_dict[datatag].append(metal[0])
    return data_dict, meta_dict
        
if __name__ == '__main__':
    main()