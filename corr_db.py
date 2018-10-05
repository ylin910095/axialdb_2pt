from sqlalchemy import *
import math, bz2
import numpy as np
import yaml
from DB import *

########################################################################
"""
From utilies.py
"""

def readin_stream(input_stream):
    """
    Utility to either read in yaml file or dictionary. Return dictionary in both
    cases. The yaml file can contain multiple documents, and if this is the case,
    it will return a list of all dictionary for each document.
    """
    if isinstance(input_stream, str):
        if input_stream.endswith('.yaml'):
            try:
                stream = file(input_stream, 'r')
                param_dict = yaml.load(stream)
                default_dict = param_dict #based default param
            except:
                try:
                    stream = file(input_stream, 'r')
                    param_dict = yaml.load_all(stream)
                    default_dict = [i for i in param_dict] # List of dictionary of each document
                except:
                    raise ValueError('Cannot open file %s '%input_stream)
        else:
            raise ValueError('Unknow suffix for file %s'%input_stream)
    else:
        default_dict = input_stream
    return default_dict

def generate_tag_baryon(src_irrep, sink_irrep, src_class, sink_class, mom, mass,
                 ensemble):
    """
    Generate data tags for particular correlator keys to identify them in
    corrfitter. The goal is that we want the tags it generated only has
    information about the irreps we use for src/sink, the classes of those
    corresponding src/sink operators, the momentum for each quark, mass for
    each quark, and the ensemble we used.

    Possible list of src/sink_irrep:
        '8p': 8 irrep (fundamental)
        '8m': 8 prime irrep
        '16p'&'16m': 16 irrep (p/m represents plus/minus eigenvalues of R12)
                     This will not work with gather_data routine
        '16': 16 irrep. This will work with gather_data routine

    This is used for both tagging the data in corrfitter and the file names
    for the pickle files.
    """

    #assume momentum are same for all three if string provided
    if isinstance(mom, str):
        if len(mom) != 3:
            raise
        mom = (mom,mom,mom)

    if src_irrep != sink_irrep:
        #I should only give one irrep choice..
        raise ValueError('source irrep must be identical to sink irrep!')

    group_construct = '%s_s_%s_%s_s_%s_'%(
            src_irrep, src_class, sink_irrep, sink_class)
    mom_construct = 'p%s_p%s_p%s_'%(mom[0], mom[1], mom[2])
    mass_construct = 'm%s_m%s_m%s_'%(mass,mass,mass) #same masses for all

    return group_construct + mom_construct + mass_construct + ensemble

def parse_tag_baryon(datatag):
    """
    Inverse of generate_tag. Given datatag, return library with keys of
        (src_irrep, sink_irrep, src_class, sink_class, mom, mass, ensemble)
    """
    datatag_split = datatag.split('_')
    src_irrep  = datatag_split[0]
    src_class  = datatag_split[2]
    sink_irrep = datatag_split[3]
    sink_class = datatag_split[5]
    mom0       = (datatag_split[6])[1:]
    mom1       = (datatag_split[7])[1:]
    mom2       = (datatag_split[8])[1:]
    mom        = (mom0,mom1,mom2)
    mass       = float((datatag_split[9])[1:])
    ensemble   = datatag_split[12]

    data_dict = dict()
    data_dict['src_irrep']  = src_irrep
    data_dict['sink_irrep'] = sink_irrep
    data_dict['src_class']  = src_class
    data_dict['sink_class'] = sink_class
    data_dict['mom']        = mom
    data_dict['mass']       = mass
    data_dict['ensemble']   = ensemble
    return data_dict

########################################################################

class Lattice_Corrlator():
    def __init__(self, db_name, corr_name, datatag, fit_type, verbose=True):
        self.datatag = datatag
        self.fit_type = fit_type
        self.db_name = db_name
        self.corr_name = corr_name
        self.avg_tsrc = False
        self.block_no = None
        self.no_tsrc = None # Number of time sources per configuration
        self.configId = None
        self.raw_configId = None
        self.nconf = None
        self.nt = None
        self.raw_data = None
        self.output_data = None
        self._all_series = None
        self.verbose = verbose

        engine = create_engine('sqlite:///'+self.db_name)
        Session = sessionmaker(bind=engine)
        session = Session()

        # First query correlator ID
        corr_entry = session.query(Correlator).filter_by(name=corr_name).all()
        if len(corr_entry) > 1:
            raise ValueError("Error in retrieving '%s': more than more entries present in" +
                             "Correlator Column")
        else:
            corr_entry = corr_entry[0]

        if self.fit_type == 'baryon':
            data_entry = session.query(Datum).filter_by(correlator_id=corr_entry.id,
                                                        ).order_by(Datum.series,
                                                                   Datum.trajectory,
                                                                   Datum.tsrc).all()
        else:
            raise ValueError("Unknow fit type: %s" %self.fit_type)
        self.raw_configId = ['%s%s_t%s'%(idata.series,
                                   (str(idata.trajectory)).zfill(5),
                                   (str(idata.tsrc).zfill(3)))
                                   for idata in data_entry]

        # Query raw data
        self.raw_data = [np.array(bz2.decompress(idata.dataBZ2).split('\n'),dtype=np.float64)
             for idata in data_entry]
        self.nt = len(self.raw_data[0]) # Obtain T

        # Delete any identical data entries
        # Do not change self.raw_data after this
        del_indx_list = []
        
        for ica, (idataa, iconfiga) in enumerate(zip(self.raw_data, self.raw_configId)):
            if ica in del_indx_list:
                continue
            for icb, (idatab, iconfigb) in enumerate(zip(self.raw_data, self.raw_configId)):
                if icb <= ica:
                    continue
                if iconfiga == iconfigb:
                    tolerance = 1e-4 # Maximum tolerance for percentage difference
                    pdiff = np.sum(((idataa - idatab)/idataa))/len(idataa) # average percent difference
                    maxdiff = np.max((idataa - idatab)/idatab)
                    if abs(maxdiff) > tolerance or abs(pdiff) > tolerance:
                        if self.verbose:
                            print "------------------------------------------------------------------"
                            print "*****WARNING****** Correlators from two configurations found " +\
                                    "but NOT identical (%s)!" %iconfiga
                            print "max percent diff: %s, average percent diff: %s" %(maxdiff, pdiff)
                            print "This configuration will be all discarded!"
                            print ica, icb
                            savconfig = iconfiga
                            print "------------------------------------------------------------------"
                        del_indx_list.append(ica)
                        del_indx_list.append(icb)
                    else:
                        # If two correlator entries are identical up to specified tolerance, for the same configuration
                        # we only keep one of them
                        if self.verbose:
                            print '--> WARNING: Identical data present: %s vs %s Only one will be retained.' %(iconfiga, iconfigb)
                        del_indx_list.append(icb)
        
        self.raw_configId = [i for j,i in enumerate(self.raw_configId) if j not in del_indx_list]
        self.raw_data = [i for j,i in enumerate(self.raw_data) if j not in del_indx_list]


        self.raw_unique_configId = set([(i.split('_'))[0] for i in self.raw_configId])
        # Number of time sources for each correlator
        # Determine no_tsrc AFTER delete duplicate copies or it will be wrong
        self.no_tsrc = int(math.floor(float(len(self.raw_configId))/float(len(self.raw_unique_configId))))

        # Check if all data have same number of timeslices
        for idata in self.raw_data:
            if len(idata) != self.nt:
                raise ValueError('Some data have inconsistent self.nt! (%s != %s)'
                                 %(len(idata),self.nt))

        self.output_data = self.raw_data
        self.configId = self.raw_configId
        self.nconf = len(self.output_data)

    def block(self, block_no, avg_tsrc):
        """
        Perform blocking to output data for both trajectories and time source.
        This will automatically perform tsrcavg. It will overwrite the
        previous blocking results if called.
        """
        self.block_no = block_no
        if avg_tsrc:
            self.avg_tsrc = True
            _tblock_no = self.block_no * self.no_tsrc
        else:
            self.avg_tsrc = False
            _tblock_no = self.block_no
        _temp_out_data = []
        _temp_out_configId = []
        counter = 0
        for idata,iconfigid in zip(self.raw_data, self.raw_configId):
            current_series = iconfigid[0] # First character in the configId string
            current_series_traj = (iconfigid.split('_t'))[0]
            current_tsrc = int((iconfigid.split('_t'))[1])

            if counter%_tblock_no == 0: # Initilization for every block
                last_series = iconfigid[0]
                last_series_traj = None
                last_tsrc = None
                _hold_data = np.zeros([_tblock_no, self.nt])
                _hold_configId = []

            # Skip particular trajectory if we cannot find enough tsrc
            if ((counter%self.no_tsrc != self.no_tsrc-1)
                and (current_series_traj != last_series_traj)
                and (counter%self.no_tsrc != 0)):
                # Revert counter to erase this particular trajectory
                for revert_indx, _tmp in enumerate(_hold_configId):
                    if _tmp.startswith(last_series_traj):
                        break
                _hold_configId = _hold_configId[:revert_indx]
                counter -= revert_indx

            if current_series != last_series: # Restart for new series
                _hold_data = np.zeros([_tblock_no, self.nt]) # Dump the unused configId
                _hold_configId = []
                counter = 0

            if current_tsrc == last_series and current_series_traj == current_series_traj:
                raise ValueError("Trying to block identical correlator!")

            _hold_data[counter%_tblock_no,:] = idata
            _hold_configId.append(iconfigid)

            if counter%_tblock_no == _tblock_no-1: # Block data once we have enough configurations
                for dk in range(np.shape(_hold_data)[0]):
                    for uk in range(np.shape(_hold_data)[0]):
                        if dk != uk:
                            if np.array_equal(_hold_data[dk,:],_hold_data[uk,:]):
                                # Another safety check
                                raise ValueError("Two correlators have same data!")
                if len(_hold_configId) == 0:
                    raise ValueError("Nothing to block? Code error?")
                _temp_out_data.append(np.sum(_hold_data, axis=0)/(np.shape(_hold_data)[0]))
                _temp_out_configId.append('+'.join(_hold_configId))

                # Zero out after averaging
                _hold_data = np.zeros([_tblock_no, self.nt])
                _hold_configId = []

            # Ready for next iteration
            last_series = current_series
            last_seriest_traj = current_series_traj
            last_tsrc = current_tsrc
            counter += 1

        self.output_data = _temp_out_data
        self.configId = _temp_out_configId
        self.nconf = len(self.output_data)

    def get_data(self):
        return self.output_data

