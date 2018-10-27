def _tsrc_generator(trajc, t0, dim_t):
    """
     By Aaron -- one of the two functions that achieve the same thing
     Approximately randomize over tstart, want something reproduceable
     Note that this generator only cares about trajectory, not series!
     And the final tsrc is given by tstart + t0.
     For 0.12fm ensembles, we use t0 = 0 and t0 = 17
    """
    series = 0 
    def prn_timeslice(traj, series,dim):
      return hash(hash(str(traj)+str(series)+"x")) % dim 
    # prn sometimes gets stuck in infinite loop
    # altering with rnchg gets out of loop
    # constructed in a way to prevent breaking those that did work after few iters
    rniter = 0 
    rnchg = ''
    tstart = -1
    while (tstart % 2 == 1): 
      tstart = prn_timeslice(str(int(trajc)+max(tstart,0))+rnchg,series,int(dim_t))
      rniter += 1
      if rniter % 10 == 9:
        rnchg = rnchg + str(tstart)
    return (int(tstart) + int(t0))%int(dim_t)

"""
For testing
"""
if __name__ == "__main__":
    config = "l4864f211b600m001907m05252m6382h-Coul.660.ildg"
    traj = 660
    t0 = 17
    tdim = 64
    tsrc = _tsrc_generator(traj, t0, tdim)
    print "config: %s"%config
    print "t0:     %s"%t0
    print "tsrc:   %s"%tsrc

