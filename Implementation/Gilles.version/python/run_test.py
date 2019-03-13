import numpy as np
import time, os
from flat_spike import flat_spike
from computeb_V2 import computeb_V2
from computeG_V2 import computeG_V2

def run_test(M, K, delta, Tmin, Tmax, ficname, debneur, outpath, sim_name, nmax=600):
    print()
    print('Test {0:s} - {1:d} neurons - computation of b and G with V2 method in python'.format(sim_name, M))
    print()
    print()
    if(len(outpath)):
        outpath += '/'
    t0 = time.time()
    DS = flat_spike(ficname, M, debneur=debneur, nmax=nmax)
    t1 = time.time()
    print('Elapsed time for flat_spike - to convert data {0:g}'.format(t1-t0))

    t0 = time.time()
    bV2 = computeb_V2(M, K, Tmin, Tmax, DS, delta)
    t1 = time.time()
    print('Elapsed time for to compute b with computeb_V2 {0:g}'.format(t1-t0))
    
    t0 = time.time()
#    print(bV2[1])
    GV2 = computeG_V2(M, K, Tmin, Tmax, DS, delta, bV2[1])
    t1 = time.time()
    print('Elapsed time for to compute G with computeG_V2 {0:g}'.format(t1-t0))
    
    os.makedirs(outpath, exist_ok=True)
    # ecriture des fichiers de sortie dans outpath
    np.savetxt(outpath+'b_'+sim_name+'.txt', bV2[0])
    np.savetxt(outpath+'G_'+sim_name+'.txt', GV2)

if __name__ == "__main__":
    M = 1
    delta = .2
    K = 5
    Tmin = 0.
    Tmax = 2.
    run_test(M, K, delta, Tmin, Tmax, '../data/n', debneur=1, outpath='/tmp/resu/', sim_name="3", nmax=100)
    M = 10
    delta = .02
    K = 5
    Tmin = 0.
    Tmax = 8.
    run_test(M, K, delta, Tmin, Tmax, '../data/M=1000/N', debneur=500, outpath='/tmp/resu/', sim_name="12", nmax=500)
#    M = 1000
#    delta = .02
#    K = 5
#    Tmin = 0.
#    Tmax = 8.
#    run_test(M, K, delta, Tmin, Tmax, '../data/M=1000/N', debneur=0, outpath='/tmp/resu/', sim_name="1", nmax=500)

    
