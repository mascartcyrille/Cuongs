# -*- coding: latin-1 -*-
# Algo proposed by PRB
import numpy as np
from flat_spike import flat_spike
from misc import get_k

def computeb_V2(M, K, Tmin, Tmax, DS, delta):

    if(np.size(DS)):
        T = DS[0]; neur = np.array(DS[1], dtype='int32')
    else:
        T = np.array([])
    Ntot = len(T) # total number of spikes
    A = K*delta  # scope
    b = np.zeros(((1+M*K), M))
    low = np.zeros(Ntot, dtype='int32')
    cnt = np.zeros(M, dtype='int32')    
    ilow  =  0
    eps = 1.e-12
    for i in np.arange(Ntot):
        t = T[i]
        r = neur[i] # neurone courant
        while(abs(T[ilow]-t)>A): # numerical precision pb!!!
            ilow = ilow +1
        low[i] = ilow

        if((Tmin<t) & (t<=Tmax)):  # precision numerique !!
#            print(r)
            cnt[r-1] += 1 
            for j in np.arange(ilow, i):
                if(abs(t-T[j])<eps):
                    continue # break is better (we have a sorted array)
            #            k = int((t-T[j])/delta)  # numerical precision
                k = get_k(t-T[j], delta)  # numerical precision
                l = neur[j]
                b[(l-1)*K+k, r-1] += 1
    b[0] = cnt
    return b, low, cnt

if __name__ == "__main__":
    M = 8
#    spike, tn = flat_spike('../data/n', M)
    spike = flat_spike('../data/n', M)
    # b est de taille M*(1+ M*K)
    delta = 0.2
    K = 5
    Tmin, Tmax = 0 , 2
#    b, mu, b3 = computeb_V2(M, K, Tmin, Tmax, spike, tn, delta)
    b, mu, b3 = computeb_V2(M, K, Tmin, Tmax, spike, delta)
    #print(tn)

    b1 = computeb_V2(1, 10, 0, 2., np.vstack([np.array([0.1,0.4, 0.45, 0.5, 0.6, 0.66]), np.ones(6, dtype='int8')]), delta=0.2)

#    b2, mu2, b23 = computeb_V2(M, K, Tmin, 1, spike, tn, delta)
    b2, mu2, b23 = computeb_V2(M, K, Tmin, 1, spike, delta)
#    print(b2)


    
