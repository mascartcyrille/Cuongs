# -*- coding: latin-1 -*-
# Algo proposed by PRB
import numpy as np
from misc import get_low_index, get_up_index

#def computeG_V2(M, K, Tmin, Tmax, T, neur, delta, depart, beta, low):
def computeG_V2(M, K, Tmin, Tmax, DS, delta, low):
    if(np.size(DS)):
        T = DS[0]; neur = np.array(DS[1], dtype='int32')
    else:
        T = np.array([])
    G = np.zeros((1+M*K, 1+M*K))
    G[0,0] = Tmax-Tmin
    A = K*delta
    #
    #   depart = integer part of Tmin-A
    # beta = integer part of Tmax
    #
    depart = get_low_index((Tmin-A), T)
    ##    beta = get_up_index(Tmax, T)
    beta = get_low_index(Tmax, T) - 1
    ##    print(depart, Tmin-A)
    ##    print(beta, Tmax)
    for i in np.arange(depart, beta+1):
        ti = T[i]
        l1 = neur[i]
        for k in np.arange(1,K+1):
            x1 = min(Tmax, ti+k*delta)
            x2 = max(Tmin, ti+(k-1)*delta)
            dx = x1 - x2
            if(dx>0):
                G[0, (l1-1)*K+k] += dx
                G[(l1-1)*K+k, 0] += dx
                
        #                print('error 1!!!!!!!!!!!')
        for j in np.arange(low[i],i):
            tj = T[j]
            l2 = neur[j]
            for k1 in np.arange(1,K+1):
                for k2 in np.arange(1,K+1):
                    x1 = min(Tmax, ti+k1*delta, tj+k2*delta)
                    x2 = max(Tmin, ti+(k1-1)*delta, tj+(k2-1)*delta)
                    dx = x1 - x2
                    if(dx>0):
                        G[(l1-1)*K+k1, (l2-1)*K+k2] += dx
                        G[(l2-1)*K+k2, (l1-1)*K+k1] += dx
                        #
                     #   print('error 2!!!!!!!!!!!')
                      #  print(Tmax,Tmin,ti+(k1+1)*delta, tj+(k2+1)*delta, ti+k1*delta, tj+k2*delta)

#diagonal part
        for k1 in np.arange(1,K+1):
            x1 = min(Tmax, ti+k1*delta)
            x2 = max(Tmin, ti+(k1-1)*delta)
            dx = x1 - x2
            if(dx >0):
                    G[(l1-1)*K+k1, (l1-1)*K+k1] += dx
            for k2 in np.arange(k1+1, K+1):
                x2 = max(Tmin, ti+(k2-1)*delta)
                dx = x1 - x2
                if(dx>0):
                    G[(l1-1)*K+k1, (l1-1)*K+k2] += dx
                    G[(l1-1)*K+k2, (l1-1)*K+k1] += dx
    return G
