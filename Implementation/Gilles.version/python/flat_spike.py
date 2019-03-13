# -*- coding: latin-1 -*-
"""On transforme ici plusieurs fichiers n*.txt contenant les spikes de chaque neurone respectif pour créer un seul tableau de sortie"""
import numpy as np

def get_digit(n):
    num = 1
    while(int(n/10)>0):
        n = n/10
        num += 1
    return num

def flat_spike(mystr, M, debneur=1, nmax=20):
    """ mystr est le début des noms de fichiers de spikes. nmax est le nb max de spikes pour 1 neurone , les extensions des fichiers d'entrée sont en .txt"""
# 8 neurones
# 20 = grande valeur pour réserver de la place
    D = np.ones((M,nmax))*np.nan
    N = get_digit(debneur+M-1)
    for i in np.arange(debneur,debneur+M):
        if(N==6):
            if(i<10):
                tmp = np.loadtxt(mystr+'00000'+str(i)+'.txt') 
            elif(i<100):
                tmp = np.loadtxt(mystr+'0000'+str(i)+'.txt')
            elif(i<1000):
                tmp = np.loadtxt(mystr+'000'+str(i)+'.txt')
            elif(i<10000):
                tmp = np.loadtxt(mystr+'00'+str(i)+'.txt')
            elif(i<100000):
                tmp = np.loadtxt(mystr+'0'+str(i)+'.txt')
            elif(i<1000000):
                tmp = np.loadtxt(mystr+str(i)+'.txt')                                                
        elif(N==5):
            if(i<10):
                tmp = np.loadtxt(mystr+'0000'+str(i)+'.txt') 
            elif(i<100):
                tmp = np.loadtxt(mystr+'000'+str(i)+'.txt')
            elif(i<1000):
                tmp = np.loadtxt(mystr+'00'+str(i)+'.txt')
            elif(i<10000):
                tmp = np.loadtxt(mystr+'0'+str(i)+'.txt')
            elif(i<100000):
                tmp = np.loadtxt(mystr+str(i)+'.txt')
        elif(N==4):
            if(i<10):
                tmp = np.loadtxt(mystr+'000'+str(i)+'.txt') 
            elif(i<100):
                tmp = np.loadtxt(mystr+'00'+str(i)+'.txt')
            elif(i<1000):
                tmp = np.loadtxt(mystr+'0'+str(i)+'.txt')
            elif(i<10000):
                tmp = np.loadtxt(mystr+str(i)+'.txt')
        elif(N==3):
            if(i<10):
                tmp = np.loadtxt(mystr+'00'+str(i)+'.txt') 
            elif(i<100):
                tmp = np.loadtxt(mystr+'0'+str(i)+'.txt')
            elif(i<1000):
                tmp = np.loadtxt(mystr+str(i)+'.txt')
        elif(N==2):
            if(i<10):
                tmp = np.loadtxt(mystr+'0'+str(i)+'.txt') 
            elif(i<100):
                tmp = np.loadtxt(mystr+str(i)+'.txt')
        else:
            tmp = np.loadtxt(mystr+str(i)+'.txt')
#        print(tmp)
#        D[i-1,0] = len(tmp)
#        D[i-1,1:(len(tmp)+1)] = tmp
        D[i-debneur,0] = len(tmp)
        D[i-debneur,1:(len(tmp)+1)] = tmp
        del tmp

#    print(D)
# D serait un DataNeur à la R? (vérifier)
#
    flat = np.ones(M*nmax)*np.nan
    tn = np.ones(M*nmax)*np.nan
    deb = int(0)
    for i in np.arange(M):
        flat[deb:(deb+int(D[i,0]))] = D[i,1:(int(D[i,0])+1)]
        tn[deb:(deb+int(D[i,0]))] = i+1 #int(i+1)
        deb += int(D[i,0])
    
    b = np.sort(flat)[0:deb]
    tn2 = tn[np.argsort(flat) ][0:deb].astype(int) # numpy >= 1.9
#    return b, tn2
    return np.vstack([b, tn2])

# DataSpike to DataNeur
def DataNeur(DS):
    if(np.size(DS)):
        T = DS[0]; neur = np.array(DS[1], dtype='int32')
    else:
        return np.array([])
    u_neur = np.unique(neur)
    irow = 0
#    ind = np.ones((len(u_neur),np.size(DS,1)))*np.nan
    ind = list()
    nmax = -1
    for r in u_neur:
        p = np.flatnonzero(neur==r)
        ind.append(list(p)) # [irow,0:len(p)] = p
        nmax = max(nmax, len(T[p]))
        irow += 1
    npind = np.array(ind)
    print(npind)
    ld = list() #DN = np.zeros((irow, nmax+1))
    for i in np.arange(irow):
#        print(i)
##        DN[irow] = np.append(np.array([len(T[indr])]), T[indr])
##        DN[irow] = 
##        (len(T[indr]), T[indr])
##        print(a)
#        DN[i,0:len(T[ind[i]])+1] = np.array([len(T[ind[i]]), T[ind[i]]])
        ld.append(list(np.append(len(T[npind[i]]), T[npind[i]])))
        #    print(ind)
    DN = np.array(ld)
    return DN
        #    return DN
    # attention à np.sort qui peut contenir des erreurs numériques

if __name__ == "__main__":
    M = 5
    spike = flat_spike('../data/n', M)
    DN = DataNeur(spike)
    print(DN)
