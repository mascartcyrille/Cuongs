""" Auxiliary functions """
import numpy as np
def get_k(x, delta, eps=1e-12):
    """ Computation of k such that x belongs to ](k-1)*delta, k*delta] """
    # si x/delta est un entier on doit retourner cet entier - 1
    k = int(x/delta)+1
    if (abs(k*delta - x)<eps):
        return k
    if (abs((k-1)*delta - x)<eps):
        return k-1
    return k


def get_low_index(x, T):
    """ Computation of ind, the first index in T such that T[ind]>x """
    """ Possible values of ind are between 0 and len(T) (included) """
    ind = 0
    n = len(T)
    while(ind <n):
        if(T[ind]<=x):
#            print(ind,T[ind])
            ind += 1
        else:
            break
    return ind

def get_up_index(x, T):
    """ Computation of ind, the last index in T such that T[ind]<=x """
    """ Possible values of ind are between -1 and (len(T)-1) (included) """
    ind = len(T)-1
    while(ind >=0):
        if(T[ind]>x):
            ind -= 1
        else:
            break
    return ind

def get_k1k2(x, delta, K, eps=1e-12):
    chi = int(x/delta)
    if(abs((chi+1)*delta - x)< eps):
        chi += 1
    M = np.array([-1, -1], dtype='int32')
    for n in np.arange(2+chi, K+1):
        M = np.vstack([M, [n, n-chi-1]])
    for n in np.arange(1+chi, K+1):
        M = np.vstack([M, [n, n-chi]])        
    M = M[1:]    
    return M


if __name__ == "__main__":
    delta = .2
    print('')
    print("get_k(delta, delta) -> should be 1")
    print(get_k(delta, delta))
    print('')
    print("get_k(1.5*delta, delta) -> should be 2")
    print(get_k(1.5*delta, delta))
    print('')
    print("get_k(3*delta, delta) -> should be 3")
    print(get_k(3*delta, delta))
    print('')
    print("get_k(2.02*delta, delta) -> should be 3")
    print(get_k(2.02*delta, delta))          
    print('')
    print('')
    print('')    
    print("get_low_index(0.2, [.1, .12, .31, .45, .76]) -> should be 2")
    print(get_low_index(0.2, [.1, .12, .31, .45, .76]))
    print('')
    print("get_low_index(0.8, [.1, .12, .31, .45, .76]) -> should be 5")
    print(get_low_index(0.8, [.1, .12, .31, .45, .76]))
    print('')
    print("get_low_index(0.04, [.1, .12, .31, .45, .76]) -> should be 0")
    print(get_low_index(0.04, [.1, .12, .31, .45, .76]))
    print('')
    print("get_low_index(0.31, [.1, .12, .31, .45, .76]) -> should be 3")
    print(get_low_index(0.31, [.1, .12, .31, .45, .76]))
    print('')
    print("get_low_index(0.3099999999999999, [.1, .12, .31, .45, .76]) -> should be 3")
    print(get_low_index(0.3099999999999999, [.1, .12, .31, .45, .76]))
    print('')
    print("get_low_index(0.1, [.1, .12, .31, .45, .76]) -> should be 1")
    print(get_low_index(0.1, [.1, .12, .31, .45, .76]))
    print('')
    print("get_low_index(0.76, [.1, .12, .31, .45, .76]) -> should be 5")
    print(get_low_index(0.76, [.1, .12, .31, .45, .76]))
    print('')
    print('')
    print("get_low_index(0.77, [.1, .12, .31, .45, .76]) -> should be 5")
    print(get_low_index(0.77, [.1, .12, .31, .45, .76]))
    print('')
    print('')
    print('')    
    print("get_up_index(0.2, [.1, .12, .31, .45, .76]) -> should be 1")
    print(get_up_index(0.2, [.1, .12, .31, .45, .76]))
    print(get_low_index(0.2, [.1, .12, .31, .45, .76])-1)
    print('')
    print("get_up_index(0.8, [.1, .12, .31, .45, .76]) -> should be 4")
    print(get_up_index(0.8, [.1, .12, .31, .45, .76]))
    print(get_low_index(0.8, [.1, .12, .31, .45, .76])-1)
    print('')
    print("get_up_index(0.04, [.1, .12, .31, .45, .76], 0.04) -> should be -1")
    print(get_up_index(0.04, [.1, .12, .31, .45, .76]))
    print(get_low_index(0.04, [.1, .12, .31, .45, .76])-1)
    print('')
    print("get_up_index(0.76, [.1, .12, .31, .45, .76], 0.04) -> should be 4")
    print(get_up_index(0.76, [.1, .12, .31, .45, .76]))
    print(get_low_index(0.76, [.1, .12, .31, .45, .76])-1)
    print('')
    print("get_up_index(0.1, [.1, .12, .31, .45, .76], 0.04) -> should be 0")
    print(get_up_index(0.1, [.1, .12, .31, .45, .76]))
    print(get_low_index(0.1, [.1, .12, .31, .45, .76])-1)
