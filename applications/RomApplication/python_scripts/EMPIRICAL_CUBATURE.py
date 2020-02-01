## Python implementation of the Empirical Cubature method
## 20-Sep-2019
import numpy as np
import pdb
from matplotlib import pyplot as plt

def EmpiricalCubatureMethod(u,s,W,DATA):
    if 'IncludeSingularValuesF' not in DATA.keys():
        DATA['IncludeSingularValuesF'] = 1
    
    if DATA['IncludeSingularValuesF'] == 1:
        #pdb.set_trace()
        #G = u[...,:] * s
        G = u[...,:] * np.ones(len(s))
        G = G.T
        # Adding the row of ones
        G = np.vstack([ G , np.ones( np.shape(G)[1] )]  )        #Perhaps I should put this outside the script, (Add the row of ones from calling script)
        #b = G @ np.sqrt(W)
        b = G @ W #The weights are 1 for the case of element selection
        bEXACT = b
    else:
        G = u.T
        b = (G @ np.sqrt(W)).T
        bEXACT = np.multiply(b,s)    
    nbEXACT = np.linalg.norm(bEXACT)

    Gnorm = np.sqrt(sum(np.multiply(G, G), 0))
    M = np.shape(G)[1]
    DATA['TOL']= 1e-4   ##Adding a tolerance for truncation
    TOL = DATA['TOL']

    # INITIALIZATIONS
    z = {}  # set of intergration points
    # Set of candidate points (those whose associated column has low norm are removed)
    y = np.arange(0,M,1)
    DATA['TOLFilterCandidatePoints'] = 1e-16
    ## Making a new norm of G, to eliminate the candidate elements with zero contribution
    GnormNOONE = np.sqrt(sum(np.multiply(G[:-1,:], G[:-1,:]), 0)) 
    if DATA['TOLFilterCandidatePoints'] > 0:
       TOL_REMOVE = DATA['TOLFilterCandidatePoints'] * np.linalg.norm(b)
       rmvpin = np.where(GnormNOONE[y] < TOL_REMOVE)
       y = np.delete(y,rmvpin)
       
    mPOS = 0 # Numbre of nonzero weights
    r = b # residual vector
    k = 1 # number of iterations   
    # Default number of points
    DATA['npoints'] = len(b)
    m = min(DATA['npoints'], len(b) )
    # END INITIALIZATIONS
    #####################################################
    normB = np.linalg.norm(b)
    nerror = np.linalg.norm(r)/normB
    H = np.array([]) # Inverse of (Gz.T @ Gz)
    nerrorACTUAL = nerror

    while nerrorACTUAL > TOL and mPOS < m and len(y) != 0: 
        #Step 1. Compute new point        
        ObjFun = G[:,y].T @ r.T
        ObjFun = ObjFun.T / Gnorm[y]
        indSORT = np.argmax(ObjFun)
        maxLOC = max(ObjFun)   
        i = y[indSORT]
        if k==1:
            #alpha = np.linalg.lstsq(G[:, [i]], b, rcond=None)[0]
            alpha = np.linalg.lstsq(G[:, [i]], b)[0]
            H = 1/(G[:,i] @ G[:,i].T)
        else:
            H, alpha = UpdateWeightsInverse(G[:,z],H,G[:,i],alpha,r)
        #Step 3. Move i from set y to set z
        if k == 1:
            z = i
        else:
            z = np.r_[z,i]
        y = np.delete(y,indSORT)
        # Step 4. Find possible negative weights
        
        if any(alpha < 0):
            #pdb.set_trace()
            print("WARNING: NEGATIVE weight found")
            indexes_neg_weight = np.where(alpha <= 0.)[0]
            y = np.append(y, (z[indexes_neg_weight]).T)
            z = np.delete(z, indexes_neg_weight)
            H = multiupdate_inverse_hermitian(H, indexes_neg_weight)
            #alpha = np.dot(H, np.dot(G[:, z].T, b))
            alpha = H @ (G[:, z].T @ b)
            alpha = alpha.reshape(len(alpha),1)
        #Step 6 Update the residual        
        if len(alpha)==1:
            r = b - (G[:,z] * alpha)
        else:
            Aux = G[:,z] @ alpha
            r = np.squeeze(b - Aux.T)     
        nerror = np.linalg.norm(r) / np.linalg.norm(b)  # Relative error (using r and b)

        if DATA['IncludeSingularValuesF']==0:
            nerrorACTUAL = s*r
            nerrorACTUAL = np.linalg.norm(nerrorACTUAL/nbEXACT)
        else:
            nerrorACTUAL = nerror        
        # STEP 7
        mPOS = np.size(z)       
        print(f'k = {k}, m = {np.size(z)}, error n(res)/n(b) (%) = {nerror*100},  Actual error % = {nerrorACTUAL*100} ')

        if k == 1:
            ERROR_GLO = np.array([nerrorACTUAL])
            NPOINTS = np.array([np.size(z)])
        else:
            ERROR_GLO = np.c_[ ERROR_GLO , nerrorACTUAL]
            NPOINTS = np.c_[ NPOINTS , np.size(z)]        
        
        k = k+1

        # if k==np.shape(u)[1]:            
        #     break
            

        

    w = alpha.T * np.sqrt(W[z])
    print(f'Total number of iterations = {k}')

    ### PLotting
    plt.plot(NPOINTS[0], ERROR_GLO[0])
    plt.title('Error Evolution')
    plt.xlabel('Number of points')
    plt.ylabel('Error %')
    plt.show()

    return G, z, w 


def UpdateWeightsInverse(A,Aast,a,xold,r):
    #No implementation for no args
    c = np.dot(A.T, a)
    d = np.dot(Aast, c).reshape(-1, 1)
    s = np.dot(a.T, a) - np.dot(c.T, d)
    aux1 = np.hstack([Aast + np.outer(d, d) / s, -d / s])
    if np.shape(-d.T / s)[1]==1:
        aux2 = np.squeeze(np.hstack([-d.T / s, 1 / s]))
    else:
        aux2 = np.hstack([np.squeeze(-d.T / s), 1 / s])
    Bast = np.vstack([aux1, aux2])
    v = np.dot(a.T, r) / s
    x = np.vstack([(xold - d * v), v])
    return Bast, x

def multiupdate_inverse_hermitian(invH, neg_indexes):
    neg_indexes = np.sort(neg_indexes)
    for i in range(np.size(neg_indexes)):
        neg_index = neg_indexes[i] - i
        invH = update_inverse_hermitian(invH, neg_index)
    return invH

def update_inverse_hermitian(invH, neg_index):
    if neg_index == np.shape(invH)[1]:
        aux = (invH[0:-1, -1] * invH[-1, 0:-1]) / invH(-1, -1)
        invH_new = invH[:-1, :-1] - aux
    else:
        aux1 = np.hstack([invH[:, 0:neg_index], invH[:, neg_index + 1:], invH[:, neg_index].reshape(-1, 1)])
        aux2 = np.vstack([aux1[0:neg_index, :], aux1[neg_index + 1:, :], aux1[neg_index, :]])
        invH_new = aux2[0:-1, 0:-1] - np.outer(aux2[0:-1, -1], aux2[-1, 0:-1]) / aux2[-1, -1]
    return invH_new

