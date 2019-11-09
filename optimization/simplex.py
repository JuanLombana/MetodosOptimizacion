import numpy as np
import numpy.linalg as npl

def optimize_simplex(c,A,b,sbfi):
    sbf = sbfi
    r=(sbf==False)
    B = np.copy(A[:,sbf])
    N = np.copy(A[:,r])
    cB=np.copy(c[sbf])
    cN=np.copy(c[r])
    BInv=npl.inv(B)
    pT=cN.T-np.dot(cB.T,np.matmul(BInv,N))
    index_sbf = np.argwhere(sbf).T[0]
    index_r = np.argwhere(r).T[0]
    print('------- ITERATION : 0 ---------')
    print('B: \r\n'+str(B))
    print('N: \r\n'+str(N))
    itera = 1
    while(len(pT[pT<0]) > 0):
        print('pT: \r\n'+str(pT))
        xb = np.dot(BInv,b)
        z=np.dot(cB.T,xb)
        print('x: \r\n'+str(index_sbf+1))
        print('xb: \r\n'+str(xb))
        print('Z: ' + str(z[0,0]))
        print('----------------')
        colN=np.argmin(pT)
        y = np.dot(BInv,N[:,colN])
        colB = np.argmin(xb.T[0]/y)
        newB = np.copy(B)
        newB[:,colB] = N[:,colN]
        newN = np.copy(N)
        newN[:,colN] = B[:,colB]
        B = np.copy(newB)
        N = np.copy(newN)
        print('------- ITERATION : '+str(itera)+' ---------')
        print('B: \r\n'+str(B))
        print('N: \r\n'+str(N))
        sbf[index_r[colN]]=True
        sbf[index_sbf[colB]]=False
        r=(sbf==False)
        new_cB = np.copy(cB)
        new_cB[colB,0]= cN[colN,0]
        new_cN = np.copy(cN)
        new_cN[colN,0]= cB[colB,0]
        cB = np.copy(new_cB)
        cN = np.copy(new_cN)
        BInv=npl.inv(B)
        pT=cN.T-np.dot(cB.T,np.matmul(BInv,N))
        index_sbf = np.argwhere(sbf).T[0]
        index_r = np.argwhere(r).T[0]
        itera=itera+1
    xb = np.dot(BInv,b)
    z=np.dot(cB.T,xb)
    index_sbf = np.argwhere(sbf).T[0]
    print('pT: \r\n'+str(pT))
    print('x: \r\n'+str(index_sbf+1))
    print('xb: \r\n'+str(xb))
    print('Z: ' + str(z[0,0]))
    print('--------END--------')