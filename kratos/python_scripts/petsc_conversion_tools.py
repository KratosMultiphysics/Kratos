# Import Kratos
import KratosMultiphysics
import numpy as np

import petsc4py
from petsc4py import PETSc

#gives a view of Kratos DistributedCSRMatrix as a petsc matrix.
#values are NOT copied. 
#Indices are copied ONLY if kratos IndexType and the one used within Petsc (by defauly int32) do not coincide

class PetscCSRWrapper:
    def __init__(self, A):  #Kratos DistributedCsrMatrix
        self.A = A 
        m = A.local_size1()
        M = A.GetRowNumbering().Size()
        n = A.GetColNumbering().LocalSize()
        N = A.GetColNumbering().Size()

        ##pointers to the kratos arrays of values
        self.a = np.array(A.GetDiagonalBlock().value_data(), copy=False)
        self.oa = np.array(A.GetOffDiagonalBlock().value_data(), copy=False)

        if(type(m) == type(PETSc.IntType)): #we can avoid doing copies
            #indices of diagonal block
            self.i = np.array(A.GetDiagonalBlock().index1_data(), copy=False)
            self.j = np.array(A.GetDiagonalBlock().index2_data(), copy=False)
            
            #indices of offdiagonal block
            self.oi = np.array(A.GetOffDiagonalBlock().index1_data(), copy=False)
            self.oj = np.array(A.GetOffDiagonalIndex2DataInGlobalNumbering(), copy=False)        

        else: #WE NEED TO COPY INDICES! since the size of the IndexType is not compatible between Kratos and PETSc
            if(A.GetComm().Rank() == 0):
                print("WARNING: a copy has been needed in the conversion to petsc since the IndexType did not match between Kratos and PETSc")
            #indices of diagonal block
            self.i = np.array(A.GetDiagonalBlock().index1_data(), dtype="int32")
            self.j = np.array(A.GetDiagonalBlock().index2_data(), dtype="int32")
            
            #indices of offdiagonal block
            self.oi = np.array(A.GetOffDiagonalBlock().index1_data(), dtype="int32")
            self.oj = np.array(A.GetOffDiagonalIndex2DataInGlobalNumbering(), dtype="int32")
        
        size = ((m, M), (n, N)) # local,global row,column sizes of the matrix
        csr = ((self.i, self.j, self.a), (self.oi, self.oj, self.oa)) # numpy arrays with the diagonal and off-diagonal arrays containing the CSR (or AIJ) structure
        self.Apetsc = PETSc.Mat().createAIJWithArrays(size, csr)

    def __del__(self):
        del self.Apetsc #we need to ensure that this is deleted before the rest

    def GetPETScMatrix(self):
        return self.Apetsc

