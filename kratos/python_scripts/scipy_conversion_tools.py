from __future__ import print_function, absolute_import, division

# Import Kratos
import KratosMultiphysics

# Import scipy
import scipy
import scipy.sparse

def to_csr(A):
    Ascipy = scipy.sparse.csr_matrix((A.value_data(), A.index2_data(), A.index1_data()), shape=(A.Size1(), A.Size2()))
    return Ascipy

