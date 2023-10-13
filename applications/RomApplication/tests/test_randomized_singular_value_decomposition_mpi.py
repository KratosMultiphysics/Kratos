#import python class test
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities

#import python packages
try:
    from mpi4py import MPI
    import numpy as np
    import numpy.matlib
    import KratosMultiphysics.RomApplication.tsqr
    from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition
    numpy_and_mpi4py_available = True
except:
    numpy_and_mpi4py_available  = False

def synthetic_matrix(degree, rows = 36,repetitions=4):
    TestMatrix = np.zeros((rows,degree))
    x = np.linspace(0,1,rows)
    for i in range(degree):
        TestMatrix[:,i] = np.power(x,i)
    return np.matlib.repmat(TestMatrix, 1, repetitions)

class TestRandomizedSVDMPI(KratosUnittest.TestCase):

    @KratosUnittest.skipUnless(numpy_and_mpi4py_available, "numpy and mpi4py are required for TSQR MPI")
    def test_radomized_svd_mpi(self):
        """According to random theory, this algorithm is being tested for a given tolerance that must always be met. 
        Random matrices are seeded ("0") within the MPI algorithm because the mpi algorithm requires it in all cases. 
        Furthermore, note that assertLess compares it to a given tolerance, which must be satisfied in all cases for random theory.
        Referring to Brunton's book "There is an extensive literature on random matrix theory, 
        where the above stereotypes (this case) are almost certainly true, meaning that they are
        true with high probability"""
        comm = MPI.COMM_WORLD      # Communications macro
        svd_truncation_tolerance = 1e-6
        rank = comm.Get_rank()
        size = comm.Get_size()
        send_data=None
        rows = None
        if rank==0:
            TestMatrix = synthetic_matrix(4) #create a matrix of known rank using polynomials
            rows = TestMatrix.shape[0]
            cols = TestMatrix.shape[1]
            # Split into sub-arrays along required axis
            arrs = np.array_split(TestMatrix, size, axis=0)
            # Flatten the sub-arrays
            raveled = [np.ravel(arr) for arr in arrs]
            rank_rows_list = [arr.shape[0] for arr in arrs]
            # Join them back up into a 1D array
            send_data = np.concatenate(raveled)
        else:
            cols = None
            rank_rows_list = None
        cols = comm.bcast(cols, root=0)
        rank_rows_list = comm.bcast(rank_rows_list, root=0)

        recvbuf = np.empty((rank_rows_list[rank], int(cols)))
        comm.Scatterv(send_data, recvbuf, root=0)

        Qi,Bi = KratosMultiphysics.RomApplication.tsqr.randomized_orthogonalization(
            recvbuf,comm,svd_truncation_tolerance,1)
        U_final,s,v = KratosMultiphysics.RomApplication.tsqr.svd_parallel(Qi, Bi, comm,
            epsilon=svd_truncation_tolerance)
        
        U_global = None
        U_global = np.array(comm.gather(U_final, root=0))

        U_local = None
        if rank==0:
            rank_rows_list = np.cumsum(rank_rows_list)
            rank_rows_list = np.insert(rank_rows_list, 0, 0)
            U_local = np.zeros((rows,U_global.shape[2]))
            for i in range(0,size):
                idxg, idyg = rank_rows_list[i], rank_rows_list[i+1]
                U_local[idxg:idyg,:] = U_global[i,:,:]

            Randomized_Reconstruction = U_local@np.diag(s)@v.T #reconstruct matrix

            #check that the difference of the reconstruction is below tolerance
            self.assertLess( np.linalg.norm(Randomized_Reconstruction - TestMatrix), svd_truncation_tolerance*np.linalg.norm(TestMatrix))

    # Cleaning
    kratos_utilities.DeleteDirectoryIfExisting("__pycache__")


if __name__=='__main__':
    KratosUnittest.main()