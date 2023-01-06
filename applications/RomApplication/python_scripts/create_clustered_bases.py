from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition
from sklearn.cluster import KMeans
import numpy as np


def ismember(A, B):
    if isinstance(A, np.int_):
        return [ np.sum(A == B) ]
    else:
        return [ np.sum(a == B) for a in A ]

def AddOverlapping(SnapshotMatrix, sub_snapshots, Number_Of_Clusters, kmeans, overlap_percentace=.1, printing=False):
    """
    This function could be implemented more efficently
    """
    if Number_Of_Clusters!=1:

        neighbors={}
        for i in range(Number_Of_Clusters):
            neighbors[i] = []

        #Add overlapping
        for i in range(np.shape(SnapshotMatrix)[1]):
            #identify the two nearest cluster centroids to state i and mark these clusters as neighbors
            this_matrix = (kmeans.cluster_centers_).T - SnapshotMatrix[:,i].reshape(np.shape(SnapshotMatrix)[0], 1)
            distance = np.zeros((Number_Of_Clusters))
            for j in range(Number_Of_Clusters):
                distance[j] = np.linalg.norm(this_matrix[:,j])
            second_nearest_idx = np.argsort(distance)[1]
            if not(sum(ismember(neighbors[kmeans.labels_[i]],second_nearest_idx))):
                neighbors[kmeans.labels_[i]].append(second_nearest_idx)

        snapshots_to_add = []
        for j in range(Number_Of_Clusters):
            N_snaps = np.shape(sub_snapshots[j])[1]
            N_neighbors = len(neighbors[j])

            N_add = int(np.ceil(N_snaps*overlap_percentace / N_neighbors))#number of snapshots to add to subset i

            for i in range(N_neighbors):
                if printing:
                    print('adding neighbors from ', neighbors[j][i], 'to cluster ', j  )
                this_matrix = sub_snapshots[neighbors[j][i]] - ((kmeans.cluster_centers_[j]).T).reshape(np.shape(SnapshotMatrix)[0], 1)
                distance = np.zeros(np.shape(sub_snapshots[neighbors[j][i]][1]))
                for k in range(len(distance)):
                    distance[k] = np.linalg.norm(this_matrix[:,k])
                indices_to_add = np.argsort(distance)
                if i==0:
                    snapshots_to_add.append(sub_snapshots[neighbors[j][i]][:,indices_to_add[:N_add]])
                else:
                    snapshots_to_add[j] =  np.c_[ snapshots_to_add[j] , sub_snapshots[neighbors[j][i]][:,indices_to_add[:N_add]] ]

        for j in range(Number_Of_Clusters):
            sub_snapshots[j] = np.c_[sub_snapshots[j], snapshots_to_add[j]]

    return sub_snapshots


def PreComputeDistances(u0, kmeans, bases):
    number_of_bases = len(kmeans.cluster_centers_)
    d = {}
    w = {}
    for i in range(number_of_bases):
        w[i]={}
        d[i]={}
        for j in range(number_of_bases):
            w[i][j]={}
            d[i][j]=-1
            for m in range(number_of_bases):
                w[i][j][m]=[0]

    regularized_centers = kmeans.cluster_centers_ - u0.T

    for m in range(number_of_bases):
        k = m+1
        for p in range(k,number_of_bases):
            d[m][p] =  (regularized_centers[m].T @ regularized_centers[m]) - (regularized_centers[p].T @ regularized_centers[p])
            d[p][m] = -1*d[m][p]
            for l in range(number_of_bases):
                w[l][m][p] = (2 * bases[l].T @ (regularized_centers[p] - regularized_centers[m] )).tolist()

    return w , d


def CreateClusteredBases( SnapshotMatrix = None, Number_Of_Clusters = 3, truncation_tolerance_svd=1e-4, add_overlapping = True):

    printing = False
    if not isinstance(SnapshotMatrix,np.ndarray):
        printing = True
        print('using a small synthetic matrix')
        #Syntetic simple matrix
        SnapshotMatrix = np.ones((10,20))
        for i in range(19):
            SnapshotMatrix[:,i+1]+=i+1

    u0 = SnapshotMatrix[:,0]

    #cluster  identification using scikit-learn's k-means
    kmeans = KMeans(n_clusters=Number_Of_Clusters).fit(SnapshotMatrix.T)

    #split snapshots into sub-sets
    sub_snapshots={}
    if printing:
        print('\nObtained clusters: ')
    for i in range(Number_Of_Clusters):
        sub_snapshots[i] = SnapshotMatrix[:,kmeans.labels_==i] #we might not neet to save this, can slice it when calling the other function...
        if printing:
            print(sub_snapshots[i])

    #add overlapping 10%
    if add_overlapping:
        sub_snapshots = AddOverlapping(SnapshotMatrix, sub_snapshots,Number_Of_Clusters, kmeans, printing=printing)
        if printing:
            print('\nClusters after overlapping ')
            for i in range(Number_Of_Clusters):
                print(sub_snapshots[i])

    #calcualte the svd of each cluster and obatins its modes
    Bases={}
    for i in range(Number_Of_Clusters):
        Bases[i],_,_,_ = RandomizedSingularValueDecomposition().Calculate(sub_snapshots[i], truncation_tolerance_svd )

    w,z0 = PreComputeDistances(u0, kmeans, Bases)

    return w,z0, Bases


if __name__=='__main__':

    w,z0, Bases = CreateClusteredBases()

    #save (Number_Of_Clusters)-basis nodal bases

    # basis_POD = {}
    # basis_POD["w"] = w
    # basis_POD["z0"] = z0

    # with open('test_send_w_z.json', 'w') as f:
    #     json.dump(basis_POD,f, indent=2)

    """
    import the nodal multiple-bases to each node before online simulation and select
    the corrent cluster based on distance
    """

