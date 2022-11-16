import numpy as np

from sklearn.cluster import KMeans
from utils.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition

""" 
This module performs the clustering given a snapshot matrix.
@author: Raul Bravo.
"""

def ismember(A, B):
    if isinstance(A, np.int_):
        return [ np.sum(A == B) ]
    else:
        return [ np.sum(a == B) for a in A ]


def add_overlapping(snapshot_matrix, sub_snapshots, number_of_clusters, kmeans, overlap_percentace=.1):
    """
    Farhat's proposed overlapping criteria. Fixed overlap percentage
    """
    if number_of_clusters!=1:

        neighbors={}
        for i in range(number_of_clusters):
            neighbors[i] = []

        #Add overlapping
        for i in range(np.shape(snapshot_matrix)[1]):

            # Identify the two nearest cluster centroids to state i and mark these clusters as neighbors
            this_matrix = (kmeans.cluster_centers_).T - snapshot_matrix[:,i].reshape(np.shape(snapshot_matrix)[0], 1)
            distance = np.zeros((number_of_clusters))

            for j in range(number_of_clusters):
                distance[j] = np.linalg.norm(this_matrix[:,j])

            second_nearest_idx = np.argsort(distance)[1]

            if not(sum(ismember(neighbors[kmeans.labels_[i]],second_nearest_idx))):
                neighbors[kmeans.labels_[i]].append(second_nearest_idx)

        snapshots_to_add = []

        for j in range(number_of_clusters):
            
            N_snaps = np.shape(sub_snapshots[j])[1]
            N_neighbors = len(neighbors[j])
            N_add = int(np.ceil(N_snaps*overlap_percentace / N_neighbors)) # Number of snapshots to add to subset i

            for i in range(N_neighbors):
                print('adding neighbors from ', neighbors[j][i], 'to cluster ', j)
                this_matrix = sub_snapshots[neighbors[j][i]] - ((kmeans.cluster_centers_[j]).T).reshape(np.shape(SnapshotMatrix)[0], 1)
                distance = np.zeros(np.shape(sub_snapshots[neighbors[j][i]][1]))

                for k in range(len(distance)):
                    distance[k] = np.linalg.norm(this_matrix[:,k])

                indices_to_add = np.argsort(distance)

                if i==0:
                    snapshots_to_add.append(sub_snapshots[neighbors[j][i]][:,indices_to_add[:N_add]])
                else:
                    snapshots_to_add[j] = np.c_[ snapshots_to_add[j], sub_snapshots[neighbors[j][i]][:,indices_to_add[:N_add]]]

        for j in range(number_of_clusters):
            sub_snapshots[j] = np.c_[sub_snapshots[j], snapshots_to_add[j]]

    return sub_snapshots


def calcualte_snapshots(snapshot_matrix, number_of_clusters=1, add_overlapping=False, save_to_file=False):

    # Clustering
    kmeans = KMeans(n_clusters=number_of_clusters).fit(snapshot_matrix.T)
    # Add random_state=0 to the KMeans initialization to get consistent results

    # Split snapshots into sub-sets
    sub_snapshots={}
    for i in range(number_of_clusters):
        sub_snapshots[i] = snapshot_matrix[:,kmeans.labels_==i]

    # Add overlapping 10 %
    if add_overlapping:
        sub_snapshots = add_overlapping(snapshot_matrix, sub_snapshots, number_of_clusters, kmeans)

    # Calcualte the svd of each cluster and obtain its modes
    bases={}
    truncation_tolerance = 1e-4

    for i in range(number_of_clusters):
        bases[i],s,_,_ = RandomizedSingularValueDecomposition().Calculate(sub_snapshots[i],truncation_tolerance)

        if save_to_file:
            np.save(f'bases_{i+1}.npy', bases[i])

            # # In case we want to control if its saved on binary/Ascii, etc...
            # with open("bases_"+str(i)+".npy", "wb") as bases_file:
            #     np.save(bases_file, Bases[i])

    return bases, sub_snapshots

def calcualte_snapshots_with_columns(snapshot_matrix, number_of_clusters=1, number_of_columns_in_the_basis=None, truncation_tolerance=1e-4):

    # Clustering
    kmeans = KMeans(n_clusters=number_of_clusters).fit(snapshot_matrix.T)
    print(f'{kmeans=}')
    # Add random_state=0 to the KMeans initialization to get consistent results

    # Split snapshots into sub-sets
    sub_snapshots={}
    for i in range(number_of_clusters):
        sub_snapshots[i] = snapshot_matrix[:,kmeans.labels_==i]

    # print(f'{sub_snapshots=}')
    print(f'{len(sub_snapshots)=}')
    # exit()

    # Calcualte the svd of each cluster and obtain its modes
    bases={}
    q = {}
    
    if number_of_columns_in_the_basis is not None:
        truncation_tolerance = 0

    for i in range(number_of_clusters):
        bases[i],_,_,_ = RandomizedSingularValueDecomposition().Calculate(sub_snapshots[i],truncation_tolerance)
        if number_of_columns_in_the_basis is not None:
            # print(f'{bases[i].shape[1]=}, {number_of_columns_in_the_basis=}')
            if bases[i].shape[1]<number_of_columns_in_the_basis:
                raise Exception(f'There are no {number_of_columns_in_the_basis} linearly independent columns spaning the range of cluster {i}')
            bases[i] = bases[i][:,:number_of_columns_in_the_basis]
        #np.save(f'bases_{i+1}.npy', bases[i])

    return bases, kmeans

# if __name__=='__main__':
#     #clustering parameters
#     SnapshotMatrix = np.load('snapshot_matrix.npy')
#     number_of_clusters = 3

#     #cluster the snaphots
#     bases, sub_snapshots = calcualte_snapshots(SnapshotMatrix, number_of_clusters)

#     #printing
#     print(f'\n\n\nThere are {number_of_clusters} clusters:\n')
#     for i in range(len(bases)):
#         print(f'cluster_{i}\'s shape is: {np.shape(sub_snapshots[i])} ')
#         print(f'Basis_{i}\'s shape is: {np.shape(bases[i])} \n')

#         with open("cluster_"+str(i)+".npy", "wb") as cluster_file:
#             np.save("cluster_"+str(i)+".npy", bases[i])