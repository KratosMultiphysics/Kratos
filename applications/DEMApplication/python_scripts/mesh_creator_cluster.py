import os

import KratosMultiphysics
import KratosMultiphysics.DEMApplication as KratosDEM

# This function takes a ".dempack" file, form the cluster mesher, and a file with the properties and creates a "DEM_Clusters.mdpa"
def WriteClusterMdpa(problem_name):

    ClusterPropertiesInfo = open(problem_name + "_cluster_properties.info" , 'r')
    ClusterMeshInfo       = open(problem_name + "_cluster_mesh.dempack", 'r')
    ClusterMdpa           = open(problem_name + "DEM_Clusters.mdpa", 'w')

    ClusterAvailableTypes  = []

    for Line in ClusterPropertiesInfo:
        ClusterMdpa.write(Line)
        if 'CLUSTER_FILE_NAME' in Line:
            Line = Line.strip('\n')# Remove the line-ending characters
            ClusterAvailableTypes.append(Line.split(' ')[1])

    ClusterType = []
    NodeNumber  = []
    NodeX       = []
    NodeY       = []
    NodeZ       = []
    ChLength    = []
    qX          = []
    qY          = []
    qZ          = []
    qW          = []

    for Line in ClusterMeshInfo:

        if Line.startswith("Begin_Cluster"):
            Line = Line.strip('\n')# Remove the line-ending characters
            for i, cluster_type in enumerate(ClusterAvailableTypes):
                if Line.split(' ')[1] + ".clu" == cluster_type:
                    ClusterType.append(i+1)

        if Line[0].isdigit():
            Line = Line.strip('\n')# Remove the line-ending characters
            NodeNumber.append(Line.split(' ')[0])
            NodeX.append(Line.split(' ')[1])
            NodeY.append(Line.split(' ')[2])
            NodeZ.append(Line.split(' ')[3])
            ChLength.append(Line.split(' ')[4])
            qX.append(Line.split(' ')[5])
            qY.append(Line.split(' ')[6])
            qZ.append(Line.split(' ')[7])
            qW.append(Line.split(' ')[8])

    ClusterMdpa.write('\nBegin Nodes\n')

    for j, node_number in enumerate(NodeNumber):
        ClusterMdpa.write(node_number + ' ' + NodeX[j] + ' ' + NodeY[j] + ' ' + NodeZ[j] + '\n')

    ClusterMdpa.write('End Nodes\n')

    ClusterMdpa.write('\nBegin Elements Cluster3D\n')

    for j, node_number in enumerate(NodeNumber):
        ClusterMdpa.write(node_number + ' ' + str(ClusterType[j]) + ' ' + node_number + '\n')

    ClusterMdpa.write('End Elements\n')

    ClusterMdpa.write('\nBegin NodalData CHARACTERISTIC_LENGTH\n')

    for j, node_number in enumerate(NodeNumber):
        ClusterMdpa.write(node_number + ' 0 ' + ChLength[j] + '\n')

    ClusterMdpa.write('End NodalData\n')

    ClusterMdpa.write('\nBegin NodalData ORIENTATION\n')

    for j, node_number in enumerate(NodeNumber):
        ClusterMdpa.write(node_number + ' 0 [4] ( ' + qX[j] + ' , ' + qY[j] + ' , ' + qZ[j] + ' , ' + qW[j] + ' )\n')

    ClusterMdpa.write('End NodalData\n')

    ClusterPropertiesInfo.close()
    ClusterMeshInfo.close()
    ClusterMdpa.close()

# This function creates a "DEM_Clusters.mdpa" from a simulation
def WriteClusterMdpaFromResults(filename_pre, pre_path, filename_post, cluster_model_part):

    ClusterMdpa_pre  = open(os.path.join(pre_path, filename_pre) + '.mdpa', 'r')
    ClusterMdpa_post = open(os.path.join(pre_path, filename_post) + '.mdpa', 'w')

    for Line in ClusterMdpa_pre:
        # The next 4 lines mean that, while we do not get to 'Begin Nodes', we simply copy the mdpa contents, in this case the Properties part
        ClusterMdpa_post.write(Line)
        if Line.startswith('Begin Nodes'):
            break

    for node in cluster_model_part.Nodes:
        ClusterMdpa_post.write(str(node.Id) + ' ' + str(node.X) + ' ' + str(node.Y) + ' ' + str(node.Z) + '\n')
    ClusterMdpa_post.write('End Nodes\n\n')

    for Line in ClusterMdpa_pre:
        # We copy the element type line from the previous mdpa
        if Line.startswith('Begin Elements'):
            ClusterMdpa_post.write(Line)
            break

    for element in cluster_model_part.Elements:
        ClusterMdpa_post.write(str(element.Id) + ' ' + str(element.Properties.Id) + ' ' + str(element.GetNode(0).Id) + '\n')
    ClusterMdpa_post.write('End Elements\n')

    ClusterMdpa_post.write('\nBegin NodalData CHARACTERISTIC_LENGTH\n')
    for node in cluster_model_part.Nodes:
        ClusterMdpa_post.write(str(node.Id) + ' 0 ' + str(node.GetSolutionStepValue(KratosDEM.CHARACTERISTIC_LENGTH)) + '\n')
    ClusterMdpa_post.write('End NodalData\n')

    ClusterMdpa_post.write('\nBegin NodalData ORIENTATION\n')
    for node in cluster_model_part.Nodes:
        ClusterMdpa_post.write(str(node.Id) + ' 0 [4] ( ' + str(node.GetSolutionStepValue(KratosMultiphysics.ORIENTATION).X) + ' , ' + str(node.GetSolutionStepValue(KratosMultiphysics.ORIENTATION).Y) + ' , ' + str(node.GetSolutionStepValue(KratosMultiphysics.ORIENTATION).Z) + ' , ' + str(node.GetSolutionStepValue(KratosMultiphysics.ORIENTATION).W) + ' )\n')
    ClusterMdpa_post.write('End NodalData\n')

    ClusterMdpa_pre.close()
    ClusterMdpa_post.close()
