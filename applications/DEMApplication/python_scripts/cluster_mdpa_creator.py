import os     #for wraps
import sys    #for wraps
import re     #for regex
import random
import math

import DEM_explicit_solver_var as DEM_parameters

def WriteClusterMdpa():

    ClusterPropertiesInfo = open(DEM_parameters["problem_name"].GetString() + "_cluster_properties.info" , 'r')
    ClusterMeshInfo       = open(DEM_parameters["problem_name"].GetString() + "_cluster_mesh.dempack", 'r')
    ClusterMdpa           = open(DEM_parameters["problem_name"].GetString() + "DEM_Clusters.mdpa", 'w')

    ClusterAvailableTypes  = []

    for Line in ClusterPropertiesInfo:
        ClusterMdpa.write(Line)
        if 'CLUSTER_FILE_NAME' in Line:
            Line = Line.strip('\n')# Remove the line-ending characters
            ClusterAvailableTypes.append(Line.split()[1])

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
        Line = Line.lstrip() # Remove the line-beginning characters

        if Line.startswith("Begin_Cluster"):
            Line = Line.strip('\n')# Remove the line-ending characters
            for i in range(len(ClusterAvailableTypes)):
                if Line.split()[1] + ".clu" == ClusterAvailableTypes[i]:
                    ClusterType.append(i+1)

        if Line[0].isdigit():
            Line = Line.strip('\n')# Remove the line-ending characters
            NodeNumber.append(Line.split()[0])
            NodeX.append(Line.split()[1])
            NodeY.append(Line.split()[2])
            NodeZ.append(Line.split()[3])
            ChLength.append(Line.split()[4])
            qX.append(Line.split()[5])
            qY.append(Line.split()[6])
            qZ.append(Line.split()[7])
            qW.append(Line.split()[8])

    ClusterMdpa.write('\n\nBegin Nodes\n')

    for j in range(len(NodeNumber)):
        ClusterMdpa.write(NodeNumber[j] + ' ' + NodeX[j] + ' ' + NodeY[j] + ' ' + NodeZ[j] + '\n')

    ClusterMdpa.write('End Nodes\n')

    ClusterMdpa.write('\nBegin Elements Cluster3D\n')

    for j in range(len(NodeNumber)):
        ClusterMdpa.write(NodeNumber[j] + ' ' + str(ClusterType[j]) + ' ' + NodeNumber[j] + '\n')

    ClusterMdpa.write('End Elements\n')

    ClusterMdpa.write('\nBegin NodalData CHARACTERISTIC_LENGTH\n')

    for j in range(len(NodeNumber)):
        ClusterMdpa.write(NodeNumber[j] + ' 0 ' + ChLength[j] + '\n')   

    ClusterMdpa.write('End NodalData\n')

    ClusterMdpa.write('\nBegin NodalData ORIENTATION\n')

    for j in range(len(NodeNumber)):
        ClusterMdpa.write(NodeNumber[j] + ' 0 [4] ( ' + qX[j] + ' , ' + qY[j] + ' , ' + qZ[j] + ' , ' + qW[j] + ' )\n')   

    ClusterMdpa.write('End NodalData\n')
 
    ClusterPropertiesInfo.close()
    ClusterMeshInfo.close()
    ClusterMdpa.close()
