# import libraries #######################################################################################################
from __future__ import print_function
import numpy as np
from itertools import islice
import sys
import os
from ctypes import *
import csv

# Definition of toolbox (Eof user modifiable part) #######################################################################

class Toolbox:

    # ********************************************************************************************
    def __init__( self,
                  SU2MeshFile,
                  SURFACE_FLOW_FILENAME,
                  SURFACE_ADJ_FILENAME_drag,
                  SURFACE_ADJ_FILENAME_lift,
                  flowHistoryFile,
                  optPatch,
                  sensitivity_scaling ):

        # Store settings
        self.SU2MeshFile = SU2MeshFile
        self.SURFACE_FLOW_FILENAME = SURFACE_FLOW_FILENAME + ".csv"
        self.SURFACE_ADJ_FILENAME_drag = SURFACE_ADJ_FILENAME_drag + ".csv"
        self.SURFACE_ADJ_FILENAME_lift = SURFACE_ADJ_FILENAME_lift + ".csv"
        self.flowHistoryFile = flowHistoryFile
        self.optPatch = optPatch
        self.sensitivity_scaling = sensitivity_scaling

        # Global variables #
        self.SU2MeshData  = {} 
            
        # Read SU2 mesh
        self.readSU2Mesh()

    # ********************************************************************************************
    def readSU2Mesh(self):
    
        print("\n> Start reading SU2 mesh file...")

        ''' imports mesh and builds python dictionary structure
            SU2MeshData            mesh data dictionary
            SU2MeshData['NDIME']   number of dimensions
            SU2MeshData['NELEM']   number of elements
            SU2MeshData['ELEM']    element array [ type, nodes, index ]
            SU2MeshData['NPOIN']   number of points
            SU2MeshData['POIN']    list of lists with all point coordinates [[x,y,z],[x,y,z]...]
            SU2MeshData['NMARK']   number of markers
            SU2MeshData['MARKS']   marker data dictionary
            SU2MeshData['MARKS']['tag_name']           marker data for 'tag_name'
            SU2MeshData['MARKS']['tag_name']['NELEM']  number of elements
            SU2MeshData['MARKS']['tag_name']['ELEM']   element array [type,nodes]
            SU2MeshData['MARKS']['tag_name']['POINTS'] stores points on marker
        '''
        
        # initialize variables
        marks  = {} 

        # open meshfile
        f_tb_read = open(self.SU2MeshFile,'r')

        # readline helper functin
        def mesh_readlines(n_lines=1):
            fileslice = islice(f_tb_read,n_lines)
            return list(fileslice)

        # scan file until end of file
        keepon = True
        while keepon:

            # read line
            line = mesh_readlines()

            # stop if line is empty
            if not line: 
                keepon = False
                break
            
            # fix white space
            line = line[0]
            line = line.replace('\t',' ')
            line = line.replace('\n',' ')

            # skip comments
            if line[0] == "%":
                pass

            # number of dimensions
            elif "NDIME=" in line:
                # save to SU2_MESH data
                self.SU2MeshData['NDIME'] = int( line.split("=")[1].strip() )
            #:if NDIME

            # elements
            elif "NELEM=" in line:
                
                # number of elements
                nelem = int( line.split("=")[1].strip() )
                # save to SU2_MESH data
                self.SU2MeshData['NELEM'] = nelem
                
                # only read nelem lines
                fileslice = islice(f_tb_read,nelem)
                
                # the data pattern
                pattern = tuple( [int] + [int]*9 )
                
                # scan next lines for element data
                elem = [ 
                    [ t(s) for t,s in zip(pattern,line.split()) ] 
                    for line in fileslice 
                ]
                
                # save to SU2_MESH data
                self.SU2MeshData['ELEM'] = elem
            #: if NELEM

            # points
            elif "NPOIN=" in line:
                
                # number of points
                npoin = int( line.split("=")[1].strip().split(' ')[0] )
                # save to SU2_MESH data
                self.SU2MeshData['NPOIN'] = npoin
                
                # only read npoin lines
                fileslice = islice(f_tb_read,npoin)
                
                # the data pattern
                pattern = tuple( [float]*3 ) # + [int] )
                
                # scan next lines for element data
                poin = [ 
                    [ t(s) for t,s in zip(pattern,line.split()) ] 
                    for line in fileslice 
                ]            

                # save to SU2_MESH data
                self.SU2MeshData['POIN'] = poin
            #:if NPOIN

            # number of markers
            elif "NMARK=" in line:
                nmark = int( line.split("=")[1].strip() )
                # save to SU2_MESH data
                self.SU2MeshData['NMARK'] = nmark
            #:if NMARK

            # a marker
            elif "MARKER_TAG=" in line:
                # marker tag
                thistag = line.split("=")[1].strip()

                # start SU2_MARK dictionary
                thismark = {} 
                # save to SU2_MARK data
                thismark['TAG'] = thistag

                # read number of marker elements
                line = mesh_readlines()[0]
                if not "MARKER_ELEMS=" in line:
                    raise Exception("Marker Specification Error")
                
                # convert string to int int
                thisnelem = int( line.split("=")[1].strip() )
                
                # save to SU2_MARK data
                thismark['NELEM'] = thisnelem
                
                # only read thisnelem lines
                fileslice = islice(f_tb_read,thisnelem)
                
                # the data pattern
                pattern = tuple( [int] + [int]*9 )
                
                # scan next lines for element data
                markelem = [ 
                    [ t(s) for t,s in zip(pattern,line.split()) ] 
                    for line in fileslice 
                ]
                
                # save to SU2_MARK data
                thismark['ELEM'] = markelem
                
                # add to marker list
                marks[thismark['TAG']] = thismark
            #:if MARKER_TAG

        #:while not end of file
        
        # save to SU2_MESH data
        self.SU2MeshData['MARKS'] = marks
        
        # Extract mesh info at optimization patch ###

        # Assign to each marker their corresponding SU2 node IDs
        for mark_tag in sorted(self.SU2MeshData['MARKS']): # !!!! Since we loop over dictionary keys, we need to define an order to always have the same order --> whe chose a simple sorting)
            this_mark = self.SU2MeshData['MARKS'][mark_tag]
            this_mark['POINTS'] = []
            for elem in this_mark['ELEM']:
                # loop over all points
                itr = 0
                for node in elem:
                    if(itr>0):
                        SU2_node_ID = int("%i " % node)
                        if SU2_node_ID not in this_mark['POINTS']:
                            this_mark['POINTS'].append(SU2_node_ID)
                    itr = itr + 1 
        
        # Close read file
        f_tb_read.close()   
        
        print("> Finished reading SU2 mesh file!\n")
      
    # ********************************************************************************************
    def convertSU2MeshToKratosMesh(self):
    
        print("> Start writing mdpa file...")
        
        # Variables in this function
        opt_mesh = "SU2ToKratos_for_optimization.mdpa"
        complete_mesh = "SU2ToKratos_for_mesh_motion.mdpa"
        KRATOSPropertyID = 1

        # First we write the file for performing the mesh-motion    
        f_tb_created = open(complete_mesh, 'w')

        # Write header
        f_tb_created.write("Begin Properties ")
        f_tb_created.write(str(KRATOSPropertyID))
        f_tb_created.write("\n")
        f_tb_created.write(" DENSITY 7870\n")
        f_tb_created.write(" YOUNG_MODULUS 200.0e9\n")
        f_tb_created.write(" THICKNESS 1.0\n")
        f_tb_created.write(" POISSON_RATIO 0.29\n")
        f_tb_created.write("End Properties\n\n")
    
        # Create elem iterator
        elemID = 0
        
        # write nodes
        f_tb_created.write("Begin Nodes\n")
        itr = 0
        for point in self.SU2MeshData['POIN']:
            f_tb_created.write(str(itr+1))
            f_tb_created.write("\t")  
            f_tb_created.write(str("%.12f"%(point[0])))
            f_tb_created.write("\t")  
            f_tb_created.write(str("%.12f"%(point[1])))
            f_tb_created.write("\t")
            if(self.SU2MeshData['NDIME'] == 2):
                f_tb_created.write(str(0.0))
            if(self.SU2MeshData['NDIME'] == 3):
                f_tb_created.write(str("%.12f"%(point[2])))
            f_tb_created.write("\n")
            itr = itr + 1
        f_tb_created.write("End Nodes\n\n")    

        # write triangular elements
        f_tb_created.write("Begin Elements ShellThinElement3D3N\n")
        for entry in self.SU2MeshData['ELEM']:
            if(entry[0] == 5):
              f_tb_created.write(str(elemID+1))
              f_tb_created.write("\t")       
              f_tb_created.write(str(KRATOSPropertyID))
              f_tb_created.write("\t")   
              f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID][1] + 1))
              f_tb_created.write("\t") 
              f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID][2] + 1))
              f_tb_created.write("\t")  
              f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID][3] + 1))          
              f_tb_created.write("\n")
              elemID = elemID + 1
        f_tb_created.write("End Elements\n\n") 

        # write tet elements (all pyramid element elements in SU2 will be modelled as two tet elements)
        f_tb_created.write("Begin Elements SmallDisplacementElement3D4N\n")
        num_pyramids = 0
        for entry in self.SU2MeshData['ELEM']:
            # Check if current element is a tet element (in SU2 identified by the prefix 10)
            if(entry[0] == 10):
                f_tb_created.write(str(elemID+1))
                f_tb_created.write("\t")       
                f_tb_created.write(str(KRATOSPropertyID))
                f_tb_created.write("\t")   
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID][1] + 1))
                f_tb_created.write("\t") 
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID][2] + 1))
                f_tb_created.write("\t")  
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID][3] + 1))
                f_tb_created.write("\t")  
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID][4] + 1))      
                f_tb_created.write("\n")
                elemID = elemID + 1
            # If its a pyramid (in SU2 identified by 14) then capture it as two tet elements by introducing another edge from a diagonal of the square base
            elif(entry[0] == 14):
                num_pyramids = num_pyramids + 1
                elemID_tet1 = elemID
                elemID_tet2 = self.SU2MeshData['NELEM'] + num_pyramids
                
                f_tb_created.write(str(elemID_tet1 +1))
                f_tb_created.write("\t")       
                f_tb_created.write(str(KRATOSPropertyID))
                f_tb_created.write("\t")   
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID_prism][2] + 1))
                f_tb_created.write("\t") 
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID_prism][3] + 1))
                f_tb_created.write("\t")  
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID_prism][1] + 1))
                f_tb_created.write("\t")  
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID_prism][5] + 1))      
                f_tb_created.write("\n")  

                f_tb_created.write(str(elemID_tet2 +1))
                f_tb_created.write("\t")       
                f_tb_created.write(str(KRATOSPropertyID))
                f_tb_created.write("\t")   
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID_prism][4] + 1))
                f_tb_created.write("\t") 
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID_prism][1] + 1))
                f_tb_created.write("\t")  
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID_prism][3] + 1))
                f_tb_created.write("\t")  
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID_prism][5] + 1))      
                f_tb_created.write("\n")         
                
        f_tb_created.write("End Elements\n\n")
        
        # write prism elements
        f_tb_created.write("Begin Elements SmallDisplacementElement3D6N\n")
        for entry in self.SU2MeshData['ELEM']:
            # Check if current element is a prism element (in SU2 identified by the prefix 13)
            if(entry[0] == 13):
                f_tb_created.write(str(elemID+1))
                f_tb_created.write("\t")       
                f_tb_created.write(str(KRATOSPropertyID))
                f_tb_created.write("\t")   
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID][4] + 1))
                f_tb_created.write("\t") 
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID][5] + 1))
                f_tb_created.write("\t")  
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID][6] + 1))
                f_tb_created.write("\t")  
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID][1] + 1)) 
                f_tb_created.write("\t")  
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID][2] + 1))
                f_tb_created.write("\t")  
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID][3] + 1))      
                f_tb_created.write("\n")
                elemID = elemID + 1
        f_tb_created.write("End Elements\n\n") 
    
        # write hex elements
        f_tb_created.write("Begin Elements SmallDisplacementElement3D8N\n")
        for entry in self.SU2MeshData['ELEM']:
            # Check if current element is a prism element (in SU2 identified by the prefix 13)
            if(entry[0] == 12):
                f_tb_created.write(str(elemID+1))
                f_tb_created.write("\t")       
                f_tb_created.write(str(KRATOSPropertyID))
                f_tb_created.write("\t")   
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID][1] + 1))
                f_tb_created.write("\t") 
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID][2] + 1))
                f_tb_created.write("\t")  
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID][3] + 1))
                f_tb_created.write("\t")  
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID][4] + 1)) 
                f_tb_created.write("\t")  
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID][5] + 1))
                f_tb_created.write("\t")  
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID][6] + 1))   
                f_tb_created.write("\t")        
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID][7] + 1))
                f_tb_created.write("\t")  
                f_tb_created.write(str(self.SU2MeshData['ELEM'][elemID][8] + 1))            
                f_tb_created.write("\n")
                elemID = elemID + 1
        f_tb_created.write("End Elements\n\n")         
        
        # write all marked patches in SU2 as separate meshes in kratos
        mark_itr = 0
        allPatchNodeIDs = []
        for current_mark in sorted(self.SU2MeshData['MARKS']): # !!! Since we loop over dictionary keys, we need to define an order to always have the same order --> whe chose a simple sorting)
            patchNodeIDs = []
            f_tb_created.write("Begin Mesh ")
            f_tb_created.write(str(mark_itr + 1))
            f_tb_created.write("\t// ")
            f_tb_created.write(current_mark)
            f_tb_created.write("\n")
            f_tb_created.write("Begin MeshNodes\n")
            for entry in self.SU2MeshData['MARKS'][current_mark]['ELEM']:
                itr = 0
                for num in entry:
                    if(itr!= 0): # The first entry represents the surface element. E.g. 5 for a triangle
                        currentNodeID = int(entry[itr])
                        if currentNodeID not in patchNodeIDs:
                            patchNodeIDs.append(currentNodeID)
                            allPatchNodeIDs.append(currentNodeID)
                            f_tb_created.write("\t") 
                            f_tb_created.write(str(currentNodeID + 1))   
                            f_tb_created.write("\n")
                    itr = itr + 1  
            f_tb_created.write("End MeshNodes\n")   
            f_tb_created.write("End Mesh\n\n")
            mark_itr = mark_itr + 1

        f_tb_created.close()
    
        # After writing the mesh-motion file we write the optimization patch as kratos mesh-file    
        f_tb_created = open(opt_mesh, 'w')

        # Write header
        f_tb_created.write("Begin Properties ")
        f_tb_created.write(str(KRATOSPropertyID))
        f_tb_created.write("\n")
        f_tb_created.write(" DENSITY 7870\n")
        f_tb_created.write(" YOUNG_MODULUS 200.0e9\n")
        f_tb_created.write(" THICKNESS 1.0\n")
        f_tb_created.write(" POISSON_RATIO 0.29\n")
        f_tb_created.write("End Properties\n\n")
    
        # write nodes
        f_tb_created.write("Begin Nodes\n")
        for point_Id in sorted(self.SU2MeshData['MARKS'][optPatch]['POINTS']): # the sorted command generates an ascending order of the points ids
            f_tb_created.write(str(point_Id+1))
            f_tb_created.write("\t")  
            f_tb_created.write(str("%.12f"%(self.SU2MeshData['POIN'][point_Id][0])))
            f_tb_created.write("\t")  
            f_tb_created.write(str("%.12f"%(self.SU2MeshData['POIN'][point_Id][1]))) 
            f_tb_created.write("\t")
            if(self.SU2MeshData['NDIME'] == 2):
                f_tb_created.write(str(0.0))
            if(self.SU2MeshData['NDIME'] == 3):
                f_tb_created.write(str("%.12f"%(self.SU2MeshData['POIN'][point_Id][2])))
            f_tb_created.write("\n")
        f_tb_created.write("End Nodes\n\n")  

        # write line elements on optimization patch
        f_tb_created.write("Begin Conditions ShapeOptimizationCondition2D2N\n")
        itr = 0
        lineElementIDs = []
        for entry in self.SU2MeshData['MARKS'][optPatch]['ELEM']:
            if(entry[0] == 3):
                lineElementIDs.append(itr)
                f_tb_created.write(str(lineElementIDs[itr] + 1))
                f_tb_created.write("\t")       
                f_tb_created.write(str(KRATOSPropertyID))
                f_tb_created.write("\t")   
                f_tb_created.write(str(entry[1] + 1))
                f_tb_created.write("\t") 
                f_tb_created.write(str(entry[2] + 1))     
                f_tb_created.write("\n")
            itr = itr + 1
        f_tb_created.write("End Conditions\n\n")   

        # write triangular elements on optimization patch
        f_tb_created.write("Begin Conditions ShapeOptimizationCondition3D3N\n")
        itr = 0
        triangleElementIDs = []
        for entry in self.SU2MeshData['MARKS'][optPatch]['ELEM']:
            if(entry[0] == 5):
                triangleElementIDs.append(itr)
                f_tb_created.write(str(triangleElementIDs[itr] + 1))
                f_tb_created.write("\t")       
                f_tb_created.write(str(KRATOSPropertyID))
                f_tb_created.write("\t")   
                f_tb_created.write(str(entry[1] + 1))
                f_tb_created.write("\t") 
                f_tb_created.write(str(entry[2] + 1))
                f_tb_created.write("\t")  
                f_tb_created.write(str(entry[3] + 1))          
                f_tb_created.write("\n")
            itr = itr + 1
        f_tb_created.write("End Conditions\n\n")       
         
        f_tb_created.close()


        print("> Finished writing mdpa file!\n")  
             
    # ********************************************************************************************
    def updateSU2Mesh(self,new_kratos_mesh):
        
        # Read in results of mesh-motion solver and store new results in the "self.SU2MeshData" object
        for kratos_node_Id in new_kratos_mesh:
            nodeIDSU2 = kratos_node_Id - 1
            self.SU2MeshData['POIN'][nodeIDSU2][0] = new_kratos_mesh[kratos_node_Id][0]
            self.SU2MeshData['POIN'][nodeIDSU2][1] = new_kratos_mesh[kratos_node_Id][1]
            self.SU2MeshData['POIN'][nodeIDSU2][2] = new_kratos_mesh[kratos_node_Id][2] 
                    
        # open file for writing
        outputfile = open(self.SU2MeshFile,'w')

        # numbers
        ndime = self.SU2MeshData['NDIME']

        # write dimension
        outputfile.write("% \n% Problem Dimension \n% \n")
        outputfile.write("NDIME= %i\n" % self.SU2MeshData['NDIME'])

        # write elements
        outputfile.write("% \n% Inner element connectivity \n% \n")
        outputfile.write("NELEM= %i\n" % self.SU2MeshData['NELEM'])
        for elem in self.SU2MeshData['ELEM']:
            for num in elem:
                outputfile.write("%i " % num)
            outputfile.write("\n")

        # write nodes
        outputfile.write("% \n% Node coordinates \n% \n")
        outputfile.write("NPOIN= %i\n" % self.SU2MeshData['NPOIN'])
        poin_itr = 0
        for poin in self.SU2MeshData['POIN']:
            for inum in range(ndime):
                outputfile.write("%#18.10e " % (poin[inum]))
            outputfile.write( "%i\n" % (int(poin_itr)) )
            poin_itr =  poin_itr + 1

        # write markers 
        outputfile.write("% \n% Boundary elements \n% \n")
        outputfile.write( "NMARK= %i\n" % self.SU2MeshData['NMARK'] )
        for mark_tag in sorted(self.SU2MeshData['MARKS'].keys()): # !!!! Since we loop over dictionary keys, we need to define an order to always have the same order --> whe chose a simple sorting)
            this_mark = self.SU2MeshData['MARKS'][mark_tag]
            outputfile.write( "MARKER_TAG= %s\n" % this_mark['TAG'] )
            outputfile.write( "MARKER_ELEMS= %i\n" % this_mark['NELEM'] )
            for elem in this_mark['ELEM']:
                for num in elem:
                    outputfile.write("%i " % num)
                outputfile.write("\n")

        # close file
        outputfile.close() 

    # ********************************************************************************************
    def getLastDragAndLiftValues(self):

        # Read last row of direct flow history file (which contains lift & drag information)
        lastRow = 'wrongLastRow'
        with open(self.flowHistoryFile, 'rt') as csvfile:
            historyReader = csv.reader(csvfile, delimiter=',', quotechar='|')
            for row in historyReader:
                lastRow = row  

        # lift is stored as second csv entry while drag is the third entry in "last row"
        drag = lastRow[2] 
        lift = lastRow[1]

        # Return values
        return float(drag),float(lift)        

    # ******************************************************************************************** 
    def computeVectorizedNormalDragSensitivities(self):
            
        # read scalar objective sensitivities from result file
        nodal_sens = {}
        with open(self.SURFACE_ADJ_FILENAME_drag, 'rt') as csvfile:
            mycsvfile = csv.reader(csvfile, delimiter=',')
            row_itr = 1
            # read values and create map
            for row in mycsvfile:
                # Ignore first row with value description and any blank line
                if(row_itr>1 and len(row)>0):
                    SU2_node_ID = int(row[0])
                    nodal_sens[SU2_node_ID] = float(row[1])
                row_itr = row_itr + 1

        # transfer nodal sensitivities to vector quantity using normal information
        
        # Compute node normals
        nodeNormals = self.computeNodalNormals(True)
        
        # Compute vectorial sensitivities and store it in dictionary using numbering of the design variables (Kratos)
        vectorizedSensitivities = {}
        for SU2_node_id in self.SU2MeshData['MARKS'][self.optPatch]['POINTS']:
            kratos_node_Id = SU2_node_id + 1
            vectorizedSensitivities[kratos_node_Id] = nodeNormals[SU2_node_id] * self.sensitivity_scaling * nodal_sens[SU2_node_id]

        # return sensitivity values according to Kratos numbering
        return vectorizedSensitivities

    # ******************************************************************************************** 
    def computeVectorizedNormalLiftSensitivities(self):
            
        # read scalar objective sensitivities from result file
        nodal_sens = {}
        with open(self.SURFACE_ADJ_FILENAME_lift, 'rt') as csvfile:
            mycsvfile = csv.reader(csvfile, delimiter=',')
            row_itr = 1
            # read values and create map
            for row in mycsvfile:
                # Ignore first row with value description and any blank line
                if(row_itr>1 and len(row)>0):
                    SU2_node_ID = int(row[0])
                    nodal_sens[SU2_node_ID] = float(row[1])
                row_itr = row_itr + 1

        # transfer nodal sensitivities to vector quantity using normal information
        
        # Compute node normals
        nodeNormals = self.computeNodalNormals(True)
        
        # Compute vectorial sensitivities and store it in dictionary using numbering of the design variables (Kratos)
        vectorizedSensitivities = {}
        for SU2_node_id in self.SU2MeshData['MARKS'][self.optPatch]['POINTS']:
            kratos_node_Id = SU2_node_id + 1
            vectorizedSensitivities[kratos_node_Id] = nodeNormals[SU2_node_id] * self.sensitivity_scaling * nodal_sens[SU2_node_id]

        # return sensitivity values according to Kratos numbering
        return vectorizedSensitivities

    # ********************************************************************************************
    def computeNodalNormals(self,normalize):

        # Note that the normal we compute has to be consistent with the optimization mesh in kratos
        
        num_elems_on_opt_patch = self.SU2MeshData['MARKS'][self.optPatch]['NELEM']
        elemNormal = {}
        cummulatedElemNormalsAtNode = {}
        nodalArea = {}
        nodeNormals = {}
        
        # initialize list of element normals
        for elem_itr in range(0,num_elems_on_opt_patch):
            elemNormal[elem_itr] = np.array([0,0,0])
        
        # initialize list of node normals
        for node in self.SU2MeshData['MARKS'][self.optPatch]['POINTS']:
            cummulatedElemNormalsAtNode[node] = []
            nodeNormals[node] = np.array([0,0,0])
            nodalArea[node] = 0.0
            
        # compute element normals including nodal area information
        elem_itr = 0 
        totalArea = 0.0
        for entry in self.SU2MeshData['MARKS'][self.optPatch]['ELEM']:
                         
            # Computes normalized normal of current element and assign the corresponding fraction of the element area to the element nodes 
            if(entry[0] == 5): # In case it is a triangle

                currentNodeID = int(entry[1])   
                lastNodeID = int(entry[len(entry)-1])
                nextNodeID = int(entry[2])

                baseVector1 = np.array([self.SU2MeshData['POIN'][nextNodeID][0] - self.SU2MeshData['POIN'][currentNodeID][0],self.SU2MeshData['POIN'][nextNodeID][1] - self.SU2MeshData['POIN'][currentNodeID][1], self.SU2MeshData['POIN'][nextNodeID][2] - self.SU2MeshData['POIN'][currentNodeID][2]])
                baseVector2 = np.array([self.SU2MeshData['POIN'][lastNodeID][0] - self.SU2MeshData['POIN'][currentNodeID][0],self.SU2MeshData['POIN'][lastNodeID][1] - self.SU2MeshData['POIN'][currentNodeID][1], self.SU2MeshData['POIN'][lastNodeID][2] - self.SU2MeshData['POIN'][currentNodeID][2]]) 
           
                areaNormal = 0.5 * np.cross(baseVector1,baseVector2)
                elemArea = np.linalg.norm(areaNormal)
                elemNormal[elem_itr] = areaNormal/elemArea

                nodalArea[int(entry[1])] = nodalArea[int(entry[1])] + elemArea / 3
                nodalArea[int(entry[2])] = nodalArea[int(entry[2])] + elemArea / 3 
                nodalArea[int(entry[3])] = nodalArea[int(entry[3])] + elemArea / 3

            elif(entry[0] == 3): # In case of lines

                node1 = int(entry[1])
                node2 = int(entry[2])

                dx =   self.SU2MeshData['POIN'][node2][1] - self.SU2MeshData['POIN'][node1][1]
                dy = -(self.SU2MeshData['POIN'][node2][0] - self.SU2MeshData['POIN'][node1][0])
                dz = 0.0

                areaNormal = np.array([dx,dy,dz])
                elemArea = np.linalg.norm(areaNormal)
                elemNormal[elem_itr] = areaNormal/elemArea

                nodalArea[int(entry[1])] = nodalArea[int(entry[1])] + elemArea / 2
                nodalArea[int(entry[2])] = nodalArea[int(entry[2])] + elemArea / 2               

            else:
                sys.exit("Given surface element in SU2 not implemented yet in the computation of the surface normal!")                
            
            elem_itr = elem_itr + 1
        
        # At each surface node collect all normalized normals from neighbor elements 
        elem_itr = 0 
        for elems in self.SU2MeshData['MARKS'][self.optPatch]['ELEM']:
            # loop over all points
            node_itr = 0
            for nodes in elems:
                # The first entry is the element type number
                if(node_itr > 0): 
                    currentNodeID = int(nodes)
                    cummulatedElemNormalsAtNode[currentNodeID].append(elemNormal[elem_itr])
                node_itr = node_itr + 1  
            elem_itr = elem_itr + 1
        
        # At each surface node form averaged sum of all neighbor element normals collected previosly
        for node in self.SU2MeshData['MARKS'][self.optPatch]['POINTS']:
            summedNodeNormal = np.array([0,0,0])
            numNormalsAtNode = len(cummulatedElemNormalsAtNode[node]) 
        
            for i in range(0,numNormalsAtNode):
                summedNodeNormal = summedNodeNormal + cummulatedElemNormalsAtNode[node][i]
            
            # Normalize the final averaged normal
            summedNodeNormal = summedNodeNormal / np.linalg.norm(summedNodeNormal)

            # In case area normal is wanted (normalization == False) then multiply by nodal area again
            if(normalize == True):
              nodeNormals[node] = summedNodeNormal
            elif(normalize == False):
              nodeNormals[node] = nodalArea[node] * summedNodeNormal            
        
        return nodeNormals 

    # ********************************************************************************************
    def compute_consistent_nodal_forces(self,reference_pressure,free_stream_static_pressure):  
          
        # read nodal pressure from file
        nodalPress = {}
        with open(self.SURFACE_FLOW_FILENAME, 'rt') as csvfile:
            mycsvfile = csv.reader(csvfile, delimiter=',')
            row_itr = 1
            # read values and create map
            for row in mycsvfile:
              # Ignore first row with value description and any blank line
              if(row_itr>1 and len(row)>0):
                SU2_ID = int(row[0])
                nodalPress[SU2_ID] = float(row[4])
              row_itr = row_itr + 1

        # transfer nodal pressures to force vectors using normal information
        normalize = False
        nodeNormals = self.computeNodalNormals(normalize)
           
        nodalForces = {}
        for SU2_node_ID in nodalPress:
          # Note that the normals are computed assuming they point away from the structure, the pressure force however just acts in direction of the structure --> multiplication with -1 in the following.
          Kratos_node_ID = SU2_node_ID+1
          nodalForces[Kratos_node_ID] = - nodeNormals[SU2_node_ID] * (nodalPress[SU2_node_ID]*reference_pressure-free_stream_static_pressure)
          # nodalForces[Kratos_node_ID][0] = - nodalPress[SU2_node_ID]*reference_pressure
          # nodalForces[Kratos_node_ID][1] = - nodalPress[SU2_node_ID]*reference_pressure
          # nodalForces[Kratos_node_ID][2] = - nodalPress[SU2_node_ID]*reference_pressure

        return nodalForces

##########################################################################################################################
