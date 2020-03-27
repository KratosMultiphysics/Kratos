import numpy as np
import os
import linecache
import sys
import time
import copy
import subprocess

import KratosMultiphysics as KM
from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
from KratosMultiphysics.CoSimulationApplication.co_simulation_interface import CoSimulationInterface
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return SolverWrapperOpenFOAM_41(parameters)


class SolverWrapperOpenFOAM_41(CoSimulationComponent):
    def __init__(self, parameters):
        super().__init__()
        
        #Settings
        self.settings = parameters["settings"]
        self.working_directory = self.settings["working_directory"].GetString()
#         input_file = self.settings["input_file"].GetString()
#         settings_file_name = os.path.join(working_directory, input_file)
#         with open(settings_file_name, 'r') as settings_file:
#             self.settings.AddParameters(cs_data_structure.Parameters(settings_file.read()))
        self.moduleType="OpenFOAM" #hard-coded, cannot be altered by naive user
        self.moduleVersion="4.1" #hard-coded, cannot be altered by naive user
        self.module=self.moduleType+"/"+self.moduleVersion
        self.application=self.settings["application"].GetString() #What type of OF-solver to be used - solver requires adaptation before running with OpenFOAM - name of adapted solver starts with 'CoCoNuT_'
        self.dimensions=self.settings["dimensions"].GetInt()
        if (self.dimensions != 2) and (self.dimensions != 3):
            sys.exit("OpenFOAM-case should be 2D or 3D.")
        self.dt = self.settings["dt"].GetDouble()  # Time step size
        self.start_time=self.settings["start_time"].GetDouble() # Start time - also the name of folder containing data at this time
        self.end_time=self.settings["end_time"].GetDouble() # End time
        self.cores=self.settings["cores"].GetInt() # Number of cores to be used in the OpenFOAM-calculation
        self.decomposeMethod=self.settings["decomposeMethod"].GetString() #Decomposition-method, can be "simple", "scotch" 
        self.newtonmax = self.settings["newtonmax"].GetInt()  # Maximal number of Newton iterations
        self.newtontol = self.settings["newtontol"].GetDouble()  # Tolerance of Newton iterations
        self.write_interval = self.settings["write_interval"].GetInt() # Number of time steps between consecutive saves performed by OpenFOAM 
        self.write_precision = self.settings["write_precision"].GetInt() # writePrecision-parameter in OpenFOAM
        self.time_precision = self.settings["time_precision"].GetInt() # timePrecision-parameter in OpenFOAM
        self.boundary_names = [_.GetString() for _ in self.settings['boundary_names'].list()] # boundary_names is the set of boundaries where the moving interface is located (will be used to go through OF-files)
        
        #Check that the boundary_names and the interface_input and interface_output are defined consistently in the JSON-file
        #For every boundary_name element, there should be one interface_input (boundary_name+"_input") element and one interface_output (boundary_name+"_output") element.
        #Make the distinction between both: interface_input/output are names of the pyKratos ModelPart - boundary_names is the name of the boundary as defined in OpenFOAM!
        if len(self.boundary_names) != len(self.settings['interface_input'].keys()):
            sys.exit("Interface_input and boundary_names should have the same length; one interface_input element corresponds with one boundary_name element!")
        if len(self.boundary_names) != len(self.settings['interface_output'].keys()):
            sys.exit("Interface_output and boundary_names should have the same length; one interface_output element corresponds with one boundary_name element!")
        for boundary in self.boundary_names:
            index=0
            key_input = self.settings['interface_input'].keys()[index]
            if not("_input" in key_input):
                sys.exit('Please make that all interface_input elements correspond to a boundary_names element followed by "_input".')
            else:
                keyBoundary=key_input.replace("_input","")
                if keyBoundary != boundary:
                    sys.exit('For OpenFOAM, please make sure that every boundary_names element is linked to an interface_input element with the following name: <boundary_names element>"_input". The corresponding elements in boundary_names and interface_input should have the same index in their respective arrays!')
            key_output = self.settings['interface_output'].keys()[index]
            if not("_output" in key_output):
                sys.exit('Please make that all interface_output elements correspond to a boundary_names element followed by "_output".')
            else:
                keyBoundary=key_output.replace("_output","")
                if keyBoundary != boundary:
                    sys.exit('For OpenFOAM, please make sure that every boundary_names element is linked to an interface_output element with the following name: <boundary_names element>"_output". The corresponding elements in boundary_names and interface_output should have the same index in their respective arrays!')
                if key_output.replace("_output","_input") != key_input:
                    sys.exit("Please make sure that the interface_input and interface_output elements occur in such a way that the corresponding elements have the same index.")
            index+=1
 
        #Check that the correct modules have been loaded
        self.check_software()  
        
        #Remove possible CoCoNuT-message from previous interrupt
        self.remove_all_messages()
        
        # Creating OpenFOAM-files - raw dictionary files are predefined in the solver_wrapper folder (and should not be moved)
        # DecomposeParDict: replace raw settings by actual settings defined by user in json-file
        if self.cores > 1: #Only if calculating in parallel
            decomposeParDict_raw_name=os.path.join(os.path.realpath(os.path.dirname(__file__)),"decomposeParDict_raw")
            decomposeParDict_name=os.path.join(self.working_directory,"system/decomposeParDict")
            with open(decomposeParDict_raw_name,'r') as rawFile:
                with open(decomposeParDict_name,'w') as newFile:
                    for line in rawFile:
                        line=line.replace('|CORES|',str(self.cores))
                        line=line.replace('|DECOMPOSEMETHOD|',str(self.decomposeMethod))
                        newFile.write(line)
            rawFile.close()
            newFile.close()
            self.write_footer(decomposeParDict_name) 
            # OpenFOAM-fields are decomposed automatically if you work in parallel
            os.system("cd " + self.working_directory + "; decomposePar -force -time "+ str(self.start_time) + " &> log.decomposePar;")  
        # ControlDict: replace raw settings by actual settings defined by user in json-file AND add function objects to write pressure and wall shear stress
        controlDict_raw_name=os.path.join(os.path.realpath(os.path.dirname(__file__)),"controlDict_raw")
        controlDict_name=os.path.join(self.working_directory,"system/controlDict")
        with open(controlDict_raw_name,'r') as rawFile:
            with open(controlDict_name,'w') as newFile:
                for line in rawFile:
                    line=line.replace('|APPLICATION|',str(self.application))
                    line=line.replace('|START_TIME|',str(self.start_time))
                    line=line.replace('|END_TIME|',str(self.end_time))
                    line=line.replace('|DT|',str(self.dt))
                    line=line.replace('|WRITE_INTERVAL|',str(self.write_interval))
                    line=line.replace('|WRITE_PRECISION|',str(self.write_precision))
                    line=line.replace('|TIME_PRECISION|',str(self.time_precision))
                    newFile.write(line)
        rawFile.close()
        newFile.close()
        nKey=0
        for key in self.boundary_names:
            if nKey == 0:
                self.write_controlDict_function(controlDict_name,"wallShearStress","libfieldFunctionObjects.so",key,True,False)
            else: 
                self.write_controlDict_function(controlDict_name,"wallShearStress","libfieldFunctionObjects.so",key,False,False)
            if nKey == (len(self.boundary_names)-1):
                self.write_controlDict_function(controlDict_name,"pressure","libfieldFunctionObjects.so",key,False,True)
            else:
                self.write_controlDict_function(controlDict_name,"pressure","libfieldFunctionObjects.so",key,False,False)
            nKey += 1
        self.write_footer(controlDict_name)
    
        # Creating Model
        self.model = cs_data_structure.Model()
        print("The model for OpenFOAM will be created. Please make sure all patch names given under the 'interface' setting are also found in the mesh used in OpenFOAM (see 'constant/polyMesh') \n")
        
        # Creating ModelParts and adding variables to these ModelParts - should happen before node addition
        for key,value in (self.settings['interface_input'].items()+self.settings['interface_output'].items()): # OpenFOAM is a completely collocated, so all variables are stored in cell centres
            self.model.CreateModelPart(key)
            mp=self.model[key]
            for var_name in value.list():
                var=vars(KM)[var_name.GetString()]
                mp.AddNodalSolutionStepVariable(var)
            
        # Adding nodes to ModelParts - should happen after variable definition
        os.system("cd "+ self.working_directory + "; writeCellCentres -time " + str(self.start_time) + " &> log.writeCellCentres;")
        self.nNodes_proc=np.zeros([len(self.settings['interface_input'].keys()),self.cores],dtype=np.int32)
        self.nNodes_tot=np.zeros(len(self.settings['interface_input'].keys()),dtype=np.int32)
        nKey=0
        for boundary in self.boundary_names:
            if self.cores == 1: # If you work in serial
                #To avoid confusion, remove all present processor-files (from a possible previous parallel calculation)
                os.system("cd " + self.working_directory + "; rm -r processor*;") 
                
                # First look up the interface wall in OpenFOAM's constant-folder - the coordinates are found using the OF-utility 'writeCellCentres'
                for i in np.arange(3):
                    if i == 0:
                        source_file=self.working_directory + "/" + str(self.start_time) + "/ccx"
                    elif i == 1:
                        source_file=self.working_directory + "/" + str(self.start_time) + "/ccy"
                    elif i == 2:
                        source_file=self.working_directory + "/" + str(self.start_time) + "/ccz"
                    try:
                        lineNameNr=self.find_string_in_file(boundary, source_file)
                        lineStartNr=lineNameNr+7 # In case of non-uniform list, this is where the list of values in the sourceFile starts
                        nNodesIndex=lineNameNr+5 # On this line, the number of cell centers on the inlet is stated
                        os.system("awk NR==" + str(nNodesIndex) + " " + source_file + " > nNodes")
                        nNodes_file=open("nNodes",'r')
                        nNodes=int(nNodes_file.readline())
                        nNodes_file.close()
                        os.system("rm nNodes")
                        tempCoordFile=np.ones([nNodes,1])*float("inf")
                        for j in np.arange(nNodes):
                            tempCoordFile[j,0]=float(linecache.getline(source_file,lineStartNr+j))
                    except ValueError: # If ValueError is triggered, it means that the source-file has a uniform coordinate in the axis you are currently looking
                        if not 'nNodes' in locals(): #If the first coordinate-file has a uniform value, the variable 'nNodes' does not exist, so you should check whether this variable exists
                            check_file=self.working_directory + "/" + str(self.start_time) + "/ccy" # if 'nNodes' does not exist, read the second sourceFile to know the number of rows
                            lineNameNr=self.find_string_in_file(boundary, check_file)
                            nNodesIndex=lineNameNr+5 # On this line, the number of cell centers on the inlet is stated
                            os.system("awk NR==" + str(nNodesIndex) + " " + check_file + " > nNodes")
                            nNodes_file=open("nNodes",'r')
                            nNodes=int(nNodes_file.readline())
                            nNodes_file.close()
                            os.system("rm nNodes")
                        indexUV=lineNameNr+4
                        os.system("awk NR==" + str(indexUV) + " " + source_file + " > unifValue")
                        unifValue_file=open("unifValue",'r')
                        unifValue=float(unifValue_file.readline().split()[-1][0:-1]) #First '-1' makes sure the value is read, but this still contains a semi-colon, so this should be removed with second index '[0:-1]'.
                        unifValue_file.close()
                        os.system("rm unifValue")
                        tempCoordFile=np.ones([nNodes,1])*float("inf")
                        for j in np.arange(nNodes):
                            tempCoordFile[j,0]=unifValue
                    if i == 0:
                        coordList_tot=np.ones([nNodes,4])*float("inf") # ID - X - Y - Z - Area
                        coordList_tot[:,0]=np.arange(nNodes) #CellID = numbers from 0 till (nNodes - 1)
                        self.nNodes_proc[nKey,0]=nNodes
                        self.nNodes_tot[nKey]=self.nNodes_proc[nKey,0]
                    coordList_tot[:,(i+1)]=tempCoordFile[:,0]
            else: # If you work in parallel
                foundProcWithInterfaceValues=False
                for p in np.arange(self.cores): # Check in each processor-folder
                    for i in np.arange(3):
                        if i == 0:
                            source_file=self.working_directory + "/processor" + str(p) + "/" + str(self.start_time) + "/ccx"
                        elif i == 1:
                            source_file=self.working_directory + "/processor" + str(p) + "/" + str(self.start_time) + "/ccy"
                        elif i == 2:
                            source_file=self.working_directory + "/processor" + str(p) + "/" + str(self.start_time) + "/ccz"
                        lineNameNr=self.find_string_in_file(boundary, source_file)
                        procContainsBoundaryNr=lineNameNr+4
                        os.system("awk NR==" + str(procContainsBoundaryNr) + " " + source_file + " > PCB" )
                        PCB_file=open("PCB",'r')
                        procContainsBoundary=str(PCB_file.readline()) #  whether there are cells adjacent to the interface in the processor-folder 
                        PCB_file.close() 
                        os.system("rm PCB")
                        if "nonuniform 0();" not in procContainsBoundary: # If this processor folder contains the interface you are looking for
                            try:
                                lineNameNr=self.find_string_in_file(boundary, source_file)
                                lineStartNr=lineNameNr+7 # In case of non-uniform list, this is where the list of values in the sourceFile starts
                                nNodesIndex=lineNameNr+5 # On this line, the number of cell centers on the inlet is stated
                                os.system("awk NR==" + str(nNodesIndex) + " " + source_file + " > nNodes")
                                nNodes_file=open("nNodes",'r')
                                nNodes=int(nNodes_file.readline())
                                nNodes_file.close()
                                os.system("rm nNodes")
                                tempCoordFile=np.ones([nNodes,1])*float("inf")
                                for j in np.arange(nNodes):
                                    tempCoordFile[j,0]=float(linecache.getline(source_file,lineStartNr+j))
                            except ValueError: # If ValueError is triggered, it means that the source-file has a uniform coordinate in the axis you are currently looking OR that the nonuniform list contains less than 11 elements - because apparently OF prints the entire array on a single line in that case... 
                                if "nonuniform" not in procContainsBoundary: # If the coordinate is a uniform value 
                                    if not 'nNodes' in locals(): #If the first coordinate-file has a uniform value, the variable 'nNodes' does not exist, so you should check whether this variable exists
                                        check_file=self.working_directory + "/processor" + str(p) + "/" + str(self.start_time) + "/ccy"
                                        lineNameNr=self.find_string_in_file(boundary, check_file)
                                        typeCcyNr=lineNameNr+5
                                        os.system("awk NR==" + str(typeCcyNr) + " " + check_file + " > typeCcy" )
                                        typeCcy_file=open("typeCcy",'r')
                                        typeCcy=str(typeCcy_file.readline())
                                        typeCcy_file.close()
                                        if "(" in typeCcy: # 10 nodes or less in this list
                                            listCcy=((typeCcy.split("(")[-1]).split(")")[0]).split(" ")
                                            nNodes=int(len(listCcy))     
                                        else : # More than 10 nodes in this list or ccy is also uniform
                                            if "nonuniform" not in typeCcy: #ccy is also uniform, check ccz
                                                check_file=self.working_directory + "/processor" + str(p) + "/" + str(self.start_time) + "/ccz"
                                                lineNameNr=lineNameNr=self.find_string_in_file(boundary, check_file)
                                                typeCczNr=lineNameNr+4
                                                os.system("awk NR==" + str(typeCczNr) + " " + check_file + " > typeCcz" )
                                                typeCcz_file=open("typeCcz",'r')
                                                typeCcz=str(typeCcz_file.readline())
                                                typeCcz_file.close()
                                                if "(" in typeCcz: # 10 nodes or less in this list 
                                                    listCcz=((typeCcz.split("(")[-1]).split(")")[0]).split(" ")
                                                    nNodes=int(len(listCcz))
                                                else: #More than 10 nodes in this list or ccz is also uniform
                                                    if "nonuniform" not in typeCcz: #ccz is also uniform - you only have 1 point in this processor (for real?)
                                                        nNodes=1
                                                    else : #there are more than 10 nodes in this list
                                                        nNodesIndex=lineNameNr+5 # On this line, the number of cell centers on the inlet is stated
                                                        os.system("awk NR==" + str(nNodesIndex) + " " + check_file + " > nNodes")
                                                        nNodes_file=open("nNodes",'r')
                                                        nNodes=int(nNodes_file.readline())
                                                        nNodes_file.close()              
                                            else: # ccy is not uniform
                                                nNodesIndex=lineNameNr+5 # On this line, the number of cell centers on the inlet is stated
                                                os.system("awk NR==" + str(nNodesIndex) + " " + check_file + " > nNodes")
                                                nNodes_file=open("nNodes",'r')
                                                nNodes=int(nNodes_file.readline())
                                                nNodes_file.close() 
                                        os.system("rm nNodes typeCcy")
                                    indexUV=lineNameNr+4
                                    os.system("awk NR==" + str(indexUV) + " " + source_file + " > unifValue")
                                    unifValue_file=open("unifValue",'r')
                                    unifValue=float(unifValue_file.readline().split()[-1][0:-1]) #First '-1' makes sure the value is read, but this still contains a semi-colon, so this should be removed with second index '[0:-1]'.
                                    unifValue_file.close()
                                    os.system("rm unifValue")
                                    tempCoordFile=np.ones([nNodes,1])*float("inf")
                                    for j in np.arange(nNodes):
                                        tempCoordFile[j,0]=unifValue
                                else : # Coordinates are a nonuniform list containing less than 11 elements and are all written on one line
                                    listToRead=((procContainsBoundary.split("(")[-1]).split(")")[0]).split(" ")
                                    nNodes=int(len(listToRead))
                                    tempCoordFile=np.ones([nNodes,1])*float("inf")
                                    for j in np.arange(nNodes):
                                            tempCoordFile[j,0]=float(listToRead[j])                                
                            if i == 0:
                                coordList=np.ones([nNodes,4])*float("inf") # ID - X - Y - Z - Area
                                coordList[:,0]=self.nNodes_tot[nKey]*np.ones(nNodes)+np.arange(nNodes) #CellID = numbers from 0 till (nNodes - 1)
                                self.nNodes_proc[nKey,p]=nNodes
                            coordList[:,(i+1)]=tempCoordFile[:,0]
                    if "nonuniform 0();" not in procContainsBoundary: # If this processor folder contains the interface you are looking for
                        if not(foundProcWithInterfaceValues) :
                            coordList_tot= coordList
                            foundProcWithInterfaceValues=True
                        else :
                            coordList_tot=np.concatenate((coordList_tot,coordList), axis=0)
                        del nNodes
                    self.nNodes_tot[nKey]=int(np.sum(self.nNodes_proc[nKey,:]))
            
            # Subsequently create each node separately in the ModelPart- both input and output!
            mp_input=self.model[boundary+"_input"]
            mp_output=self.model[boundary+"_output"]
            for i in np.arange(int(self.nNodes_tot[nKey])):
                mp_input.CreateNewNode(coordList_tot[i,0],coordList_tot[i,1],coordList_tot[i,2],coordList_tot[i,3])
                mp_output.CreateNewNode(coordList_tot[i,0],coordList_tot[i,1],coordList_tot[i,2],coordList_tot[i,3])
            nKey+=1    
        
        # Print node distribution among processor folders
        self.write_node_distribution()
    
        # Create CoSimulationInterfaces
        self.interface_input = CoSimulationInterface(self.model, self.settings["interface_input"])
        self.interface_output = CoSimulationInterface(self.model, self.settings["interface_output"])

        # Create Variables
        self.pressure=vars(KM)['PRESSURE']
        self.shear=KM.KratosGlobals.GetVariable('WALLSHEARSTRESS') #used longer definition here because wallshearstress might be overwritten in future updates - this method has more clear error message if variable is not found
        self.displacement=vars(KM)['DISPLACEMENT']
        
             
    def Initialize(self):
        super().Initialize()
                
        # Define timestep and physical time
        self.timestep=0
        self.physical_time=self.start_time
        
        # If no pointDisplacement file is defined yet, initialize a pointDisplacement file in the start time folder
        # Normally, after restart or from the second iteration onwards, a pointDisplacement-file already exists. In that case, that pointDisplacement-file will be used (and is NOT overwritten)
        pointDisp_raw_name=os.path.join(os.path.realpath(os.path.dirname(__file__)),"pointDisplacement_raw")
        if self.cores == 1:
            pointDisp_name=os.path.join(self.working_directory,str(self.physical_time),"pointDisplacement")
            if not(os.path.isfile(pointDisp_name)):
                self.write_pointDisplacement_file(pointDisp_raw_name,pointDisp_name,0)
        else:
            for p in np.arange(self.cores):
                pointDisp_name=os.path.join(self.working_directory,"processor"+str(p),str(self.physical_time),"pointDisplacement")
                if not(os.path.isfile(pointDisp_name)):
                    self.write_pointDisplacement_file(pointDisp_raw_name,pointDisp_name,p)
        
        # Don't forget to start OpenFOAM-loop!
        if self.cores == 1:
            cmd = self.application + "&> log." + self.application
        else:
            cmd = "mpirun -np " + str(self.cores) + " " + self.application + " -parallel &> log." + self.application
        self.openfoam_process = subprocess.Popen(cmd, cwd=self.working_directory, shell=True) 
        
                                        

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        
        # Prepare new time step folder and reset the number of iterations
        self.timestep += 1
        self.iteration = 0
        self.physical_time += self.dt
        newPath=os.path.join(self.working_directory, str(self.physical_time))
        if os.path.isdir(newPath):
            print("\n\n\n Warning! In 5s, CoCoNuT will overwrite existing time step folder: "+str(newPath) + ". \n\n\n")
            time.sleep(5)
            os.system("rm -rf " + newPath)
        os.system("mkdir "+ newPath)
        print('\t Time step '+str(self.timestep))
        
        self.send_message('next') # Let OpenFOAM go to next time step
        self.wait_message('next_ready') # Let OpenFOAM wait for input data
    

    def SolveSolutionStep(self, interface_input): # NOT CHANGED YET! PURELY COPIED FROM FLUENT WRAPPER!!!!!!
        self.iteration += 1
        print(f'\t\tIteration {self.iteration}')

        # store incoming displacements
        self.interface_input.SetPythonList(interface_input.GetPythonList())

        # update X,Y,Z in interface
        for key in [_[0] for _ in self.interface_input.model_parts_variables]:
            for node in self.model[key].Nodes:
                disp = node.GetSolutionStepValue(self.displacement)
                node.X = node.X0 + disp[0]
                node.Y = node.Y0 + disp[1]
                node.Z = node.Z0 + disp[2]

        # write interface data to OpenFOAM-file
        self.write_node_input()
            
        # let Fluent run, wait for data
        self.send_message('continue')
        self.wait_message('continue_ready')

        # read data from OpenFOAM
        self.read_node_output()
        
        # return interface_output object
        return self.interface_output.deepcopy()


    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        # Let OpenFOAM check whether it needs to save this timestep (in OF-solver: runTime.write())
        
        if not(self.timestep % self.write_interval):
            self.send_message('save')
            self.wait_message('save_ready')
#         else:
#             # This is done to remove the OpenFOAM-subfolder containing the wallShearStress, pressure and displacement for that particular time step
#             os.system("rm -r " + os.path.join(self.working_directory, self.physical_time))
#             pass     
            
            
    def Finalize(self):
        super().Finalize()
        
        self.send_message('stop')
        self.wait_message('stop_ready')
        
        self.openfoam_process.kill()
                
        print("OpenFOAM was stopped with the Finalize() method defined in CoCoNuT.")


    def GetInterfaceInput(self):
        return self.interface_input.deepcopy()


    def SetInterfaceInput(self):
        Exception("This solver interface provides no mapping.")


    def GetInterfaceOutput(self):
        return self.interface_output.deepcopy()


    def SetInterfaceOutput(self):
        Exception("This solver interface provides no mapping.")

    
    def write_header(self,fileLoc,className,objectName):
        f=open(fileLoc,'w')
        f.write(r'/*--------------------------------*- C++ -*----------------------------------*\\'+"\n")
        f.write(r'| =========                 |                                                 |'+"\n")
        f.write(r'| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |'+"\n")
        f.write(r'|  \\    /   O peration     | Version:  4.x                                   |'+"\n")
        f.write(r'|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |'+"\n")
        f.write(r'|    \\/     M anipulation  |                                                 |'+"\n")
        f.write(r'\*---------------------------------------------------------------------------*/'+"\n")
        f.write(r'FoamFile'+"\n")
        f.write(r'{'+"\n")
        f.write('\t version \t\t 4.1;'+"\n")
        f.write('\t format \t\t ascii;'+"\n")
        f.write('\t class \t\t ' + className + ';'+"\n")
        f.write('\t object \t\t ' + objectName + ';'+"\n")
        f.write('}'+"\n")
        f.write(r'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //'+"\n")
        f.write("\n")  
        f.close()
        
    def read_node_output(self):
        nKey=0
        maxIt=100
        for boundary in self.boundary_names:
            wss_tmp=np.zeros([self.nNodes_tot[nKey],3])
            pres_tmp=np.zeros([self.nNodes_tot[nKey],3])
            mp = self.model[boundary+"_output"]
            it=0
            if self.cores == 1:
                wss_file= os.path.join(self.working_directory, str(self.physical_time), "wallShearStress")
                pres_file= os.path.join(self.working_directory, str(self.physical_time), "static(p)")
                if self.nNodes_proc[nKey,0] == 1: #Single node on interface (Really??)
                    #Read wall shear stress
                    while True:
                        try:
                            lineNameNr_wss=self.find_string_in_file(boundary, wss_file)
                            indexUV=lineNameNr_wss+4
                            os.system("awk NR==" + str(indexUV) + " " + wss_file + " > unifValue")
                            unifValue_file=open("unifValue",'r')
                            unifValue=unifValue_file.readline().split("\t")[-1][0:-1] #First '-1' makes sure the value is read, but this still contains a semi-colon, so this should be removed with second index '[0:-1]'.
                            unifValue_file.close()
                            os.system("rm unifValue ")
                            for i in np.arange(self.nNodes[nKey,0]):
                                    wss_tmp[i,0]=float(unifValue[0][1:])
                                    wss_tmp[i,1]=float(unifValue[1])
                                    wss_tmp[i,2]=float(unifValue[2][0:-2])
                            break
                        except ValueError:
                            time.sleep(1)
                            it+=1
                        if it > maxIt:
                            os.system("pkill " + self.application)
                            sys.exit("CoCoNUT stopped in read_node_output - please check the output-files of OpenFOAM")      
                    #Read pressure
                    while True:
                        try:
                            lineNameNr_pres=self.find_string_in_file(boundary, pres_file)
                            indexUV=lineNameNr_pres+4
                            os.system("awk NR==" + str(indexUV) + " " + pres_file + " > unifValue")
                            unifValue_file=open("unifValue",'r')
                            unifValue=unifValue_file.readline().split()[-1][0:-1] #First '-1' makes sure the value is read, but this still contains a semi-colon, so this should be removed with second index '[0:-1]'.
                            unifValue_file.close()
                            os.system("rm unifValue ")
                            for i in np.arange(self.nNodes[nKey,0]):
                                pres_tmp[i]=float(unifValue)
                            break
                        except ValueError:
                            time.sleep(1)
                            it+=1
                        if it > maxIt:
                            os.system("pkill " + self.application)
                            sys.exit("CoCoNUT stopped in read_node_output - please check the output-files of OpenFOAM")    
                elif self.nNodes_proc[nKey,0] < 11: # 10 or less elements on the interface
                    # read wall shear stress
                    while True:
                        try:
                            lineNameNr_wss=self.find_string_in_file(boundary, wss_file)
                            listNr=lineNameNr_wss+4
                            os.system("awk NR==" + str(listNr) + " " + wss_file + " > listNr" )
                            listNr_file=open("listNr",'r')
                            listToRead=str(listNr_file.readline()) #  whether there are cells adjacent to the interface in the processor-folder  
                            listNr_file.close()
                            os.system("rm listNr")
                            listToRead=((listToRead.split("(")[-1]).split(")")[0]).split(" ")
                            for i in np.arange(self.nNodes_proc[nKey,0]):
                                wss_tmp[i,0]=float(listToRead[i])
                                wss_tmp[i,1]=float(listToRead[i])
                                wss_tmp[i,2]=float(listToRead[i])
                            break
                        except ValueError:
                            time.sleep(1)
                            it+=1
                        if it > maxIt:
                            os.system("pkill " + self.application)
                            sys.exit("CoCoNUT stopped in read_node_output - please check the output-files of OpenFOAM")  
                    # read pressure
                    while True:
                        try:
                            lineNameNr_pres=self.find_string_in_file(boundary, pres_file)
                            listNr=lineNameNr_pres+4
                            os.system("awk NR==" + str(listNr) + " " + pres_file + " > listNr" )
                            listNr_file=open("listNr",'r')
                            listToRead=str(listNr_file.readline()) #  whether there are cells adjacent to the interface in the processor-folder  
                            listNr_file.close()
                            os.system("rm listNr")
                            listToRead=((listToRead.split("(")[-1]).split(")")[0]).split(" ")
                            for i in np.arange(self.nNodes_proc[nKey,0]):
                                pres_tmp[i]=float(listToRead[i])
                            break
                        except ValueError:
                            time.sleep(1)
                            it+=1
                        if it > maxIt:
                            os.system("pkill " + self.application)
                            sys.exit("CoCoNUT stopped in read_node_output - please check the output-files of OpenFOAM")  
                else: # More than 10 elements on interface
                    # read wall shear stress
                    for i in np.arange(self.nNodes_proc[nKey,0]):
                        while True:
                            try:
                                lineNameNr_wss=self.find_string_in_file(boundary, wss_file)
                                lineStartNr=lineNameNr_wss+7 # In case of non-uniform list, this is where the list of values in the sourceFile starts
                                str_line=linecache.getline(wss_file,lineStartNr+i).split(" ")
                                wss_tmp[i,0]=float(str_line[0][1:])
                                wss_tmp[i,1]=float(str_line[1])
                                wss_tmp[i,2]=float(str_line[2][0:-2])   
                                break
                            except ValueError:
                                time.sleep(1)
                                it+=1
                            if it > maxIt:
                                os.system("pkill " + self.application)
                                sys.exit("CoCoNUT stopped in read_node_output - please check the output-files of OpenFOAM")  
                    # read pressure
                    for i in np.arange(self.nNodes_proc[nKey,0]):
                        while True:
                            try:
                                lineNameNr_pres=self.find_string_in_file(boundary, pres_file)
                                lineStartNr=lineNameNr_pres+7 # In case of non-uniform list, this is where the list of values in the sourceFile starts
                                str_line=linecache.getline(pres_file,lineStartNr+i)
                                pres_tmp[i]=float(str_line[0:-2]) # remove '\n'    
                                break
                            except ValueError:
                                time.sleep(1)
                                it+=1
                            if it > maxIt:
                                os.system("pkill " + self.application)
                                sys.exit("CoCoNUT stopped in read_node_output - please check the output-files of OpenFOAM")                      
            else:
                nodeCount=0
                for p in np.arange(self.cores):
                    wss_file= os.path.join(self.working_directory, ("processor"+str(p)), str(self.physical_time), "wallShearStress")
                    pres_file= os.path.join(self.working_directory, ("processor"+str(p)), str(self.physical_time), "static(p)")
                    if self.nNodes_proc[nKey,p] == 1: #Single node of the interface on this processor
                        # read wall shear stress
                        while True:
                            try:
                                lineNameNr_wss=self.find_string_in_file(boundary, wss_file)
                                indexUV=lineNameNr_wss+4
                                os.system("awk NR==" + str(indexUV) + " " + wss_file + " > unifValue")
                                unifValue_file=open("unifValue",'r')
                                unifValue=float(unifValue_file.readline().split()[-1][0:-1]) #First '-1' makes sure the value is read, but this still contains a semi-colon, so this should be removed with second index '[0:-1]'.
                                unifValue_file.close()
                                for i in np.arange(self.nNodes[nKey,p]):
                                    wss_tmp[nodeCount+i,0]=float(unifValue[0][1:])
                                    wss_tmp[nodeCount+i,1]=float(unifValue[1])
                                    wss_tmp[nodeCount+i,2]=float(unifValue[2][0:-2])
                                break
                            except ValueError:
                                time.sleep(1)
                                it+=1
                            if it > maxIt:
                                os.system("pkill " + self.application)
                                sys.exit("CoCoNUT stopped in read_node_output - please check the output-files of OpenFOAM")  
                        # read pressure
                        while True:
                            try:
                                lineNameNr_pres=self.find_string_in_file(boundary, pres_file)
                                indexUV=lineNameNr_pres+4
                                os.system("awk NR==" + str(indexUV) + " " + pres_file + " > unifValue")
                                unifValue_file=open("unifValue",'r')
                                unifValue=float(unifValue_file.readline().split()[-1][0:-1]) #First '-1' makes sure the value is read, but this still contains a semi-colon, so this should be removed with second index '[0:-1]'.
                                unifValue_file.close()
                                for i in np.arange(self.nNodes[nKey,p]):
                                    pres_tmp[nodeCount+i]=float(unifValue)
                                break
                            except ValueError:
                                time.sleep(1)
                                it+=1
                            if it > maxIt:
                                os.system("pkill " + self.application)
                                sys.exit("CoCoNUT stopped in read_node_output - please check the output-files of OpenFOAM")  
                    elif self.nNodes_proc[nKey,p] < 11: # 10 or less elements of the interface on this processor, nonuniform list recorded on a single line
                        # read wall shear stress
                        while True:
                            try:
                                lineNameNr_wss=self.find_string_in_file(boundary, wss_file)
                                listNr=lineNameNr_wss+4
                                os.system("awk NR==" + str(listNr) + " " + wss_file + " > listNr" )
                                listNr_file=open("listNr",'r')
                                listToRead=str(listNr_file.readline()) #  whether there are cells adjacent to the interface in the processor-folder  
                                listNr_file.close()
                                os.system("rm listNr")
                                listToRead=((listToRead.split("(")[-1]).split(")")[0]).split(" ")
                                for i in np.arange(self.nNodes_proc[nKey,p]):
                                    wss_tmp[nodeCount+i,0]=float(listToRead[0])
                                    wss_tmp[nodeCount+i,1]=float(listToRead[1])
                                    wss_tmp[nodeCount+i,2]=float(listToRead[2])
                                break
                            except ValueError:
                                time.sleep(1)
                                it+=1
                            if it > maxIt:
                                os.system("pkill " + self.application)
                                sys.exit("CoCoNUT stopped in read_node_output - please check the output-files of OpenFOAM")  
                        # read pressure
                        while True:
                            try:
                                lineNameNr_pres=self.find_string_in_file(boundary, pres_file)
                                listNr=lineNameNr_pres+4
                                os.system("awk NR==" + str(listNr) + " " + pres_file + " > listNr" )
                                listNr_file=open("listNr",'r')
                                listToRead=str(listNr_file.readline()) #  whether there are cells adjacent to the interface in the processor-folder  
                                listNr_file.close()
                                os.system("rm listNr")
                                listToRead=((listToRead.split("(")[-1]).split(")")[0]).split(" ")
                                for i in np.arange(self.nNodes_proc[nKey,p]):
                                    pres_tmp[nodeCount+i]=float(listToRead[i])
                                break
                            except ValueError:
                                time.sleep(1)
                                it+=1
                            if it > maxIt:
                                os.system("pkill " + self.application)
                                sys.exit("CoCoNUT stopped in read_node_output - please check the output-files of OpenFOAM")  
                    else: # More than 10 elements on interface
                        # read wall shear stress
                        for i in np.arange(self.nNodes_proc[nKey,p]):
                            while True:
                                try:
                                    lineNameNr_wss=self.find_string_in_file(boundary, wss_file)
                                    lineStartNr=lineNameNr_wss+7 # In case of non-uniform list, this is where the list of values in the sourceFile starts
                                    str_line=linecache.getline(pres_file,lineStartNr+i).split(" ")
                                    wss_tmp[nodeCount+i,0]=float(str_line[0][1:])
                                    wss_tmp[nodeCount+i,1]=float(str_line[1])
                                    wss_tmp[nodeCount+i,2]=float(str_line[2][0:-2])
                                    break
                                except ValueError:
                                    time.sleep(1)
                                    it+=1
                                if it > maxIt:
                                    os.system("pkill " + self.application)
                                    sys.exit("CoCoNUT stopped in read_node_output - please check the output-files of OpenFOAM")  
                        # read pressure
                        for i in np.arange(self.nNodes_proc[nKey,p]):
                            while True:
                                try:
                                    lineNameNr_pres=self.find_string_in_file(boundary, pres_file)
                                    lineStartNr=lineNameNr_pres+7 # In case of non-uniform list, this is where the list of values in the sourceFile starts
                                    str_line=linecache.getline(pres_file,lineStartNr+i)
                                    pres_tmp[nodeCount+i]=float(str_line[0:-2])
                                    break
                                except ValueError:
                                    time.sleep(1)
                                    it+=1
                                if it > maxIt:
                                    os.system("pkill " + self.application)
                                    sys.exit("CoCoNUT stopped in read_node_output - please check the output-files of OpenFOAM")  
                    nodeCount += self.nNodes_proc[nKey,p]
            
            # store pressure and traction in Nodes
            index=0
            for node in mp.Nodes: # Easier than in Fluent because the sequence of nodes stays the same
                node.SetSolutionStepValue(self.shear, 0, wss_tmp[index])
                node.SetSolutionStepValue(self.pressure, 0, pres_tmp[index])
                index += 1
            
            # go to next interface
            nKey += 1
      
        
    #writeFooter: to write OpenFOAM-footer in file at location 'fileLoc'
    def write_footer(self,fileLoc):
        f=open(fileLoc,'a+')
        f.write("\n")
        f.write(r'// ************************************************************************* //'+"\n")
        f.close()        
    
    
    def write_node_input(self):
        # The values in the pointDisplacement-file at the interface need to be overwritten with the newest values for the interface_input
        nKey=0
        for boundary in self.boundary_names:
            mp = self.model[boundary+"_input"]
            if self.cores == 1: #Working in serial operation
                disp_file = os.path.join(self.working_directory, str(self.physical_time), 'pointDisplacement')
                if self.iteration == 1: #first iteration of new time step: disp_file does not exist yet
                    pointDisp_raw_name=os.path.join(os.path.realpath(os.path.dirname(__file__)),"pointDisplacement_raw")
                    self.write_pointDisplacement_file(pointDisp_raw_name, disp_file,0)
                startNr=self.find_string_in_file(boundary, disp_file)
                os.system("head -n " + str(startNr+1) + " " + disp_file + " > tempDisp")
                if self.nNodes_proc[nKey,0] == 1: #Only 1 node on the interface (Really??)
                    with open('tempDisp', 'a+') as file:
                        file.write("\t { \n")
                        file.write("\t\t type  \t fixedValue; \n")
                        with mp.Nodes[0] as node:
                            dispX=node.X-node.X0
                            dispY=node.Y-node.Y0
                            dispZ=node.Z-node.Z0
                            file.write('\t\t value \t uniform (' + f'{dispX:27.17e} {dispY:27.17e} {dispZ:27.17e}'+'); \n')
                    file.close()
                elif self.nNodes_proc[nKey,0] < 11: # 10 or less elements on interface
                    with open('tempDisp', 'a+') as file:
                        file.write("\t { \n")
                        file.write("\t\t type  \t fixedValue; \n")
                        file.write('\t\t value \t nonuniform List<vector> (')
                        for node in mp.Nodes:
                            dispX=node.X-node.X0
                            dispY=node.Y-node.Y0
                            dispZ=node.Z-node.Z0
                            file.write(' (' + f'{dispX:27.17e} {dispY:27.17e} {dispZ:27.17e}'+') ')
                        file.write(');\n')
                    file.close()
                else :    
                    with open('tempDisp', 'a+') as file:
                        file.write("\t { \n")
                        file.write("\t\t type  \t fixedValue; \n")
                        file.write('\t\t value \t nonuniform List<vector> ( \n')
                        for node in mp.Nodes:
                            dispX=node.X-node.X0
                            dispY=node.Y-node.Y0
                            dispZ=node.Z-node.Z0
                            file.write(' (' + f'{dispX:27.17e} {dispY:27.17e} {dispZ:27.17e}'+') \n')
                        file.write(');\n')
                    file.close()
                os.system("wc -l " + disp_file + " > lengthDisp")
                lengthDisp_file=open("lengthDisp",'r')
                length_disp=int(lengthDisp_file.readline().split(" ")[0])
                lengthDisp_file.close()
                os.system("tail -n " + str(length_disp-(startNr+1)) + " " + disp_file + " > tempDisp2")
                startToEndNr=self.find_string_in_file("}", "tempDisp2")
                os.system("tail -n " + str(length_disp-(startNr+1)-startToEndNr) + " " + disp_file + " > tempDisp3")
                os.system("cat tempDisp tempDisp3 > "+ disp_file)
                os.system("rm tempDisp* lengthDisp")
            else : # Working in parallel mode
                nNodes_previousProc=0
                for p in np.arange(self.cores): 
                    disp_file = os.path.join(self.working_directory, ("processor"+str(p)),str(self.physical_time), 'pointDisplacement')
                    if self.iteration == 1: #first iteration of new time step: disp_file does not exist yet
                        pointDisp_raw_name=os.path.join(os.path.realpath(os.path.dirname(__file__)),"pointDisplacement_raw")
                        self.write_pointDisplacement_file(pointDisp_raw_name, disp_file,p)
                    startNr=self.find_string_in_file(boundary, disp_file)
                    os.system("head -n " + str(startNr+1) + " " + disp_file + " > tempDisp")
                    if self.nNodes_proc[nKey,p] == 1: #Only 1 node on the interface in this processor-folder
                        with open('tempDisp', 'a+') as file:
                            file.write("\t { \n")
                            file.write("\t\t type  \t fixedValue; \n")
                            with mp.Nodes[nNodes_previousProc] as node:
                                dispX=node.X-node.X0
                                dispY=node.Y-node.Y0
                                dispZ=node.Z-node.Z0
                                file.write('\t\t value \t uniform (' + f'{dispX:27.17e} {dispY:27.17e} {dispZ:27.17e}'+'); \n')
                        file.close()
                    elif self.nNodes_proc[nKey,p] < 11: # 10 or less elements on the interface in this processor-folder
                        with open('tempDisp', 'a+') as file:
                            file.write("\t { \n")
                            file.write("\t\t type  \t fixedValue; \n")
                            file.write('\t\t value \t nonuniform List<vector> (')
                            for i in np.arange(self.nNodes_proc[nKey,p]):
                                with mp.Nodes[nNodes_previousProc+i] as node:
                                    dispX=node.X-node.X0
                                    dispY=node.Y-node.Y0
                                    dispZ=node.Z-node.Z0
                                    file.write(' (' + f'{dispX:27.17e} {dispY:27.17e} {dispZ:27.17e}'+') ')
                            file.write(');\n')
                        file.close()
                    else :    
                        with open('tempDisp', 'a+') as file:
                            file.write("\t { \n")
                            file.write("\t\t type  \t fixedValue; \n")
                            file.write('\t\t value \t nonuniform List<vector> ( \n')
                            for i in np.arange(self.nNodes_proc[nKey,p]):
                                with mp.Nodes[nNodes_previousProc+i] as node:
                                    dispX=node.X-node.X0
                                    dispY=node.Y-node.Y0
                                    dispZ=node.Z-node.Z0
                                    file.write(' (' + f'{dispX:27.17e} {dispY:27.17e} {dispZ:27.17e}'+') \n')
                            file.write(');\n')
                        file.close()
                    os.system("wc -l " + disp_file + " > lengthDisp")
                    lengthDisp_file=open("lengthDisp",'r')
                    length_disp=int(lengthDisp_file.readline())
                    lengthDisp_file.close()
                    os.system("tail -n " + str(length_disp-(startNr+1)) + " " + disp_file + " > tempDisp2")
                    startToEndNr=self.find_string_in_file("}", "tempDisp2")
                    os.system("tail -n " + str(length_disp-(startNr+1)-startToEndNr) + " " + disp_file + " > tempDisp3")
                    os.system("cat tempDisp tempDisp3 > "+ disp_file)
                    nNodes_previousProc+=self.nNodes_proc[nKey,p]
                    os.system("rm tempDisp* lengthDisp")
            nKey += 1     
    
    
    def write_node_distribution(self):
        nKey=0
        for boundary in self.boundary_names:
            tmp= 'CoCoNuT_nNodes_processor_' + boundary
            file_name= os.path.join(self.working_directory, tmp)
            with open(file_name, 'w') as file:
                if self.cores > 1:
                    for i in np.arange(len(self.nNodes_proc[nKey,:])):
                        file.write("Processor " + str(i) + ": " + str(int(self.nNodes_proc[nKey,i]))+'\n')
                else :
                        file.write("Serial: " + str(int(self.nNodes_proc[nKey,0]))+'\n')
            file.close()
            nKey += 1
     
            
    def write_controlDict_function(self, filename, funcname, libOFname, patchname, writeStart, writeEnd):
        with open(filename,'a+') as file:
            if writeStart:
                file.write("functions \n")
                file.write("{ \n ")
            file.write(" \n \t " + funcname + "_" + patchname +" \n")
            file.write("\t { \n")
            file.write("\t\t type  \t " + funcname + "; \n")
            file.write('\t\t libs \t ("' + libOFname + '"); \n')
            file.write('\t\t patches ( ' + patchname + ' ); \n')
            file.write('\t\t writeControl \t timeStep; \n')
            file.write('\t\t writeInterval \t 1; \n')
            file.write('\t\t log \t false; \n')
            if funcname == "pressure":
                file.write('\t\t calcTotal \t no; \n')
                file.write('\t\t calcCoeff \t no; \n')
                file.write('\t\t rho \t rho; \n')
                print("\n\n Please check the 'rho' option in the static pressure definition in controlDict! This might vary from OF-solver to OF-solver.\n\n")
            file.write("\t } \n")
            if writeEnd:
                file.write("} \n ")
        file.close()
            
    
    def write_pointDisplacement_file(self,pointDisp_raw_name,pointDisp_name,procNr):
        with open(pointDisp_raw_name,'r') as rawFile: 
                with open(pointDisp_name,'w') as newFile:
                    for line in rawFile:
                        newFile.write(line)
                    nKey=0
                    for boundary in self.boundary_names: 
                        newFile.write(" \n \t "  + boundary +" \n")
                        newFile.write("\t { \n")
                        if self.nNodes_proc[nKey,procNr] > 0:
                            newFile.write("\t\t type  \t fixedValue; \n")
                            newFile.write("\t\t value \t uniform (0 0 0); \n")
                        else:
                            newFile.write("\t\t type  \t calculated; \n")
                            newFile.write("\t\t value \t nonuniform 0(); \n")
                        newFile.write("\t } \n")
                    newFile.write("} \n")
                    newFile.close()
                    self.write_footer(pointDisp_name)
        rawFile.close()
    
    def find_string_in_file(self,string_name,file_name):
        index=-1
        with open(file_name) as f:
            for num, line in enumerate(f):
                if string_name in line:
                    index=num
                    break
        f.close()
        return index
    
    def send_message(self, message):
        file = os.path.join(self.working_directory, message + ".coco")
        open(file, 'w').close()
        return


    def wait_message(self, message):
        waitTimeLimit=0.5*60 # 10 minutes maximum waiting time for a single flow solver iteration
        cumulTime=0
        file = os.path.join(self.working_directory, message + ".coco")
        while not os.path.isfile(file):
            time.sleep(0.01)
            cumulTime += 0.01
            if cumulTime > waitTimeLimit:
                os.system("pkill " + self.application)
                sys.exit("CoCoNuT timed out in the OpenFOAM solver_wrapper, waiting for message: "+ message + ".coco.")
        os.remove(file)
        return


    def check_message(self, message):
        file = os.path.join(self.working_directory, message + ".coco")
        if os.path.isfile(file):
            os.remove(file)
            return True
        return False


    def remove_all_messages(self):
        for file_name in os.listdir(self.working_directory):
            if file_name.endswith('.coco'):
                file = os.path.join(self.working_directory, file_name)
                os.remove(file)
    
    
    def check_software(self):
        if os.system(self.application+' -help &> checkSoftware') != 0:
            sys.exit("You either did not load the module for OpenFOAM/4.1, did not compile the solver you are trying to use or did not source $FOAM_BASH prior to execution of CoCoNuT.")

        # The statement above could work for other version of OpenFOAM is the solver is also compiled for that version. Therefore, the version number is checked explicity (don't forget to remove the end-of-line variable at the end of the String versionNr
        with open('checkSoftware','r') as f:
            lastLine=f.readlines()[-2] # Second last line contains 'Build: XX' with XX the version number 
        f.close()
        os.system('rm checkSoftware')
        versionNr=lastLine.split(' ')[-1]
        if versionNr[:-1] != self.moduleVersion :
                sys.exit("OpenFOAM 4.1 should be loaded! Currently, another version of OpenFOAM is loaded")
        
        #Check that a 'dynamicMeshDict' is present in the constant-folder of the case
        dynamicMeshDict_file=os.path.join(self.working_directory,"constant/dynamicMeshDict")
        if not(os.path.isfile(dynamicMeshDict_file)):
            sys.exit("Please make sure that a 'dynamicMeshDict' file is defined in the constant-folder of your working-directory.")       
    