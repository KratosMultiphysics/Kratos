from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.IncompressibleFluidApplication 
import KratosMultiphysics.ChimeraApplication as chimeraApp

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver():
    return Partitioned_chimera_solver()

class Partitioned_chimera_solver:
    
    
    ##constructor. the constructor shall only take care of storing the settings 
    ##and the pointer to the main_model part. This is needed since at the point of constructing the 
    ##model part is still not filled and the variables are not yet allocated
    ##
    ##real construction shall be delayed to the function "Initialize" which 
    ##will be called once the model is already filled
    def __init__(self): 
        
        self.patchNameDistanceCalcMap = {} 
        self.patchNameBinPointLocaterMap = {}
        self.patchNameModelPartMap = {}
        self.patchNameSolverMap = {}
        self.overLapDistance = -1
        #### Uitility for extracting the surface meshes on both background and patches
        #### These surface meshes are used to apply the boundary conditions
        self.extractorUtil = chimeraApp.CustomHoleCuttingProcess()
        self.varableExtractor = chimeraApp.CustomExtractVariablesProcess()
        
        print("Construction of patitioned chimera solver finished !")
        
        
    def GetMinimumBufferSize(self):
        return 3;

    def AddBackgroundModelPart(self,background_model_part, background_solver):
        self.background_model_part = background_model_part
        self.background_solver     = background_solver
        self.background_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.background_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        self.background_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)  
        self.background_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        
    def AddChimeraPatchModelPart(self,patchName, patchModelPart, patchSolver):
        patchModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        patchModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        patchModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        patchModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)
        self.patchNameModelPartMap[patchName] = patchModelPart
        self.patchNameSolverMap[patchName]    = patchSolver

    def AddVariables(self):
        pass

    def ImportModelPart(self):
        pass

    def AddDofs(self):
        pass
    
    #### Here I am for now passing only one patchSurfaceName, But there should be a mechanism to pass multiple as
    #### there can be many of them.
    def Initialize(self, patchSurfaceName):
        #### Even here the name "Inlet3D_interface" is specific. There should be a mechanism to generalize this name. 
        #### TODO :: This should be changed. 
        #self.CalculatePatchDistance(patchSurfaceName)

        print ("Partitioned Chimera solver initialization finished.")
        

    #### Here I am for now passing only one patchSurfaceName, But there should be a mechanism to pass multiple as
    #### there can be many of them.
    def CalculatePatchDistance(self, patchSurfaceName):
        #### Even here the name "Inlet3D_interface" is specific. There should be a mechanism to generalize this name. 
        #### TODO :: This should be changed. 
        for patchName, patch_model_part in self.patchNameModelPartMap.items():
            self.background_distance_calculator = KratosMultiphysics.CalculateSignedDistanceTo3DConditionSkinProcess(patch_model_part.GetSubModelPart('Inlet3D_interface'), 
                                                                                         self.background_model_part)
            self.background_distance_calculator.Execute()
            self.Bin = KratosMultiphysics.BinBasedFastPointLocator3D(self.background_model_part) # Locator of fluid element that contains specific point in space
            self.Bin.UpdateSearchDatabase()            
        
        #### Even here the name "Inlet3D_interface" is specific. There should be a mechanism to generalize this name. 
        #### TODO :: This should be changed. 
        for patchName, patch_model_part in self.patchNameModelPartMap.items():
            patch_distance_calc = KratosMultiphysics.CalculateSignedDistanceTo3DConditionSkinProcess(
                                                            patch_model_part.GetSubModelPart('Inlet3D_interface'), 
                                                            patch_model_part)
            patch_distance_calc.Execute()
            self.patchNameDistanceCalcMap[patchName]  = patch_distance_calc
            patchBin = KratosMultiphysics.BinBasedFastPointLocator3D(patch_model_part) # Locator of fluid element that contains specific point in space
            patchBin.UpdateSearchDatabase()             
            self.patchNameBinPointLocaterMap[patchName] = patchBin            

        
    def GetOutputVariables(self):
        pass
    
    def SetOverlapDistance(self, overlapDistance):
        self.overLapDistance = overlapDistance
        
    def ComputeDeltaTime(self):
        pass
        
    def SaveRestart(self):
        pass #one should write the restart file here
        
    def Solve(self):
        self.background_solver.Solve()       
        for patchName, patch_solver in self.patchNameSolverMap.items():
            patch_solver.Solve()   
        
    def SetEchoLevel(self, level):
        self.solver.SetEchoLevel(level)

    def Clear(self):
        self.solver.Clear()
        
    def Check(self):
        self.solver.Check()
        
    def GetPressureOnChimeraParts(self):
        patchBin = KratosMultiphysics.BinBasedFastPointLocator3D(self.background_model_part) # Locator of fluid element that contains specific point in space
        patchBin.UpdateSearchDatabase()
        self.background_distance_calculator.MappingPressureToStructure(patchBin)

    def UpdateOverlapBoundaries(self):
        self.CalculatePatchDistance('')
        
        self.background_cut_volume_modelPart = KratosMultiphysics.ModelPart('test')
        testMPInner = KratosMultiphysics.ModelPart('testInner')
        testMPOuterSurface = KratosMultiphysics.ModelPart('testOuterSurf')
        if(self.overLapDistance < 0):
            print('\n\nWARNING :: Setting default overlap distance of 0.2. \n If this is not desired please set the overlap distance before calling this function !\n\n')
            self.overLapDistance = 0.2
        self.extractorUtil.ExtractMeshBetweenLimits(self.background_model_part, self.background_cut_volume_modelPart, -self.overLapDistance, self.overLapDistance)
        self.extractorUtil.ExtractMeshAtCentroidDistance(self.background_model_part, testMPOuterSurface, -1*self.overLapDistance)
        for patchName, patch_model_part in self.patchNameModelPartMap.items():
            self.extractorUtil.ExtractMeshBetweenLimits(patch_model_part, testMPInner, -1*(self.overLapDistance + 0.750*self.overLapDistance), -1*(self.overLapDistance - 0.750*self.overLapDistance))
                    
        #self.varableExtractor.ExtractVariable(self.background_cut_volume_modelPart, self.patchNameModelPartMap['square'].GetSubModelPart('Inlet3D_interface'), KratosMultiphysics.PRESSURE)        
        return self.background_cut_volume_modelPart, testMPInner
        
