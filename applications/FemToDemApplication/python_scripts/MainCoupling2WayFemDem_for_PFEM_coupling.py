

import KratosMultiphysics

import KratosMultiphysics.FemToDemApplication.MainDEM_for_coupling as DEM
import KratosMultiphysics.FemToDemApplication.MainFEM_for_PFEM_coupling as FEM
import KratosMultiphysics.MeshingApplication.mmg_process as MMG
import KratosMultiphysics.FemToDemApplication.MainCoupling2WayFemDem as MainCoupling2WayFemDem

def Wait():
    input("Press Something")

#============================================================================================================================
class MainCoupled2WayFemDem_Solution_for_PFEM_coupling(MainCoupling2WayFemDem.MainCoupled2WayFemDem_Solution):
#============================================================================================================================
    def __init__(self, Model, path = ""):
        self.model = Model
        # Initialize solutions

        if path == "":
            DEMProjectParametersFile = open("ProjectParametersDEM.json", 'r')
        else:
            DEMProjectParametersFile = open(os.path.join(path, "ProjectParametersDEM.json"), 'r')
        DEM_project_parameters = KratosMultiphysics.Parameters(DEMProjectParametersFile.read())
        self.DEM_Solution = DEM.DEM_for_coupling_Solution(Model, DEM_project_parameters)

        self.DEM_Solution.Initialize()
        self.FEM_Solution = FEM.FEM_for_PFEM_coupling_Solution(Model, path, self.DEM_Solution._GetSolver())

        # Initialize Remeshing files
        self.DoRemeshing = self.FEM_Solution.ProjectParameters["AMR_data"]["activate_AMR"].GetBool()
        if self.DoRemeshing:
            self.mmg_parameter_file = open("MMGParameters.json",'r')
            self.mmg_parameters = KratosMultiphysics.Parameters(self.mmg_parameter_file.read())
            self.RemeshingProcessMMG = MMG.MmgProcess(Model, self.mmg_parameters)
        self.domain_size = self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        self.InitializePlotsFiles()
        self.echo_level = 0
        self.is_slave = False

