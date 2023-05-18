# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA 

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation


def Create(*args):
    return ComputeNodalCouplingForce(*args)

class ComputeNodalCouplingForce(CoSimulationCouplingOperation):

    def __init__(self, settings, solver_wrappers, process_info, data_communicator):
        super().__init__(settings, process_info, data_communicator)
        self.model = solver_wrappers[self.settings["solver"].GetString()].model
        self.model_part_name = self.settings["model_part_name"].GetString()
        self.model_part = self.model[self.model_part_name]
        self.penalty = self.settings["penalty"].GetDouble()  
        self.dt= self.settings["timestep"].GetDouble() 


    def InitializeCouplingIteration(self):

        for node in self.model_part.Nodes: 
                total_mass = node.GetSolutionStepValue(KM.NODAL_MAUX)
                if(total_mass!=0): 
                    velocity_dem = node.GetSolutionStepValue(KM.LINEAR_MOMENTUM) / total_mass
                    node.SetSolutionStepValue(KM.LAGRANGE_DISPLACEMENT, self.dt * (velocity_dem - node.GetSolutionStepValue(KM.VELOCITY)))
                    body_force_density= node.GetSolutionStepValue(KM.LAGRANGE_DISPLACEMENT)* self.penalty* -1
                    node.SetSolutionStepValue(KM.VOLUME_ACCELERATION, node.GetSolutionStepValue(KM.VOLUME_ACCELERATION) + body_force_density)

        for elem in self.model_part.Elements:
             if(elem.GetNodes()[0].GetSolutionStepValue(KM.NODAL_MAUX))!=0:  
                for i in range(elem.GetGeometry().IntegrationPointsNumber()): #gauss quadrature
                    w = 1 #w = (elem.GetGeometry().IntegrationPoints[i].Weight())?
                    J = (elem.GetGeometry().DeterminantOfJacobian(i))
                    shape_functions = elem.GetGeometry().ShapeFunctionsValues()
                    for n in range(len(elem.GetNodes())):
                        for m in range(len(elem.GetNodes())):
                            vol = self.penalty * w * J * shape_functions[(i,n)] * shape_functions[(i,m)]
                            elem.GetNodes()[n].SetSolutionStepValue(SMA.POINT_LOAD, elem.GetNodes()[n].GetSolutionStepValue(SMA.POINT_LOAD) + vol * elem.GetNodes()[n].GetSolutionStepValue(KM.LAGRANGE_DISPLACEMENT))  


    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "solver"                : "UNSPECIFIED",
            "model_part_name"       : "",
            "penalty"               : 1e5,
            "timestep"              : 1e-5
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults
