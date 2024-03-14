import KratosMultiphysics
from KratosMultiphysics.deprecation_management import DeprecationManager
import KratosMultiphysics.MPMApplication as KratosMPM
from KratosMultiphysics.MPMApplication.apply_mpm_particle_neumann_condition_process import ApplyMPMParticleNeumannConditionProcess

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyMPMCouplingInterfaceNeumannConditionProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyMPMCouplingInterfaceNeumannConditionProcess(ApplyMPMParticleNeumannConditionProcess):
    def __init__(self, Model, settings ):

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name"           : "PLEASE_SPECIFY_MODEL_PART_NAME",
                "material_points_per_condition"   : 0,
                "variable_name"             : "POINT_LOAD",
                "constrained"               : "fixed",
                "option"                    : ""
            }  """ )

        context_string = type(self).__name__
        old_name = 'particles_per_condition'
        new_name = 'material_points_per_condition'
        if DeprecationManager.HasDeprecatedVariable(context_string, settings, old_name, new_name):
            DeprecationManager.ReplaceDeprecatedVariableName(settings, old_name, new_name)

        settings.ValidateAndAssignDefaults(default_parameters)
        self.model = Model
        self.model_part_name = settings["model_part_name"].GetString()

        # Initiate base class - Neumann condition
        super(ApplyMPMCouplingInterfaceNeumannConditionProcess, self).__init__(Model, settings)

        # Set INTERFACE flag active
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INTERFACE, True, self.model_part.Conditions)


    def ExecuteBeforeSolutionLoop(self):
        # Get updated model_part
        if (self.model_part_name.startswith('Background_Grid.')):
            self.model_part_name = self.model_part_name.replace('Background_Grid.','')
        mpm_material_model_part_name = "MPM_Material." + self.model_part_name
        self.model_part = self.model[mpm_material_model_part_name]

        ### Translate conditions with INTERFACE flag into a new model part "MPM_Coupling_Neumann_Interface" responsible for coupling 
        # Create coupling model part
        if not self.model.HasModelPart("MPM_Coupling_Neumann_Interface"):
            self.model.CreateModelPart("MPM_Coupling_Neumann_Interface")
        self.coupling_model_part = self.model.GetModelPart("MPM_Coupling_Neumann_Interface").CreateSubModelPart(self.model_part_name)

        # Prepare coupling model part
        self._prepare_coupling_model_part(self.coupling_model_part)

        # Create nodes and fill coupling model part
        for mpc in self.model_part.Conditions:
            if (mpc.Is(KratosMultiphysics.INTERFACE)):
                node_id         = mpc.Id
                node_coordinate = mpc.CalculateOnIntegrationPoints(KratosMPM.MPC_COORD, self.model_part.ProcessInfo)[0]
                self.coupling_model_part.CreateNewNode(node_id, node_coordinate[0], node_coordinate[1], node_coordinate[2])


    def ExecuteInitializeSolutionStep(self):
        ### Clone delta time
        self.coupling_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME] = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

        # ### Send Point Load from coupling node to condition

        for coupling_node in self.coupling_model_part.Nodes:
            coupling_id  = coupling_node.Id
            point_load = coupling_node.GetSolutionStepValue(KratosMultiphysics.CONTACT_FORCE)

            self.model_part.GetCondition(coupling_id).SetValuesOnIntegrationPoints(KratosMPM.POINT_LOAD, [point_load], self.model_part.ProcessInfo)



    def ExecuteFinalizeSolutionStep(self):
        ### Get kinematic variables from mp to coupling_mp
        for mpc in self.model_part.Conditions:
            if (mpc.Is(KratosMultiphysics.INTERFACE)):
                coupling_id   = mpc.Id

                displacement = mpc.CalculateOnIntegrationPoints(KratosMPM.MPC_DISPLACEMENT, self.model_part.ProcessInfo)[0]
                velocity = mpc.CalculateOnIntegrationPoints(KratosMPM.MPC_VELOCITY, self.model_part.ProcessInfo)[0]
                
                self.coupling_model_part.GetNode(coupling_id).SetSolutionStepValue(KratosMultiphysics.VELOCITY,0,velocity)
                self.coupling_model_part.GetNode(coupling_id).SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,0,displacement)

    # Local functions
    def _prepare_coupling_model_part(self, coupling_model_part):
        # Define domain size
        domain_size = self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        coupling_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = domain_size

        # Add variables
        coupling_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        coupling_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        coupling_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONTACT_FORCE)
        
