import KratosMultiphysics
from KratosMultiphysics.deprecation_management import DeprecationManager
import KratosMultiphysics.MPMApplication as KratosMPM
from KratosMultiphysics.MPMApplication.apply_mpm_particle_dirichlet_condition_process import ApplyMPMParticleDirichletConditionProcess

from math import sqrt

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyMPMCouplingInterfaceDirichletConditionProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyMPMCouplingInterfaceDirichletConditionProcess(ApplyMPMParticleDirichletConditionProcess):
    def __init__(self, Model, settings ):

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name"           : "PLEASE_SPECIFY_MODEL_PART_NAME",
                "material_points_per_condition"   : 0,
                "imposition_type"           : "penalty",
                "penalty_factor"            : 0,
                "constrained"               : "fixed",
                "option"                    : "",
                "is_equal_distributed"      : false
            }  """ )

        context_string = type(self).__name__
        old_name = 'particles_per_condition'
        new_name = 'material_points_per_condition'
        if DeprecationManager.HasDeprecatedVariable(context_string, settings, old_name, new_name):
            DeprecationManager.ReplaceDeprecatedVariableName(settings, old_name, new_name)

        settings.ValidateAndAssignDefaults(default_parameters)
        self.model = Model
        self.model_part_name = settings["model_part_name"].GetString()

        # Initiate base class - Dirichlet condition
        super(ApplyMPMCouplingInterfaceDirichletConditionProcess, self).__init__(Model, settings)

        # Set INTERFACE flag active
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INTERFACE, True, self.model_part.Conditions)


    def ExecuteBeforeSolutionLoop(self):
        # Get updated model_part
        if (self.model_part_name.startswith('Background_Grid.')):
            self.model_part_name = self.model_part_name.replace('Background_Grid.','')
        mpm_material_model_part_name = "MPM_Material." + self.model_part_name
        self.model_part = self.model[mpm_material_model_part_name]

        ### Translate conditions with INTERFACE flag into a new model part "MPM_Coupling_Dirichlet_Interface" responsible for coupling with structure
        # Create coupling model part
        if not self.model.HasModelPart("MPM_Coupling_Dirichlet_Interface"):
            self.model.CreateModelPart("MPM_Coupling_Dirichlet_Interface")
        self.coupling_model_part = self.model.GetModelPart("MPM_Coupling_Dirichlet_Interface").CreateSubModelPart(self.model_part_name)

        # Prepare coupling model part
        self._prepare_coupling_model_part(self.coupling_model_part)

        # Create nodes and fill coupling model part
        for mpc in self.model_part.Conditions:
            if (mpc.Is(KratosMultiphysics.INTERFACE)):
                node_id         = mpc.Id
                node_coordinate = mpc.CalculateOnIntegrationPoints(KratosMPM.MPC_COORD, self.model_part.ProcessInfo)[0]
                coupling_node   = self.coupling_model_part.CreateNewNode(node_id, node_coordinate[0], node_coordinate[1], node_coordinate[2])

                ## Set Displacement and Normal
                normal = mpc.CalculateOnIntegrationPoints(KratosMPM.MPC_NORMAL, self.model_part.ProcessInfo)[0]
                coupling_node.SetSolutionStepValue(KratosMultiphysics.NORMAL,0,normal)


    def ExecuteInitializeSolutionStep(self):
        ### Clone delta time
        self.coupling_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME] = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

        ### Send displacement from coupling_mp to mp
        for coupling_node in self.coupling_model_part.Nodes:
            coupling_id  = coupling_node.Id

            ## IMPOSED DISPLACEMENT
            total_displacement = coupling_node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,0)
            old_displacement = self.model_part.GetCondition(coupling_id).CalculateOnIntegrationPoints(KratosMPM.MPC_DISPLACEMENT, self.model_part.ProcessInfo)[0]
            incremental_displacement = total_displacement - old_displacement
            self.model_part.GetCondition(coupling_id).SetValuesOnIntegrationPoints(KratosMPM.MPC_IMPOSED_DISPLACEMENT, [incremental_displacement], self.model_part.ProcessInfo)

            ## ADD VELOCITY
            current_velocity = coupling_node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,0)
            self.model_part.GetCondition(coupling_id).SetValuesOnIntegrationPoints(KratosMPM.MPC_VELOCITY, [current_velocity], self.model_part.ProcessInfo)

            ## ADD NORMAL
            normal = coupling_node.GetSolutionStepValue(KratosMultiphysics.NORMAL,0)
            # Check and see whether the normal is not zero
            norm_normal = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2])
            if norm_normal > 1.e-10:
                self.model_part.GetCondition(coupling_id).SetValuesOnIntegrationPoints(KratosMPM.MPC_NORMAL, [normal], self.model_part.ProcessInfo)


    def ExecuteFinalizeSolutionStep(self):
        ### Get contact force from mp to coupling_mp
        for mpc in self.model_part.Conditions:
            if (mpc.Is(KratosMultiphysics.INTERFACE)):
                coupling_id   = mpc.Id
                contact_force = mpc.CalculateOnIntegrationPoints(KratosMPM.MPC_CONTACT_FORCE, self.model_part.ProcessInfo)[0]
                self.coupling_model_part.GetNode(coupling_id).SetSolutionStepValue(KratosMultiphysics.CONTACT_FORCE,0,contact_force)


    # Local functions
    def _prepare_coupling_model_part(self, coupling_model_part):
        # Define domain size
        domain_size = self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        coupling_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = domain_size

        # Add variables
        coupling_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        coupling_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        coupling_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONTACT_FORCE)
        coupling_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
