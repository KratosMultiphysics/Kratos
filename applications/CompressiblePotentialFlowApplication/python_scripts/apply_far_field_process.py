import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp

def DotProduct(A,B):
    result = 0
    for i,j in zip(A,B):
        result += i*j
    return result

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyFarFieldProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyFarFieldProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "inlet_potential": 1.0,
                "free_stream_velocity": [1.0,0.0,0]
            }  """ )
        settings.ValidateAndAssignDefaults(default_parameters);

        self.far_field_model_part = Model[settings["model_part_name"].GetString()]
        self.free_stream_velocity = settings["free_stream_velocity"].GetVector()
        self.inlet_potential_0 = settings["inlet_potential"].GetDouble()
        self.far_field_model_part.ProcessInfo.SetValue(CPFApp.VELOCITY_INFINITY,self.free_stream_velocity)

    def ExecuteInitializeSolutionStep(self):
        self.Execute()

    def Execute(self):
        reference_inlet_node = self._FindFarthestUpstreamBoundaryNode()
        self._AssignFarFieldBoundaryConditions(reference_inlet_node)

    def _FindFarthestUpstreamBoundaryNode(self):
        # The farthes upstream boundary node is the node with smallest
        # projection of its distante to a given reference node onto
        # the free stream velocity.

        # Select a random node in the mesh as reference (e.g. the first node)
        for node in self.far_field_model_part.Nodes:
            reference_node = node
            break

        # Find the farthest upstream boundary node
        temporal_smallest_projection = 1e30
        for node in self.far_field_model_part.Nodes:
            distance_to_reference = node - reference_node
            # Projecting the distance to the reference onto the free stream velocity
            distance_projection = DotProduct(distance_to_reference, self.free_stream_velocity)

            if(distance_projection < temporal_smallest_projection):
                temporal_smallest_projection = distance_projection
                reference_inlet_node = node

        return reference_inlet_node

    def _AssignFarFieldBoundaryConditions(self, reference_inlet_node):
        # A Dirichlet condition is applyed at the inlet nodes and
        # a Neumann condition is applyed at the outlet nodes

        for cond in self.far_field_model_part.Conditions:
            normal = cond.GetValue(KratosMultiphysics.NORMAL)

            # Computing the projection of the free stream velocity onto the normal
            velocity_projection = DotProduct(normal, self.free_stream_velocity)

            if( velocity_projection < 0):
                # A negative projection means inflow (i.e. inlet condition)
                self._AssignDirichletFarFieldBoundaryCondition(reference_inlet_node, cond)
            else:
                # A positive projection means outlow (i.e. outlet condition)
                self._AssignNeumannFarFieldBoundaryCondition(cond)

    def _AssignDirichletFarFieldBoundaryCondition(self, reference_inlet_node, cond):
        for node in cond.GetNodes():
            # Computing the value of the potential at the inlet
            inlet_potential = DotProduct( node - reference_inlet_node, self.free_stream_velocity)

            # Fixing the potential at the inlet nodes
            node.Fix(CPFApp.VELOCITY_POTENTIAL)
            node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL,0,inlet_potential + self.inlet_potential_0)

            # Applying Dirichlet condition in the adjoint problem
            if self.far_field_model_part.HasNodalSolutionStepVariable(CPFApp.ADJOINT_VELOCITY_POTENTIAL):
                node.Fix(CPFApp.ADJOINT_VELOCITY_POTENTIAL)
                node.SetSolutionStepValue(CPFApp.ADJOINT_VELOCITY_POTENTIAL,0,inlet_potential)

    def _AssignNeumannFarFieldBoundaryCondition(self, cond):
        cond.SetValue(CPFApp.VELOCITY_INFINITY, self.free_stream_velocity)




