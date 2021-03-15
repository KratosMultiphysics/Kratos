import KratosMultiphysics
from KratosMultiphysics import Parameters, Logger
import KratosMultiphysics.CompressiblePotentialFlowApplication as KCPFApp
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface
import KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis as potential_flow_analysis
import time as timer
import math

def DotProduct(A,B):
    result = 0
    for i,j in zip(A,B):
        result += i*j
    return result

def _GetModelPart(model, solver_settings):
    #TODO can be removed once model is fully available
    model_part_name = solver_settings["model_part_name"].GetString()
    if not model.HasModelPart(model_part_name):
        model_part = model.CreateModelPart(model_part_name, 2)
        domain_size = solver_settings["domain_size"].GetInt()
        if domain_size < 0:
            raise Exception('Please specify a "domain_size" >= 0!')
        model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)
    else:
        model_part = model.GetModelPart(model_part_name)

    return model_part

# ==============================================================================
class AdjointResponseFunction(ResponseFunctionInterface):
    """Linear static adjoint strain energy response function.
    - runs the primal analysis (writes the primal results to an .h5 file)
    - reads the primal results from the .h5 file into the adjoint model part
    - uses primal results to calculate value
    - uses primal results to calculate gradient by running the adjoint analysis

    Attributes
    ----------
    primal_analysis : Primal analysis object of the response function
    adjoint_analysis : Adjoint analysis object of the response function
    """
    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier
        self.response_settings = response_settings
        # Create the primal solver
        with open(self.response_settings["primal_settings"].GetString(),'r') as parameter_file:
            primal_parameters = Parameters( parameter_file.read() )

        if not primal_parameters.Has("reform_dofs_at_each_step") or not primal_parameters["reform_dofs_at_each_step"].GetBool():
            if not primal_parameters.Has("reform_dofs_at_each_step"):
                primal_parameters.AddEmptyValue("reform_dofs_at_each_step")
            primal_parameters["reform_dofs_at_each_step"].SetBool(True)
            wrn_msg = 'This solver requires the setting reform the dofs at each step in optimization.'
            wrn_msg += 'The solver setting has been set to True')
            Logger.PrintWarning(self._GetLabel(), wrn_msg)

        self.primal_model_part = _GetModelPart(model, primal_parameters["solver_settings"])

        self.primal_analysis = potential_flow_analysis.PotentialFlowAnalysis(model, primal_parameters)

        self.primal_data_transfer_with_python = self.response_settings["primal_data_transfer_with_python"].GetBool()

        # Create the adjoint solver
        adjoint_parameters = self._GetAdjointParameters()
        adjoint_model = KratosMultiphysics.Model()
        self.adjoint_model_part = _GetModelPart(adjoint_model, adjoint_parameters["solver_settings"])

        # TODO find out why it is not possible to use the same model_part
        self.adjoint_analysis = potential_flow_analysis.PotentialFlowAnalysis(adjoint_model, adjoint_parameters)

        self.primal_state_variables = [KCPFApp.VELOCITY_POTENTIAL, KCPFApp.AUXILIARY_VELOCITY_POTENTIAL]

    def Initialize(self):
        self.primal_analysis.Initialize()
        self.adjoint_analysis.Initialize()

    def InitializeSolutionStep(self):
        # Run the primal analysis.
        # TODO if primal_analysis.status==solved: return
        Logger.PrintInfo(self._GetLabel(), "Starting primal analysis for response:", self.identifier)
        Logger.PrintInfo("\n> Starting primal analysis for response:", self.identifier)
        startTime = timer.time()
        if not self.primal_analysis.time < self.primal_analysis.end_time:
            self.primal_analysis.end_time += 1
        self.primal_analysis.RunSolutionLoop()
        Logger.PrintInfo(self._GetLabel(), "Time needed for solving the primal analysis = ",round(timer.time() - startTime,2),"s")

    def CalculateValue(self):
        startTime = timer.time()
        self._GetResponseFunctionUtility().InitializeSolutionStep()
        value = self._GetResponseFunctionUtility().CalculateValue(self.primal_model_part)
        Logger.PrintInfo(self._GetLabel(), "Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")
        print("LIFT VALUE", value)
        self._value = value

    def CalculateGradient(self):
        # synchronize the modelparts
        self._SynchronizeAdjointFromPrimal()
        startTime = timer.time()
        Logger.PrintInfo("\n> Starting adjoint analysis for response:", self.identifier)
        if not self.adjoint_analysis.time < self.adjoint_analysis.end_time:
            self.adjoint_analysis.end_time += 1
        self.adjoint_analysis.RunSolutionLoop()
        Logger.PrintInfo("> Time needed for solving the adjoint analysis = ",round(timer.time() - startTime,2),"s")

    def GetValue(self):
        #switching to negative to minimize the negative of lift
        return -self._value

    def GetNodalGradient(self, variable):
        if variable != KratosMultiphysics.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))
        gradient = {}
        for node in self.adjoint_model_part.Nodes:
            #switching to negative to minimize the negative of lift
            gradient[node.Id] = -1*node.GetSolutionStepValue(variable)

        return gradient

    def Finalize(self):
        self.primal_analysis.Finalize()
        self.adjoint_analysis.Finalize()

    def _GetResponseFunctionUtility(self):
        return self.adjoint_analysis._GetSolver()._GetResponseFunction()

    def _SynchronizeAdjointFromPrimal(self):
        Logger.PrintInfo(self._GetLabel(), "Synchronize primal and adjoint modelpart for response:", self.identifier)

        if len(self.primal_model_part.Nodes) != len(self.adjoint_model_part.Nodes):
            raise RuntimeError("_SynchronizeAdjointFromPrimal: Model parts have a different number of nodes!")

        # TODO this should happen automatically
        for primal_node, adjoint_node in zip(self.primal_model_part.Nodes, self.adjoint_model_part.Nodes):
            adjoint_node.X0 = primal_node.X0
            adjoint_node.Y0 = primal_node.Y0
            adjoint_node.Z0 = primal_node.Z0
            adjoint_node.X = primal_node.X
            adjoint_node.Y = primal_node.Y
            adjoint_node.Z = primal_node.Z

        # Put primal solution on adjoint model
        if self.primal_data_transfer_with_python:
            Logger.PrintInfo(self._GetLabel(), "Transfer primal state to adjoint model part.")
            variable_utils = KratosMultiphysics.VariableUtils()
            for variable in self.primal_state_variables:
                variable_utils.CopyModelPartNodalVar(variable, self.primal_model_part, self.adjoint_model_part, 0)


    def _GetAdjointParameters(self):
        adjoint_settings = self.response_settings["adjoint_settings"].GetString()
        with open(self.response_settings["adjoint_settings"].GetString(),'r') as parameter_file:
            adjoint_parameters = Parameters( parameter_file.read() )

        return adjoint_parameters

    def _GetLabel(self):
        type_labels = {
            "adjoint_lift_potential_jump" : "LiftPotentialJump"
        }
        response_type = self.response_settings["response_type"].GetString()
        return "Adjoint" + type_labels[response_type]  +"Response"



class AngleOfAttackResponseFunction(ResponseFunctionInterface):
    def __init__(self, response_id, response_settings, model):
        self.model = model

    def Initialize(self):
        self.main_model_part = self.model["MainModelPart"]
        self.free_stream_velocity = KratosMultiphysics.Vector(3,0.0)
        self.free_stream_velocity[0] = 1.0
        for node in self.main_model_part.GetSubModelPart("TrailingEdgeNode").Nodes:
            self.te_node = node
            break
        for node in self.main_model_part.GetSubModelPart("LeadingEdgeNode").Nodes:
            self.le_node = node
            break

    def CalculateValue(self):
        self.aoa = self._CalculateAOA(self.te_node.X,self.te_node.Y, self.le_node.X, self.le_node.Y)
        print("AOA:", self.aoa)
    def CalculateGradient(self):
        epsilon = 1e-9
        self.te_x_gradient = (self._CalculateAOA(self.te_node.X + epsilon,self.te_node.Y, self.le_node.X, self.le_node.Y) - self.aoa) / epsilon
        self.te_y_gradient = (self._CalculateAOA(self.te_node.X ,self.te_node.Y+ epsilon, self.le_node.X, self.le_node.Y) - self.aoa) / epsilon
        self.le_x_gradient = (self._CalculateAOA(self.te_node.X ,self.te_node.Y, self.le_node.X+ epsilon, self.le_node.Y) - self.aoa) / epsilon
        self.le_y_gradient = (self._CalculateAOA(self.te_node.X ,self.te_node.Y, self.le_node.X, self.le_node.Y+ epsilon) - self.aoa) / epsilon

    def GetValue(self):
        return self.aoa

    def GetNodalGradient(self, variable):
        zero_vector = KratosMultiphysics.Vector(3, 0.0)
        gradient ={node.Id : zero_vector for node in self.main_model_part.Nodes}

        shape_gradient = KratosMultiphysics.Vector(3, 0.0)
        shape_gradient[0] = self.te_x_gradient
        shape_gradient[1] = self.te_y_gradient
        gradient[self.te_node.Id] = shape_gradient

        shape_gradient = KratosMultiphysics.Vector(3, 0.0)
        shape_gradient[0] = self.le_x_gradient
        shape_gradient[1] = self.le_y_gradient
        gradient[self.le_node.Id] = shape_gradient

        return gradient

    def _CalculateAOA(self, te_x, te_y, le_x, le_y):

        velocity_norm = self.free_stream_velocity.norm_2()
        chord_vector = KratosMultiphysics.Vector(3, 0.0)
        chord_vector[0] = le_x-te_x
        chord_vector[1] = le_y-te_y
        chord_norm = chord_vector.norm_2()
        aoa = math.acos(DotProduct(self.free_stream_velocity, -1*chord_vector)/(velocity_norm*chord_norm))

        return aoa

class ChordLengthResponseFunction(ResponseFunctionInterface):


    def __init__(self, response_id, response_settings, model):
        self.model = model

    def Initialize(self):
        self.main_model_part = self.model["MainModelPart"]
        for node in self.main_model_part.GetSubModelPart("TrailingEdgeNode").Nodes:
            self.te_node = node
            break
        for node in self.main_model_part.GetSubModelPart("LeadingEdgeNode").Nodes:
            self.le_node = node
            break

    def _ComputeChord(self, te_x, te_y, le_x, le_y):
        chord = math.sqrt((te_x-le_x)**2+(te_y-le_y)**2)

        return chord

    def CalculateValue(self):
        self.chord = self._ComputeChord(self.te_node.X,self.te_node.Y, self.le_node.X,self.le_node.Y)
        print("CHORD:", self.chord )
    def CalculateGradient(self):
        epsilon = 1e-9
        self.te_x_gradient = (self._ComputeChord(self.te_node.X + epsilon,self.te_node.Y, self.le_node.X, self.le_node.Y) - self.chord) / epsilon
        self.te_y_gradient = (self._ComputeChord(self.te_node.X ,self.te_node.Y+ epsilon, self.le_node.X, self.le_node.Y) - self.chord) / epsilon
        self.le_x_gradient = (self._ComputeChord(self.te_node.X ,self.te_node.Y, self.le_node.X+ epsilon, self.le_node.Y) - self.chord) / epsilon
        self.le_y_gradient = (self._ComputeChord(self.te_node.X ,self.te_node.Y, self.le_node.X, self.le_node.Y+ epsilon) - self.chord) / epsilon

    def GetValue(self):
        return self.chord

    def GetNodalGradient(self, variable):
        zero_vector = KratosMultiphysics.Vector(3, 0.0)
        gradient ={node.Id : zero_vector for node in self.main_model_part.Nodes}

        shape_gradient = KratosMultiphysics.Vector(3, 0.0)
        shape_gradient[0] = self.te_x_gradient
        shape_gradient[1] = self.te_y_gradient
        gradient[self.te_node.Id] = shape_gradient

        shape_gradient = KratosMultiphysics.Vector(3, 0.0)
        shape_gradient[0] = self.le_x_gradient
        shape_gradient[1] = self.le_y_gradient
        gradient[self.le_node.Id] = shape_gradient

        return gradient



class PerimeterResponseFunction(ResponseFunctionInterface):

    def __init__(self, response_id, response_settings, model):
        self.model = model

    def Initialize(self):
        self.main_model_part = self.model["MainModelPart"]
        self.body_model_part = self.main_model_part.GetSubModelPart("Body2D_Body")

    def _ComputePerimeter(self,  model_part):
        return KCPFApp.PotentialFlowUtilities.CalculateArea(model_part.Conditions)

    def CalculateValue(self):
        pass

    def CalculateGradient(self):
        pass

    def GetValue(self):

        return self._ComputePerimeter(self.body_model_part)

    def GetNodalGradient(self, variable):
        gradient = {}
        epsilon = 1e-6
        initial_perimeter =  self._ComputePerimeter(self.body_model_part)
        print("PERIMETER", initial_perimeter)
        for node in self.body_model_part.Nodes:
            shape_gradient = KratosMultiphysics.Vector(3, 0.0)


            node.X += epsilon
            current_perimeter =  self._ComputePerimeter(self.body_model_part)
            shape_gradient[0] = (current_perimeter-initial_perimeter)/epsilon
            node.X -= epsilon

            node.Y += epsilon
            current_perimeter =  self._ComputePerimeter(self.body_model_part)
            shape_gradient[1] = (current_perimeter-initial_perimeter)/epsilon
            node.Y -= epsilon

            gradient[node.Id] = shape_gradient

        return gradient

