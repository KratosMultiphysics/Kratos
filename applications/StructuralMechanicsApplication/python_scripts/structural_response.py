# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics as KM
from KratosMultiphysics import Parameters, Logger
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

try:
    from KratosMultiphysics.MeshMovingApplication.mesh_moving_analysis import MeshMovingAnalysis
except ImportError as err:
    print("MeshMovingAnalysis not available ")


import time as timer

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
class StrainEnergyResponseFunction(ResponseFunctionInterface):
    """Linear strain energy response function. It triggers the primal analysis and
    uses the primal analysis results to calculate response value and gradient.

    Attributes
    ----------
    primal_model_part : Model part of the primal analysis object
    primal_analysis : Primal analysis object of the response function
    response_function_utility: Cpp utilities object doing the actual computation of response value and gradient.
    """

    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier
        self.model = model

        with open(response_settings["primal_settings"].GetString()) as parameters_file:
            ProjectParametersPrimal = Parameters(parameters_file.read())

        self.primal_model_part = _GetModelPart(model, ProjectParametersPrimal["solver_settings"])

        self.primal_analysis = StructuralMechanicsAnalysis(model, ProjectParametersPrimal)
        self.primal_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)

        self.response_function_utility = StructuralMechanicsApplication.StrainEnergyResponseFunctionUtility(self.primal_model_part, response_settings)

    def Initialize(self):
        self.primal_analysis.Initialize()
        self.response_function_utility.Initialize()

    def InitializeSolutionStep(self):
        self.primal_analysis.time = self.primal_analysis._GetSolver().AdvanceInTime(self.primal_analysis.time)
        self.primal_analysis.InitializeSolutionStep()

    def CalculateValue(self):
        Logger.PrintInfo("StrainEnergyResponse", "Starting primal analysis for response", self.identifier)

        startTime = timer.time()
        self.primal_analysis._GetSolver().Predict()
        self.primal_analysis._GetSolver().SolveSolutionStep()
        Logger.PrintInfo("StrainEnergyResponse", "Time needed for solving the primal analysis",round(timer.time() - startTime,2),"s")

        startTime = timer.time()
        value = self.response_function_utility.CalculateValue()
        self.primal_model_part.ProcessInfo[StructuralMechanicsApplication.RESPONSE_VALUE] = value
        Logger.PrintInfo("StrainEnergyResponse", "Time needed for calculating the response value",round(timer.time() - startTime,2),"s")

    def CalculateGradient(self):
        Logger.PrintInfo("StrainEnergyResponse", "Starting gradient calculation for response", self.identifier)

        startTime = timer.time()
        self.response_function_utility.CalculateGradient(False)
        Logger.PrintInfo("StrainEnergyResponse", "Time needed for calculating gradients",round(timer.time() - startTime,2),"s")

    def FinalizeSolutionStep(self):
        self.primal_analysis.FinalizeSolutionStep()
        self.primal_analysis.OutputSolutionStep()

    def Finalize(self):
        self.primal_analysis.Finalize()

    def GetValue(self):
        return self.primal_model_part.ProcessInfo[StructuralMechanicsApplication.RESPONSE_VALUE]

    def GetNodalGradient(self, variable):
        if variable != KratosMultiphysics.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))
        gradient = {}
        for node in self.primal_model_part.Nodes:
            gradient[node.Id] = node.GetSolutionStepValue(variable)
        return gradient


class StrainEnergyResponseFunctionWithStages(StrainEnergyResponseFunction):
    def CalculateValue(self):
        Logger.PrintInfo("StrainEnergyResponse", "Starting primal and gradient calculation for response", self.identifier)
        value = 0.0

        ## This order should be passed from the input JSON file.

        #modelpart_order = ["Parts_Solid_vol_1","Parts_Solid_vol_2","Parts_Solid_vol_3","Parts_Solid_vol_4"]
        #load_modelpart_order = ["SurfaceLoad3D_force_vol_1","SurfaceLoad3D_force_vol_2","SurfaceLoad3D_force_vol_3","SurfaceLoad3D_force_vol_4"]

        modelpart_order = ["Parts_Shell_part1","Parts_Shell_part3","Parts_Shell_part2","Parts_Shell_roof1", "Parts_Shell_roof3", "Parts_Shell_roof2"]
        #modelpart_order = ["Parts_Shell_part1","Parts_Shell_part3","Parts_Shell_part2"]
        load_modelpart_order = []

        startTime = timer.time()
        to_super_impose = False
        i = 1
        for model_part_name in modelpart_order:
            print()
            print("######")
            print("###### Activating the sub modelpart : ", model_part_name)
            print("######")
            KM.VariableUtils().SetFlag(KM.ACTIVE, False, self.primal_model_part.Elements)
            KM.VariableUtils().SetFlag(KM.ACTIVE, False, self.primal_model_part.Conditions)
            KM.VariableUtils().SetHistoricalVariableToZero(KM.MESH_DISPLACEMENT, self.primal_model_part.Nodes)
            for node in self.primal_model_part.Nodes:
                node.Free(KM.MESH_DISPLACEMENT_X)
                node.Free(KM.MESH_DISPLACEMENT_Y)
                node.Free(KM.MESH_DISPLACEMENT_Z)

            for j in range(0,i):
                elem_sub_model_part = self.primal_model_part.GetSubModelPart(modelpart_order[j])
                KM.VariableUtils().SetFlag(KM.ACTIVE, True, elem_sub_model_part.Elements)

                # cond_sub_model_part = self.primal_model_part.GetSubModelPart(load_modelpart_order[j])
                # KM.VariableUtils().SetFlag(KM.ACTIVE, True, cond_sub_model_part.Conditions)

                self.__Reset(elem_sub_model_part, KM.DISPLACEMENT)
                #self.__Reset(elem_sub_model_part, KM.ROTATION)

            self.primal_analysis.InitializeSolutionStep()
            self.primal_analysis._GetSolver().Predict()
            self.primal_analysis._GetSolver().SolveSolutionStep()
            self.primal_analysis._GetSolver().get_mechanical_solution_strategy().MoveMesh()
            self.primal_analysis.FinalizeSolutionStep()

            ## Use mesh motion solver to update the mesh.
            for k_mp in range(0,i):
                mp = self.primal_model_part.GetSubModelPart(modelpart_order[k_mp])
                for node in mp.Nodes:
                    node.SetSolutionStepValue(KM.MESH_DISPLACEMENT, node.GetSolutionStepValue(KM.DISPLACEMENT))
                for node in mp.Nodes:
                    node.Fix(KM.MESH_DISPLACEMENT_X)
                    node.Fix(KM.MESH_DISPLACEMENT_Y)
                    node.Fix(KM.MESH_DISPLACEMENT_Z)

            #self.__MoveMesh()

            value += self.response_function_utility.CalculateValue()

            self.response_function_utility.CalculateGradient(to_super_impose)
            to_super_impose = True
            i += 1
            self.primal_analysis.OutputSolutionStep()
            #input("Press Enter to continue...")
        #self.__UpdateInitialConfiguration()
        self.primal_model_part.ProcessInfo[StructuralMechanicsApplication.RESPONSE_VALUE] = value
        for node in self.primal_model_part.Nodes:
            node.Free(KM.MESH_DISPLACEMENT_X)
            node.Free(KM.MESH_DISPLACEMENT_Y)
            node.Free(KM.MESH_DISPLACEMENT_Z)

        Logger.PrintInfo("StrainEnergyResponse", "Time needed for primal and gradient calculation",round(timer.time() - startTime,2),"s")

    def CalculateGradient(self):
        pass

    def __Reset(self, mp, variable):
        for node in mp.Nodes:
            node.SetSolutionStepValue(variable, [0.0,0.0,0.0])

    def __MoveMesh(self):
        print("################ : StartMeshMovement")
        default_mesh_solver_settings = KM.Parameters("""
        {
            "apply_mesh_solver" : true,
            "problem_data":{
                "echo_level"    : 0,
                "start_time"    : 0.0,
                "end_time"      : 1.0,
                "parallel_type" : "OpenMP"
            },
            "solver_settings" : {
                "domain_size"     : 3,
                "echo_level"      : 0,
                "solver_type"     : "structural_similarity",
                "model_part_name" : "NONE",
                "model_import_settings"              : {
                    "input_type"     : "use_input_model_part"
                },
                "time_stepping" : {
                    "time_step"       : 1.0
                },
                "linear_solver_settings" : {
                    "solver_type" : "amgcl",
                    "smoother_type":"ilu0",
                    "krylov_type": "gmres",
                    "coarsening_type": "aggregation",
                    "max_iteration": 200,
                    "verbosity" : 0,
                    "tolerance": 1e-7
                },
                "compute_reactions"                : false,
                "calculate_mesh_velocity"          : false
            },
            "processes" : {
                "boundary_conditions_process_list" : []
            }
        }""")
        default_mesh_solver_settings["solver_settings"]["model_part_name"].SetString(self.primal_model_part.Name+"_MeshPart")
        mesh_moving_analysis = MeshMovingAnalysis(self.model, default_mesh_solver_settings)
        mesh_moving_analysis.Initialize()
        mesh_moving_analysis.RunSolutionLoop()
        mesh_moving_analysis.Finalize()
        print("################ : EndMeshMovement")

    def __UpdateInitialConfiguration(self):
        for node in self.primal_model_part.Nodes:
            node.X = node.X0
            #pass

# ==============================================================================
class EigenFrequencyResponseFunction(StrainEnergyResponseFunction):
    """Eigenfrequency response function. The internal procedure is the same as
    for the StrainEnergyResponseFunction. It triggers the primal analysis and
    uses the primal analysis results to calculate response value and gradient.
    Only the response_function_utility is a different object.

    Attributes
    ----------
    primal_model_part : Model part of the primal analysis object
    primal_analysis : Primal analysis object of the response function
    response_function_utility: Cpp utilities object doing the actual computation of response value and gradient.
    """

    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier

        with open(response_settings["primal_settings"].GetString()) as parameters_file:
            ProjectParametersPrimal = Parameters(parameters_file.read())

        eigen_solver_settings = ProjectParametersPrimal["solver_settings"]["eigensolver_settings"]

        max_required_eigenfrequency = int(max(response_settings["traced_eigenfrequencies"].GetVector()))
        if max_required_eigenfrequency is not eigen_solver_settings["number_of_eigenvalues"].GetInt():
            Logger.PrintWarning("EigenFrequencyResponse", "Specified number of eigenvalues in the primal analysis and the max required eigenvalue according the response settings do not match!!!")
            Logger.PrintWarning("EigenFrequencyResponse", "Primal parameters were adjusted accordingly!\n")
            eigen_solver_settings["number_of_eigenvalues"].SetInt(max_required_eigenfrequency)

        if not eigen_solver_settings.Has("normalize_eigenvectors"):
            eigen_solver_settings.AddEmptyValue("normalize_eigenvectors")
            eigen_solver_settings["normalize_eigenvectors"].SetBool(True)
            Logger.PrintWarning("EigenFrequencyResponse", "Eigenfrequency response function requires mass normalization of eigenvectors!")
            Logger.PrintWarning("EigenFrequencyResponse", "Primal parameters were adjusted accordingly!\n")

        if not eigen_solver_settings["normalize_eigenvectors"].GetBool():
            eigen_solver_settings["normalize_eigenvectors"].SetBool(True)
            Logger.PrintWarning("EigenFrequencyResponse", "Eigenfrequency response function requires mass normalization of eigenvectors!")
            Logger.PrintWarning("EigenFrequencyResponse", "Primal parameters were adjusted accordingly!\n")

        self.primal_model_part = _GetModelPart(model, ProjectParametersPrimal["solver_settings"])

        self.primal_analysis = StructuralMechanicsAnalysis(model, ProjectParametersPrimal)
        self.primal_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)

        self.response_function_utility = StructuralMechanicsApplication.EigenfrequencyResponseFunctionUtility(self.primal_model_part, response_settings)

# ==============================================================================
class MassResponseFunction(ResponseFunctionInterface):
    """Mass response function. It reads the materials for the model part and
    calculates response value and gradient.

    Attributes
    ----------
    model_part : Model part object of the response function
    response_function_utility: Cpp utilities object doing the actual computation of response value and gradient.
    """

    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier

        self.response_settings = response_settings
        self.model = model
        self.model_part_needs_to_be_imported = False

        model_part_name = response_settings["model_part_name"].GetString()
        input_type = response_settings["model_import_settings"]["input_type"].GetString()
        if input_type == "mdpa":
            self.model_part = self.model.CreateModelPart(model_part_name, 2)
            domain_size = response_settings["domain_size"].GetInt()
            if domain_size not in [2, 3]:
                raise Exception("MassResponseFunction: Invalid 'domain_size': {}".format(domain_size))
            self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)
            self.model_part_needs_to_be_imported = True
        elif input_type == "use_input_model_part":
            self.model_part = self.model.GetModelPart(model_part_name)
        else:
            raise Exception("Other model part input options are not yet implemented.")

        self.response_function_utility = StructuralMechanicsApplication.MassResponseFunctionUtility(self.model_part, response_settings)

        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)

    def Initialize(self):
        if self.model_part_needs_to_be_imported:
            # import model part
            model_part_io = KratosMultiphysics.ModelPartIO(self.response_settings["model_import_settings"]["input_filename"].GetString())
            model_part_io.ReadModelPart(self.model_part)

        # Add constitutive laws and material properties from json file to model parts.
        material_settings = KratosMultiphysics.Parameters("""{"Parameters": {} }""")
        materials_file_name = self.response_settings["material_import_settings"]["materials_filename"]
        material_settings["Parameters"].AddValue("materials_filename", materials_file_name)
        KratosMultiphysics.ReadMaterialsUtility(material_settings, self.model)
        self.response_function_utility.Initialize()

    def CalculateValue(self):
        Logger.PrintInfo("MassResponse", "Starting primal analysis for response", self.identifier)

        startTime = timer.time()
        value = self.response_function_utility.CalculateValue()
        self.model_part.ProcessInfo[StructuralMechanicsApplication.RESPONSE_VALUE] = value
        Logger.PrintInfo("MassResponse", "Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")

    def CalculateGradient(self):
        Logger.PrintInfo("MassResponse", "Starting gradient calculation for response", self.identifier)

        startTime = timer.time()
        self.response_function_utility.CalculateGradient()
        Logger.PrintInfo("MassResponse", "Time needed for calculating gradients",round(timer.time() - startTime,2),"s")

    def GetValue(self):
        return self.model_part.ProcessInfo[StructuralMechanicsApplication.RESPONSE_VALUE]

    def GetNodalGradient(self, variable):
        if variable != KratosMultiphysics.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))
        gradient = {}
        for node in self.model_part.Nodes:
            gradient[node.Id] = node.GetSolutionStepValue(variable)
        return gradient

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

        self.primal_model_part = _GetModelPart(model, primal_parameters["solver_settings"])

        self.primal_analysis = StructuralMechanicsAnalysis(model, primal_parameters)

        self.primal_data_transfer_with_python = self.response_settings["primal_data_transfer_with_python"].GetBool()

        # Create the adjoint solver
        adjoint_parameters = self._GetAdjointParameters()
        adjoint_model = KratosMultiphysics.Model()
        self.adjoint_model_part = _GetModelPart(adjoint_model, adjoint_parameters["solver_settings"])

        # TODO find out why it is not possible to use the same model_part
        self.adjoint_analysis = StructuralMechanicsAnalysis(adjoint_model, adjoint_parameters)

        self.primal_state_variables = [KratosMultiphysics.DISPLACEMENT]
        if primal_parameters["solver_settings"].Has("rotation_dofs"):
            if primal_parameters["solver_settings"]["rotation_dofs"].GetBool():
                self.primal_state_variables.append(KratosMultiphysics.ROTATION)

    def Initialize(self):
        self.primal_analysis.Initialize()
        self.adjoint_analysis.Initialize()

    def InitializeSolutionStep(self):
        # Run the primal analysis.
        # TODO if primal_analysis.status==solved: return
        Logger.PrintInfo(self._GetLabel(), "Starting primal analysis for response:", self.identifier)
        startTime = timer.time()
        if not self.primal_analysis.time < self.primal_analysis.end_time:
            self.primal_analysis.end_time += 1
        self.primal_analysis.RunSolutionLoop()
        Logger.PrintInfo(self._GetLabel(), "Time needed for solving the primal analysis = ",round(timer.time() - startTime,2),"s")

    def CalculateValue(self):
        startTime = timer.time()
        value = self._GetResponseFunctionUtility().CalculateValue(self.primal_model_part)
        Logger.PrintInfo(self._GetLabel(), "Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")

        self.primal_model_part.ProcessInfo[StructuralMechanicsApplication.RESPONSE_VALUE] = value

    def CalculateGradient(self):
        # synchronize the modelparts
        self._SynchronizeAdjointFromPrimal()
        startTime = timer.time()
        Logger.PrintInfo(self._GetLabel(), "Starting adjoint analysis for response:", self.identifier)
        if not self.adjoint_analysis.time < self.adjoint_analysis.end_time:
            self.adjoint_analysis.end_time += 1
        self.adjoint_analysis.RunSolutionLoop()
        Logger.PrintInfo(self._GetLabel(), "Time needed for solving the adjoint analysis = ",round(timer.time() - startTime,2),"s")

    def GetValue(self):
        return self.primal_model_part.ProcessInfo[StructuralMechanicsApplication.RESPONSE_VALUE]

    def GetNodalGradient(self, variable):
        if variable != KratosMultiphysics.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))
        gradient = {}
        for node in self.adjoint_model_part.Nodes:
            gradient[node.Id] = node.GetSolutionStepValue(variable)
        return gradient

    def Finalize(self):
        self.primal_analysis.Finalize()
        self.adjoint_analysis.Finalize()

    def _GetResponseFunctionUtility(self):
        return self.adjoint_analysis._GetSolver().response_function

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

        if adjoint_settings == "auto":
            Logger.PrintInfo(self._GetLabel(), "Automatic set up adjoint parameters for response:", self.identifier)

            if not self.primal_data_transfer_with_python:
                raise Exception("Auto setup of adjoint parameters does only support primal data transfer with python.")

            with open(self.response_settings["primal_settings"].GetString(),'r') as parameter_file:
                primal_parameters = Parameters( parameter_file.read() )

            # check that HDF5 process is not there
            if primal_parameters["processes"].Has("list_other_processes"):
                for i in range(0,primal_parameters["processes"]["list_other_processes"].size()):
                    process = primal_parameters["processes"]["list_other_processes"][i]
                    raise Exception("Auto setup of adjoint parameters does not support {} in list_other_processes".format(process["python_module"].GetString()))

            # clone primal settings as base for adjoint
            adjoint_parameters = primal_parameters.Clone()

            # analysis settings
            solver_settings = adjoint_parameters["solver_settings"]
            primal_solver_type = solver_settings["solver_type"].GetString()
            if primal_solver_type != "static":
                raise Exception("Auto setup of adjoint parameters does not support {} solver_type. Only available for 'static'".format(primal_solver_type))
            solver_settings["solver_type"].SetString("adjoint_"+primal_solver_type)

            if not solver_settings.Has("compute_reactions"):
                solver_settings.AddEmptyValue("compute_reactions")
            solver_settings["compute_reactions"].SetBool(False)

            if not solver_settings.Has("move_mesh_flag"):
                solver_settings.AddEmptyValue("move_mesh_flag")
            solver_settings["move_mesh_flag"].SetBool(False)

            if solver_settings.Has("scheme_settings"):
                depr_msg = '\nDEPRECATION-WARNING: "scheme_settings" is deprecated, please remove it from your json parameters.\n'
                Logger.PrintWarning(__name__, depr_msg)
                solver_settings.RemoveValue("scheme_settings")

            if solver_settings["model_import_settings"]["input_type"].GetString() == "use_input_model_part":
                solver_settings["model_import_settings"]["input_type"].SetString("mdpa")
                if solver_settings["model_import_settings"].Has("input_filename"):
                    file_name = solver_settings["model_import_settings"]["input_filename"].GetString()
                else:
                    Logger.PrintWarning(self._GetLabel(), "Automatic adjoint settings creator assumes the model_part_name as input_filename.")
                    solver_settings["model_import_settings"].AddEmptyValue("input_filename")
                    file_name = solver_settings["model_part_name"].GetString()
                solver_settings["model_import_settings"]["input_filename"].SetString(file_name)

            # Dirichlet conditions: change variables
            for i in range(0,primal_parameters["processes"]["constraints_process_list"].size()):
                process = adjoint_parameters["processes"]["constraints_process_list"][i]
                variable_name = process["Parameters"]["variable_name"].GetString()
                process["Parameters"]["variable_name"].SetString("ADJOINT_"+variable_name)

            # Neumann conditions - do not modify to read the same load values as in primal:

            # Output process:
            # TODO how to add the output process? How find out about the variables?
            if adjoint_parameters.Has("output_processes"):
                Logger.PrintInfo(self._GetLabel(), "Output process is removed for adjoint analysis. To enable it define adjoint_parameters yourself.")
                adjoint_parameters.RemoveValue("output_processes")

            # sensitivity settings
            adjoint_parameters["solver_settings"].AddValue("sensitivity_settings", self.response_settings["sensitivity_settings"])

            # response settings
            adjoint_parameters["solver_settings"].AddValue("response_function_settings", self.response_settings)

        else: # adjoint parameters file is explicitely given - do not change it.
            with open(self.response_settings["adjoint_settings"].GetString(),'r') as parameter_file:
                adjoint_parameters = Parameters( parameter_file.read() )

        return adjoint_parameters


    def _GetLabel(self):
        type_labels = {
            "adjoint_nodal_displacement" : "NodalDisplacement",
            "adjoint_linear_strain_energy" : "StrainEnergy",
            "adjoint_local_stress" : "LocalStress",
            "adjoint_max_stress" : "MaxStress",
            "adjoint_nodal_reaction" : "NodalReaction"
        }
        response_type = self.response_settings["response_type"].GetString()
        return "Adjoint" + type_labels[response_type] + "Response"
