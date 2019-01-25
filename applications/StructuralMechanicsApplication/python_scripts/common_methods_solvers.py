from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

def AuxiliarDampingSettings():
    damping_settings = KratosMultiphysics.Parameters("""
    {
        "damping_settings" :
        {
            "determine_rayleigh_damping" : false,
            "determine_rayleigh_damping_settings" : {
                "echo_level"          : 0,
                "write_on_properties" : true,
                "damping_ratio_0"     : 0.0,
                "damping_ratio_1"     : -1.0,
                "eigen_values_vector" : [0.0],
                "eigen_system_settings" : {
                    "solver_type"                : "eigen_eigensystem"
                }
            }
        }
    }
    """)

    return damping_settings

def AuxiliarEigenSettings(solver_type):
    if solver_type == "feast" or solver_type == "FEAST":
        eigen_system_settings = KratosMultiphysics.Parameters("""
        {
            "eigen_system_settings" : {
                "solver_type"                : "feast",
                "print_feast_output"         : false,
                "perform_stochastic_estimate": true,
                "solve_eigenvalue_problem"   : true,
                "lambda_min"                 : 0.0,
                "lambda_max"                 : 4.0e5,
                "number_of_eigenvalues"      : 2,
                "search_dimension"           : 15,
                "linear_solver_settings": {
                    "solver_type": "SkylineLUComplexSolver"
                }
            }
        }
        """)
    else:
        eigen_system_settings = KratosMultiphysics.Parameters("""
        {
            "eigen_system_settings" : {
                "solver_type"       : "eigen_eigensystem"
            }
        }
        """)
        eigen_system_settings["eigen_system_settings"]["solver_type"].SetString(solver_type)
    return eigen_system_settings

def ComputeDampingCoefficients(model, settings, damping_settings):
    import KratosMultiphysics.kratos_utilities as kratos_utils
    if kratos_utils.AreApplicationsAvailable(["ExternalSolversApplication", "EigenSolversApplication"]):
        if kratos_utils.IsApplicationAvailable("ExternalSolversApplication"):
            from KratosMultiphysics import ExternalSolversApplication
        if kratos_utils.IsApplicationAvailable("EigenSolversApplication"):
            from KratosMultiphysics import EigenSolversApplication
    else:
        raise Exception("ExternalSolversApplication or EigenSolversApplication not available")

    # The general damping ratios
    damping_ratio_0 = damping_settings["determine_rayleigh_damping_settings"]["damping_ratio_0"].GetDouble()
    damping_ratio_1 = damping_settings["determine_rayleigh_damping_settings"]["damping_ratio_1"].GetDouble()

    # We get the model parts which divide the problem
    main_model_part_name = settings["model_part_name"].GetString()
    current_process_info = model[main_model_part_name].ProcessInfo
    existing_computation = current_process_info.Has(StructuralMechanicsApplication.EIGENVALUE_VECTOR)
    structural_parts = ExtractModelParts(model, settings)
    for part in structural_parts:
        # Create auxiliar parameters
        compute_damping_coefficients_settings = KratosMultiphysics.Parameters("""
        {
            "echo_level"          : 0,
            "damping_ratio_0"     : 0.0,
            "damping_ratio_1"     : -1.0,
            "eigen_values_vector" : [0.0]
        }
        """)

        # Setting custom parameters
        compute_damping_coefficients_settings["echo_level"].SetInt(damping_settings["determine_rayleigh_damping_settings"]["echo_level"].GetInt())
        compute_damping_coefficients_settings["damping_ratio_0"].SetDouble(damping_ratio_0)
        compute_damping_coefficients_settings["damping_ratio_1"].SetDouble(damping_ratio_1)

        # We check if the values are previously defined
        properties = part.GetProperties()
        for prop in properties:
            if prop.Has(StructuralMechanicsApplication.SYSTEM_DAMPING_RATIO):
                damping_settings["determine_rayleigh_damping_settings"]["damping_ratio_0"].SetDouble(prop.GetValue(StructuralMechanicsApplication.SYSTEM_DAMPING_RATIO))
                break
        for prop in properties:
            if prop.Has(StructuralMechanicsApplication.SECOND_SYSTEM_DAMPING_RATIO):
                damping_settings["determine_rayleigh_damping_settings"]["damping_ratio_1"].SetDouble(prop.GetValue(StructuralMechanicsApplication.SECOND_SYSTEM_DAMPING_RATIO))
                break

        # We have computed already the eigen values
        current_process_info = part.ProcessInfo
        precomputed_eigen_values = damping_settings["determine_rayleigh_damping_settings"]["eigen_values_vector"].GetVector()
        if len(precomputed_eigen_values) > 1:
            compute_damping_coefficients_settings["eigen_values_vector"].SetVector(precomputed_eigen_values)
        else:
            # We check if BC defined in some node
            has_bc = False
            for node in part.Nodes:
                if node.IsFixed(KratosMultiphysics.DISPLACEMENT_X):
                    has_bc = True
                    break
                if node.IsFixed(KratosMultiphysics.DISPLACEMENT_Y):
                    has_bc = True
                    break
                if node.IsFixed(KratosMultiphysics.DISPLACEMENT_Z):
                    has_bc = True
                    break
            # If not computed eigen values already
            if not existing_computation and has_bc:
                KratosMultiphysics.Logger.PrintInfo("::[MechanicalSolver]::", "EIGENVALUE_VECTOR not previously computed. Computing automatically, take care")
                from KratosMultiphysics import eigen_solver_factory
                eigen_linear_solver = eigen_solver_factory.ConstructSolver(damping_settings["determine_rayleigh_damping_settings"]["eigen_system_settings"])
                builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolverWithConstraints(eigen_linear_solver)
                eigen_scheme = StructuralMechanicsApplication.EigensolverDynamicScheme()
                eigen_solver = StructuralMechanicsApplication.EigensolverStrategy(part, eigen_scheme, builder_and_solver)
                eigen_solver.Solve()
            
            if has_bc:
                eigenvalue_vector = current_process_info.GetValue(StructuralMechanicsApplication.EIGENVALUE_VECTOR)
                compute_damping_coefficients_settings["eigen_values_vector"].SetVector(eigenvalue_vector)

        # We compute the coefficients
        coefficients_vector = StructuralMechanicsApplication.ComputeDampingCoefficients(compute_damping_coefficients_settings)

        # We set the values
        if damping_settings["determine_rayleigh_damping_settings"]["write_on_properties"].GetBool():
            for prop in part.Properties:
                prop.SetValue(StructuralMechanicsApplication.RAYLEIGH_ALPHA, coefficients_vector[0])
                prop.SetValue(StructuralMechanicsApplication.RAYLEIGH_BETA, coefficients_vector[1])
        else:
            current_process_info.SetValue(StructuralMechanicsApplication.RAYLEIGH_ALPHA, coefficients_vector[0])
            current_process_info.SetValue(StructuralMechanicsApplication.RAYLEIGH_BETA, coefficients_vector[1])

def ExtractModelParts(model, settings):
    """Extracting a list of SubModelParts to be added to the ComputingModelPart
    If the MainModelPart is to be added to the ComputingModelPart, then this is
    done directly, no need to also add the SubModelParts, since they are contained
    in the MainModelpart
    """
    params_list = settings["problem_domain_sub_model_part_list"]
    list_model_part_names = [params_list[i].GetString() for i in range(params_list.size())]
    main_model_part_name = settings["model_part_name"].GetString()
    model_parts_list = []
    for model_part_name  in list_model_part_names:
        # only SubModelParts of the MainModelPart can be used!
        if model_part_name != main_model_part_name:
            full_model_part_name = main_model_part_name + "." + model_part_name
        else:
            full_model_part_name = main_model_part_name
        model_parts_list.append(model[full_model_part_name])
    return model_parts_list
