# Import necessary Kratos libraries
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *
import optimization_settings as opt_settings

print("\n> ==============================================================================================================")
print("> Initializing optimizer")
print("> ==============================================================================================================\n")

# Import shape optimizer factory from ShapeOptimizationApplication
import optimizer_factory

# Define an analyzer
def analyzer(X,controls,iterator,response):
    
    # ------------------------------------------------------------------------------------------------------
    # Input format:
    # ------------------------------------------------------------------------------------------------------
    #
    # X = { unique_node_id: [x,y,z],
    #       unique_node_id: [x,y,z],
    #       ... }
    #
    # controls = { "unique_response_id": {"calc_value": 1/0, "calc_gradient": 1/0},
    #              "unique_response_id": {"calc_value": 1/0, "calc_gradient": 1/0},
    #              ... }
    #
    # iterator = "mayor_itr.sub_itr" (Initial design iterator = "1.0" or "1.1")
    #
    # ------------------------------------------------------------------------------------------------------
    # Output format:
    # ------------------------------------------------------------------------------------------------------
    #
    # response = { "unique_response_id": {"value": xxx, "grad": {unique_node_id: [dx,dy,dz], ... }},
    #              "unique_response_id": {"value": xxx, "grad": {unique_node_id: [dx,dy,dz], ... }},
    #              ... }
    #
    # response must contain only those entries which were called through the controls
    #
    # ------------------------------------------------------------------------------------------------------

    if( iterator=="1.0" or iterator=="1.1" ):

        # Some initialization for the disciplinary solvers might go here
        
    # ------------------------------------------------------------------------------------------------------

    # Compute some primal field 
    # - Note that objective and constraint values may be computed within the same solve
    # - Note also that through the controls, we can make sure that whenever an adjoint gradient is to be 
    #   calculated, the primal field is computed first. At the same time , additional unnecessary primal 
    #   function calls may be avoided in case different adjoint solutions need to be caluclated based on the 
    #   same primal field. This is possible since through the controls we know from the beginning all 
    #   the function calls that the optimizer asks for in this step, so we may optimize the workflow for each
    #   combination.
    if( controls["some_objective_name"]["calc_value"] or 
        controls["some_objective_name"]["calc_gradient"] or 
        controls["some_constraint_name"]["calc_value"] or
        controls["some_constraint_name"]["calc_gradient"] ):
        
        # Run e.g. some CFD or CSM code
        # ...-> objective_value, constraint_value
      
        # Store objective value in response container
        response["some_objective_name"]["value"] = objective_value

        # Compute and store constraint
        response["some_constraint_name"]["value"] = constraint_value

    # ------------------------------------------------------------------------------------------------------

    # Compute adjoint field and ojbective gradient if required
    if( controls["some_objective_name"]["calc_gradient"] ):

        # Call disciplinary solver to compute objective gradients
        # ... -> dFdX

        # Store gradients in response container
        response["some_objective_name"]["gradient"] = dFdX

    # ------------------------------------------------------------------------------------------------------   
    
    # Compute adjoint field and constraint gradient if required
    if( controls["some_constraint_name"]["calc_gradient"] ):

        # Call disciplinary solver to compute constraint gradients

        # Store gradients in response container
        response["some_constraint_name"]["gradient"] = dFdX

    # ------------------------------------------------------------------------------------------------------ 

# Initalize model_part here to have it available for further use in this main script
design_surface = ModelPart("design_surface")

# Create an optimizer object 
optimizer = optimizer_factory.CreateOptimizer(design_surface,opt_settings,analyzer)

print("\n> ==============================================================================================================")
print("> Starting optimization")
print("> ==============================================================================================================\n")

# Call the optimize function of the optimizer
optimizer.optimize()

print("\n> ==============================================================================================================")
print("> Finished shape optimization!")
print("> ==============================================================================================================\n")