## This script collects the available convergence criteria to be used in the SolidMechanicsApplication

from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

try:
  import KratosMultiphysics.MeshingApplication as MeshingApplication
  missing_meshing_dependencies = False
  missing_application = ''
except ImportError as e:
    missing_meshing_dependencies = True
    # extract name of the missing application from the error message
    import re
    missing_application = re.search(r'''.*'KratosMultiphysics\.(.*)'.*''','{0}'.format(e)).group(1)

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Convergence criteria class
class convergence_criterion:
    def __init__(self, convergence_criterion_parameters, model_part = None, adaptative_parameters = None, processes_list = None):
        # Note that all the convergence settings are introduced via a Kratos parameters object.
        
        D_RT = convergence_criterion_parameters["displacement_relative_tolerance"].GetDouble()
        D_AT = convergence_criterion_parameters["displacement_absolute_tolerance"].GetDouble()
        R_RT = convergence_criterion_parameters["residual_relative_tolerance"].GetDouble()
        R_AT = convergence_criterion_parameters["residual_absolute_tolerance"].GetDouble()
        
        echo_level = convergence_criterion_parameters["echo_level"].GetInt()
        
        if(echo_level >= 1):
            print("::[Mechanical Solver]:: CONVERGENCE CRITERION : ", convergence_criterion_parameters["convergence_criterion"].GetString())

        if (missing_meshing_dependencies == True):
            if (convergence_criterion_parameters["convergence_criterion"].GetString() == "AdaptativeErrorCriteria"):
                raise NameError('The AdaptativeErrorCriteria can not be used without compiling the MeshingApplication')

        rotation_dofs = False
        if(convergence_criterion_parameters.Has("rotation_dofs")):
            if(convergence_criterion_parameters["rotation_dofs"].GetBool()):
                rotation_dofs = True
        
        # Convergence criteria if there are rotation DOFs in the problem
        if(rotation_dofs == True):
            if(convergence_criterion_parameters["convergence_criterion"].GetString() == "Displacement_criterion"):
                self.mechanical_convergence_criterion = StructuralMechanicsApplication.DisplacementAndOtherDoFCriteria(D_RT, D_AT)
                self.mechanical_convergence_criterion.SetEchoLevel(echo_level)
                
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "Residual_criterion"):
                self.mechanical_convergence_criterion = StructuralMechanicsApplication.ResidualDisplacementAndOtherDoFCriteria(R_RT, R_AT)
                self.mechanical_convergence_criterion.SetEchoLevel(echo_level)
                
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "And_criterion"):
                Displacement = StructuralMechanicsApplication.DisplacementAndOtherDoFCriteria(D_RT, D_AT)
                Displacement.SetEchoLevel(echo_level)
                Residual = StructuralMechanicsApplication.ResidualDisplacementAndOtherDoFCriteria(R_RT, R_AT)
                Residual.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria(Residual, Displacement)
                
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "Or_criterion"):
                Displacement = StructuralMechanicsApplication.DisplacementAndOtherDoFCriteria(D_RT, D_AT)
                Displacement.SetEchoLevel(echo_level)
                Residual = StructuralMechanicsApplication.ResidualDisplacementAndOtherDoFCriteria(R_RT, R_AT)
                Residual.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KratosMultiphysics.OrCriteria(Residual, Displacement)
            
        # Convergence criteria without rotation DOFs        
        else:
            if(convergence_criterion_parameters["convergence_criterion"].GetString() == "Displacement_criterion"):
                self.mechanical_convergence_criterion = KratosMultiphysics.DisplacementCriteria(D_RT, D_AT)
                self.mechanical_convergence_criterion.SetEchoLevel(echo_level)
                
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "Residual_criterion"):
                self.mechanical_convergence_criterion = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
                self.mechanical_convergence_criterion.SetEchoLevel(echo_level)
                    
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "And_criterion"):
                Displacement = KratosMultiphysics.DisplacementCriteria(D_RT, D_AT)
                Displacement.SetEchoLevel(echo_level)
                Residual = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
                Residual.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria(Residual, Displacement)
                
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "Or_criterion"):
                Displacement = KratosMultiphysics.DisplacementCriteria(D_RT, D_AT)
                Displacement.SetEchoLevel(echo_level)
                Residual = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
                Residual.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KratosMultiphysics.OrCriteria(Residual, Displacement)
                
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "AdaptativeErrorCriteria"):
                self.mechanical_convergence_criterion = MeshingApplication.ErrorMeshCriteria(model_part, adaptative_parameters, processes_list)
                
        

