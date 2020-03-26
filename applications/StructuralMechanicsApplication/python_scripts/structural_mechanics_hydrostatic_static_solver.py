from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
import structural_mechanics_solver
import numpy as np
import time

def CreateSolver(model, custom_settings):
    return HydroStaticMechanicalSolver(model, custom_settings)


class HydroStaticMechanicalSolver(structural_mechanics_solver.MechanicalSolver):
    """The structural mechanics hydrostatic static solver.

    This class creates the mechanical solvers for hydrostatic analysis. It currently
    supports partition with linear free surface update, monolithic with volume conservation, monolithic
    strategies.

    Public member variables:
    arc_length_settings -- settings for the arc length method.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, model, custom_settings):
        # Set defaults and validate custom settings.
        default_settings = KratosMultiphysics.Parameters("""
        {
            "coupling_iterations_settings" : {
                "residual_absolute_tolerance": 1e-9,
                "max_iteration"     : 100,
                "coupling_scheme"   : "gauss_seidel",
                "relaxation"        :  1.0,
                "filtering_epsilon" : 1e-15,
                "skip_first_time_step" : false
            }
        }
        """)
        self.validate_and_transfer_matching_settings(custom_settings, default_settings)
        self.coupling_itertions_settings = default_settings["coupling_iterations_settings"]
        
        # Validate the remaining settings in the base class.

        # Construct the base solver.
        super(HydroStaticMechanicalSolver, self).__init__(model, custom_settings)
        self.analysis_type = self.settings["analysis_type"].GetString()
        self.print_on_rank_zero("::[HydroStaticMechanicalSolver]:: ", "Construction finished")

    def Initialize(self):
        self.print_on_rank_zero("::[HydroStaticMechanicalSolver]:: ", "Initializing ...")
        super(HydroStaticMechanicalSolver, self).Initialize() # The mechanical solver is created here.
        self.skip_first_time_step = self.coupling_itertions_settings["skip_first_time_step"].GetBool()         
        self.print_on_rank_zero("::[HydroStaticMechanicalSolver]:: ", "Finished initialization.")
        
        if "partitioned" in self.analysis_type :
            with open("coupling_iterations.csv", "w") as f:
                f.write("iteration,computing_time\n")
        else:
            with open("computing_time.csv", "w") as f:
                f.write("computing_time\n")

    
            
    def SolveSolutionStep(self):
        start_time = time.process_time()
        ########## For conserving the volume with structural solver########
        
        if "partitioned" in self.analysis_type :
            self.coupling_scheme_type = self.coupling_itertions_settings["coupling_scheme"].GetString()
            self.IntializeCouplingScheme()
            if(self.skip_first_time_step):
                super(HydroStaticMechanicalSolver, self).SolveSolutionStep()
                super(HydroStaticMechanicalSolver, self).FinalizeSolutionStep()
            else:
                converged_signal = False
                max_iteration = self.coupling_itertions_settings["max_iteration"].GetInt()
                tolerance = self.coupling_itertions_settings["residual_absolute_tolerance"].GetDouble()
                displacements = self.GetDisplacements()
                
                for outer_iteration in range(0, max_iteration):
                    self.ApplyDisplacements(displacements)
                    #Conservation of volume happens in Solver.InitializeSolutionStep()
                    if (outer_iteration > 0): # IntializeSolutionStep is already called in the first iteration
                        super(HydroStaticMechanicalSolver, self).InitializeSolutionStep() 

                    super(HydroStaticMechanicalSolver, self).SolveSolutionStep()
                    super(HydroStaticMechanicalSolver, self).FinalizeSolutionStep()
                    displacements_ = self.GetDisplacements()
                    res = displacements_ - displacements
                    mag_res = np.linalg.norm(res)
                    KratosMultiphysics.Logger.Print("____________________________________________________"+"\n")            
                    KratosMultiphysics.Logger.PrintInfo(" Coupling iteration Nr ", outer_iteration+1)
                    KratosMultiphysics.Logger.PrintInfo(" Max residual magnitude for coupling iteration ", tolerance)
                    KratosMultiphysics.Logger.PrintInfo(" Residual magnitude for FSI iteration ", mag_res)
                    KratosMultiphysics.Logger.Print("_____________________________________________________"+"\n")
                    if (mag_res < tolerance):
                        converged_signal = True
                        break
                    else:
                        displacements += self.coupling_scheme.ComputeUpdate(res,displacements)
                
                if (converged_signal):
                    KratosMultiphysics.Logger.PrintInfo("Solver msg ", "Coupling iterations convergence achieved")
                    end_time = time.process_time()
                    computing_time = end_time-start_time
                    with open("coupling_iterations.csv", "a") as f:
                        f.write(str(outer_iteration+1)+" , "+ str(computing_time)+"\n")
                else:
                    KratosMultiphysics.Logger.PrintInfo("Solver msg ", "Max iterations reached")
                    raise Exception ("Partitioned solver didn't converge")
            self.skip_first_time_step = False                
        elif "monolithic" in self.analysis_type:
            super(HydroStaticMechanicalSolver, self).SolveSolutionStep()
            KratosMultiphysics.Logger.PrintInfo("Solver msg ","Monolithic Hydrostatic Solver End")
            end_time = time.process_time()
            computing_time = end_time-start_time
            with open("computing_time.csv", "a") as f:
                f.write(str(computing_time)+"\n")
        end_time = time.process_time()
        KratosMultiphysics.Logger.PrintInfo("Computing Time ", end_time-start_time)
        

    def FinalizeSolutionStep(self):
        super(HydroStaticMechanicalSolver, self).FinalizeSolutionStep()
        if "partitioned" in self.analysis_type :
            self.coupling_scheme.AdvanceTimeStep()

    def IntializeCouplingScheme(self):
        import sys
        sys.path.append("/lusers/temp/navaneeth/Software/coupling_schemes")
        from aitken import Aitken
        from iqnils import IQNILS 
        from iqnils_with_filtering import IQNILSWithFiltering 
        from gauss_seidel import GaussSeidel
        scheme_type = self.coupling_itertions_settings["coupling_scheme"].GetString()
        relaxation = self.coupling_itertions_settings["relaxation"].GetDouble()
        
        filtering_epsilon = self.coupling_itertions_settings["filtering_epsilon"].GetDouble()
        
        if scheme_type == "aitken":
            self.coupling_scheme = Aitken(relaxation,1.0)
        elif scheme_type == "iqnils":      
            self.coupling_scheme = IQNILS(100,1,relaxation)
        elif scheme_type == "iqnils_with_filtering":
            self.coupling_scheme = IQNILSWithFiltering(0,relaxation, filtering_epsilon)
        elif scheme_type == "gauss_seidel":
            self.coupling_scheme = GaussSeidel(relaxation)
        else : 
            err_msg =  "The requested coupling scheme \"" + scheme_type + "\" is not available!\n"
            err_msg += "Available options are: \"gauss_seidel\", \"aitken\", \"iqnils\", \"iqnils_with_filtering\""
            raise Exception(err_msg)
     

        

    def GetDisplacements(self):
        displacements = []
        computing_model_part = self.GetComputingModelPart() 
        for node in computing_model_part.Nodes:
            displacement_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,0)
            displacement_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0)
            displacement_z = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,0)
            displacements.append(displacement_x)
            displacements.append(displacement_y)
            displacements.append(displacement_z)
        return np.array(displacements)

    def ApplyDisplacements(self, displacements):
        i=0
        computing_model_part = self.GetComputingModelPart() 
        for node in computing_model_part.Nodes:       
            displacement_x = displacements[i]
            displacement_y = displacements[i+1]
            displacement_z = displacements[i+2]
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,displacement_x)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,displacement_y)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,displacement_z)

            node.X = node.X0 + displacement_x
            node.Y = node.Y0 + displacement_y
            node.Z = node.Z0 + displacement_z
            i = i+3
    

    
    def _create_solution_scheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

    def _create_mechanical_solution_strategy(self):
    
        if self.analysis_type == "partitioned":
            do_vol_cons = False
            do_line_search = False
            add_fluid_linearization  = False
            if (self.settings["line_search"].GetBool()):
                input("Warning::  Line search will only work for monolithic hydrostatic solvers, press key to continue")

            mechanical_solution_strategy = self._create_hydrostatic_strategy(do_line_search, do_vol_cons, add_fluid_linearization)
       
        elif self.analysis_type == "partitioned_with_linearization":
            do_vol_cons = False
            do_line_search = False
            add_fluid_linearization  = True
            if (self.settings["line_search"].GetBool()):
                input("Warning::  Line search will only work for monolithic hydrostatic solvers, press key to continue")
 
            mechanical_solution_strategy = self._create_hydrostatic_strategy(do_line_search, do_vol_cons, add_fluid_linearization) 

        elif self.analysis_type == "monolithic_vol_cons":
            do_vol_cons = True
            do_line_search = self.settings["line_search"].GetBool()
            add_fluid_linearization  = True
            mechanical_solution_strategy = self._create_hydrostatic_strategy(do_line_search,do_vol_cons, add_fluid_linearization)
        elif self.analysis_type == "monolithic_vol_cons_constraint": 
            do_vol_cons = False
            do_line_search = self.settings["line_search"].GetBool()
            add_fluid_linearization  = True 
            mechanical_solution_strategy = self._create_hydrostatic_strategy(do_line_search, do_vol_cons, add_fluid_linearization)
            self.prop_id_list = []
            for prop in self.GetComputingModelPart().Properties:
                if prop.Has(StructuralMechanicsApplication.FLUID_VOLUME):
                    self.prop_id_list.append(prop.Id)
            for prop_id in self.prop_id_list: 
                self.GetComputingModelPart().GetProperties()[prop_id].SetValue(StructuralMechanicsApplication.DO_UPDATE_FROM_DEL_VOL, True)   
        elif self.analysis_type ==  "monolithic_vol_cons_without_linearization":
            do_vol_cons = True
            do_line_search = self.settings["line_search"].GetBool()
            add_fluid_linearization  = False
            mechanical_solution_strategy = self._create_hydrostatic_strategy(do_line_search,do_vol_cons, add_fluid_linearization)
        else:
            err_msg =  "The requested analysis type \"" + self.analysis_type + "\" is not available!\n"
            err_msg += "Available options are: \"partitioned\", \"partitioned_with_linearization\", \"monolithic_vol_cons_constraint\", \"monolithic_vol_cons\", \"monolithic_vol_cons_without_linearization\""
            raise Exception(err_msg)
        return mechanical_solution_strategy

   
    def _create_hydrostatic_strategy(self, line_search, conserve_volume, add_fluid_linearization ):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self.get_solution_scheme()
        linear_solver = self.get_linear_solver()
        mechanical_convergence_criterion = self.get_convergence_criterion()
        builder_and_solver = self.get_builder_and_solver()

        return StructuralMechanicsApplication.NewtonRaphsonWithHydrostaticLoadStrategy(
                                                                computing_model_part,
                                                                mechanical_scheme,
                                                                linear_solver,
                                                                mechanical_convergence_criterion,
                                                                builder_and_solver,
                                                                self.settings["max_iteration"].GetInt(),
                                                                self.settings["compute_reactions"].GetBool(),
                                                                self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                self.settings["move_mesh_flag"].GetBool(),line_search, conserve_volume, add_fluid_linearization) 
