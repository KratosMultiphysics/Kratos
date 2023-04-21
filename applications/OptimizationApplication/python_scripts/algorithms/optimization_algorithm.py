# ==============================================================================
#  KratosOptimizationApplication
#
#  License:         BSD License
#                   license: OptimizationApplication/license.txt
#
#  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
#
# ===============================================================================

import KratosMultiphysics as KM
from KratosMultiphysics import Parameters, Logger
from KratosMultiphysics.vtk_output_process import VtkOutputProcess

# Additional imports
import csv
# ==============================================================================
class OptimizationAlgorithm:

    @abstractmethod

    def InitializeOptimizationLoop( self ):

        # create vtkIOs
        self.root_model_parts_vtkIOs = {}
        for root_model_part,data_fields in self.root_model_part_data_field_names.items():
            vtk_parameters = Parameters()
            root_controlling_model_part = self.model_parts_controller.GetRootModelPart(root_model_part)
            vtk_parameters.AddString("model_part_name",root_controlling_model_part.Name)
            vtk_parameters.AddBool("write_ids",False)
            vtk_parameters.AddString("file_format","ascii")
            vtk_parameters.AddBool("output_sub_model_parts",False)
            vtk_parameters.AddString("output_path","Optimization_Results")
            vtk_parameters.AddEmptyArray("nodal_solution_step_data_variables")
            for nodal_result in data_fields:
                vtk_parameters["nodal_solution_step_data_variables"].Append(nodal_result)
            
            controlling_model_part_vtkIO = VtkOutputProcess(self.model, vtk_parameters)
            controlling_model_part_vtkIO.ExecuteInitialize()
            controlling_model_part_vtkIO.ExecuteBeforeSolutionLoop()
            self.root_model_parts_vtkIOs[root_model_part] = controlling_model_part_vtkIO

    # --------------------------------------------------------------------------
    def RunOptimizationLoop( self ):
        raise RuntimeError("Algorithm base class is called. Please check your implementation of the function >> RunOptimizationLoop << .")

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop( self ):
        for vtkIO in self.root_model_parts_vtkIOs.values():
            vtkIO.ExecuteFinalize() 

    # --------------------------------------------------------------------------
    def SetResponseValue( self,response,value ):

        for objective in self.opt_parameters["objectives"]:
            if objective["name"].GetString() == response:
                objective["prev_itr_value"].SetDouble(objective["value"].GetDouble())
                objective["value"].SetDouble(value)
                if self.opt_parameters["opt_itr"].GetInt()<2:
                   objective["init_value"].SetDouble(value)   
                
                return  True

        for constraint in self.opt_parameters["constraints"]:
            if constraint["name"].GetString() == response:
                constraint["prev_itr_value"].SetDouble(constraint["value"].GetDouble())
                constraint["value"].SetDouble(value) 
                if self.opt_parameters["opt_itr"].GetInt()<2:
                    constraint["init_value"].SetDouble(value)
                    valid_types = ["initial_value_equality","smaller_than_initial_value","bigger_than_initial_value"]
                    if constraint["type"].GetString() in valid_types:
                        constraint["ref_value"].SetDouble(value)

                ref_value = constraint["ref_value"].GetDouble()

                type = constraint["type"].GetString()
                is_active = False
                if type == "equality" or type == "initial_value_equality":
                    is_active = True
                elif (type == "smaller_than" or type == "smaller_than_initial_value") and value>ref_value:
                    is_active = True
                elif (type == "bigger_than" or type == "bigger_than_initial_value")  and value<ref_value:
                    is_active = True

                if is_active:
                    self.num_active_consts += 1

                constraint["prev_itr_is_active"].SetBool(constraint["is_active"].GetBool())
                constraint["is_active"].SetBool(is_active)

                return is_active       

    # --------------------------------------------------------------------------
    def _InitializeCSVLogger(self):
        self.complete_log_file_name = "optimization_log.csv"
        with open(self.complete_log_file_name, 'w') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("{:>4s}".format("itr"))
            for objective in self.opt_parameters["objectives"]:
                row.append("{:>4s}".format("f: "+str(objective["name"].GetString())))
                row.append("{:>4s}".format("abs[%]"))
                row.append("{:>4s}".format("rel[%]"))
                row.append("{:>4s}".format("weight"))

            for constraint in self.opt_parameters["constraints"]:
                row.append("{:>4s}".format("g: "+str(constraint["name"].GetString())))
                row.append("{:>4s}".format("ref_val "))
                row.append("{:>4s}".format("ref_diff[%]"))
                row.append("{:>4s}".format("weight"))
                row.append("{:>4s}".format("is_active"))


            row.append("{:>4s}".format("projection_step_size"))
            if self.opt_parameters["constraints"].size() > 0:
                row.append("{:>4s}".format("correction_step_size"))
                row.append("{:>4s}".format("sin_alpha"))

            historyWriter.writerow(row)

    # --------------------------------------------------------------------------
    def _WriteCurrentOptItrToCSVFile( self ):
        
        with open(self.complete_log_file_name, 'a') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("{:>4d}".format(self.optimization_iteration))

            Logger.Print("")

            for objective in self.opt_parameters["objectives"]:
                objectivs_current_value = objective["value"].GetDouble()
                objectivs_weight = objective["weight"].GetDouble()
                objective_initial_value = objectivs_current_value
                objective_previous_value = objectivs_current_value

                if self.optimization_iteration>1:                
                    objective_initial_value = objective["init_value"].GetDouble() 
                    objective_previous_value = objective["prev_itr_value"].GetDouble() 

                rel_change = 100 * (objectivs_current_value-objective_previous_value)/objective_previous_value
                abs_change = 100 * (objectivs_current_value-objective_initial_value)/objective_initial_value

                Logger.Print("  ===== objective: ",objective["name"].GetString())
                Logger.Print("                   current value: ",objectivs_current_value)
                Logger.Print("                   rel_change: ",rel_change)
                Logger.Print("                   abs_change: ",abs_change)
                Logger.Print("                   weight: ",objectivs_weight)
                row.append(" {:> .5E}".format(objectivs_current_value))
                row.append(" {:> .5E}".format(abs_change))
                row.append(" {:> .5E}".format(rel_change)) 
                row.append(" {:> .5E}".format(objectivs_weight))

            for constraint in self.opt_parameters["constraints"]:
                constraint_current_value = constraint["value"].GetDouble() 
                constraint_ref_val = constraint["ref_value"].GetDouble()
                constraint_weight = constraint["weight"].GetDouble()
                is_active = constraint["is_active"].GetBool()
                abs_change = 100 * abs(constraint_current_value-constraint_ref_val)/abs(constraint_ref_val)

                Logger.Print("  ===== constraint: ",constraint["name"].GetString())
                Logger.Print("                   current value: ",constraint_current_value)
                Logger.Print("                   ref value: ",constraint_ref_val)
                Logger.Print("                   change: ",abs_change)
                Logger.Print("                   weight: ",constraint_weight)
                Logger.Print("                   is_active: ",is_active)

                row.append(" {:> .5E}".format(constraint_current_value))
                row.append(" {:> .5E}".format(constraint_ref_val))
                row.append(" {:> .5E}".format(abs_change))
                row.append(" {:> .5E}".format(constraint_weight))
                row.append(" {:> .5E}".format(is_active))
                

            Logger.Print("  ===== projection_step_size: ",self.opt_parameters["projection_step_size"].GetDouble())            
            row.append(" {:> .5E}".format(self.opt_parameters["projection_step_size"].GetDouble()))

            if self.opt_parameters["constraints"].size() > 0 and self.opt_parameters["num_active_consts"].GetInt()>0:
                Logger.Print("  ===== correction_step_size: ",self.opt_parameters["correction_step_size"].GetDouble())
                row.append(" {:> .5E}".format(self.opt_parameters["correction_step_size"].GetDouble()))
                Logger.Print("  ===== sin alpha: ",self.opt_parameters["sin_alpha"].GetDouble())
                row.append(" {:> .5E}".format(self.opt_parameters["sin_alpha"].GetDouble()))

            historyWriter.writerow(row)

# ==============================================================================
