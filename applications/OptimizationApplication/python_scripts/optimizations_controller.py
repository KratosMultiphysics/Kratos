# ==============================================================================
#  KratosOptimizationApplication
#
#  License:         BSD License
#                   license: OptimizationApplication/license.txt
#
#  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
#
# ==============================================================================

# Kratos Core and Apps
import KratosMultiphysics as KM

# Additional imports
import KratosMultiphysics.OptimizationApplication.algorithms.algorithm_steepest_descent as steepest_descent
import KratosMultiphysics.OptimizationApplication.algorithms.algorithm_gradient_projection as gradient_projection
import time as timer

# ==============================================================================
def CreateController(optimizations_settings,model,model_parts_controller,analyses_controller,responses_controller,controls_controller):
    return OptimizationsController(optimizations_settings,model,model_parts_controller,analyses_controller,responses_controller,controls_controller)

# ==============================================================================
class OptimizationsController:
    # --------------------------------------------------------------------------
    def __init__(self,optimizations_settings,model,model_parts_controller,analyses_controller,responses_controller,controls_controller):
        
        self.optimizations_settings = optimizations_settings
        self.controls_controller = controls_controller
        self.responses_controller = responses_controller
        self.model = model

        default_settings = KM.Parameters("""
        {
            "name"                : "OPT_NAME",
            "type"                : "OPT_TYPE",
            "settings"                : {
                "objectives": [],
                "objectives_improvements": [],
                "constraints": [],
                "constraints_types": [],
                "constraints_ref_values": [],
                "controls": [],
                "objectives_controls_weights": [],
                "constraints_controls_weights": [],
                "controls_lower_bounds": [],
                "controls_lower_bounds_values": [],
                "controls_upper_bounds": [],
                "controls_upper_bounds_values": [],                
                "algorithm" : "ALG_NAME",
                "algorithm_settings" : {
                    "max_iterations"     : 100,
                    "step_size"  : 0.1,
                    "relative_tolerance" : 1e-3
                }
            }
        }""")

        if not self.optimizations_settings.size()>0:
            raise RuntimeError("OptimizationsController: optimizations list in optimizer's parameters can not be empty !") 


        self.optimizations = {}
        self.optimizations_types={}


        self.optimizations_objectives={}
        self.optimizations_objectives_improvements={}


        self.optimizations_constraints={}
        self.optimizations_constraints_types={}
        self.optimizations_constraints_ref_values={}
        
        
        self.optimizations_controls={}
        self.optimizations_objectives_controls_weights={}
        self.optimizations_constraints_controls_weights={}


        self.optimizations_algorithm={}
        self.supported_opt_types = ["gradient_based"]
        self.supported_algorithms = ["gradient_projection"]



        for itr in range(self.optimizations_settings.size()):
            opt_settings = self.optimizations_settings[itr]
            opt_name = opt_settings["name"].GetString()
            opt_type = opt_settings["type"].GetString()   
            # check for name
            if opt_name in self.optimizations.keys():  
                raise RuntimeError("OptimizationsController: Optimization name '{}' already exists.".format(opt_name))    


            # check for type
            if not opt_type in self.supported_opt_types:  
                raise RuntimeError("OptimizationsController: Optimization type '{}' is not supported, supported types {}.".format(opt_type,self.supported_opt_types))                  
            self.optimizations_types[opt_name]=opt_type


            # checks for objectives
            objectives_names = opt_settings["settings"]["objectives"].GetStringArray()
            if not len(objectives_names)>0:  
                raise RuntimeError("OptimizationsController: Objectives list of optimization '{}' can not be empty.".format(opt_name))   
            if len(set(objectives_names)) != len(objectives_names):
                raise RuntimeError("OptimizationsController: Objectives list of optimization '{}' has duplicate response names '{}'.".format(opt_name,objectives_names))
            self.responses_controller.CheckIfResponsesExist(objectives_names)
            self.optimizations_objectives[opt_name]=objectives_names


            # check for objectives improvements
            if not opt_settings["settings"]["objectives_improvements"].IsVector():  
                raise RuntimeError("OptimizationsController:'objectives_improvements' of optimization '{}' should be vector of doubles.".format(opt_name)) 
            objectives_improvements = opt_settings["settings"]["objectives_improvements"].GetVector()     
            if not len(objectives_improvements) == len(objectives_names):
                raise RuntimeError("OptimizationsController:'objectives_improvements' of optimization '{}' should be of the same size of objectives list.".format(opt_name))
            self.optimizations_objectives_improvements[opt_name]=objectives_improvements


            # checks for constraints
            if opt_settings["settings"]["constraints"].size()>0:
                constraints_names = opt_settings["settings"]["constraints"].GetStringArray()
                if len(set(constraints_names)) != len(constraints_names):
                    raise RuntimeError("OptimizationsController: Constraint list of optimization '{}' has duplicate response names .".format(opt_name))
                self.responses_controller.CheckIfResponsesExist(constraints_names)

                for constraints_name in constraints_names:
                    if constraints_name in objectives_names:
                        raise RuntimeError("OptimizationsController: Response {} in optimization {} is used as both objective and constraint.".format(constraints_name,opt_name))
                self.optimizations_constraints[opt_name]=constraints_names

                constraints_types = opt_settings["settings"]["constraints_types"].GetStringArray()  
                if not len(constraints_types) == len(constraints_names):
                    raise RuntimeError("OptimizationsController:'constraints_types' of optimization '{}' should be of the same size of constraint list.".format(opt_name))
                for index, type in enumerate(constraints_types):
                    if not type in ["equality","smaller_than","bigger_than","initial_value_equality","smaller_than_initial_value","bigger_than_initial_value"]: 
                        raise RuntimeError("OptimizationsController: constraint type {} of constraint {} of optimization '{}' should be either 'equality' or 'inequality'.".format(type,constraints_names[index],opt_name,))
                self.optimizations_constraints_types[opt_name]=constraints_types                

                if not opt_settings["settings"]["constraints_ref_values"].IsVector():  
                    raise RuntimeError("OptimizationsController:'constraints_ref_values' of optimization '{}' should be vector of doubles.".format(opt_name)) 
                constraints_ref_values = opt_settings["settings"]["constraints_ref_values"].GetVector() 
                if not len(constraints_ref_values) == len(constraints_names):
                    raise RuntimeError("OptimizationsController:'constraints_ref_values' of optimization '{}' should be of the same size of constraint list.".format(opt_name))
                self.optimizations_constraints_ref_values[opt_name]=constraints_ref_values
                
            
            # checks for controls
            controls_names = opt_settings["settings"]["controls"].GetStringArray()
            if not len(controls_names)>0:  
                raise RuntimeError("OptimizationsController: Controls list of optimization '{}' can not be empty.".format(opt_name))   
            if len(set(controls_names)) != len(controls_names):
                raise RuntimeError("OptimizationsController: Controls list of optimization '{}' has duplicate control names .".format(opt_name))
            self.controls_controller.CheckIfControlsExist(controls_names)

            for key,values in self.optimizations_controls.items():
                compare_results = set(controls_names) & set(values)
                if len(compare_results):
                    raise RuntimeError("OptimizationsController: Controls {} can not be shared between optimizations '{}' and '{}' .".format(compare_results,key,opt_name))
            
            self.optimizations_controls[opt_name]=controls_names            


            # checks for objectives controls weights
            if opt_settings["settings"]["objectives_controls_weights"].size() != len(controls_names):
                raise RuntimeError("OptimizationsController: 'objectives_controls_weights' of optimization '{}' should be of the same size of control list .".format(opt_name))            

            objectives_controls_weights=[]
            for i in range(opt_settings["settings"]["objectives_controls_weights"].size()):
                if not opt_settings["settings"]["objectives_controls_weights"][i].IsNumber():
                    raise RuntimeError("OptimizationsController: entry {} of 'objectives_controls_weights' of optimization '{}' should be a number .".format(i+1,opt_name))
                else:
                    objectives_controls_weights.append(opt_settings["settings"]["objectives_controls_weights"][i].GetDouble())

            self.optimizations_objectives_controls_weights[opt_name] = objectives_controls_weights


            # checks for constraints controls weights
            if opt_settings["settings"]["constraints_controls_weights"].size() != len(controls_names):
                raise RuntimeError("OptimizationsController: 'constraints_controls_weights' of optimization '{}' should be of the same size of control list .".format(opt_name))            

            constraints_controls_weights=[]
            for i in range(opt_settings["settings"]["constraints_controls_weights"].size()):
                if not opt_settings["settings"]["constraints_controls_weights"][i].IsNumber():
                    raise RuntimeError("OptimizationsController: entry {} of 'constraints_controls_weights' of optimization '{}' should be a number .".format(i+1,opt_name))
                else:
                    constraints_controls_weights.append(opt_settings["settings"]["constraints_controls_weights"][i].GetDouble())

            self.optimizations_constraints_controls_weights[opt_name] = constraints_controls_weights

            
            # checks for algorithms settings
            algorithm = opt_settings["settings"]["algorithm"].GetString()
            if not algorithm in self.supported_algorithms:
                raise RuntimeError("OptimizationsController: Optimization algorithm '{}' is not supported, supported types {}.".format(algorithm,self.supported_algorithms))                  

            if algorithm == "gradient_projection":
                self.optimizations[opt_name] = gradient_projection.AlgorithmGradientProjection(opt_name,opt_settings["settings"],model,model_parts_controller,analyses_controller,responses_controller,controls_controller)
    
        
    # --------------------------------------------------------------------------
    def Initialize(self):
        for opt in self.optimizations.values():
            opt.InitializeOptimizationLoop()

    # --------------------------------------------------------------------------
    def Optimize(self,opt_name):
        if not opt_name in self.optimizations.keys():
            raise RuntimeError("OptimizationsController:Optimize: Optimization {} doesn not exist !.".format(opt_name))
        
        self.optimizations[opt_name].RunOptimizationLoop()
    # --------------------------------------------------------------------------
    def OptimizeAll(self):
        for opt in self.optimizations.values():
            opt.RunOptimizationLoop()



            

               

