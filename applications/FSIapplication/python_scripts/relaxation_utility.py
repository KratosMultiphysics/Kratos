from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics

def CreateInterfaceUtility(coupling_strategy_parameters):
    return RelaxationUtility(coupling_strategy_parameters)


class RelaxationUtility:
    
    def __init__(self,coupling_strategy_parameters):
        
        # Common initialization for all strategies
        self.vector_space = KratosMultiphysics.UblasSparseSpace()
        
        # Relaxation method initialization
        self.acceleration_type = coupling_strategy_parameters["acceleration_type"].GetString()  # Acceleration types available: "Aitken", "Fixed", "Full"
        self.w_0 = coupling_strategy_parameters["w_0"].GetDouble()                              # Relaxation parameter initialization
        
        
    def ExecuteInitializeSolutionStep(self):
        pass
        
        
    def InterfaceSolutionUpdate(self,step,nl_it,iteration_guess_value,interface_residual):
        # Inputs:
        # - step: step id.
        # - nl_it: non-linear interface iteration id.
        # - iteration_guess_value: solution guess at the beginning of the iteration (array type)
        # - interface_residual: interface residual coming form the "iteration_guess_value" (array type)
                     
        if self.acceleration_type == "Aitken":
            if nl_it==1:
                self.res_1 = interface_residual

                val_correct = iteration_guess_value + self.w_0*interface_residual
                                    
            elif nl_it>1:
                
                self.res_2 = interface_residual

                aux1 = self.vector_space.Dot(self.res_1,self.res_2-self.res_1)
                aux2 = self.vector_space.Dot(self.res_2-self.res_1,self.res_2-self.res_1)

                w_1 = -self.w_0*(aux1/aux2)
                
                val_correct = iteration_guess_value + w_1*self.res_2               
                    
                # Update values
                self.res_1 = interface_residual
                self.w_0 = w_1
                
        elif self.acceleration_type == "Fixed":
            val_correct = iteration_guess_value + self.w_0*interface_residual
            
        elif self.acceleration_type == "Full":
            val_correct = iteration_guess_value + interface_residual
            
        return val_correct
        
        
    def ExecuteFinalizeSolutionStep(self):
        pass
    
    
