## Dictionary of the constant parameters used in the Variational Formulation
from sympy import *

params = {
            "dim": 2,			                                    # Dimension
            "mu": Symbol('mu', positive = True),			        # Dynamic viscosity 
            "h" : Symbol('h', positive = True),	        	        # Element size
            "lambda" : Symbol('lambda', positive = True),           # Thermal Conductivity of the fluid
            "c_v" : Symbol('c_v', positive = True),			        # Specific Heat at Constant volume
            "gamma": Symbol('gamma',positive = True),			    # Gamma (Cp/Cv) 
            "stab_c1" : Symbol('stab_c1', positive = True),			# Algorithm constant
            "stab_c2" : Symbol('stab_c2', positive = True),			# Algorithm constant
        
	}

'''
## 2D test AIR 20º
params = {
            "dim": 2,			                            # Dimension
            "mu": 1.846*0.00005,			# Dynamic viscosity 
            "nu" :1.568*0.00005,			# Kinematic viscosity 
            "h" : 1,	        	# Element size
            "lambda" : 0.0257,         # Thermal Conductivity of the fluid
            "c_v" : 718,			# Specific Heat at Constant volume
            "c_p" : 1005,			# Specific Heat at Constant Pressure
            "gamma": 1.4,			# Gamma (Cp/Cv) 
            "stab_c1" : 4,			# Algorithm constant
            "stab_c2" : 2,			# Algorithm constant
            
	}

'''
