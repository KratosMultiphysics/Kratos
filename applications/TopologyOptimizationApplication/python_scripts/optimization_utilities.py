# ======================================================================================================================================
# IMPORTS
# ======================================================================================================================================
# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division 

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.TopologyOptimizationApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# Further necessary imports
import math
import time

# ======================================================================================================================================
class OptimizationUtilities:

	# ----------------------------------------------------------------------------------------------------------------------------------
        def __init__(self, opt_model_part, analyzer, controller, config):

            self.opt_model_part = opt_model_part
            self.analyzer = analyzer
            self.controller = controller
            self.config = config

            # active constraints
            self.only_C_id = None

	# ----------------------------------------------------------------------------------------------------------------------------------
        def set_active_constraints(self,only_C_id):

            self.only_C_id = only_C_id

	# ----------------------------------------------------------------------------------------------------------------------------------
        def update_design_oc(self, opt_itr):

            start_time = time.time()
            print("  Optimality Criterion Method (OC) chosen to solve the optimization problem")

            # Check if Grey Scale Filter should be used
            q = 1
            if (self.config.grey_scale_filter == 1):
                if (opt_itr < 15):
                    q = 1
                else:
                    q = min(self.config.q_max, 1.01*q)

                print("  Grey Scale Filter activated, q = ", q)
            else:
                print("  Grey Scale Filter deactivated, q = ", q)

            # Update Densities procedure
            l1      = 0.0
            l2      = 1000000000.0
            move    = 0.2

            while((l2-l1)/(l1+l2) > 0.001):
                lmid       = 0.5*(l2+l1)
                nele       = 0
                x_new      = 0.0
                sum_X_Phys = 0.0

                for element_i in self.opt_model_part.Elements:
                    dcdx       = element_i.GetValue(DCDX)
                    dvdx       = element_i.GetValue(DVDX)
                    x_old      = element_i.GetValue(X_PHYS_OLD)
                    solid_void = element_i.GetValue(SOLID_VOID)

                    # Update Density Method
                    # When q = 1, Grey Scale Filter is not activated, i.e., the results are in the classical OC update method

                    # NORMAL elements
                    if (solid_void==0):
                        x_new = max(0.0, max(x_old - move, min(1.0, pow(min(x_old + move, x_old * math.sqrt(-dcdx/dvdx/lmid)),q) )))

                    # ACTIVE elements (solid elements)
                    elif (solid_void==1):
                        x_new = 1;

                    # PASSIVE elements (void elements)
                    elif (solid_void==2):
                        x_new = 0;

                    # If no element identification was found
                    else:
                        print("This value for SOLID_VOID does not exist.")

                    # Update of the calculated X_PHYS for the next iteration
                    element_i.SetValue(X_PHYS, x_new)

                    # Updating additional quantities to determine the correct Lagrange Multiplier (lmid)
                    sum_X_Phys = sum_X_Phys + x_new
                    nele = nele + 1

                if (sum_X_Phys > (self.config.volume_fraction * nele)):
                    l1 = lmid
                else:
                    l2 = lmid

            end_time = time.time()
            print("  Updating of values performed               [ spent time = ", round(end_time - start_time, 6), "] ")

# ======================================================================================================================================
