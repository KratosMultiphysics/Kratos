# ==============================================================================
'''
 KratosShapeOptimizationApplication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 (Released on march 05, 2007).

 Copyright (c) 2016: Daniel Baumgaertner
                     daniel.baumgaertner@tum.de
                     Chair of Structural Analysis
                     Technische Universitaet Muenchen
                     Arcisstrasse 21 80333 Munich, Germany

 Permission is hereby granted, free  of charge, to any person obtaining
 a  copy  of this  software  and  associated  documentation files  (the
 "Software"), to  deal in  the Software without  restriction, including
 without limitation  the rights to  use, copy, modify,  merge, publish,
 distribute,  sublicense and/or  sell copies  of the  Software,  and to
 permit persons to whom the Software  is furnished to do so, subject to
 the following condition:

 Distribution of this code for  any  commercial purpose  is permissible
 ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

 The  above  copyright  notice  and  this permission  notice  shall  be
 included in all copies or substantial portions of the Software.

 THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
 EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
 CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
 TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''
#==============================================================================
#
#   Project Name:        KratosShape                            $
#   Created by:          $Author:    daniel.baumgaertner@tum.de $
#                        $Author:           armin.geiser@tum.de $
#   Date:                $Date:                   December 2016 $
#   Revision:            $Revision:                         0.0 $
#
# ==============================================================================

# ------------------------------------------------------------------------------
# Imports
# ------------------------------------------------------------------------------
# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# ==============================================================================
def CreateSolver(model_part, opt_settings):

    # Dictionary to store solvers of all response functions defined in the optimization settings
    solver = {}

    # Collect all responses
    specified_responses = {}
    for response_id in opt_settings.objectives:
        specified_responses[response_id] = opt_settings.objectives[response_id]

    if not specified_responses:
        raise ValueError("No objective function specified!")

    for response_id in opt_settings.constraints:
        specified_responses[response_id] = opt_settings.constraints[response_id]

    # Creat response function solver according to specified settings and add relevant variables
    if "strain_energy" in specified_responses.keys():
        model_part.AddNodalSolutionStepVariable(STRAIN_ENERGY_SHAPE_GRADIENT)
        solver["strain_energy"] = StrainEnergyResponseFunction(model_part, specified_responses["strain_energy"])
    else:
        raise ValueError("Specified response function not implemented. Implemented response functions are: \"strain_energy\"")

    return solver

# ==============================================================================
