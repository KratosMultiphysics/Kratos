from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication  as KratosDelaunay


def AddVariables(main_model_part):
  
    main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL);
    main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H);

    main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONTACT_FORCE)
    main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONTACT_NORMAL);
    
    main_model_part.AddNodalSolutionStepVariable(KratosDelaunay.OFFSET);
    main_model_part.AddNodalSolutionStepVariable(KratosDelaunay.SHRINK_FACTOR);
    main_model_part.AddNodalSolutionStepVariable(KratosDelaunay.MEAN_ERROR);
    main_model_part.AddNodalSolutionStepVariable(KratosDelaunay.RIGID_WALL);
            
    print("::[Pfem Extra Vars]:: Variables ADDED")

