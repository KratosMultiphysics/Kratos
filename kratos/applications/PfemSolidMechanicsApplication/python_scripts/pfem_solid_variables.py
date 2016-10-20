from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication     as KratosSolid
import KratosMultiphysics.PfemBaseApplication           as KratosPfemBase
import KratosMultiphysics.ContactMechanicsApplication   as KratosContact
import KratosMultiphysics.PfemSolidMechanicsApplication as KratosPfemSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()


def AddVariables(self, main_model_part):

    main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL);
    main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H);
    
    main_model_part.AddNodalSolutionStepVariable(KratosPfemBase.OFFSET);
    main_model_part.AddNodalSolutionStepVariable(KratosPfemBase.SHRINK_FACTOR);
    main_model_part.AddNodalSolutionStepVariable(KratosPfemBase.MEAN_ERROR);
    main_model_part.AddNodalSolutionStepVariable(KratosPfemBase.RIGID_WALL);
    
    main_model_part.AddNodalSolutionStepVariable(KratosSolid.DETERMINANT_F);
    
    main_model_part.AddNodalSolutionStepVariable(KratosPfemSolid.WALL_TIP_RADIUS);
    main_model_part.AddNodalSolutionStepVariable(KratosPfemSolid.WALL_REFERENCE_POINT);
    
    main_model_part.AddNodalSolutionStepVariable(KratosContact.CONTACT_STRESS);
    
            
    print("::[Pfem Solid Vars]:: Variables ADDED")

