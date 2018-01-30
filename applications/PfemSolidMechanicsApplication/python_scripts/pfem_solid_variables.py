from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication     as KratosSolid
import KratosMultiphysics.PfemApplication           as KratosPfem
import KratosMultiphysics.ContactMechanicsApplication   as KratosContact
import KratosMultiphysics.PfemSolidMechanicsApplication as KratosPfemSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()


def AddVariables(main_model_part):
  
    main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL);
    main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H);

    main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONTACT_FORCE)
    main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONTACT_NORMAL);
    
    main_model_part.AddNodalSolutionStepVariable(KratosPfem.OFFSET);
    main_model_part.AddNodalSolutionStepVariable(KratosPfem.SHRINK_FACTOR);
    main_model_part.AddNodalSolutionStepVariable(KratosPfem.MEAN_ERROR);
    main_model_part.AddNodalSolutionStepVariable(KratosPfem.RIGID_WALL);
       
    main_model_part.AddNodalSolutionStepVariable(KratosContact.CONTACT_STRESS);
    main_model_part.AddNodalSolutionStepVariable(KratosContact.EFFECTIVE_CONTACT_STRESS);
    main_model_part.AddNodalSolutionStepVariable(KratosContact.EFFECTIVE_CONTACT_FORCE);
    
    #transfer variables
    main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DETERMINANT_F);
    main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CAUCHY_STRESS_VECTOR)
    main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DEFORMATION_GRADIENT)
            
    print("::[Pfem Solid Vars]:: Variables ADDED")

