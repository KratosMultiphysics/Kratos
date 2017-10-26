//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher 
//
	        

// System includes


// External includes 


// Project includes
#include "shell_utilities.hpp"


namespace Kratos
{
    typedef Properties PropertiesType;

    void ShellUtilities::CheckVariables()
    {
        // verify that the variables are correctly initialized
        if(DISPLACEMENT.Key() == 0)
            KRATOS_ERROR << "DISPLACEMENT has Key zero! (check if the application is correctly registered" << std::endl;
        if(ROTATION.Key() == 0)
            KRATOS_ERROR << "ROTATION has Key zero! (check if the application is correctly registered" << std::endl;
        if(VELOCITY.Key() == 0)
            KRATOS_ERROR << "VELOCITY has Key zero! (check if the application is correctly registered" << std::endl;
        if(ACCELERATION.Key() == 0)
            KRATOS_ERROR << "ACCELERATION has Key zero! (check if the application is correctly registered" << std::endl;
        if(DENSITY.Key() == 0)
            KRATOS_ERROR << "DENSITY has Key zero! (check if the application is correctly registered" << std::endl;
        if(SHELL_CROSS_SECTION.Key() == 0)
            KRATOS_ERROR << "SHELL_CROSS_SECTION has Key zero! (check if the application is correctly registered" << std::endl;
        if(THICKNESS.Key() == 0)
            KRATOS_ERROR << "THICKNESS has Key zero! (check if the application is correctly registered" << std::endl;
        if(CONSTITUTIVE_LAW.Key() == 0)
            KRATOS_ERROR << "CONSTITUTIVE_LAW has Key zero! (check if the application is correctly registered" << std::endl;
    }

    void ShellUtilities::CheckDofs(GeometryType& rGeom)
    {
        // verify that the dofs exist
        for (unsigned int i = 0; i < rGeom.size(); i++)
        {
            if (rGeom[i].SolutionStepsDataHas(DISPLACEMENT) == false)
                KRATOS_ERROR << "missing variable DISPLACEMENT on node " << rGeom[i].Id() << std::endl;

            if (rGeom[i].HasDofFor(DISPLACEMENT_X) == false ||
                rGeom[i].HasDofFor(DISPLACEMENT_Y) == false ||
                rGeom[i].HasDofFor(DISPLACEMENT_Z) == false)
                KRATOS_ERROR << "missing one of the dofs for the variable DISPLACEMENT on node " << rGeom[i].Id() << std::endl;

            if (rGeom[i].SolutionStepsDataHas(ROTATION) == false)
                KRATOS_ERROR << "missing variable ROTATION on node " << rGeom[i].Id() << std::endl;

            if (rGeom[i].HasDofFor(ROTATION_X) == false ||
                rGeom[i].HasDofFor(ROTATION_Y) == false ||
                rGeom[i].HasDofFor(ROTATION_Z) == false)
                KRATOS_ERROR << "missing one of the dofs for the variable ROTATION on node " << rGeom[i].Id() << std::endl;

            if (rGeom[i].GetBufferSize() < 2)
                KRATOS_ERROR << "This Element needs at least a buffer size = 2" << std::endl;
        }
    }

    void ShellUtilities::CheckProperties(const Element* pTheElement, const ProcessInfo& rCurrentProcessInfo, const bool IsThickShell)
    {
        // check properties
        if(pTheElement->pGetProperties() == nullptr)
            KRATOS_ERROR << "Properties not provided for element " << pTheElement->Id() << std::endl;

        const PropertiesType & props = pTheElement->GetProperties();

        const GeometryType& geom = pTheElement->GetGeometry(); // TODO check if this can be const

        if(props.Has(SHELL_CROSS_SECTION)) // if the user specified a cross section ...
        {
            const ShellCrossSection::Pointer & section = props[SHELL_CROSS_SECTION];
            if(section == nullptr)
                KRATOS_ERROR << "SHELL_CROSS_SECTION not provided for element " << pTheElement->Id() << std::endl;
    
            section->Check(props, geom, rCurrentProcessInfo);
        }
        else if (props.Has(SHELL_ORTHOTROPIC_LAYERS))
        {
            CheckSpecificProperties(pTheElement, props, IsThickShell);
    
            // perform detailed orthotropic check later in shell_cross_section
        }
        else // ... allow the automatic creation of a homogeneous section from a material and a thickness
        {
            CheckSpecificProperties(pTheElement, props, IsThickShell);
    
            ShellCrossSection::Pointer dummySection = ShellCrossSection::Pointer(new ShellCrossSection());
            dummySection->BeginStack();
            dummySection->AddPly(props[THICKNESS], 0.0, 5, pTheElement->pGetProperties());
            dummySection->EndStack();
            dummySection->SetSectionBehavior(ShellCrossSection::Thick);
            dummySection->Check(props, geom, rCurrentProcessInfo);
        }

    }

    void ShellUtilities::CheckSpecificProperties(const Element* pTheElement, const PropertiesType & rProps, const bool IsThickShell)
    {
        if (!rProps.Has(CONSTITUTIVE_LAW))
        KRATOS_ERROR << "CONSTITUTIVE_LAW not provided for element " << pTheElement->Id() << std::endl;
        const ConstitutiveLaw::Pointer& claw = rProps[CONSTITUTIVE_LAW];
        if (claw == nullptr)
            KRATOS_ERROR << "CONSTITUTIVE_LAW not provided for element " << pTheElement->Id() << std::endl;

        if(!rProps.Has(THICKNESS))
            KRATOS_ERROR << "THICKNESS not provided for element " << pTheElement->Id() << std::endl;
        if(rProps[THICKNESS] <= 0.0)
            KRATOS_ERROR << "wrong THICKNESS value provided for element " << pTheElement->Id() << std::endl;
            
        if(!rProps.Has(DENSITY))
            KRATOS_ERROR << "DENSITY not provided for element " << pTheElement->Id() << std::endl;
        if(rProps[DENSITY] < 0.0)
            KRATOS_ERROR << "wrong DENSITY value provided for element " << pTheElement->Id() << std::endl;
        
        if(IsThickShell)
        {
            //Check constitutive law has been verified with Stenberg stabilization
            //applicable for 5-parameter shells only.
            ConstitutiveLaw::Features pclFeatures;
            claw->GetLawFeatures(pclFeatures);
            if (pclFeatures.mOptions.IsNot(ConstitutiveLaw::STENBERG_STABILIZATION_SUITABLE))
            {
                std::cout << "\nWARNING:\nThe current constitutive law has not been checked with Stenberg shear stabilization."
                    << "\nPlease check results carefully."
                    << std::endl;
            }
        }
    }  
  
}  // namespace Kratos.


