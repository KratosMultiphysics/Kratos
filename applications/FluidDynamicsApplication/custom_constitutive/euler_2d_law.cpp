//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/euler_2d_law.h"

#include "fluid_dynamics_application_variables.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

    //******************************CONSTRUCTOR*******************************************
    //************************************************************************************

    Euler2DLaw::Euler2DLaw()
        : ConstitutiveLaw() {}

    //******************************COPY CONSTRUCTOR**************************************
    //************************************************************************************

    Euler2DLaw::Euler2DLaw(const Euler2DLaw& rOther)
        : ConstitutiveLaw(rOther) {}

    //********************************CLONE***********************************************
    //************************************************************************************

    ConstitutiveLaw::Pointer Euler2DLaw::Clone() const {
        Euler2DLaw::Pointer p_clone(new Euler2DLaw(*this));
        return p_clone;
    }

    //*******************************DESTRUCTOR*******************************************
    //************************************************************************************

    Euler2DLaw::~Euler2DLaw() {}

    void  Euler2DLaw::CalculateMaterialResponseCauchy (Parameters& rValues) {
        // Get values to compute the constitutive law:
        Flags &Options = rValues.GetOptions();

        Vector& StressVector = rValues.GetStressVector();

        //computation of stress
        StressVector = ZeroVector(3);

        if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
            Matrix& C = rValues.GetConstitutiveMatrix();
            noalias(C) = ZeroMatrix(3,3);
        }

    }


    //*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
    //************************************************************************************

    void Euler2DLaw::GetLawFeatures(Features& rFeatures) {
        //Set the type of law
        rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
        rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
        rFeatures.mOptions.Set( ISOTROPIC );

        //Set strain measure required by the consitutive law
        rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
        rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

        //Set the strain size
        rFeatures.mStrainSize = 3;

        //Set the spacedimension
        rFeatures.mSpaceDimension = 2;

    }

} // Namespace Kratos
