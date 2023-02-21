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
        : FluidConstitutiveLaw() {}

    //******************************COPY CONSTRUCTOR**************************************
    //************************************************************************************

    Euler2DLaw::Euler2DLaw(const Euler2DLaw& rOther)
        : FluidConstitutiveLaw(rOther) {}

    //********************************CLONE***********************************************
    //************************************************************************************

    ConstitutiveLaw::Pointer Euler2DLaw::Clone() const {
        return Kratos::make_shared<Euler2DLaw>(*this);
    }

    //*******************************DESTRUCTOR*******************************************
    //************************************************************************************

    Euler2DLaw::~Euler2DLaw() {}

    ConstitutiveLaw::SizeType Euler2DLaw::WorkingSpaceDimension() {
        return 2;
    }

    ConstitutiveLaw::SizeType Euler2DLaw::GetStrainSize() const {
        return 3;
    }

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


    int Euler2DLaw::Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo) const
    {

        return 0;
    }

    double Euler2DLaw::GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const {
        return 0.0;
    }

    void Euler2DLaw::save(Serializer& rSerializer) const {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FluidConstitutiveLaw )
    }

    void Euler2DLaw::load(Serializer& rSerializer) {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FluidConstitutiveLaw )
    }


} // Namespace Kratos
