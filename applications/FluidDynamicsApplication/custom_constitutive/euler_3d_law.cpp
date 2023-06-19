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
#include "custom_constitutive/euler_3d_law.h"

#include "fluid_dynamics_application_variables.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

    //******************************CONSTRUCTOR*******************************************
    //************************************************************************************

    Euler3DLaw::Euler3DLaw()
        : FluidConstitutiveLaw() {}

    //******************************COPY CONSTRUCTOR**************************************
    //************************************************************************************

    Euler3DLaw::Euler3DLaw(const Euler3DLaw& rOther)
        : FluidConstitutiveLaw(rOther) {}

    //********************************CLONE***********************************************
    //************************************************************************************

    ConstitutiveLaw::Pointer Euler3DLaw::Clone() const {
        return Kratos::make_shared<Euler3DLaw>(*this);
    }

    //*******************************DESTRUCTOR*******************************************
    //************************************************************************************

    Euler3DLaw::~Euler3DLaw() {}

    ConstitutiveLaw::SizeType Euler3DLaw::WorkingSpaceDimension() {
        return 3;
    }

    ConstitutiveLaw::SizeType Euler3DLaw::GetStrainSize() const {
        return 6;
    }

    void  Euler3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues) {
        // Get values to compute the constitutive law:
        Flags &Options = rValues.GetOptions();

        Vector& StressVector = rValues.GetStressVector();

        //computation of stress
        StressVector = ZeroVector(6);

        if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
            Matrix& C = rValues.GetConstitutiveMatrix();
            noalias(C) = ZeroMatrix(6,6);
        }

    }


    int Euler3DLaw::Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo) const
    {

        return 0;
    }

    double Euler3DLaw::GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const {
        return 0.0;
    }

    void Euler3DLaw::save(Serializer& rSerializer) const {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FluidConstitutiveLaw )
    }

    void Euler3DLaw::load(Serializer& rSerializer) {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FluidConstitutiveLaw )
    }


} // Namespace Kratos
