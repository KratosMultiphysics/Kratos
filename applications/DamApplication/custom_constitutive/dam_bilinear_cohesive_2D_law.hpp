//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Joaquín Irazábal González
//

#pragma once

// System includes

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/bilinear_cohesive_2D_law.hpp"
#include "dam_application_variables.h"

namespace Kratos
{

    class KRATOS_API(DAM_APPLICATION) DamBilinearCohesive2DLaw : public BilinearCohesive2DLaw
    {

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DamBilinearCohesive2DLaw);

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        // Default Constructor
        DamBilinearCohesive2DLaw() : BilinearCohesive2DLaw() {}

        ConstitutiveLaw::Pointer Clone() const override
        {
            return Kratos::make_shared<DamBilinearCohesive2DLaw>(*this);
        }

        // Copy Constructor
        DamBilinearCohesive2DLaw (DamBilinearCohesive2DLaw const& rOther) : BilinearCohesive2DLaw(rOther) {}

        // Destructor
        ~DamBilinearCohesive2DLaw() override {}

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        double& GetValue( const Variable<double>& rThisVariable, double& rValue ) override;

        void SetValue( const Variable<double>& rVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo ) override;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    protected:

        // Member Variables
        double mUpliftPressure;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        void ComputeStressVector(Vector& rStressVector,
                                 ConstitutiveLawVariables& rVariables,
                                 Parameters& rValues) override;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    private:

        // Serialization

        friend class Serializer;

        void save(Serializer& rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BilinearCohesive2DLaw )
        }

        void load(Serializer& rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BilinearCohesive2DLaw )
        }

    }; // Class DamBilinearCohesive2DLaw
}  // namespace Kratos.