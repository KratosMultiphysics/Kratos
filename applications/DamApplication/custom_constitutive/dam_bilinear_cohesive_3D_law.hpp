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
#include "custom_constitutive/bilinear_cohesive_3D_law.hpp"
#include "dam_application_variables.h"

namespace Kratos
{

    class KRATOS_API(DAM_APPLICATION) DamBilinearCohesive3DLaw : public BilinearCohesive3DLaw
    {

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DamBilinearCohesive3DLaw);

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        // Default Constructor
        DamBilinearCohesive3DLaw() : BilinearCohesive3DLaw() {}

        ConstitutiveLaw::Pointer Clone() const override
        {
            return Kratos::make_shared<DamBilinearCohesive3DLaw>(*this);
        }

        // Copy Constructor
        DamBilinearCohesive3DLaw (DamBilinearCohesive3DLaw const& rOther) : BilinearCohesive3DLaw(rOther) {}

        // Destructor
        ~DamBilinearCohesive3DLaw() override {}

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
            KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BilinearCohesive3DLaw )
        }

        void load(Serializer& rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BilinearCohesive3DLaw )
        }

    }; // Class DamBilinearCohesive3DLaw
}  // namespace Kratos.