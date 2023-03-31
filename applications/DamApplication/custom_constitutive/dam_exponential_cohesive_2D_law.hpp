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
#include "custom_constitutive/exponential_cohesive_2D_law.hpp"
#include "dam_application_variables.h"

namespace Kratos
{

    class KRATOS_API(DAM_APPLICATION) DamExponentialCohesive2DLaw : public ExponentialCohesive2DLaw
    {

    public:

        /// Definition of the base class
        typedef DamExponentialCohesive2DLaw BaseType;

        KRATOS_CLASS_POINTER_DEFINITION(DamExponentialCohesive2DLaw);

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        // Default Constructor
        DamExponentialCohesive2DLaw()
        {
        }

        ConstitutiveLaw::Pointer Clone() const override
        {
            return Kratos::make_shared<DamExponentialCohesive2DLaw>(DamExponentialCohesive2DLaw(*this));
        }

        // Copy Constructor
        DamExponentialCohesive2DLaw (const DamExponentialCohesive2DLaw& rOther) : ExponentialCohesive2DLaw(rOther)
        {
        }

        // Destructor
        ~DamExponentialCohesive2DLaw() override
        {
        }

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
            KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ExponentialCohesive2DLaw )
        }

        void load(Serializer& rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ExponentialCohesive2DLaw )
        }

    }; // Class DamExponentialCohesive2DLaw
}  // namespace Kratos.
