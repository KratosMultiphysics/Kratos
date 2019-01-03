//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

#if !defined (KRATOS_EXPONENTIAL_COHESIVE_2D_LAW_H_INCLUDED)
#define  KRATOS_EXPONENTIAL_COHESIVE_2D_LAW_H_INCLUDED

// System includes

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/exponential_cohesive_3D_law.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) ExponentialCohesive2DLaw : public ExponentialCohesive3DLaw
{

public:

    /// Definition of the base class
    typedef ExponentialCohesive3DLaw BaseType;

    KRATOS_CLASS_POINTER_DEFINITION(ExponentialCohesive2DLaw);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    ExponentialCohesive2DLaw()
    {
    }

    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<ExponentialCohesive2DLaw>(ExponentialCohesive2DLaw(*this));
    }

    // Copy Constructor
    ExponentialCohesive2DLaw (const ExponentialCohesive2DLaw& rOther) : ExponentialCohesive3DLaw(rOther)
    {
    }

    // Destructor
    ~ExponentialCohesive2DLaw() override
    {
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeConstitutiveLawVariables(ConstitutiveLawVariables& rVariables, Parameters& rValues) override;

    void ComputeEquivalentStrain(ConstitutiveLawVariables& rVariables, Parameters& rValues) override;

    void ComputeCriticalDisplacement(ConstitutiveLawVariables& rVariables, Parameters& rValues) override;

    void ComputeConstitutiveMatrix(Matrix& rConstitutiveMatrix,
                                    ConstitutiveLawVariables& rVariables,
                                    Parameters& rValues) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

}; // Class ExponentialCohesive2DLaw
}  // namespace Kratos.
#endif // KRATOS_EXPONENTIAL_COHESIVE_2D_LAW_H_INCLUDED  defined
