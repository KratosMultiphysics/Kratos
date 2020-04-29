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

#if !defined (KRATOS_SIMPLIFIED_BILINEAR_2D_LAW_H_INCLUDED)
#define  KRATOS_SIMPLIFIED_BILINEAR_2D_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/simplified_bilinear_3D_law.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) SimplifiedBilinear2DLaw : public SimplifiedBilinear3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(SimplifiedBilinear2DLaw);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    SimplifiedBilinear2DLaw()
    {
    }

    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<SimplifiedBilinear2DLaw>(SimplifiedBilinear2DLaw(*this));
    }

    // Copy Constructor
    SimplifiedBilinear2DLaw (const SimplifiedBilinear2DLaw& rOther) : SimplifiedBilinear3DLaw(rOther)
    {
    }

    // Destructor
    ~SimplifiedBilinear2DLaw() override
    {
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ComputeEquivalentStrain(ConstitutiveLawVariables& rVariables, Parameters& rValues) override;

    void ComputeConstitutiveMatrix(Matrix& rConstitutiveMatrix,
                                   ConstitutiveLawVariables& rVariables,
                                   Parameters& rValues) override;

    void ComputeStressVector(Vector& rStressVector,
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

}; // Class SimplifiedBilinear2DLaw
}  // namespace Kratos.
#endif // KRATOS_SIMPLIFIED_BILINEAR_2D_LAW_H_INCLUDED  defined
