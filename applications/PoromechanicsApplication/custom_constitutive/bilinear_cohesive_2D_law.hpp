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

#if !defined (KRATOS_BILINEAR_COHESIVE_2D_LAW_H_INCLUDED)
#define  KRATOS_BILINEAR_COHESIVE_2D_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/bilinear_cohesive_3D_law.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) BilinearCohesive2DLaw : public BilinearCohesive3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(BilinearCohesive2DLaw);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    BilinearCohesive2DLaw()
    {
    }

    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<BilinearCohesive2DLaw>(BilinearCohesive2DLaw(*this));
    }

    // Copy Constructor
    BilinearCohesive2DLaw (const BilinearCohesive2DLaw& rOther) : BilinearCohesive3DLaw(rOther)
    {
    }

    // Destructor
    ~BilinearCohesive2DLaw() override
    {
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void GetLawFeatures(Features& rFeatures) override;

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

}; // Class BilinearCohesive2DLaw
}  // namespace Kratos.
#endif // KRATOS_BILINEAR_COHESIVE_2D_LAW_H_INCLUDED  defined
