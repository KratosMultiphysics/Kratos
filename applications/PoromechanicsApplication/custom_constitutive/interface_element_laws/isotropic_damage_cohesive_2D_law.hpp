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

#if !defined (KRATOS_ISOTROPIC_DAMAGE_COHESIVE_2D_LAW_H_INCLUDED)
#define  KRATOS_ISOTROPIC_DAMAGE_COHESIVE_2D_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/interface_element_laws/isotropic_damage_cohesive_3D_law.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) IsotropicDamageCohesive2DLaw : public IsotropicDamageCohesive3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(IsotropicDamageCohesive2DLaw);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    IsotropicDamageCohesive2DLaw()
    {
    }

    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<IsotropicDamageCohesive2DLaw>(IsotropicDamageCohesive2DLaw(*this));
    }

    // Copy Constructor
    IsotropicDamageCohesive2DLaw (const IsotropicDamageCohesive2DLaw& rOther) : IsotropicDamageCohesive3DLaw(rOther)
    {
    }

    // Destructor
    ~IsotropicDamageCohesive2DLaw() override
    {
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void GetLawFeatures(Features& rFeatures) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ComputeEquivalentStrain(ConstitutiveLawVariables& rVariables, Parameters& rValues) override;

    void GetElasticConstitutiveMatrix(Matrix& rElasticConstitutiveMatrix,
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

}; // Class IsotropicDamageCohesive2DLaw
}  // namespace Kratos.
#endif // KRATOS_ISOTROPIC_DAMAGE_COHESIVE_2D_LAW_H_INCLUDED  defined
