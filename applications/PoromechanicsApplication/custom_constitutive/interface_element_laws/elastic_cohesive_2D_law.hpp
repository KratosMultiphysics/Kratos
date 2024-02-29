//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana and Danilo Cavalcanti
//

#if !defined (KRATOS_ELASTIC_COHESIVE_2D_LAW_H_INCLUDED)
#define  KRATOS_ELASTIC_COHESIVE_2D_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/interface_element_laws/elastic_cohesive_3D_law.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) ElasticCohesive2DLaw : public ElasticCohesive3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ElasticCohesive2DLaw);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    ElasticCohesive2DLaw()
    {
    }

    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<ElasticCohesive2DLaw>(ElasticCohesive2DLaw(*this));
    }

    // Copy Constructor
    ElasticCohesive2DLaw (const ElasticCohesive2DLaw& rOther) : ElasticCohesive3DLaw(rOther)
    {
    }

    // Destructor
    ~ElasticCohesive2DLaw() override
    {
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void GetLawFeatures(Features& rFeatures) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

}; // Class ElasticCohesive2DLaw
}  // namespace Kratos.
#endif // KRATOS_ELASTIC_COHESIVE_2D_LAW_H_INCLUDED  defined
