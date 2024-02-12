//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Danilo Cavalcanti and Ignasi de Pouplana
//

#if !defined (KRATOS_ELASTOPLASTIC_MODIFIED_MOHRCOULOMB_COHESIVE_2D_LAW_H_INCLUDED)
#define  KRATOS_ELASTOPLASTIC_MODIFIED_MOHRCOULOMB_COHESIVE_2D_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/interface_element_laws/elastoplastic_mod_mohr_coulomb_cohesive_3D_law.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) ElastoPlasticModMohrCoulombCohesive2DLaw : public ElastoPlasticModMohrCoulombCohesive3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ElastoPlasticModMohrCoulombCohesive2DLaw);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    ElastoPlasticModMohrCoulombCohesive2DLaw()
    {
    }

    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<ElastoPlasticModMohrCoulombCohesive2DLaw>(ElastoPlasticModMohrCoulombCohesive2DLaw(*this));
    }

    // Copy Constructor
    ElastoPlasticModMohrCoulombCohesive2DLaw (const ElastoPlasticModMohrCoulombCohesive2DLaw& rOther) : ElastoPlasticModMohrCoulombCohesive3DLaw(rOther)
    {
    }

    // Destructor
    ~ElastoPlasticModMohrCoulombCohesive2DLaw() override
    {
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void GetLawFeatures(Features& rFeatures) override;

    void InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues ) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    double GetShearResultantStressVector(Vector& StressVector) override;

    void GetElasticConstitutiveMatrix(Matrix& rElasticConstitutiveMatrix, ConstitutiveLawVariables& rVariables, Parameters& rValues) override;

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

}; // Class ElastoPlasticModMohrCoulombCohesive2DLaw
}  // namespace Kratos.
#endif // KRATOS_ELASTOPLASTIC_MODIFIED_MOHRCOULOMB_COHESIVE_2D_LAW_H_INCLUDED  defined
