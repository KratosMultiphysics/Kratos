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

#if !defined (KRATOS_ELASTOPLASTIC_MOHRCOULOMB_COHESIVE_2D_LAW_H_INCLUDED)
#define  KRATOS_ELASTOPLASTIC_MOHRCOULOMB_COHESIVE_2D_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/interface_element_laws/elastoplastic_mohr_coulomb_cohesive_3D_law.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) ElastoPlasticMohrCoulombCohesive2DLaw : public ElastoPlasticMohrCoulombCohesive3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ElastoPlasticMohrCoulombCohesive2DLaw);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    ElastoPlasticMohrCoulombCohesive2DLaw()
    {
    }

    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<ElastoPlasticMohrCoulombCohesive2DLaw>(ElastoPlasticMohrCoulombCohesive2DLaw(*this));
    }

    // Copy Constructor
    ElastoPlasticMohrCoulombCohesive2DLaw (const ElastoPlasticMohrCoulombCohesive2DLaw& rOther) : ElastoPlasticMohrCoulombCohesive3DLaw(rOther)
    {
    }

    // Destructor
    ~ElastoPlasticMohrCoulombCohesive2DLaw() override
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

    void StressVectorInstersectionYieldSurfaces(Vector& rStressVector, const double ts, const double ts_intersection, const double ft) override;

    void ConstitutiveMatrixInstersectionYieldSurfaces(Vector& StressVector, Matrix& rConstitutiveMatrix, ConstitutiveLawVariables& rVariables) override;
    
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

}; // Class ElastoPlasticMohrCoulombCohesive2DLaw
}  // namespace Kratos.
#endif // KRATOS_ELASTOPLASTIC_MOHRCOULOMB_COHESIVE_2D_LAW_H_INCLUDED  defined
