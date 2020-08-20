//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Javier San Mauro Saiz
//                   Joaquin Irazabal Gonzalez
//

#if !defined (KRATOS_SIMPLIFIED_BILINEAR_3D_LAW_H_INCLUDED)
#define  KRATOS_SIMPLIFIED_BILINEAR_3D_LAW_H_INCLUDED

// System includes

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/bilinear_cohesive_3D_law.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) SimplifiedBilinear3DLaw : public BilinearCohesive3DLaw
{

public:

    /// Definition of the base class
    typedef BilinearCohesive3DLaw BaseType;

    KRATOS_CLASS_POINTER_DEFINITION(SimplifiedBilinear3DLaw);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    SimplifiedBilinear3DLaw()
    {
    }

    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<SimplifiedBilinear3DLaw>(SimplifiedBilinear3DLaw(*this));
    }

    // Copy Constructor
    SimplifiedBilinear3DLaw (const SimplifiedBilinear3DLaw& rOther) : BilinearCohesive3DLaw(rOther)
    {
    }

    // Destructor
    ~SimplifiedBilinear3DLaw() override
    {
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) override;

    void InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues ) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void FinalizeMaterialResponseCauchy (Parameters & rValues) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeConstitutiveLawVariables(ConstitutiveLawVariables& rVariables, Parameters& rValues) override;

    void ComputeEquivalentStrain(ConstitutiveLawVariables& rVariables, Parameters& rValues) override;

    void CheckLoadingFunction(ConstitutiveLawVariables& rVariables, Parameters& rValues) override;

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

}; // Class SimplifiedBilinear3DLaw
}  // namespace Kratos.
#endif // KRATOS_SIMPLIFIED_BILINEAR_3D_LAW_H_INCLUDED  defined
