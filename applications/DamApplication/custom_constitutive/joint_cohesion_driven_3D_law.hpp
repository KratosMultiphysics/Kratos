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

#if !defined (KRATOS_JOINT_COHESION_DRIVEN_3D_LAW_H_INCLUDED)
#define  KRATOS_JOINT_COHESION_DRIVEN_3D_LAW_H_INCLUDED

// System includes

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/bilinear_cohesive_3D_law.hpp"
#include "dam_application_variables.h"

namespace Kratos
{

class KRATOS_API(DAM_APPLICATION) JointCohesionDriven3DLaw : public BilinearCohesive3DLaw
{

public:

    /// Definition of the base class
    typedef BilinearCohesive3DLaw BaseType;

    KRATOS_CLASS_POINTER_DEFINITION(JointCohesionDriven3DLaw);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    JointCohesionDriven3DLaw()
    {
    }

    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<JointCohesionDriven3DLaw>(JointCohesionDriven3DLaw(*this));
    }

    // Copy Constructor
    JointCohesionDriven3DLaw (const JointCohesionDriven3DLaw& rOther) : BilinearCohesive3DLaw(rOther)
    {
    }

    // Destructor
    ~JointCohesionDriven3DLaw() override
    {
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) const override;

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

}; // Class JointCohesionDriven3DLaw
}  // namespace Kratos.
#endif // KRATOS_JOINT_COHESION_DRIVEN_3D_LAW_H_INCLUDED  defined
