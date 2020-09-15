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

#if !defined (KRATOS_JOINT_COHESION_DRIVEN_2D_LAW_H_INCLUDED)
#define  KRATOS_JOINT_COHESION_DRIVEN_2D_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/joint_cohesion_driven_3D_law.hpp"
#include "dam_application_variables.h"

namespace Kratos
{

class KRATOS_API(DAM_APPLICATION) JointCohesionDriven2DLaw : public JointCohesionDriven3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(JointCohesionDriven2DLaw);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    JointCohesionDriven2DLaw()
    {
    }

    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<JointCohesionDriven2DLaw>(JointCohesionDriven2DLaw(*this));
    }

    // Copy Constructor
    JointCohesionDriven2DLaw (const JointCohesionDriven2DLaw& rOther) : JointCohesionDriven3DLaw(rOther)
    {
    }

    // Destructor
    ~JointCohesionDriven2DLaw() override
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

}; // Class JointCohesionDriven2DLaw
}  // namespace Kratos.
#endif // KRATOS_JOINT_COHESION_DRIVEN_2D_LAW_H_INCLUDED  defined
