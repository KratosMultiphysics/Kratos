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

#if !defined (KRATOS_JOINT_STRESS_DRIVEN_2D_LAW_H_INCLUDED)
#define  KRATOS_JOINT_STRESS_DRIVEN_2D_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/joint_stress_driven_3D_law.hpp"
#include "dam_application_variables.h"

namespace Kratos
{

class KRATOS_API(DAM_APPLICATION) JointStressDriven2DLaw : public JointStressDriven3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(JointStressDriven2DLaw);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    JointStressDriven2DLaw()
    {
    }

    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<JointStressDriven2DLaw>(JointStressDriven2DLaw(*this));
    }

    // Copy Constructor
    JointStressDriven2DLaw (const JointStressDriven2DLaw& rOther) : JointStressDriven3DLaw(rOther)
    {
    }

    // Destructor
    ~JointStressDriven2DLaw() override
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

}; // Class JointStressDriven2DLaw
}  // namespace Kratos.
#endif // KRATOS_JOINT_STRESS_DRIVEN_2D_LAW_H_INCLUDED  defined
