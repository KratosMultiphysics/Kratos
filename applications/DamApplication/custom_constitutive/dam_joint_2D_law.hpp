//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Joaquín Irazábal González
//                   Ignasi de Pouplana
//
//

#pragma once

// Application includes
#include "custom_constitutive/dam_joint_3D_law.hpp"

namespace Kratos
{

    class KRATOS_API(DAM_APPLICATION) DamJoint2DLaw : public DamJoint3DLaw
    {

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DamJoint2DLaw);

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        // Default Constructor
        DamJoint2DLaw()
        {
        }

        ConstitutiveLaw::Pointer Clone() const override
        {
            return Kratos::make_shared<DamJoint2DLaw>(DamJoint2DLaw(*this));
        }

        // Copy Constructor
        DamJoint2DLaw (const DamJoint2DLaw& rOther) : DamJoint3DLaw(rOther)
        {
        }

        // Destructor
        ~DamJoint2DLaw() override
        {
        }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        void GetLawFeatures(Features& rFeatures) override;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    protected:

        // Member Variables

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        void ComputeEquivalentStrain(ConstitutiveLawVariables& rVariables,
                                     Parameters& rValues) override;

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

    }; // Class DamJoint2DLaw
}  // namespace Kratos.
