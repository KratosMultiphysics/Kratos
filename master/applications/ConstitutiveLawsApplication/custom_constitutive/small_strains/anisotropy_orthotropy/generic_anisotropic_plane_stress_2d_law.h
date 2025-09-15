// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "generic_anisotropic_law.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// The size type definition
using SizeType = std::size_t;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/*
 * @class GenericAnisotropicLaw
 * @ingroup ConstitutiveLawsApplication
 * @brief This CL takes into account the material anisotropy in terms of young modulus, poisson ratio, orientation and strengths. 3D CLs and plane strain.
 * @details See "Nonlinear behavior of existing pre-tensioned concrete beams: Experimental study and finite element modeling with the constitutive Serial-Parallel rule of mixtures", DOI: https://doi.org/10.1016/j.istruc.2024.106990
 * @author Alejandro Cornejo
 */
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) GenericAnisotropicPlaneStress2DLaw
    : public GenericAnisotropicLaw<2>
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = GenericAnisotropicLaw<2>;

    /// Counted pointer of GenericAnisotropicLaw
    KRATOS_CLASS_POINTER_DEFINITION(GenericAnisotropicPlaneStress2DLaw);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * Constructor.
    */
    GenericAnisotropicPlaneStress2DLaw()
    {
    }

    /**
    * Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<GenericAnisotropicPlaneStress2DLaw>(*this);
    }

    // Copy constructor
    GenericAnisotropicPlaneStress2DLaw(GenericAnisotropicPlaneStress2DLaw const& rOther)
        : BaseType(rOther)
    {
    }

    /**
    * Destructor.
    */
    ~GenericAnisotropicPlaneStress2DLaw() override
    {
    }

    /**
     * creates a new constitutive law pointer
     * @param NewParameters The configuration parameters of the new constitutive law
     * @return a Pointer to the new constitutive law
     */
    ConstitutiveLaw::Pointer Create(Kratos::Parameters NewParameters) const override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{


    /**
     * @brief This method computes the orthotropic elastic constitutive matrix
     * This method is overriden in the derived class for plane stress
     */
    void CalculateOrthotropicElasticMatrix(
        BoundedMatrixVoigtType &rElasticityTensor,
        const Properties &rMaterialProperties) override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
    }

    ///@}

}; // Class GenericAnisotropicLaw

} // namespace Kratos
