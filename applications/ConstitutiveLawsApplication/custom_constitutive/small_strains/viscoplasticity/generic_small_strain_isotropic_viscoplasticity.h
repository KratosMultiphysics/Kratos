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

#pragma once

// System includes

// External includes

// Project includes
#include "custom_constitutive/small_strains/plasticity/generic_small_strain_isotropic_plasticity.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class GenericSmallStrainIsotropicViscoPlasticity
 * @ingroup ConstitutiveLawsApp
 * @brief 
 * @author Alejandro Cornejo
 */
template <class TConstLawIntegratorType>
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) GenericSmallStrainIsotropicViscoPlasticity
    : public GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>
{
public:
    ///@name Type Definitions
    ///@{

    // The size type definition
    using SizeType = std::size_t;

    using ConstLawIntegratorType = TConstLawIntegratorType;

    /// The define the working dimension size, already defined in the integrator
    static constexpr SizeType Dimension = TConstLawIntegratorType::Dimension;

    /// The define the Voigt size, already defined in the  integrator
    static constexpr SizeType VoigtSize = TConstLawIntegratorType::VoigtSize;

    /// Definition of the base class
    using BaseType = GenericSmallStrainIsotropicPlasticity<TConstLawIntegratorType>;

    /// Counted pointer of GenericSmallStrainIsotropicPlasticity
    KRATOS_CLASS_POINTER_DEFINITION(GenericSmallStrainIsotropicViscoPlasticity);

    // /// The node definition
    // typedef Node NodeType;

    // /// The geometry definition
    // typedef Geometry<NodeType> GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * Default constructor.
    */
   GenericSmallStrainIsotropicViscoPlasticity()
    {
    }

    /**
    * Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<GenericSmallStrainIsotropicViscoPlasticity<TConstLawIntegratorType>>(*this);
    }

    /**
    * Copy constructor.
    */
   GenericSmallStrainIsotropicViscoPlasticity(const GenericSmallStrainIsotropicViscoPlasticity &rOther)
        : BaseType(rOther)
    {
    }

    /**
    * Destructor.
    */
    ~GenericSmallStrainIsotropicViscoPlasticity() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{


    /**
     * @brief Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * Finalize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues) override;

    // /**
    //  * @brief Returns whether this constitutive Law has specified variable (double)
    //  * @param rThisVariable the variable to be checked for
    //  * @return true if the variable is defined in the constitutive law
    //  */
    // bool Has(const Variable<double> &rThisVariable) override;

    // /**
    //  * @brief Returns whether this constitutive Law has specified variable (Vector)
    //  * @param rThisVariable the variable to be checked for
    //  * @return true if the variable is defined in the constitutive law
    //  */
    // bool Has(const Variable<Vector> &rThisVariable) override;

    // /**
    //  * @brief Returns whether this constitutive Law has specified variable (Matrix)
    //  * @param rThisVariable the variable to be checked for
    //  * @return true if the variable is defined in the constitutive law
    //  */
    // bool Has(const Variable<Matrix> &rThisVariable) override;

    // /**
    //  * @brief Sets the value of a specified variable (double)
    //  * @param rVariable the variable to be returned
    //  * @param rValue new value of the specified variable
    //  * @param rCurrentProcessInfo the process info
    //  */
    // void SetValue(
    //     const Variable<double> &rThisVariable,
    //     const double& rValue,
    //     const ProcessInfo& rCurrentProcessInfo
    //     ) override;

    // /**
    //  * @brief Sets the value of a specified variable (Vector)
    //  * @param rThisVariable the variable to be returned
    //  * @param rValue new value of the specified variable
    //  * @param rCurrentProcessInfo the process info
    //  */
    // void SetValue(
    //     const Variable<Vector> &rThisVariable,
    //     const Vector& rValue,
    //     const ProcessInfo& rCurrentProcessInfo
    //     ) override;

    // /**
    //  * @brief Returns the value of a specified variable (double)
    //  * @param rThisVariable the variable to be returned
    //  * @param rValue a reference to the returned value
    //  * @return rValue output: the value of the specified variable
    //  */
    // double& GetValue(
    //     const Variable<double> &rThisVariable,
    //     double& rValue
    //     ) override;

    // /**
    //  * @brief Returns the value of a specified variable (Vector)
    //  * @param rThisVariable the variable to be returned
    //  * @param rValue a reference to the returned value
    //  * @return rValue output: the value of the specified variable
    //  */
    // Vector& GetValue(
    //     const Variable<Vector> &rThisVariable,
    //     Vector& rValue
    //     ) override;

    // /**
    //  * @brief Returns the value of a specified variable (matrix)
    //  * @param rThisVariable the variable to be returned
    //  * @param rValue a reference to the returned value
    //  * @return rValue output: the value of the specified variable
    //  */
    // Matrix& GetValue(
    //     const Variable<Matrix>& rThisVariable,
    //     Matrix& rValue
    //     ) override;

    // /**
    //  * @brief Returns the value of a specified variable (double)
    //  * @param rParameterValues the needed parameters for the CL calculation
    //  * @param rThisVariable the variable to be returned
    //  * @param rValue a reference to the returned value
    //  * @param rValue output: the value of the specified variable
    //  */
    // double& CalculateValue(
    //     ConstitutiveLaw::Parameters& rParameterValues,
    //     const Variable<double>& rThisVariable,
    //     double& rValue) override;

    // /**
    //  * @brief Returns the value of a specified variable (vector)
    //  * @param rParameterValues the needed parameters for the CL calculation
    //  * @param rThisVariable the variable to be returned
    //  * @param rValue a reference to the returned value
    //  * @param rValue output: the value of the specified variable
    //  */
    // Vector& CalculateValue(
    //     ConstitutiveLaw::Parameters& rParameterValues,
    //     const Variable<Vector>& rThisVariable,
    //     Vector& rValue
    //     ) override;

    // /**
    //  * @brief Returns the value of a specified variable (matrix)
    //  * @param rParameterValues the needed parameters for the CL calculation
    //  * @param rThisVariable the variable to be returned
    //  * @param rValue a reference to the returned value
    //  * @param rValue output: the value of the specified variable
    //  */
    // Matrix& CalculateValue(
    //     ConstitutiveLaw::Parameters& rParameterValues,
    //     const Variable<Matrix>& rThisVariable,
    //     Matrix& rValue
    //     ) override;

    // /**
    //  * @brief This function provides the place to perform checks on the completeness of the input.
    //  * @details It is designed to be called only once (or anyway, not often) typically at the beginning
    //  * of the calculations, so to verify that nothing is missing from the input or that no common error is found.
    //  * @param rMaterialProperties The properties of the material
    //  * @param rElementGeometry The geometry of the element
    //  * @param rCurrentProcessInfo The current process info instance
    //  * @return 0 if OK, 1 otherwise
    //  */
    // int Check(
    //     const Properties& rMaterialProperties,
    //     const GeometryType& rElementGeometry,
    //     const ProcessInfo& rCurrentProcessInfo
    //     ) const override;

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

    // double& GetThreshold() { return mThreshold; }
    // double& GetPlasticDissipation() { return mPlasticDissipation; }
    // Vector& GetPlasticStrain() { return mPlasticStrain; }

    // void SetThreshold(const double Threshold) { mThreshold = Threshold; }
    // void SetPlasticDissipation(const double PlasticDissipation) { mPlasticDissipation = PlasticDissipation; }
    // void SetPlasticStrain(const array_1d<double, VoigtSize>& rPlasticStrain) { mPlasticStrain = rPlasticStrain; }

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

    void save(Serializer &rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
    }

    void load(Serializer &rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
    }

    ///@}

}; // Class GenericYieldSurface

} // namespace Kratos
