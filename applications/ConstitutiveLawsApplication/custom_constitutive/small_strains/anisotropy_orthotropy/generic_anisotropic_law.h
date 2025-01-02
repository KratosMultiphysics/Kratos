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
//  Collaborator:    Lucia Barbu
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

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
template <SizeType TDim>
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) GenericAnisotropicLaw
    : public ConstitutiveLaw
{
public:
    ///@name Type Definitions
    ///@{

    /// The node definition
    typedef Node NodeType;

    /// The geometry definition
    typedef Geometry<NodeType> GeometryType;

    /// Definition of the machine precision tolerance
    static constexpr double machine_tolerance = std::numeric_limits<double>::epsilon();

    /// Static definition of the dimension
    static constexpr SizeType Dimension = TDim;

    /// Static definition of the VoigtSize
    static constexpr SizeType VoigtSize = (Dimension == 3) ? 6 : 3;

    /// The definition of the bounded matrix type
    using BoundedMatrixType = BoundedMatrix<double, Dimension, Dimension>;

    /// The definition of the bounded matrix type
    using BoundedMatrixVoigtType = BoundedMatrix<double, VoigtSize, VoigtSize>;

    /// Counted pointer of GenericAnisotropicLaw
    KRATOS_CLASS_POINTER_DEFINITION(GenericAnisotropicLaw);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * Constructor.
    */
    GenericAnisotropicLaw()
    {
    }

    /**
    * Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<GenericAnisotropicLaw<TDim>>(*this);
    }

    // Copy constructor
    GenericAnisotropicLaw(GenericAnisotropicLaw const& rOther)
        : ConstitutiveLaw(rOther),
        mpIsotropicCL(rOther.mpIsotropicCL)
    {
    }

    /**
    * Destructor.
    */
    ~GenericAnisotropicLaw() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * creates a new constitutive law pointer
     * @param NewParameters The configuration parameters of the new constitutive law
     * @return a Pointer to the new constitutive law
     */
    ConstitutiveLaw::Pointer Create(Kratos::Parameters NewParameters) const override;

    /**
     * @brief Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override
    {
        return Dimension;
    };

    /**
     * @brief Voigt tensor size:
     */
    SizeType GetStrainSize() const override
    {
        return VoigtSize;
    };

    /**
     * Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Computes the material response in terms of Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Finalize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Finalize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Finalize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Finalize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Returns the value of a specified variable (Vector)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue) override;

    /**
     * @brief Returns the value of a specified variable (double)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override;

    /**
     * @brief Returns the value of a specified variable (Matrix)
     * @param rThisVariable the variable to be returned
     * @return rValue output: the value of the specified variable
     */
    Matrix& GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (boolean)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<bool>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (double)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<double>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (Vector)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Vector>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (Matrix)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Matrix>& rThisVariable) override;

    /**
     * @brief Calculates the value of a specified variable (double)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    double& CalculateValue(
        Parameters& rParameterValues,
        const Variable<double>& rThisVariable,
        double& rValue) override;

    /**
     * @brief Calculates the value of a specified variable (Vector)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    Vector& CalculateValue(
        Parameters& rParameterValues,
        const Variable<Vector>& rThisVariable,
        Vector& rValue) override;

    /**
     * @brief Calculates the value of a specified variable (Matrix)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    Matrix& CalculateValue(
        Parameters& rParameterValues,
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue) override;

    /**
     * This is to be called at the very beginning of the calculation
     * (e.g. from InitializeElement) in order to initialize all relevant
     * attributes of the constitutive law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    void InitializeMaterial(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues) override;


    /**
     * @brief This computes the mapper operator between the stresses in the isotropic
     * "ficticious" space and the real anisotropic space. S_iso = As*Sa_niso
     * @param rValues The values of the constitutive la
     * @param rAs The stress mapper operator
     * @param rAs The stress mapper operator inverse
     * @note Eq.(2.39) S. Oller book: Comportamiento mecánico de los materiales compuestos
     */
    void CalculateAnisotropicStressMapperMatrix(
        const Properties& rProperties,
        BoundedMatrixVoigtType &rAs,
        BoundedMatrixVoigtType& rAsInv
        );

    /**
     * @brief This computes the rotation matrix, from global to local
     */
    void CalculateRotationMatrixVoigt(
        const Properties& rProperties,
        BoundedMatrixVoigtType &rT);

    /**
     * @brief This computes the mapper operator between the strain in the isotropic
     * "ficticious" space and the real anisotropic space
     * @param rValues The values of the constitutive la
     * @param rAs The mapper operator
     * @param rAs The mapper operator inverse
     * @note Eq.(2.35) S. Oller book: Comportamiento mecánico de los materiales compuestos
     */
    void CalculateAnisotropicStrainMapperMatrix(
        const BoundedMatrixVoigtType& rAnisotropicElasticMatrix,
        const BoundedMatrixVoigtType& rIsotropicElasticMatrix,
        const BoundedMatrixVoigtType &rAs,
        BoundedMatrixVoigtType& rAe);

    /**
     * Initialize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponsePK2(Parameters& rValues) override;

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresInitializeMaterialResponse() override
    {
        return mpIsotropicCL->RequiresInitializeMaterialResponse();
    }

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresFinalizeMaterialResponse() override
    {
        return mpIsotropicCL->RequiresFinalizeMaterialResponse();
    }

    /**
     * @brief This method computes the orthotropic elastic constitutive matrix
     * This method is overriden in the derived class for plane stress
     */
    virtual void CalculateOrthotropicElasticMatrix(
        BoundedMatrixVoigtType &rElasticityTensor,
        const Properties &rMaterialProperties);

    /**
     * @brief This method checks the properties in the nested CL
     */
    int Check(const Properties &rMaterialProperties,
              const GeometryType &rElementGeometry,
              const ProcessInfo &rCurrentProcessInfo) const override;

    /**
     * @brief This method computes an estimation of the tangent constitutive matrix, which relates the
     * stress increment with the strain increment.
     */
    void CalculateTangentTensor(ConstitutiveLaw::Parameters &rValues);

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

    /**
     * @brief This method sets the constitutive law of the isotropic space
     */
    ConstitutiveLaw::Pointer GetIsotropicConstitutiveLaw()
    {
        return mpIsotropicCL;
    }

    /**
     * @brief This method sets the constitutive law of the isotropic space
     */
    void SetIsotropicConstitutiveLaw(ConstitutiveLaw::Pointer pIsotropicConstitutiveLaw)
    {
        mpIsotropicCL = pIsotropicConstitutiveLaw;
    }

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

    ConstitutiveLaw::Pointer mpIsotropicCL;

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.save("IsotropicCL", mpIsotropicCL);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.load("IsotropicCL", mpIsotropicCL);
    }

    ///@}

}; // Class GenericAnisotropicLaw

} // namespace Kratos
