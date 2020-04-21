// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo 
//  Collaborator:    Lucia Barbu
//

#if !defined(KRATOS_GENERIC_ANISOTROPIC_LAW_H_INCLUDED)
#define KRATOS_GENERIC_ANISOTROPIC_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"
#include "custom_utilities/constitutive_law_utilities.h"


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
 * @class GenericAnisotropicLaw3D
 * @ingroup StructuralMechanicsApplication
 * @brief This CL takes into account the material anisotropy in terms of 
 * young modulus, poisson ratio, orientation and strengths. 
 * @details See "Comportamiento mecánico de los materiales compuestos" by S. Oller,
 * Chapter "Anisotropía mecánica"
 * @author Alejandro Cornejo
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) GenericAnisotropicLaw3D
    : public ConstitutiveLaw
{
  public:
    ///@name Type Definitions
    ///@{

    /// The node definition
    typedef Node<3> NodeType;

    /// The geometry definition
    typedef Geometry<NodeType> GeometryType;

    /// Definition of the machine precision tolerance
    static constexpr double machine_tolerance = std::numeric_limits<double>::epsilon();

    /// The size type definition
    typedef std::size_t SizeType;
    
    /// Static definition of the dimension
    static constexpr SizeType Dimension = 3;
    
    /// Static definition of the VoigtSize
    static constexpr SizeType VoigtSize = 6;

    /// Counted pointer of GenericAnisotropicLaw3D
    KRATOS_CLASS_POINTER_DEFINITION(GenericAnisotropicLaw3D);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * Constructor.
    */
    GenericAnisotropicLaw3D()
    {
    }

    /**
    * Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<GenericAnisotropicLaw3D>(*this);
    }

    // Copy constructor
    GenericAnisotropicLaw3D(GenericAnisotropicLaw3D const& rOther)
        : ConstitutiveLaw(rOther), 
        mpIsotropicCL(rOther.mpIsotropicCL)
    {
    }

    /**
    * Destructor.
    */
    ~GenericAnisotropicLaw3D() override
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
        return mpIsotropicCL->WorkingSpaceDimension();
    };

    /**
     * @brief Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return mpIsotropicCL->GetStrainSize();
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
     * Initialize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponsePK2(Parameters& rValues) override;  


    /**
     * @brief This method computes the tangent tensor
     * @param rValues The constitutive law parameters and flags
     */
    void CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues);

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

    void CalculateOrthotropicElasticMatrix(
        Matrix &rElasticityTensor,
        const Properties &rMaterialProperties);

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

    /**
     * @brief This method computes the elastic tensor
     * @param rElasticityTensor The elastic tensor
     * @param rMaterialProperties The material properties
     */
    void CalculateElasticMatrix(
        Matrix& rElasticityTensor,
        const Properties& rMaterialProperties);

    void CalculateCauchyGreenStrain(
        ConstitutiveLaw::Parameters &rValues,
        Vector &rStrainVector);

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

}; // Class GenericAnisotropicLaw3D

} // namespace Kratos

#endif
