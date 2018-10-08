//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Tescheamacher
//                   Riccardo Rossi
//

#if !defined(KRATOS_PLANE_STRESS_2D_KINEMATICALLY_ENRICHED_LAW_H_INCLUDED)
#define KRATOS_PLANE_STRESS_2D_KINEMATICALLY_ENRICHED_LAW_H_INCLUDED

// System includes

// Project includes
#include "includes/constitutive_law.h"

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

class KRATOS_API(IGA_STRUCTURAL_MECHANICS_APPLICATION) PlaneStress2dKinematicallyEnrichedLaw
    : public ConstitutiveLaw
{
  public:
    ///@name Type Definitions
    ///@{

    /// Definition of the base class
    typedef ConstitutiveLaw BaseType;

    /// The index definition
    typedef std::size_t IndexType;

    /// The size definition
    typedef std::size_t SizeType;

    /// Counted pointer of GenericYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(PlaneStress2dKinematicallyEnrichedLaw);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * Default constructor.
    */
    PlaneStress2dKinematicallyEnrichedLaw()
    {
    }

    /**
    * Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<PlaneStress2dKinematicallyEnrichedLaw>(*this);
    }

    /**
    * Copy constructor.
    */
    PlaneStress2dKinematicallyEnrichedLaw(const PlaneStress2dKinematicallyEnrichedLaw& rOther)
        : ConstitutiveLaw(rOther)
    {
    }
    /**
    * Destructor.
    */
    ~PlaneStress2dKinematicallyEnrichedLaw() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override
    {
        return 3;
    };

    /**
     * @brief Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return 3;
    };


    /**
    * returns the value of a specified variable
    * @param rThisVariable the variable to be returned
    * @param rValue a reference to the returned value
    * @param rValue output: the value of the specified variable
    */
    virtual double& GetValue(
        const Variable<double>& rThisVariable,
        double& rValue) override;
    /**
    * @brief Returns the value of a specified variable (Vector)
    * @param rThisVariable the variable to be returned
    * @param rValue a reference to the returned value
    * @return rValue output: the value of the specified variable
    */
    virtual Vector& GetValue(
        const Variable<Vector>& rThisVariable,
        Vector& rValue) override;
    /**
    * @brief Returns the value of a specified variable (Matrix)
    * @param rThisVariable the variable to be returned
    * @return rValue output: the value of the specified variable
    */
    virtual Matrix& GetValue(
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue) override;

    void InitializeMaterial(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues);

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
     * to be called at the end of each solution step
     * (e.g. from Element::FinalizeSolutionStep)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    //void FinalizeSolutionStep(
    //    const Properties& rMaterialProperties,
    //    const GeometryType& rElementGeometry,
    //    const Vector& rShapeFunctionsValues,
    //    const ProcessInfo& rCurrentProcessInfo) override;

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
     * @brief Calculates the value of a specified variable (Matrix)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    //Matrix& CalculateValue(
    //    ConstitutiveLaw::Parameters& rParameterValues,
    //    const Variable<Matrix>& rThisVariable,
    //    Matrix& rValue
    //    ) override;

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

  private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

      Matrix m_D0;

      double m_compression_parameter_A;
      double m_compression_parameter_B;

      double m_tension_parameter_A;

      double m_compressive_strength;
      double m_tensile_strength;

      double m_beta;
      double m_Gf;

      double m_E;
      double m_nu;

      double m_K;

      Matrix m_eigen_vectors;
      Vector m_eigen_values;

      double m_damage_treshold_compression_initial;
      double m_damage_treshold_tension_initial;

      double m_damage_t;
      double m_damage_c;

    //// Converged values
    //Vector mPrevStressVector = ZeroVector(6);
    //Vector mPrevInelasticStrainVector = ZeroVector(6);

    //// Non Converged values
    //Vector mNonConvPrevStressVector = ZeroVector(6);
    //Vector mNonConvPrevInelasticStrainVector = ZeroVector(6);

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

      void CalculateMaterialResponseInternal(
          const Vector& rStrainVector,
          Vector& rStressVector,
          Matrix& rConstitutiveLaw);

    /**
     * @brief This method computes the elastic tensor
     * @param rElasticityTensor The elastic tensor
     * @param rMaterialProperties The material properties
     */
    void CalculateElasticityMatrix(
        Matrix& rElasticityMatrix);

    void CalculateDamageTension(
        const double& rDamageTreshholdTension,
        double& rDamageTension);

    void CalculateDamageCriterionTension(
        const Vector& rStressVector,
        double& rDamageTreshholdTension);

    void CalculateDamageCompression(
        const double& rDamageTreshholdCompression,
        double& rDamageCompression);

    void CalculateDamageCriterionCompression(
        const Vector& rStressVector,
        double& rDamageTreshholdCompression);

    void SpectralDecompositionStrain(
        const Vector& rStrainVector,
        Matrix& rQ_CW_Tension);

    void SpectralDecomposition(
        const Vector& StrainVector,
        Vector& StressVectorTension,
        Vector& StressVectorCompression,
        Matrix& PVectorTension,
        Matrix& PVectorCompression);
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
        //rSerializer.save("PrevStressVector", mPrevStressVector);
        //rSerializer.save("PrevInelasticStrainVector", mPrevInelasticStrainVector);
        //rSerializer.save("NonConvPrevStressVector", mNonConvPrevStressVector);
        //rSerializer.save("NonConvPrevInelasticStrainVector", mNonConvPrevInelasticStrainVector);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
        //rSerializer.load("PrevStressVector", mPrevStressVector);
        //rSerializer.load("PrevInelasticStrainVector", mPrevInelasticStrainVector);
        //rSerializer.load("NonConvPrevStressVector", mNonConvPrevStressVector);
        //rSerializer.load("NonConvPrevInelasticStrainVector", mNonConvPrevInelasticStrainVector);
    }

    ///@}

}; // Class GenericYieldSurface

} // namespace Kratos
#endif
