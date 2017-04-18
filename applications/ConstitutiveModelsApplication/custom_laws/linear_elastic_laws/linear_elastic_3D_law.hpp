//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                December 2016 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined (KRATOS_LINEAR_ELASTIC_3D_LAW_H_INCLUDED)
#define  KRATOS_LINEAR_ELASTIC_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_laws/elastic_3D_law.hpp"

namespace Kratos
{
  /**
   * Defines a linear isotropic constitutive law
   * This material law is defined by the parameters:
   * 1) YOUNG MODULUS
   * 2) POISSON RATIO 
   * As there are no further parameters the functionality is valid
   * for small and large displacements elasticity.
   */

  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) LinearElastic3DLaw : public Elastic3DLaw
  {
  public:
  
    ///@name Type Definitions
    ///@{
	
    /// Pointer definition of LinearElastic3DLaw
    KRATOS_CLASS_POINTER_DEFINITION(LinearElastic3DLaw);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LinearElastic3DLaw();

    /// Copy constructor.
    LinearElastic3DLaw (const LinearElastic3DLaw& rOther);

    /// Clone.
    ConstitutiveLaw::Pointer Clone() const;

    /// Assignment operator.
    LinearElastic3DLaw& operator=(const LinearElastic3DLaw& rOther);

    /// Destructor.
    virtual ~LinearElastic3DLaw();

    ///@}
    ///@name Operators
    ///@{

    
    ///@}
    ///@name Operations
    ///@{
    

    /**
     * Computes the material response:
     * PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    void CalculateMaterialResponsePK2 (Parameters & rValues);

    /**
     * Computes the material response:
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff (Parameters & rValues);

    
    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures);

    ///@}
    ///@name Access
    ///@{

    /**
     * Set Values
     */
    
    void SetValue(const Variable<Vector>& rThisVariable,
                  const Vector& rValue,
                  const ProcessInfo& rCurrentProcessInfo) override;
    
    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "LinearElastic3DLaw";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "LinearElastic3DLaw";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "LinearElastic3DLaw Data";
    }


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

    VectorType  mInitialStrainVector;
    
    ///@}
    ///@name Protected Operators
    ///@{
    
    ///@}
    ///@name Protected Operations
    ///@{
    ///@}

    
    /**
     * Calculates the stresses for given strain state
     * @param rStrainVector
     * @param rConstitutiveMatrix
     * @param rStressVector the stress vector corresponding to the deformation
     */
    virtual void CalculateStress(const Vector &rStrainVector,
				 const Matrix &rConstitutiveMatrix,
				 Vector& rStressVector);


    /**
     * calculates the linear elastic constitutive matrix in terms of Young's modulus and
     * Poisson ratio
     * @param E the Young's modulus
     * @param NU the Poisson ratio
     * @return the linear elastic constitutive matrix
     */


    virtual void CalculateLinearElasticMatrix(Matrix& rConstitutiveMatrix,
					      const double &rYoungModulus,
					      const double &rPoissonCoefficient );


    /**
     * Adds Initial Strain to Strain Vector
     * @param rStrainVector
     */
    void AddInitialStrainVector(Vector& rStrainVector);

    
    
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
    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Elastic3DLaw )
    }

    virtual void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Elastic3DLaw )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
  }; // Class LinearElastic3DLaw

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  ///@}

  ///@} addtogroup block
    
}  // namespace Kratos.
#endif // KRATOS_LINEAR_ELASTIC_3D_LAW_H_INCLUDED  defined 
