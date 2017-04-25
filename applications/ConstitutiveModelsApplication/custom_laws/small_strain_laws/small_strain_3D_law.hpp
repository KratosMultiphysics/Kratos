//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined (KRATOS_SMALL_STRAIN_3D_LAW_H_INCLUDED)
#define  KRATOS_SMALL_STRAIN_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_laws/constitutive_3D_law.hpp"
#include "custom_models/constitutive_model.hpp"

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

  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) SmallStrain3DLaw : public Constitutive3DLaw
  {
  public:
  
    ///@name Type Definitions
    ///@{

    typedef ConstitutiveModel                                     ModelType; //small_strain model
    typedef typename ModelType::Pointer                    ModelTypePointer;
	
    /// Pointer definition of SmallStrain3DLaw
    KRATOS_CLASS_POINTER_DEFINITION(SmallStrain3DLaw);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SmallStrain3DLaw();
    
    /// Constructor.
    SmallStrain3DLaw(ModelTypePointer pModel);

    /// Copy constructor.
    SmallStrain3DLaw (const SmallStrain3DLaw& rOther);

    /// Clone.
    ConstitutiveLaw::Pointer Clone() const;

    /// Assignment operator.
    SmallStrain3DLaw& operator=(const SmallStrain3DLaw& rOther);

    /// Destructor.
    virtual ~SmallStrain3DLaw();

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
        buffer << "SmallStrain3DLaw";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SmallStrain3DLaw";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "SmallStrain3DLaw Data";
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

    //constitutive model
    ModelTypePointer mpModel;

    //internal elastic variables
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


    virtual void CalculateConstitutiveMatrix(Matrix& rConstitutiveMatrix,
					     const Properties& rMaterialProperties);


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
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Constitutive3DLaw )
    }

    virtual void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Constitutive3DLaw )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
  }; // Class SmallStrain3DLaw

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  ///@}

  ///@} addtogroup block
    
}  // namespace Kratos.
#endif // KRATOS_SMALL_STRAIN_3D_LAW_H_INCLUDED  defined 
