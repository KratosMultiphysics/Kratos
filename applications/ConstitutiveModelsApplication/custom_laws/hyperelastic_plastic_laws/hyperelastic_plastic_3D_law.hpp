//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                December 2016 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined (KRATOS_HYPERELASTIC_PLASTIC_3D_LAW_H_INCLUDED)
#define  KRATOS_HYPERELASTIC_PLASTIC_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_laws/hyperelastic_laws/hyperelastic_3d_law.hpp"
#include "custom_laws/plasticity_laws/plasticity_model.hpp"

namespace Kratos
{
  /**
   * Defines a hyperelastic-plastic isotropic constitutive law in 3D 
   * the functionality is limited to large displacements plasticity.
   */
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) HyperElasticPlastic3DLaw : public HyperElastic3DLaw
  {
  public:

    ///@name Type Definitions
    ///@{      
    typedef HyperElasticityModel                                               ElasticModelType;
    typedef HardeningLaw                                                       HardeningLawType; 
    typedef YieldCriterion<HardeningLawType>                                 YieldCriterionType;
    typedef PlasticityModel<ElasticModelType,YieldCriterionType>                      ModelType;

       
    /// Pointer definition of HyperElasticPlastic3DLaw
    KRATOS_CLASS_POINTER_DEFINITION( HyperElasticPlastic3DLaw );
    
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    HyperElasticPlastic3DLaw();

    /// Constructor.
    HyperElasticPlastic3DLaw(ModelType::Pointer pModel); 

    /// Copy constructor.
    HyperElasticPlastic3DLaw (const HyperElasticPlastic3DLaw& rOther);

    /// Clone.
    ConstitutiveLaw::Pointer Clone() const override;
    
    /// Assignment operator.
    HyperElasticPlastic3DLaw& operator=(const HyperElasticPlastic3DLaw& rOther);
    
    /// Destructor.
    virtual ~HyperElasticPlastic3DLaw();

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
    virtual void CalculateMaterialResponsePK2(Parameters & rValues) override;


    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures) override;

    
    ///@}
    ///@name Access
    ///@{

    
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
      buffer << "HyperElasticPlastic3DLaw";
      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "HyperElasticPlastic3DLaw";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "HyperElasticPlastic3DLaw Data";
      rOStream << "InverseDeformationGradientF0 "<< mInverseDeformationGradientF0;
      rOStream << "CauchyGreenVector "<< mCauchyGreenVector;
      rOStream << "DeterminantF0 "<< mDeterminantF0;
      mpModel->PrintData(rOStream);
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
    	
    ///@}
    ///@name Protected Operators
    ///@{
    
    ///@}
    ///@name Protected Operations
    ///@{


    /**
     * Initialize ModelData type:
     */
    virtual void InitializeModelData(Parameters& rValues, ModelDataType& rModelValues) override;	

    /**
     * Finalize ModelData type:
     */
    virtual void FinalizeModelData(Parameters& rValues, ModelDataType& rModelValues) override;
    

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
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HyperElastic3DLaw )
    }

    virtual void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HyperElastic3DLaw )    	
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
  }; // Class HyperElasticPlastic3DLaw
  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.
#endif // KRATOS_HYPERELASTIC_PLASTIC_3D_LAW_H_INCLUDED defined
