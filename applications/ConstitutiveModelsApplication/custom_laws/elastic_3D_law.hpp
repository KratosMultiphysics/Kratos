//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined (KRATOS_ELASTIC_3D_LAW_H_INCLUDED)
#define  KRATOS_ELASTIC_3D_LAW_H_INCLUDED

// System includes
#include <iostream>
#include <cmath>

// External includes

// Project includes
#include "includes/constitutive_law.h"
#include "custom_utilities/constitutive_law_utilities.hpp"
#include "custom_models/constitutive_model_data.hpp"

namespace Kratos
{
  /**
   * Defines a elastic constitutive law in 3D
   */
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) Elastic3DLaw : public ConstitutiveLaw
  {
  protected:

    using VoigtIndexType = const unsigned int(*)[2];

  public:

    ///@name Type Definitions
    ///@{
    typedef ProcessInfo                                           ProcessInfoType;
    typedef ConstitutiveLaw                                              BaseType;
    typedef ConstitutiveModelData::SizeType                              SizeType;

    typedef ConstitutiveModelData::VectorType                          VectorType;
    typedef ConstitutiveModelData::MatrixType                          MatrixType;
    typedef ConstitutiveModelData::ModelData                        ModelDataType;
    typedef ConstitutiveModelData::ConstitutiveLawData    ConstitutiveLawDataType;
       
    
    /// Pointer definition of Elastic3DLaw
    KRATOS_CLASS_POINTER_DEFINITION( Elastic3DLaw );
	
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Elastic3DLaw();
   
    /// Copy constructor.
    Elastic3DLaw(const Elastic3DLaw& rOther);

    /// Clone.
    ConstitutiveLaw::Pointer Clone() const override;

    /// Assignment operator.
    Elastic3DLaw& operator=(const Elastic3DLaw& rOther);

    /// Destructor.
    virtual ~Elastic3DLaw();

    ///@}
    ///@name Operators
    ///@{

    
    ///@}
    ///@name Operations
    ///@{

    /**
     * Material parameters are inizialized
     */
    void InitializeMaterial(const Properties& rMaterialProperties,
			    const GeometryType& rElementGeometry,
			    const Vector& rShapeFunctionsValues ) override;

    /**
     * Step Initialize
     */
    void InitializeSolutionStep(const Properties& rMaterialProperties,
                                const GeometryType& rElementGeometry, //this is just to give the array of nodes
				const Vector& rShapeFunctionsValues ,
				const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Step Finalize
     */
    void FinalizeSolutionStep(const Properties& rMaterialProperties,
			      const GeometryType& rElementGeometry, //this is just to give the array of nodes
			      const Vector& rShapeFunctionsValues ,
			      const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Computes the material response:
     * PK1 stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    virtual void CalculateMaterialResponsePK1(Parameters & rValues) override;

    /**
     * Computes the material response:
     * PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    virtual void CalculateMaterialResponsePK2(Parameters & rValues) override;

    /**
     * Computes the material response:
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    virtual void CalculateMaterialResponseKirchhoff (Parameters & rValues) override;


    /**
     * Computes the material response:
     * Cauchy stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    virtual void CalculateMaterialResponseCauchy (Parameters & rValues) override;


    /**
     * Updates the material response:
     * Cauchy stresses and Internal Variables
     * @param rValues
     * @see   Parameters
     */
    virtual void FinalizeMaterialResponsePK1(Parameters & rValues) override;

    /**
     * Updates the material response:
     * Cauchy stresses and Internal Variables
     * @param rValues
     * @see   Parameters
     */
    virtual void FinalizeMaterialResponsePK2(Parameters & rValues) override;

    /**
     * Updates the material response:
     * Cauchy stresses and Internal Variables
     * @param rValues
     * @see   Parameters
     */
    virtual void FinalizeMaterialResponseKirchhoff(Parameters & rValues) override;

    /**
     * Updates the material response:
     * Cauchy stresses and Internal Variables
     * @param rValues
     * @see   Parameters
     */
    virtual void FinalizeMaterialResponseCauchy(Parameters & rValues) override;


    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures) override;
    

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{
        
    /**
     * Has Values
     */   
    bool Has(const Variable<double>& rThisVariable) override;
    
    bool Has(const Variable<Vector>& rThisVariable) override;
    
    bool Has(const Variable<Matrix>& rThisVariable) override;

    /**
     * Set Values
     */
    void SetValue(const Variable<double>& rVariable,
                  const double& rValue,
                  const ProcessInfo& rCurrentProcessInfo) override;
    
    void SetValue(const Variable<Vector>& rThisVariable,
                  const Vector& rValue,
                  const ProcessInfo& rCurrentProcessInfo) override;
    
    void SetValue(const Variable<Matrix>& rThisVariable,
                  const Matrix& rValue,
                  const ProcessInfo& rCurrentProcessInfo) override;
    
    /**
     * Get Values
     */
    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override;
    
    Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue) override;

    Matrix& GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue) override;

    ///@}
    ///@name Inquiry
    ///@{
    
    /**
     * Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override
    {
      return 3;
    };

    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
      return 6;
    };

    
    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Elastic3DLaw";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Elastic3DLaw";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "Elastic3DLaw Data";
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
     * Get voigt index tensor:
     */
    virtual VoigtIndexType GetVoigtIndexTensor()
    {
      return this->msIndexVoigt3D6C;
    }

    /**
     * Initialize ModelData type:
     */
    virtual void InitializeModelData(Parameters& rValues, ModelDataType& rModelValues);	


    /**
     * Finalize ModelData type:
     */
    virtual void FinalizeModelData(Parameters& rValues, ModelDataType& rModelValues);
    
    /**
     * Calculates the variables of the domain (element)
     * @param rValues
     * @see   Parameters
     * @param rDataValues
     * @see   ConstitutiveLawData
     */
    virtual void CalculateDomainVariables(Parameters& rValues, ModelDataType& rModelValues);
    
    /**
     * Calculates the Pressure of the domain (element)
     * @param rValues
     * @see   Parameters
     * @param rPressure the calculated pressure to be returned
     */
    virtual double& CalculateDomainVariable(Parameters& rValues, const Variable<double>& rThisVariable, double& rVariable);

    
    /**
     * Calculates the Temperature of the domain (element)
     * @param rValues
     * @see   Parameters
     * @param rTemperature the calculated temperature to be returned
     */
    virtual double& CalculateDomainTemperature(Parameters& rValues, double& rTemperature);

    /**
     * Calculates the Pressure of the domain (element)
     * @param rValues
     * @see   Parameters
     * @param rPressure the calculated pressure to be returned
     */
    virtual double& CalculateDomainPressure (Parameters& rValues, double& rPressure);


    /**
     * This function is designed to be called when before the material response
     * to check if all needed parameters for the constitutive are initialized
     * @param Parameters
     * @return
     */
    virtual bool CheckParameters(Parameters& rValues);


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
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

    virtual void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
  }; // Class Elastic3DLaw

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.
#endif // KRATOS_ELASTIC_3D_LAW_H_INCLUDED  defined 
