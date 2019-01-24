//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_CONSTITUTIVE_MODEL_H_INCLUDED )
#define  KRATOS_CONSTITUTIVE_MODEL_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "custom_utilities/constitutive_model_utilities.hpp"

#include "custom_models/constitutive_model_data.hpp"

namespace Kratos
{
  ///@addtogroup ConstitutiveModelsApplication
  ///@{

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

  /// Short class definition.
  /** Detail class definition.
   */
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) ConstitutiveModel
  {
  protected:

    using VoigtIndexType = const unsigned int(*)[2];

  public:

    //state flags
    KRATOS_DEFINE_LOCAL_FLAG( ADD_HISTORY_VECTOR );
    KRATOS_DEFINE_LOCAL_FLAG( HISTORY_STRAIN_MEASURE );
    KRATOS_DEFINE_LOCAL_FLAG( HISTORY_STRESS_MEASURE );

    ///@name Type Definitions
    ///@{
    typedef ConstitutiveModelData::SizeType                    SizeType;
    typedef ConstitutiveModelData::VectorType                VectorType;
    typedef ConstitutiveModelData::MatrixType                MatrixType;
    typedef ConstitutiveModelData::ModelData              ModelDataType;
    typedef ConstitutiveModelData::MaterialData        MaterialDataType;

    typedef ConstitutiveModelData::StrainMeasureType  StrainMeasureType;
    typedef ConstitutiveModelData::StressMeasureType  StressMeasureType;

    /// Pointer definition of ConstitutiveModel
    KRATOS_CLASS_POINTER_DEFINITION( ConstitutiveModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ConstitutiveModel();

    /// Copy constructor.
    ConstitutiveModel(ConstitutiveModel const& rOther);

    /// Clone.
    virtual ConstitutiveModel::Pointer Clone() const;

    /// Assignment operator.
    ConstitutiveModel& operator=(ConstitutiveModel const& rOther);

    /// Destructor.
    virtual ~ConstitutiveModel();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * Initialize member data
     */
    virtual void InitializeMaterial(const Properties& rProperties);


    /**
     * Initialize member data
     */
    virtual void InitializeModel(ModelDataType& rValues);

    /**
     * Finalize member data
     */
    virtual void FinalizeModel(ModelDataType& rValues);


    /**
     * Calculate Strain Energy Density Functions
     */
    virtual void CalculateStrainEnergy(ModelDataType& rValues, double& rDensityFunction);


    /**
     * Calculate Stresses
     */
    virtual void CalculateStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix);

    virtual void CalculateIsochoricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix);

    virtual void CalculateVolumetricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix);


    /**
     * Calculate Constitutive Tensor
     */
    virtual void CalculateConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutive);

    virtual void CalculateIsochoricConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutive);

    virtual void CalculateVolumetricConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutive);


    /**
     * Calculate Stress and Constitutive Tensor
     */
    virtual void CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutive);

    virtual void CalculateIsochoricStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutive);

    virtual void CalculateVolumetricStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutive);


    /**
     * Check
     */
    virtual int Check(const Properties& rProperties, const ProcessInfo& rCurrentProcessInfo);

    ///@}
    ///@name Access
    ///@{

    /**
     * Has Values
     */
    virtual bool Has(const Variable<double>& rThisVariable);

    /**
     * Set Values
     */
    virtual void SetValue(const Variable<double>& rVariable, const double& rValue,
			  const ProcessInfo& rCurrentProcessInfo);

    virtual void SetValue(const Variable<Vector>& rThisVariable, const Vector& rValue,
			  const ProcessInfo& rCurrentProcessInfo);

    virtual void SetValue(const Variable<Matrix>& rThisVariable, const Matrix& rValue,
			  const ProcessInfo& rCurrentProcessInfo);

    /**
     * Get Values
     */
    virtual double& GetValue(const Variable<double>& rThisVariable, double& rValue);

    /**
     * method to ask the constituitve model the list of variables (dofs) needed from the domain
     * @param rScalarVariables : list of scalar dofs
     * @param rComponentVariables :  list of vector dofs
     */
    virtual void GetDomainVariablesList(std::vector<Variable<double> >& rScalarVariables,
					std::vector<Variable<array_1d<double,3> > >& rComponentVariables);

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "ConstitutiveModel";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ConstitutiveModel";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
      rOStream << "ConstitutiveModel Data";
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

    Flags mOptions;

    //initial or historical strains/stresses
    VectorType mHistoryVector;

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
    ///@name Serialization
    ///@{
    friend class Serializer;


    virtual void save(Serializer& rSerializer) const
    {
      rSerializer.save("mOptions",mOptions);
      rSerializer.save("mHistoryVector",mHistoryVector);
    }

    virtual void load(Serializer& rSerializer)
    {
      rSerializer.load("mOptions",mOptions);
      rSerializer.load("mHistoryVector",mHistoryVector);
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class ConstitutiveModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CONSTITUTIVE_MODEL_H_INCLUDED  defined
