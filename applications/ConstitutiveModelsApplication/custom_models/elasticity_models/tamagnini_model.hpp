//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                MCiantia $
//   Date:                $Date:                    JULY 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_TAMAGNINI_MODEL_H_INCLUDED )
#define  KRATOS_TAMAGNINI_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/elasticity_models/borja_model.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) TamagniniModel : public BorjaModel
  {
  public:
    
    ///@name Type Definitions
    ///@{
    
    /// Pointer definition of TamagniniModel
    KRATOS_CLASS_POINTER_DEFINITION( TamagniniModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.    
    TamagniniModel(); 

    /// Copy constructor.
    TamagniniModel(TamagniniModel const& rOther);

    /// Assignment operator.
    TamagniniModel& operator=(TamagniniModel const& rOther);

    /// Clone.
    virtual ConstitutiveModel::Pointer Clone() const override;


    /// Destructor.
    virtual ~TamagniniModel();


    ///@}
    ///@name Operators
    ///@{

    
    ///@}
    ///@name Operations
    ///@{
    
    
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
        buffer << "TamagniniModel";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "TamagniniModel";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "TamagniniModel Data";
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
     * Calculate Stresses
     */
    virtual void CalculateAndAddStressTensor(HyperElasticDataType& rVariables, MatrixType& rStressMatrix) override;

    /**
     * Calculate Constitutive Tensor
     */
    virtual void CalculateAndAddConstitutiveTensor(HyperElasticDataType& rVariables, Matrix& rConstitutiveMatrix) override;


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


    virtual void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BorjaModel )
    }

    virtual void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BorjaModel )      
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class TamagniniModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_BORJA_MODEL_H_INCLUDED  defined 


