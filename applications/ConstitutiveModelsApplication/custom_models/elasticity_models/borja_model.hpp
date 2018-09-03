//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_BORJA_MODEL_H_INCLUDED )
#define  KRATOS_BORJA_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/elasticity_models/hencky_hyper_elastic_model.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) BorjaModel : public HenckyHyperElasticModel
  {
  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of BorjaModel
    KRATOS_CLASS_POINTER_DEFINITION( BorjaModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BorjaModel();

    /// Copy constructor.
    BorjaModel(BorjaModel const& rOther);

    /// Assignment operator.
    BorjaModel& operator=(BorjaModel const& rOther);

    /// Clone.
    ConstitutiveModel::Pointer Clone() const override;


    /// Destructor.
    ~BorjaModel() override;


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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "BorjaModel";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "BorjaModel";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "BorjaModel Data";
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
    void CalculateAndAddStressTensor(HyperElasticDataType& rVariables, MatrixType& rStressMatrix) override;

    /**
     * Calculate Constitutive Tensor
     */
    void CalculateAndAddConstitutiveTensor(HyperElasticDataType& rVariables, Matrix& rConstitutiveMatrix) override;

    /**
      * Calculate some strain invariants
      */
    void SeparateVolumetricAndDeviatoricPart( const MatrixType& rA, double & rVolumetric, MatrixType& rDev, double & devNorm);

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


    void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HenckyHyperElasticModel )
    }

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HenckyHyperElasticModel )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class BorjaModel

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


