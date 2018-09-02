//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  IPouplana $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_MODIFIED_MISES_YIELD_SURFACE_H_INCLUDED )
#define  KRATOS_MODIFIED_MISES_YIELD_SURFACE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/yield_surfaces/yield_surface.hpp"

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
  template<class THardeningRule>
  class ModifiedMisesYieldSurface : public YieldSurface<THardeningRule>
  {
  public:
    ///@name Type Definitions
    ///@{

    typedef ConstitutiveModelData::MatrixType                          MatrixType;
    typedef ConstitutiveModelData::VectorType                          VectorType;
    typedef ConstitutiveModelData::ModelData                        ModelDataType;
    typedef ConstitutiveModelData::MaterialData                  MaterialDataType;

    typedef YieldSurface<THardeningRule>                                 BaseType;
    typedef typename BaseType::Pointer                            BaseTypePointer;
    typedef typename BaseType::PlasticDataType                    PlasticDataType;

    /// Pointer definition of ModifiedMisesYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION( ModifiedMisesYieldSurface );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ModifiedMisesYieldSurface() : BaseType() {}

    /// Copy constructor.
    ModifiedMisesYieldSurface(ModifiedMisesYieldSurface const& rOther) : BaseType(rOther) {}

    /// Assignment operator.
    ModifiedMisesYieldSurface& operator=(ModifiedMisesYieldSurface const& rOther)
    {
      BaseType::operator=(rOther);
      return *this;
    }

    /// Clone.
    BaseTypePointer Clone() const override
    {
      return Kratos::make_shared<ModifiedMisesYieldSurface>(*this);
    }

    /// Destructor.
    ~ModifiedMisesYieldSurface() override {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * Calculate Yield Condition
     */

    double& CalculateYieldCondition(const PlasticDataType& rVariables, double & rYieldCondition) override
    {
      KRATOS_TRY

      const ModelDataType& rModelData = rVariables.GetModelData();

      // Compute I1
      const MatrixType& rStrainMatrix = rVariables.GetStrainMatrix();
      double I1 = 0.0;

      for(unsigned int i=0; i<3; i++)
	{
	  I1 += rStrainMatrix(i,i);
	}

      // Compute J2
      MatrixType DeviatoricStrain;
      noalias(DeviatoricStrain) = rStrainMatrix;

      for(unsigned int i=0; i<3; i++)
	{
	  DeviatoricStrain(i,i) -= I1/3.0;
	}

      MatrixType Auxiliar;
      noalias(Auxiliar) = prod(DeviatoricStrain,DeviatoricStrain);
      double J2 = 0.0;

      for(unsigned int i=0; i<3; i++)
	{
	  J2 += Auxiliar(i,i);
	}

      J2 *= 0.5;

      // Compute Equivalent Strain (rYieldCondition)
      const Properties& rMaterialProperties = rModelData.GetMaterialProperties();
      const double& StrengthRatio = rMaterialProperties[STRENGTH_RATIO];
      const double& PoissonRatio  = rMaterialProperties[POISSON_RATIO];

      rYieldCondition  =  I1*(StrengthRatio-1.0)/(2.0*StrengthRatio*(1.0-2.0*PoissonRatio));
      rYieldCondition +=  sqrt( I1*I1*(StrengthRatio-1.0)*(StrengthRatio-1.0)/((1.0-2.0*PoissonRatio)*(1.0-2.0*PoissonRatio)) + J2*12.0*StrengthRatio/((1.0+PoissonRatio)*(1.0+PoissonRatio)) )/(2.0*StrengthRatio);

      return rYieldCondition;

      KRATOS_CATCH(" ")

    }
    /**
     * Calculate State Function
     */

    double& CalculateStateFunction(const PlasticDataType& rVariables, double & rStateFunction) override
    {
      KRATOS_TRY

      rStateFunction = this->mHardeningRule.CalculateHardening(rVariables,rStateFunction);

      return rStateFunction;

      KRATOS_CATCH(" ")
    }

    /**
     * Calculate State Function derivative
     */

    double& CalculateDeltaStateFunction(const PlasticDataType& rVariables, double & rDeltaStateFunction) override
    {
      KRATOS_TRY

      rDeltaStateFunction = this->mHardeningRule.CalculateDeltaHardening(rVariables,rDeltaStateFunction);

      return rDeltaStateFunction;

      KRATOS_CATCH(" ")
    }

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
      buffer << "ModifiedMissesYieldSurface" ;
      return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "ModifiedMissesYieldSurface";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "ModifiedMissesYieldSurface Data";
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
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
    }

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class ModifiedMisesYieldSurface

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MODIFIED_MISES_YIELD_SURFACE_H_INCLUDED  defined
