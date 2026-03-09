//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_NON_LINEAR_RATE_DEPENDENT_PLASTICITY_MODEL_H_INCLUDED )
#define  KRATOS_NON_LINEAR_RATE_DEPENDENT_PLASTICITY_MODEL_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/non_linear_associative_plasticity_model.hpp"

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
  template<class TElasticityModel, class TYieldSurface>
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) NonLinearRateDependentPlasticityModel : public NonLinearAssociativePlasticityModel<TElasticityModel,TYieldSurface>
  {
  public:

    ///@name Type Definitions
    ///@{

    //elasticity model
    typedef TElasticityModel                               ElasticityModelType;
    typedef typename  ElasticityModelType::Pointer      ElasticityModelPointer;

    //yield surface
    typedef TYieldSurface                                     YieldSurfaceType;
    typedef typename YieldSurfaceType::Pointer             YieldSurfacePointer;

    //derived type
    typedef NonLinearAssociativePlasticityModel<ElasticityModelType,YieldSurfaceType>   DerivedType;

    //base type
    typedef PlasticityModel<ElasticityModelType,YieldSurfaceType>     BaseType;

    //common types
    typedef typename BaseType::Pointer                         BaseTypePointer;
    typedef typename BaseType::SizeType                               SizeType;
    typedef typename BaseType::VoigtIndexType                   VoigtIndexType;
    typedef typename BaseType::MatrixType                           MatrixType;
    typedef typename BaseType::ModelDataType                     ModelDataType;
    typedef typename BaseType::MaterialDataType               MaterialDataType;
    typedef typename BaseType::PlasticDataType                 PlasticDataType;
    typedef typename BaseType::InternalVariablesType     InternalVariablesType;

    /// Pointer definition of NonLinearRateDependentPlasticityModel
    KRATOS_CLASS_POINTER_DEFINITION( NonLinearRateDependentPlasticityModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NonLinearRateDependentPlasticityModel() : DerivedType() {}

    /// Constructor.
    NonLinearRateDependentPlasticityModel(ElasticityModelPointer pElasticityModel, YieldSurfacePointer pYieldSurface) : DerivedType(pElasticityModel, pYieldSurface) {}

    /// Copy constructor.
    NonLinearRateDependentPlasticityModel(NonLinearRateDependentPlasticityModel const& rOther) : DerivedType(rOther) {}

    /// Assignment operator.
    NonLinearRateDependentPlasticityModel& operator=(NonLinearRateDependentPlasticityModel const& rOther)
    {
      DerivedType::operator=(rOther);
      return *this;
    }

    /// Clone.
    ConstitutiveModel::Pointer Clone() const override
    {
      return Kratos::make_shared<NonLinearRateDependentPlasticityModel>(*this);
    }

    /// Destructor.
    ~NonLinearRateDependentPlasticityModel() override {}


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
      buffer << "NonLinearRateDependentPlasticityModel" ;
      return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "NonLinearRateDependentPlasticityModel";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "NonLinearRateDependentPlasticityModel Data";
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


    // calculate return mapping
    bool CalculateReturnMapping(PlasticDataType& rVariables, MatrixType& rStressMatrix) override
    {
      KRATOS_TRY

      bool converged    = false;

      //Start 1rst Newton Raphson iteration
      rVariables.State().Set(ConstitutiveModelData::PLASTIC_RATE_REGION,true);
      rVariables.RateFactor=1; //plastic rate region on
      converged = CalculateRateDependentReturnMapping(rVariables,rStressMatrix);

      // if(!converged)
      //   std::cout<<" ConstitutiveLaw did not converge on the rate dependent return mapping"<<std::endl;

      const ModelDataType& rModelData   = rVariables.GetModelData();
      const double& rDeltaTime          = rModelData.GetProcessInfo()[DELTA_TIME];

      const Properties& rProperties     = rModelData.GetProperties();
      const double& rPlasticStrainRate  = rProperties[PLASTIC_STRAIN_RATE];


      double MaterialDeltaPlasticStrain = rPlasticStrainRate * rDeltaTime;

      //std::cout<<" DeltaPlasticStrain: "<<rPlasticVariables.DeltaPlasticStrain<<" MaterialDeltaPlasticStrain: "<<MaterialDeltaPlasticStrain<<std::endl;

      double& rDeltaGamma = rVariables.DeltaInternal.Variables[0];

      if(  sqrt(2.0/3.0) * rDeltaGamma < MaterialDeltaPlasticStrain ){

	//std::cout<<" DeltaPlasticStrain: "<<rPlasticVariables.DeltaPlasticStrain<<" MaterialDeltaPlasticStrain: "<<MaterialDeltaPlasticStrain<<std::endl;

	//Start 2nd Newton Raphson iteration
	rVariables.State().Set(ConstitutiveModelData::PLASTIC_RATE_REGION,false);
	rVariables.RateFactor=0; //plastic rate region on

	converged = CalculateRateIndependentReturnMapping(rVariables,rStressMatrix);

	// if(!converged)
	//   std::cout<<" ConstitutiveLaw did not converge on the rate independent return mapping"<<std::endl;

      }

      return converged;

      KRATOS_CATCH(" ")
    }


    bool CalculateRateDependentReturnMapping(PlasticDataType& rVariables, MatrixType& rStressMatrix)
    {
      KRATOS_TRY

      //Set convergence parameters
      unsigned int iter    = 0;
      double Tolerance     = 1e-5;
      double MaxIterations = 50;

      //start
      double DeltaDeltaGamma    = 0;
      double DeltaStateFunction = 0;
      double DeltaPlasticStrain = 0;

      double& rEquivalentPlasticStrainOld  = this->mPreviousInternal.Variables[0];
      double& rEquivalentPlasticStrain     = rVariables.Internal.Variables[0];
      double& rDeltaGamma                  = rVariables.DeltaInternal.Variables[0];


      const ModelDataType& rModelData      = rVariables.GetModelData();
      const double& rDeltaTime             = rModelData.GetProcessInfo()[DELTA_TIME];

      const Properties& rProperties        = rModelData.GetProperties();
      const double& rPlasticStrainRate     = rProperties[PLASTIC_STRAIN_RATE];


      rEquivalentPlasticStrain = 0;
      rDeltaGamma = sqrt(3.0*0.5) * rPlasticStrainRate * rDeltaTime;

      DeltaPlasticStrain       = sqrt(2.0/3.0) * rDeltaGamma;
      rEquivalentPlasticStrain = rEquivalentPlasticStrainOld + DeltaPlasticStrain;

      double StateFunction     =  this->mYieldSurface.CalculateStateFunction( rVariables, StateFunction );

      double alpha = 1;
      while ( fabs(StateFunction)>=Tolerance && iter<=MaxIterations)
	{
	  //Calculate Delta State Function:
	  DeltaStateFunction = this->mYieldSurface.CalculateDeltaStateFunction( rVariables, DeltaStateFunction );

	  //Calculate DeltaGamma:
	  DeltaDeltaGamma  = StateFunction/DeltaStateFunction;
	  rDeltaGamma += DeltaDeltaGamma;

	  //Update Equivalent Plastic Strain:
	  DeltaPlasticStrain       = sqrt(2.0/3.0) * rDeltaGamma;

	  //alpha = CalculateLineSearch( rVariables, alpha );

	  rEquivalentPlasticStrain = rEquivalentPlasticStrainOld + alpha * DeltaPlasticStrain;

	  //Calculate State Function:
	  StateFunction = this->mYieldSurface.CalculateStateFunction( rVariables, StateFunction );

	  iter++;
	}


      if(iter>MaxIterations)
	return false;


      return true;

      KRATOS_CATCH(" ")
    }


    bool CalculateRateIndependentReturnMapping(PlasticDataType& rVariables, MatrixType& rStressMatrix)
    {
      KRATOS_TRY

      //Set convergence parameters
      unsigned int iter    = 0;
      double Tolerance     = 1e-5;
      double MaxIterations = 50;

      //start
      double DeltaDeltaGamma    = 0;
      double DeltaStateFunction = 0;
      double DeltaPlasticStrain = 0;

      double& rEquivalentPlasticStrainOld  = this->mPreviousInternal.Variables[0];
      double& rEquivalentPlasticStrain     = rVariables.Internal.Variables[0];
      double& rDeltaGamma                  = rVariables.DeltaInternal.Variables[0];

      rEquivalentPlasticStrain = 0;
      rDeltaGamma = 1e-40;  //this can not be zero (zig-zag in the iterative loop if is zero)

      DeltaPlasticStrain       = sqrt(2.0/3.0) * rDeltaGamma;
      rEquivalentPlasticStrain = rEquivalentPlasticStrainOld + DeltaPlasticStrain;

      double StateFunction     =  rVariables.TrialStateFunction;

      double alpha = 1;
      while ( fabs(StateFunction)>=Tolerance && iter<=MaxIterations)
	{
	  //Calculate Delta State Function:
	  DeltaStateFunction = this->mYieldSurface.CalculateDeltaStateFunction( rVariables, DeltaStateFunction );

	  //Calculate DeltaGamma:
	  DeltaDeltaGamma  = StateFunction/DeltaStateFunction;
	  rDeltaGamma += DeltaDeltaGamma;

	  //Update Equivalent Plastic Strain:
	  DeltaPlasticStrain       = sqrt(2.0/3.0) * rDeltaGamma;

	  //alpha = CalculateLineSearch( rVariables, alpha );

	  rEquivalentPlasticStrain = rEquivalentPlasticStrainOld + alpha * DeltaPlasticStrain;

	  //Calculate State Function:
	  StateFunction = this->mYieldSurface.CalculateStateFunction( rVariables, StateFunction );

	  iter++;
	}


      if(iter>MaxIterations)
	return false;


      return true;

      KRATOS_CATCH(" ")
    }


    // implex protected methods
    void CalculateImplexReturnMapping(PlasticDataType& rVariables, MatrixType& rStressMatrix) override
    {
      KRATOS_TRY

      const ModelDataType& rModelData = rVariables.GetModelData();

      const double& rEquivalentPlasticStrain     = rVariables.GetInternalVariables()[0];
      const double& rEquivalentPlasticStrainOld  = this->mPreviousInternal.Variables[0];

      double& rDeltaGamma                        = rVariables.Internal.Variables[0];
      const double& rDeltaTime                   = rModelData.GetProcessInfo()[DELTA_TIME];


      //1.-Computation of the plastic Multiplier
      rDeltaGamma = sqrt(3.0*0.5) * ( rEquivalentPlasticStrain - rEquivalentPlasticStrainOld );

      //2.- Update back stress, plastic strain and stress
      this->UpdateStressConfiguration(rVariables,rStressMatrix);

      //3.- Calculate thermal dissipation and delta thermal dissipation
      if( rDeltaGamma > 0 ){

	const Properties& rProperties  = rModelData.GetProperties();

	const double& rPlasticStrainRate =rProperties[PLASTIC_STRAIN_RATE];

	double MaterialDeltaPlasticStrain = rPlasticStrainRate * rDeltaTime;

	//plastic rate region on
	rVariables.State().Set(ConstitutiveModelData::PLASTIC_RATE_REGION,true);
	rVariables.RateFactor=1;

	if( sqrt(2.0/3.0) * rDeltaGamma < MaterialDeltaPlasticStrain ){
	  //plastic rate region off

	  rVariables.State().Set(ConstitutiveModelData::PLASTIC_RATE_REGION,false);
	  rVariables.RateFactor=0;
	}

	this->CalculateImplexThermalDissipation( rVariables );

	rVariables.State().Set(ConstitutiveModelData::PLASTIC_REGION,true);

      }
      else{

	this->mThermalVariables.PlasticDissipation = 0;
	this->mThermalVariables.DeltaPlasticDissipation = 0;
      }

      KRATOS_CATCH(" ")
    }


    double& CalculateLineSearch(PlasticDataType& rVariables, double& rAlpha)
    {
      KRATOS_TRY

      //Set convergence parameters
      unsigned int iter    = 0;
      double MaxIterations = 10;

      //start preserve initial variables
      double& rDeltaGamma  = rVariables.DeltaInternal.Variables[0];
      double DeltaGamma = sqrt(2.0/3.0) * rDeltaGamma;

      double StateFunction = this->mYieldSurface.CalculateStateFunction( rVariables, StateFunction );
      double R0 = sqrt(2.0/3.0) * rDeltaGamma * StateFunction;

      //double Residual0 = StateFunction;

      double& rEquivalentPlasticStrain     = rVariables.GetInternalVariables()[0];
      const double& rEquivalentPlasticStrainOld  = this->mPreviousInternal.Variables[0];

      rEquivalentPlasticStrain = rEquivalentPlasticStrainOld + sqrt(2.0/3.0) * rDeltaGamma;
      StateFunction = this->mYieldSurface.CalculateStateFunction( rVariables, StateFunction );

      double R1 = sqrt(2.0/3.0) * rDeltaGamma * StateFunction;

      rAlpha = 1;

      if(R0*R1<0){

	double R2 = R1;

	if(fabs(R1)<fabs(R0))
	  R2=R0;
	double R0start = R0;


	double nabla = 0;
	double delta = 1;

	//if( Residual0 < StateFunction ){

	while ( fabs(R2/R0start)>0.3 && iter<MaxIterations && (R1*R0)<0 && fabs(R1)>1e-7 && fabs(R0)>1e-7 )
	  {

	    rAlpha = 0.5*(nabla+delta);

	    rDeltaGamma *= rAlpha;

	    rEquivalentPlasticStrain  = rEquivalentPlasticStrainOld + sqrt(2.0/3.0) * rDeltaGamma;

	    StateFunction = this->mYieldSurface.CalculateStateFunction( rVariables, StateFunction );

	    R2 = sqrt(2.0/3.0) * rDeltaGamma * StateFunction;

	    rDeltaGamma /= rAlpha;


	    if(R2*R1<0){
	      nabla = rAlpha;
	      R0 = R2;
	    }
	    else if(R2*R0<0){
	      delta = rAlpha;
	      R1 = R2;
	    }
	    else{
	      break;
	    }

	    iter++;
	  }
	//}

      }

      rDeltaGamma = DeltaGamma;

      if( rAlpha != 1)
	std::cout<<" [ LINE SEARCH: (Iterations: "<<iter<<", rAlpha: "<<rAlpha<<") ] "<<std::endl;


      if(rAlpha>1 || rAlpha<=0)
	rAlpha=1;

      return rAlpha;

      KRATOS_CATCH(" ")
    }

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
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, DerivedType )
    }

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, DerivedType )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class NonLinearRateDependentPlasticityModel

  ///@}

  ///@name Type Definitions
  ///@{

  ///@}
  ///@name Input and output
  ///@{

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NON_LINEAR_RATE_DEPENDENT_PLASTICITY_MODEL_H_INCLUDED  defined
