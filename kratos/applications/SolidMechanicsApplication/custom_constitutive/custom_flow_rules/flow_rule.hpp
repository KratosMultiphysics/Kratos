//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_FLOW_RULE_H_INCLUDED )
#define  KRATOS_FLOW_RULE_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/properties.h"
#include "custom_constitutive/custom_yield_criteria/yield_criterion.hpp"
#include "custom_constitutive/custom_hardening_laws/hardening_law.hpp"

namespace Kratos
{
  ///@addtogroup ApplicationNameApplication
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
  class FlowRule
  {
  public:
    ///@name Type Definitions
    ///@{

    typedef YieldCriterion::Pointer    YieldCriterionPointer;
    typedef HardeningLaw::Pointer        HardeningLawPointer;
    typedef Properties::Pointer            PropertiesPointer;


    struct ProcessVariables
    {

      //it has to be changed to a flags variable
      bool ImplexActive;

      bool PlasticRegion;

      bool ReturnMappingComputed;

      bool PlasticStrainRateBehaviour;

    public:
      
      void initialize()
      {
	ImplexActive  = false;
	PlasticRegion = false;
	ReturnMappingComputed = false;
	PlasticStrainRateBehaviour = false;
      };

    };

    struct PlasticFactors
    {
      double Beta0;
      double Beta1;
      double Beta2;
      double Beta3;
      double Beta4;	   

      Matrix  Normal;
      Matrix  Dev_Normal;
    };


    struct RadialReturnVariables
    {
      Vector TrialIsoStressVector;

      double NormIsochoricStress;
      double TrialStateFunction;

      double DeltaGamma;
      double DeltaBeta;

      double LameMu_bar;
      double TimeStep;

      double Temperature;

      ProcessVariables Control;

    public:
      
      void clear()
        {
	  NormIsochoricStress = 0;
	  TrialStateFunction  = 0;

	  DeltaGamma  = 0;
	  DeltaBeta   = 0;

	  LameMu_bar  = 0;
	  TimeStep    = 0;
	  Temperature = 0;
	}
      

      void initialize()
        {
	  Control.initialize();
	  clear();
	}

    };

    struct InternalVariables
    {
        double EquivalentPlasticStrain;
        double DeltaPlasticStrain;

        //needed in IMPLEX calculation
        double EquivalentPlasticStrainOld;

        
    public:

        void clear()
        {
	  EquivalentPlasticStrain = 0;
	  DeltaPlasticStrain = 0;
	  EquivalentPlasticStrainOld = 0;
	}

    private:

      friend class Serializer;

      // A private default constructor necessary for serialization
      
      void save(Serializer& rSerializer) const
      {
	rSerializer.save("EquivalentPlasticStrain",EquivalentPlasticStrain);
	rSerializer.save("DeltaPlasticStrain",DeltaPlasticStrain);
	rSerializer.save("EquivalentPlasticStrainOld",EquivalentPlasticStrainOld);
      };

      void load(Serializer& rSerializer)
      {
	rSerializer.load("EquivalenPlasticStrain",EquivalentPlasticStrain);
	rSerializer.save("DeltaPlasticStrain",DeltaPlasticStrain);
	rSerializer.load("EquivalenPlasticStrainOld",EquivalentPlasticStrainOld);
      };
    };

    struct ThermalVariables
    {
      double PlasticDissipation;
      double DeltaPlasticDissipation;

    public:

        void clear()
        {
	  PlasticDissipation = 0;
	  DeltaPlasticDissipation = 0;
	}

    };

    /// Pointer definition of FlowRule
      KRATOS_CLASS_POINTER_DEFINITION(FlowRule);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FlowRule()
    {
      //KRATOS_ERROR(std::logic_error, "calling the default constructor in FlowRule ... illegal operation!!","");
    };

    /// Copy constructor.
    FlowRule(FlowRule const& rOther)
    :mInternalVariables(rOther.mInternalVariables)
    ,mpYieldCriterion(rOther.mpYieldCriterion)
    ,mpHardeningLaw(rOther.mpHardeningLaw)
    ,mpProperties(rOther.mpProperties)
    {
    };
	
    /// Assignment operator.
    FlowRule& operator=(FlowRule const& rOther) 
    { 
       mInternalVariables = rOther.mInternalVariables;
       mThermalVariables  = rOther.mThermalVariables;
       mpYieldCriterion   = rOther.mpYieldCriterion;
       mpHardeningLaw     = rOther.mpHardeningLaw;
       mpProperties       = rOther.mpProperties;

       return *this; 
    };

    /// Destructor.
    virtual ~FlowRule() {};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    void InitializeMaterial (YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw, const Properties& rMaterialProperties)
    {
      
      mpYieldCriterion = pYieldCriterion;
      mpHardeningLaw   = pHardeningLaw;
      mpProperties     = &rMaterialProperties;

      mpHardeningLaw->InitializeMaterial(rMaterialProperties);	
      mpYieldCriterion->InitializeMaterial(mpHardeningLaw);	

      mInternalVariables.clear();
      mThermalVariables.clear();
    };


    const Properties & GetProperties()
    {
      return *mpProperties;
    };
	
    const InternalVariables & GetInternalVariables()
    {
      return mInternalVariables;
    };
	

    const ThermalVariables & GetThermalVariables()
    {
      return mThermalVariables;
    };

    virtual bool CalculateReturnMapping( RadialReturnVariables& rReturnMappingVariables, Matrix& rIsoStressMatrix )
    {
	    KRATOS_ERROR(std::logic_error, "calling the base class function in FlowRule ... illegal operation!!","");
	    return 0;
    };

 
    virtual void CalculateScalingFactors(const RadialReturnVariables& rReturnMappingVariables, PlasticFactors& rScalingFactors )
    {
	    KRATOS_ERROR(std::logic_error, "calling the base class function in FlowRule ... illegal operation!!","");
    };
    

    virtual bool UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables )
    {
	    KRATOS_ERROR(std::logic_error, "calling the base class function in FlowRule ... illegal operation!!","");

	    return 0;
    };


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    // /// Turn back information as a string.
    // virtual std::string Info() const;

    // /// Print information about this object.
    // virtual void PrintInfo(std::ostream& rOStream) const;

    // /// Print object's data.
    // virtual void PrintData(std::ostream& rOStream) const;


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
    
    InternalVariables   mInternalVariables;
    ThermalVariables    mThermalVariables;

    YieldCriterionPointer   mpYieldCriterion;
    HardeningLawPointer       mpHardeningLaw;

    const Properties*           mpProperties;
	  
    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    virtual double& CalculateStressNorm ( Matrix & rStressMatrix, double& rStressNorm )
    {
            KRATOS_ERROR(std::logic_error, "calling the base class function in FlowRule ... illegal operation!!","");

	    return rStressNorm;
    };



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

    // A private default constructor necessary for serialization

    virtual void save(Serializer& rSerializer) const
    {
	    rSerializer.save("InternalVariables",mInternalVariables);
	    rSerializer.save("YieldCriterion",mpYieldCriterion);
	    rSerializer.save("HardeningLaw",mpHardeningLaw);
	    //rSerializer.save("Properties",mpProperties);
    };

    virtual void load(Serializer& rSerializer)
    {
	    rSerializer.load("InternalVariables",mInternalVariables);
	    rSerializer.load("YieldCriterion",mpYieldCriterion);
	    rSerializer.load("HardeningLaw",mpHardeningLaw);
	    //rSerializer.load("Properties",mpProperties);
    };

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

  }; // Class FlowRule

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  // /// input stream function
  // inline std::istream& operator >> (std::istream& rIStream,
  // 				    FlowRule& rThis);

  // /// output stream function
  // inline std::ostream& operator << (std::ostream& rOStream,
  // 				    const FlowRule& rThis)
  // {
  //   rThis.PrintInfo(rOStream);
  //   rOStream << std::endl;
  //   rThis.PrintData(rOStream);

  //   return rOStream;
  // }
  ///@}

  ///@} addtogroup block



  ///@}
  ///@ Template Operations
  ///@{

  ///@}


}  // namespace Kratos.

#endif // KRATOS_FLOW_RULE_H_INCLUDED  defined 


