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

    struct PlasticFactors
    {
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

      dobule DeltaGamma;
      dobule DeltaBeta;

      double LameMu_bar;
      dobule TimeStep;
    }

    struct InternalVariables
    {
        double EquivalentPlasticStrain;
        double DeltaPlasticStrain;
        
        void clear ()
        {
	  EquivalentPlasticStrain = 0;
	  DeltaPlasticStrain = 0;
	}
    }

    /// Pointer definition of FlowRule
      KRATOS_CLASS_POINTER_DEFINITION(FlowRule);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FlowRule()
    {
       KRATOS_ERROR(std::logic_error, "calling the default constructor in FlowRule ... illegal operation!!","");
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
       mpYieldCriterion   = rOther.mpYieldCriterion.Clone();
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
    
    void InitializeMaterial (YieldCriterion & rYieldCriterion, HardeningLaw & rHardeningLaw, const Properties& rProperties)
    {
      mInternalVariables.clear();

      mpYieldCriterion = &rYieldCriterion;
      mpHardeningLaw   = &rHardeningLaw;
      mpProperties     = &rProperties;

      mpHardeningLaw.InitializeMaterial(rProperties);	
      mpYieldCriterion.InitializeMaterial(rHardeningLaw);	
    };


    Properties & GetProperties()
    {
      return *mpProperties;
    };
	

    virtual void CalculateReturnMapping( RadialReturnVariables& rReturnMappingVariables, Matrix& rIsoStressMatrix )
    {
	    KRATOS_ERROR(std::logic_error, "calling the base class function in FlowRule ... illegal operation!!","");
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

    /// Turn back information as a string.
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;


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

    YieldCriterion*     mpYieldCriterion;
    HardeningLaw*       mpHardeningLaw;
    Properties*         mpProperties;
	  
    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    virtual double& CalculateNormStress ( Matrix & rStressMatrix, double& rNormStress )
    {
	    return 0;
    };

     virtual bool CalculateConstistencyCondition( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables )
    {
	    KRATOS_ERROR(std::logic_error, "calling the base class function in FlowRule ... illegal operation!!","");

	    return 0;
    };


    virtual void UpdateConfiguration( RadialReturnVariables& rReturnMappingVariables, Matrix & rIsoStressMatrix )
    {
	    KRATOS_ERROR(std::logic_error, "calling the base class function in FlowRule ... illegal operation!!","");

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
	    rSerializer.save("Properties",mpProperties);
    };

    virtual void load(Serializer& rSerializer)
    {
	    rSerializer.load("InternalVariables",mInternalVariables);
	    rSerializer.load("YieldCriterion",mpYieldCriterion);
	    rSerializer.load("HardeningLaw",mpHardeningLaw);
	    rSerializer.load("Properties",mpProperties);
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


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    FlowRule& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const FlowRule& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}

  ///@} addtogroup block



  ///@}
  ///@ Template Operations
  ///@{

  ///@}


}  // namespace Kratos.

#endif // KRATOS_FLOW_RULE_H_INCLUDED  defined 


