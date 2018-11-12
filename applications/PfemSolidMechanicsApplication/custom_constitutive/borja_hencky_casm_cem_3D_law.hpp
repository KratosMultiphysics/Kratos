//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                    LHauser $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                 October 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined (KRATOS_BORJA_HENCKY_CASM_CEM_PLASTIC_3D_LAW_H_INCLUDED)
#define       KRATOS_BORJA_HENCKY_CASM_CEM_PLASTIC_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/non_linear_hencky_plastic_3D_law.hpp"
#include "custom_constitutive/custom_flow_rules/borja_casm_cem_explicit_plastic_flow_rule.hpp"
#include "custom_constitutive/custom_yield_criteria/casm_cem_yield_criterion.hpp"
#include "custom_constitutive/custom_hardening_laws/casm_cem_hardening_law.hpp"


namespace Kratos
{
/**
 * Defines a hyperelastic-plastic isotropic constitutive law J2 in 3D 
 * With stress split in an isochoric and volumetric parts
 * This material law is defined by the parameters needed by the yield criterion:

 * The functionality is limited to large displacements 
 */



class BorjaHenckyCasmCemPlastic3DLaw 
  : public NonLinearHenckyElasticPlastic3DLaw

{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;

    typedef FlowRule::Pointer            				     FlowRulePointer;
    typedef YieldCriterion::Pointer   				 YieldCriterionPointer;
    typedef HardeningLaw::Pointer        				 HardeningLawPointer;
    typedef Properties::Pointer           				 PropertiesPointer;
    typedef HardeningLaw::PlasticVariables      PlasticVariablesType;

    /**
     * Counted pointer of BorjaHenckyCasmCemPlastic3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION( BorjaHenckyCasmCemPlastic3DLaw );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    BorjaHenckyCasmCemPlastic3DLaw();


    BorjaHenckyCasmCemPlastic3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw); 

    /**
     * Copy constructor.
     */
    BorjaHenckyCasmCemPlastic3DLaw (const BorjaHenckyCasmCemPlastic3DLaw& rOther);


    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const;

    /**
     * Destructor.
     */
    virtual ~BorjaHenckyCasmCemPlastic3DLaw();

    /**
     * Operators
     */

		virtual double& GetValue( const Variable<double>& rThisVariable, double& rValue );

		virtual void SetValue( const Variable<double>& rThisVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo );

		virtual void SetValue( const Variable<Vector>& rThisVarialbe, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo );

		void SetPlasticVariables( const double& rInitialPreconPressure, const double& rInitialBonding); 
		void GetHardeningParameters(double& rPreconPressure, double& rBonding); 
		const double GetBonding(); 
		const double GetPreconPressure();

		int Check( const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo); 
    /**
     * Operations needed by the base class:
     */


    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param props
     * @param geom
     * @param CurrentProcessInfo
     * @return
     */
    //int Check(const Properties& rProperties, const GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo);



    /**
     * Input and output
     */
    /**
     * Turn back information as a string.
     */
    //virtual String Info() const;
    /**
     * Print information about this object.
     */
    //virtual void PrintInfo(std::ostream& rOStream) const;
    /**
     * Print object's data.
     */
    //virtual void PrintData(std::ostream& rOStream) const;

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


    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, NonLinearHenckyElasticPlastic3DLaw )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, NonLinearHenckyElasticPlastic3DLaw )
    }



}; // Class BorjaHenckyCasmCemPlastic3DLaw
}  // namespace Kratos.
#endif // KRATOS_BORJA_HENCKY_CASM_CEM_PLASTIC_3D_LAW_H_INCLUDED

