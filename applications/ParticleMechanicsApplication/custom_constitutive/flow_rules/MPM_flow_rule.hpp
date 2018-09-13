//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//

#if !defined(KRATOS_MPM_FLOW_RULE_H_INCLUDED )
#define  KRATOS_MPM_FLOW_RULE_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/properties.h"
#include "utilities/math_utils.h"

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
class KRATOS_API(PARTICLE_MECHANICS_APPLICATION) MPMFlowRule
{
public:
    ///@name Type Definitions
    ///@{

    typedef YieldCriterion::Pointer    YieldCriterionPointer;
    typedef HardeningLaw::Pointer        HardeningLawPointer;
    typedef const Properties*              PropertiesPointer;


    KRATOS_DEFINE_LOCAL_FLAG( IMPLEX_ACTIVE );
    KRATOS_DEFINE_LOCAL_FLAG( PLASTIC_REGION );
    KRATOS_DEFINE_LOCAL_FLAG( PLASTIC_RATE_REGION );
    KRATOS_DEFINE_LOCAL_FLAG( RETURN_MAPPING_COMPUTED );


    struct PlasticFactors
    {
        double Beta0;
        double Beta1;
        double Beta2;
        double Beta3;
        double Beta4;
        double Beta5;

        Matrix  Normal;
        Matrix  Dev_Normal;
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

        void print()
        {
            std::cout<<" Internal Thermal Variables "<<std::endl;
            std::cout<<" PlasticDissipation: "<<PlasticDissipation<<std::endl;
            std::cout<<" DeltaPlasticDissipation: "<<DeltaPlasticDissipation<<std::endl;
        }

    private:

        friend class Serializer;

        // A private default constructor necessary for serialization

        void save(Serializer& rSerializer) const
        {
            rSerializer.save("PlasticDissipation",PlasticDissipation);
            rSerializer.save("DeltaPlasticDissipation",DeltaPlasticDissipation);
        };

        void load(Serializer& rSerializer)
        {
            rSerializer.load("PlasticDissipation",PlasticDissipation);
            rSerializer.load("DeltaPlasticDissipation",DeltaPlasticDissipation);
        };

    };


    struct RadialReturnVariables
    {
        Flags  Options;

        double NormIsochoricStress;
        double TrialStateFunction;
        double TraceStress;

        double DeltaGamma;
        double DeltaBeta;

        double LameMu_bar;
        double Friction;
        double Dilatancy;
        double BulkModulus;
        double DeltaTime;

        double IncrementalPlasticShearStrain;

        double Temperature;

        double CharacteristicSize;

        Matrix TrialIsoStressMatrix;

        Matrix StrainMatrix;
        Matrix MainDirections;

        ThermalVariables Thermal;

    public:

        void clear()
        {
            NormIsochoricStress = 0;
            TrialStateFunction  = 0;
            TraceStress = 0;

            DeltaGamma  = 0;
            DeltaBeta   = 0;

            LameMu_bar  = 0;
            Friction    = 0;
            Dilatancy   = 0;
            BulkModulus = 0;
            DeltaTime   = 1;
            Temperature = 0;
        }


        void initialize()
        {
            Options.Set(IMPLEX_ACTIVE,false);
            Options.Set(PLASTIC_REGION,false);
            Options.Set(PLASTIC_RATE_REGION,false);
            Options.Set(RETURN_MAPPING_COMPUTED,false);

            clear();
        }

    };

    struct InternalVariables
    {
        double EquivalentPlasticStrain;
        double DeltaPlasticStrain;
        double AccumulatedPlasticVolumetricStrain;
        double AccumulatedPlasticDeviatoricStrain;
        double DeltaPlasticVolumetricStrain;
        double DeltaPlasticDeviatoricStrain;
        Matrix  Normal;

        double LameMu_bar;
        double Friction;
        double Dilatancy;
        double BulkModulus;
        //needed in IMPLEX calculation
        double EquivalentPlasticStrainOld;


    public:

        void clear()
        {
            EquivalentPlasticStrain = 0;
            AccumulatedPlasticVolumetricStrain = 0;
            AccumulatedPlasticDeviatoricStrain = 0; 
            DeltaPlasticStrain = 0;
            DeltaPlasticVolumetricStrain = 0;
            DeltaPlasticDeviatoricStrain = 0;
            Normal.clear();
            LameMu_bar  = 0;
            Friction    = 0;
            Dilatancy   = 0;
            BulkModulus = 0;
            EquivalentPlasticStrainOld = 0;
        }


        void print()
        {
            std::cout<<" Internal Variables "<<std::endl;
            std::cout<<" EquivalentPlasticStrain: "<<EquivalentPlasticStrain<<std::endl;
            std::cout<<" DeltaPlasticStrain: "<<DeltaPlasticStrain<<std::endl;
            std::cout<<" AccumulatedPlasticVolumetricStrain: "<<AccumulatedPlasticVolumetricStrain<<std::endl;
            std::cout<<" DeltaPlasticVolumetricStrain: "<<DeltaPlasticVolumetricStrain<<std::endl;
            std::cout<<" AccumulatedPlasticDeviatoricStrain: "<<AccumulatedPlasticDeviatoricStrain<<std::endl;
            std::cout<<" DeltaPlasticDeviatoricStrain: "<<DeltaPlasticDeviatoricStrain<<std::endl;
            std::cout<<" EquivalentPlasticStrainOld: "<<EquivalentPlasticStrainOld<<std::endl;
        }

    private:

        friend class Serializer;

        // A private default constructor necessary for serialization

        void save(Serializer& rSerializer) const
        {
            rSerializer.save("EquivalentPlasticStrain",EquivalentPlasticStrain);
            rSerializer.save("DeltaPlasticStrain",DeltaPlasticStrain);
            rSerializer.save("AccumulatedPlasticVolumetricStrain",AccumulatedPlasticVolumetricStrain);
            rSerializer.save("DeltaPlasticVolumetricStrain",DeltaPlasticVolumetricStrain);
            rSerializer.save("AccumulatedPlasticDeviatoricStrain",AccumulatedPlasticDeviatoricStrain);
            rSerializer.save("DeltaPlasticDeviatoricStrain",DeltaPlasticDeviatoricStrain);
            rSerializer.save("EquivalentPlasticStrainOld",EquivalentPlasticStrainOld);
        };

        void load(Serializer& rSerializer)
        {
            rSerializer.load("EquivalentPlasticStrain",EquivalentPlasticStrain);
            rSerializer.load("DeltaPlasticStrain",DeltaPlasticStrain);
            rSerializer.load("AccumulatedPlasticVolumetricStrain",AccumulatedPlasticVolumetricStrain);
            rSerializer.load("DeltaPlasticVolumetricStrain",DeltaPlasticVolumetricStrain);
            rSerializer.load("AccumulatedPlasticDeviatoricStrain",AccumulatedPlasticDeviatoricStrain);
            rSerializer.load("DeltaPlasticDeviatoricStrain",DeltaPlasticDeviatoricStrain);           
            rSerializer.load("EquivalentPlasticStrainOld",EquivalentPlasticStrainOld);
        };
    };


    /// Pointer definition of FlowRule
    KRATOS_CLASS_POINTER_DEFINITION(MPMFlowRule);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MPMFlowRule()
    {
        //KRATOS_ERROR << std::logic_error, "calling the default constructor in FlowRule ... illegal operation!!", "" )
    };

    /// Initialization constructor.
    MPMFlowRule(YieldCriterionPointer pYieldCriterion)
        :mpYieldCriterion(pYieldCriterion)
    {
    };

    /// Copy constructor.
    MPMFlowRule(MPMFlowRule const& rOther)
        :mInternalVariables(rOther.mInternalVariables)
        ,mThermalVariables(rOther.mThermalVariables)
        ,mpYieldCriterion(rOther.mpYieldCriterion)
         //,mpHardeningLaw(rOther.mpHardeningLaw)
         //,mpProperties(rOther.mpProperties)
    {
    };


    /// Assignment operator.
    MPMFlowRule& operator=(MPMFlowRule const& rOther)
    {
        mInternalVariables = rOther.mInternalVariables;
        mThermalVariables  = rOther.mThermalVariables;
        mpYieldCriterion   = rOther.mpYieldCriterion;
        //mpHardeningLaw     = rOther.mpHardeningLaw;
        //mpProperties       = rOther.mpProperties;

        return *this;
    };

    /// Destructor.
    virtual ~MPMFlowRule() {};


    ///@}
    ///@name Operators
    ///@{


    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this flow rule
     */
    virtual MPMFlowRule::Pointer Clone() const
    {
        MPMFlowRule::Pointer p_clone(new MPMFlowRule(*this));
        return p_clone;
    }


    ///@}
    ///@name Operations
    ///@{

    virtual void InitializeMaterial (YieldCriterionPointer& pYieldCriterion, HardeningLawPointer& pHardeningLaw, const Properties& rMaterialProperties)
    {
        //set yield criterion
        mpYieldCriterion = pYieldCriterion;
        mpYieldCriterion->InitializeMaterial(pHardeningLaw, rMaterialProperties);

        //initialize material variables
        mInternalVariables.clear();
        mThermalVariables.clear();

    };

    virtual void InitializeMaterial (const Properties& rMaterialProperties)
    {

        mpYieldCriterion->GetHardeningLaw().InitializeMaterial(rMaterialProperties);

        //initialize material variables
        mInternalVariables.clear();
        mThermalVariables.clear();

    };

    const Properties & GetProperties()
    {
        return mpYieldCriterion->GetHardeningLaw().GetProperties();
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
        KRATOS_ERROR << "Calling the base class function (CalculateReturnMapping) in MPM FlowRule:: illegal operation!" << std::endl;
        return 0;
    };

    virtual bool CalculateReturnMapping( RadialReturnVariables& rReturnMappingVariables, const Matrix& rIncrementalDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen)
    {
        KRATOS_ERROR << "Calling the base class function (CalculateReturnMapping) in MPM FlowRule:: illegal operation!" << std::endl;
        return 0;
    };

    virtual void ComputeElastoPlasticTangentMatrix( const RadialReturnVariables& rReturnMappingVariables, const Matrix& rElasticLeftCauchyGreen, const double& rAlpha, Matrix& rElastoPlasticMatrix)
    {
        KRATOS_ERROR << "Calling the base class function (ComputeElastoPlasticTangentMatrix) in MPM FlowRule:: illegal operation!" << std::endl;
    };

    virtual void CalculateScalingFactors(const RadialReturnVariables& rReturnMappingVariables, PlasticFactors& rScalingFactors )
    {
        KRATOS_ERROR << "Calling the base class function (CalculateScalingFactors) in MPM FlowRule:: illegal operation!" << std::endl;
    };

    virtual bool UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables )
    {
        KRATOS_ERROR << "Calling the base class function (UpdateInternalVariables) in MPM FlowRule:: illegal operation!" << std::endl;

        return 0;
    };

    virtual Matrix GetElasticLeftCauchyGreen( RadialReturnVariables& rReturnMappingVariables )
    {
        KRATOS_ERROR << "Calling the base class function (GetElasticLeftCauchyGreen) in MPM FlowRule:: illegal operation!" << std::endl;
    };

    virtual unsigned int GetPlasticRegion()
    {
        KRATOS_ERROR << "Calling the base class function (GetPlasticRegion) in MPM FlowRule:: illegal operation!" << std::endl;
    };

    virtual void CalculatePrincipalStressTrial(const RadialReturnVariables& rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Matrix& rStressMatrix)
    {
        KRATOS_ERROR << "Calling the base class function (CalculatePrincipalStressTrial) in MPM FlowRule:: illegal operation!" << std::endl;
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

    //HardeningLawPointer     mpHardeningLaw;
    //PropertiesPointer       mpProperties;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    virtual double& CalculateStressNorm ( Matrix & rStressMatrix, double& rStressNorm )
    {
        KRATOS_ERROR << "Calling the base class function (CalculateStressNorm) in MPM FlowRule:: illegal operation!" << std::endl;
    
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
        rSerializer.save("ThermalVariables",mThermalVariables);
        rSerializer.save("YieldCriterion",mpYieldCriterion);
    };

    virtual void load(Serializer& rSerializer)
    {
        rSerializer.load("InternalVariables",mInternalVariables);
        rSerializer.load("ThermalVariables",mThermalVariables);
        rSerializer.load("YieldCriterion",mpYieldCriterion);
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

#endif // KRATOS_MPM_FLOW_RULE_H_INCLUDED  defined 


