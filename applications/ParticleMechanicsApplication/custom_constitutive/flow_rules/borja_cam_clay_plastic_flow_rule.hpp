//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


#if !defined(KRATOS_BORJA_CAM_CLAY_PLASTIC_FLOW_RULE_H_INCLUDED )
#define      KRATOS_BORJA_CAM_CLAY_PLASTIC_FLOW_RULE_H_INCLUDED

// System includes
#include <cmath>

// External includes

// Project includes
#include "custom_constitutive/flow_rules/MPM_flow_rule.hpp"


namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

//struct MCStressInvariants {

//double MeanStress;
//double J2InvSQ;
//double LodeAngle;

//};

//struct MCSmoothingConstants {

//double A;
//double B;

//};
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
class BorjaCamClayPlasticFlowRule
    :public MPMFlowRule
{



public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NonLinearAssociativePlasticFlowRule
    KRATOS_CLASS_POINTER_DEFINITION( BorjaCamClayPlasticFlowRule );

    struct MaterialParameters
    {
        double PreconsolidationPressure;
        double PlasticHardeningModulus;
        double ConsistencyParameter;

    public:
        void PrintInfo()
        {
            KRATOS_INFO("MPMFlowRule.MaterialParameters") << "PreconsolidationPressure = " <<  PreconsolidationPressure  << std::endl;
            KRATOS_INFO("MPMFlowRule.MaterialParameters") << "PlasticHardeningModulus  = " <<  PlasticHardeningModulus   << std::endl;
            KRATOS_INFO("MPMFlowRule.MaterialParameters") << "ConsistencyParameter     = " <<  ConsistencyParameter      << std::endl;
        }

    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BorjaCamClayPlasticFlowRule();

    /// Initialization constructor.
    BorjaCamClayPlasticFlowRule(YieldCriterionPointer pYieldCriterion);

    /// Copy constructor.
    BorjaCamClayPlasticFlowRule(BorjaCamClayPlasticFlowRule const& rOther);

    /// Assignment operator.
    BorjaCamClayPlasticFlowRule& operator=(BorjaCamClayPlasticFlowRule const& rOther);

    // CLONE
    MPMFlowRule::Pointer Clone() const override;

    /// Destructor.
    ~BorjaCamClayPlasticFlowRule() override;

    bool CalculateReturnMapping( RadialReturnVariables& rReturnMappingVariables, const Matrix& rIncrementalDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen) override;

    bool UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables ) override;

    Matrix GetElasticLeftCauchyGreen(RadialReturnVariables& rReturnMappingVariables) override;

    unsigned int GetPlasticRegion() override;

    void ComputeElastoPlasticTangentMatrix(const RadialReturnVariables& rReturnMappingVariables, const Matrix& rNewElasticLeftCauchyGreen, const double& alfa, Matrix& rConsistMatrix) override;

    void CalculatePrincipalStressTrial(const RadialReturnVariables& rReturnMappingVariables, const Matrix& rNewElasticLeftCauchyGreen, Matrix& rStressMatrix) override;

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
    Vector mElasticPrincipalStrain;
    Vector mPlasticPrincipalStrain;

    Vector mPrincipalStressUpdated;

    unsigned int mRegion;
    bool mLargeStrainBool;

    MaterialParameters mMaterialParameters;

    double mInitialVolumetricStrain;

    double mStateFunction;
    Vector mStateFunctionFirstDerivative ;
    Vector mStateFunctionSecondDerivative;

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
    void InitializeMaterial(YieldCriterionPointer& pYieldCriterion, HardeningLawPointer& pHardeningLaw, const Properties& rProp) override;

    void InitializeMaterialParameters();

    void CalculatePrincipalStressVector(const Vector& rPrincipalStrain, Vector& rPrincipalStress);

    void CalculateMeanStress(const double& rVolumetricStrain, const double& rDeviatoricStrain, double& rMeanStress);

    void CalculateDeviatoricStress(const double& rVolumetricStrain, const Vector& rDeviatoricStrainVector, Vector& rDeviatoricStress);

    void CalculatePrincipalStrainFromStrainInvariants(Vector& rPrincipalStrain, const double& rVolumetricStrain, const double& rDeviatoricStrain, const Vector& rDirectionVector);

    void CalculateStrainInvariantsFromPrincipalStrain(const Vector& rPrincipalStrain, double& rVolumetricStrain, double& rDeviatoricStrain, Vector& rDeviatoricStrainVector);

    bool CalculateConsistencyCondition(RadialReturnVariables& rReturnMappingVariables, Vector& rPrincipalStress, Vector& rPrincipalStrain, unsigned int& region, Vector& rPrincipalStressUpdated);

    void CalculateLHSMatrix(Matrix& rLHSMatrix, const Vector& rPrincipalStressVector, const Vector& rUnknownVector, const double& rK_p);

    void CalculateHessianMatrix_2x2(Matrix& rHessianMatrix);

    void ComputeElasticMatrix_2X2(const Vector& rPrincipalStressVector, const double& rVolumetricStrain, const double& rDeviatoricStrain, Matrix& rElasticMatrix);

    void ComputePlasticMatrix_2X2(const Vector& rPrincipalStressVector, const double& rVolumetricStrain, const double& rDeviatoricStrain, const Matrix& rElasticMatrix, Matrix& rPlasticMatrix);

    void ReturnStressFromPrincipalAxis(const Matrix& rEigenVectors, const Vector& rPrincipalStress, Matrix& rStressMatrix);

    void CalculateTransformationMatrix(const Matrix& rMainDirection, Matrix& rA);

    void UpdateStateVariables(const Vector rPrincipalStress, const double rAlpha = 0.0, const double rConsistencyParameter = 0.0);

    double GetPI();

    //virtual void GetPrincipalStressAndStrain(Vector& PrincipalStresses, Vector& PrincipalStrains);
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

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class NonLinearAssociativePlasticFlowRule

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
// 				    NonLinearAssociativePlasticFlowRule& rThis);

// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
// 				    const NonLinearAssociativePlasticFlowRule& rThis)
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

#endif // KRATOS_BORJA_CAM_CLAY_PLASTIC_FLOW_RULE_H_INCLUDED  defined
