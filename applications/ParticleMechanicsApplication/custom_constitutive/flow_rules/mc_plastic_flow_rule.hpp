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


#if !defined(KRATOS_MC_PLASTIC_FLOW_RULE_H_INCLUDED )
#define      KRATOS_MC_PLASTIC_FLOW_RULE_H_INCLUDED

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
class MCPlasticFlowRule
    :public MPMFlowRule
{



public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NonLinearAssociativePlasticFlowRule
    KRATOS_CLASS_POINTER_DEFINITION( MCPlasticFlowRule );

    // Variable material parameters which can change due to hardening
    struct MaterialParameters
    {
        double Cohesion;
        double FrictionAngle;
        double DilatancyAngle;

    public:
        void PrintInfo()
        {
            KRATOS_INFO("MPMFlowRule.MaterialParameters") << "Cohesion       = " << Cohesion       << std::endl;
            KRATOS_INFO("MPMFlowRule.MaterialParameters") << "FrictionAngle  = " << FrictionAngle  << std::endl;
            KRATOS_INFO("MPMFlowRule.MaterialParameters") << "DilatancyAngle = " << DilatancyAngle << std::endl;
        }

    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MCPlasticFlowRule();

    /// Initialization constructor.
    MCPlasticFlowRule(YieldCriterionPointer pYieldCriterion);

    /// Copy constructor.
    MCPlasticFlowRule(MCPlasticFlowRule const& rOther);

    /// Assignment operator.
    MCPlasticFlowRule& operator=(MCPlasticFlowRule const& rOther);

    // CLONE
    MPMFlowRule::Pointer Clone() const override;

    /// Destructor.
    ~MCPlasticFlowRule() override;

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
    Vector mElasticPreviousPrincipalStrain;
    Vector mPrincipalStressTrial;
    Vector mPrincipalStressUpdated;
    unsigned int mRegion;
    bool mLargeStrainBool;
    double mEquivalentPlasticStrain;

    MaterialParameters mMaterialParameters;

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

    virtual void ComputePlasticHardeningParameter(const Vector& rHenckyStrainVector, const double& rAlpha, double& rH);

    bool CalculateConsistencyCondition(RadialReturnVariables& rReturnMappingVariables, Vector& rPrincipalStress, Vector& rPrincipalStrain, unsigned int& region, Vector& rPrincipalStressUpdated);

    void ComputeElasticMatrix_3X3(const RadialReturnVariables& rReturnMappingVariables, Matrix& rElasticMatrix);

    void CalculateDepSurface(Matrix& rElasticMatrix, Vector& rFNorm, Vector& rGNorm, Matrix& rAuxDep);

    void CalculateDepLine(Matrix& rInvD, Vector& rFNorm, Vector& rGNorm, Matrix& rAuxDep);

    void CalculateElastoPlasticMatrix(const RadialReturnVariables& rReturnMappingVariables, unsigned int& rRegion, Vector& DiffPrincipalStress, Matrix& rDep);

    void ReturnStressFromPrincipalAxis(const Matrix& rEigenVectors, const Vector& rPrincipalStress, Matrix& rStressMatrix);


    void CalculateInverseElasticMatrix(const RadialReturnVariables& rReturnMappingVariables, Matrix& rInverseElasticMatrix);

    void CalculateElasticMatrix(const RadialReturnVariables& rReturnMappingVariables, Matrix& rElasticMatrix);

    void CalculateModificationMatrix(const RadialReturnVariables& rReturnMappingVariables, Matrix& rAuxT, Matrix& rInvAuxT);

    void CalculateTransformationMatrix(const Matrix& rMainDirection, Matrix& rA);


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

#endif // KRATOS_MATSUOKA_NAKAI_PLASTIC_FLOW_RULE_H_INCLUDED  defined
