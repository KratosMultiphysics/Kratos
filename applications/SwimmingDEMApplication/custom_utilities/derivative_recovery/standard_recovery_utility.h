//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Guillermo Casas (gcasas@cimne.upc.edu)
//
//

#ifndef KRATOS_STANDARD_RECOVERY_UTILITY_H
#define KRATOS_STANDARD_RECOVERY_UTILITY_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "derivative_recovery_utility.h"
#include "includes/kratos_parameters.h"

// Application includes


namespace Kratos
{
///@addtogroup SwimmingDEMApplication
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

class KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility : public DerivativeRecoveryUtility
{
public:
    ///@name Type Definitions
    ///@{
    using DoubleVarType = DerivativeRecoveryUtility::DoubleVarType;
    using ComponentVarType = DerivativeRecoveryUtility::ComponentVarType;
    using TensorComponentVarType = DerivativeRecoveryUtility::TensorComponentVarType;
    using ArrayVarType = DerivativeRecoveryUtility::ArrayVarType;

    /// Pointer definition of StandardRecoveryUtility
    KRATOS_CLASS_POINTER_DEFINITION(StandardRecoveryUtility);

    ///@}
    ///@name Life Cycle
    ///@{
    /// Constructor with Kratos parameters.
    StandardRecoveryUtility(
        ModelPart& rModelPart,
        Parameters rParameters);

    /// Constructor with Kratos model
    StandardRecoveryUtility(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    ~StandardRecoveryUtility() override {}

    ///@}
    ///@name Operators
    ///@{

    void Initialize() override;

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
        buffer << "StandardRecoveryUtility" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "StandardRecoveryUtility";}



    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:

void AddPartialTimeDerivative(const DoubleVarType& rVariable, const DoubleVarType& rTimeDerivativeVariable) override;

// void AddPartialTimeDerivative(const DoubleVarType& rScalarVariable, const DoubleVarType& rTimeDerivativeVariable) override;

// void AddPartialTimeDerivative(const ComponentVarType& rScalarComponent, const DoubleVarType& rTimeDerivativeVariable) override;

void AddPartialTimeDerivative(const ArrayVarType& rVectorVariable, const ArrayVarType& rTimeDerivativeVariable) override;

void CalculateGradient(const DoubleVarType& rScalarVariable, const ArrayVarType& rGradientVariable) override;

void CalculateGradient(const ArrayVarType& rVectorVariable, const TensorVarType& rGradientVariable) override;
// void CalculateGradient(const VariableData& rVariable, const VariableData& rGradientVariable) override;

// // void CalculateGradient(const ArrayVarType& rVariable
// //                        const ArrayVarType& rGradient0Variable,
// //                        const ArrayVarType& rGradient1Variable,
// //                        const ArrayVarType& rGradient2Variable,) override;

// // void CalculateGradient(const ArrayVarType& rVectorVariable,
// //                        const ArrayVarType& rComponent0GradientVariable,
// //                        const ArrayVarType& rComponent1GradientVariable,
// //                        const ArrayVarType& rComponent2GradientVariable) override;

void CalculateDivergence(const ArrayVarType& rVectorVariable, const DoubleVarType& rDivergenceVariable) override;

// void CalculateLaplacian(const VariableData& rVariable, const VariableData& rLaplacianVariable) override;

// // void CalculateLaplacian(const ComponentVarType& rScalarComponent, const DoubleVarType& rLaplacianVariable) override;

// // void CalculateLaplacian(const ArrayVarType& rVectorComponent, const ArrayVarType& rLaplacianVariable) override;

void CalculateMaterialDerivative(const DoubleVarType& rVariable, const DoubleVarType& rMaterialDerivativeVariable) override;

// // void CalculateMaterialDerivative(const DoubleVarType& rScalarVariable, const DoubleVarType& rMaterialDerivativeVariable) override;

// // void CalculateMaterialDerivative(const ComponentVarType& rScalarComponent, const DoubleVarType& rMaterialDerivativeVariable) override;

void CalculateMaterialDerivative(const ArrayVarType& rScalarComponent, const ArrayVarType& rMaterialDerivativeVariable) override;

void CalculateRotational(const ArrayVarType rVectorVariable, const ArrayVarType& rRotationalVariable) override;

// // void CalculateRotational(const ArrayVarType rVectorVariable, const ArrayVarType& rRotationalVariable) override;

// void CheckDefaultsAndSettings(Parameters rParameters) override;

void CheckDefaultsAndSettings(Parameters rParameters) override;

private:
    ///@name Static Member Variables
    ///@{
    static int GetVectorizedMatrixIndex(const unsigned int i, const unsigned int j) {
        // Following the Kratos convention
        return 3 * i + j;
    }
    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    // template<class TVariable, class TDerivedVariable>
    // void AddPartialTimeDerivative(const TVariable& rVariable, const TDerivedVariable& rTimeDerivativeVariable);

    // template<>
    // void AddPartialTimeDerivative(const DoubleVarType& rVariable, const DoubleVarType& rTimeDerivativeVariable);

    // template<ComponentVarType, DoubleVarType>
    // void AddPartialTimeDerivative(const ComponentVarType& rVariable, const DoubleVarType& rTimeDerivativeVariable);

    // template<ArrayVarType, ArrayVarType>
    // void AddPartialTimeDerivative(const ArrayVarType& rVariable, const ArrayVarType& rTimeDerivativeVariable);

    template<class TScalarVariable>
    void CalculateScalarGradient(const TScalarVariable& rScalarVariable, const ArrayVarType& rGradientVariable);


    template<class TScalarVariable>
    void AddScalarPartialTimeDerivative(const TScalarVariable& rScalarVariable, const DoubleVarType& rGradientVariable);

    // template<class TVariable, class TDerivedVariable>
    // void CalculateMaterialDerivative(const TVariable& rVariable, const TDerivedVariable& rMaterialDerivativeVariable);

    template<class TScalarVariable>
    void CalculateScalarMaterialDerivative(const TScalarVariable& rScalarVariable, const DoubleVarType& rMaterialDerivativeVariable);

    // template<ArrayVarType, ArrayVarType>
    // void CalculateMaterialDerivative(const ArrayVarType& rVectorVariable,
    //                                  const ArrayVarType& rMaterialDerivativeVariable);
    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Default constructor.
    StandardRecoveryUtility() = delete;

    /// Assignment operator.
    StandardRecoveryUtility& operator=(StandardRecoveryUtility const& rOther) = delete;

    /// Copy constructor.
    StandardRecoveryUtility(StandardRecoveryUtility const& rOther) = delete;

    ///@}

}; // Class StandardRecoveryUtility

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_STANDARD_RECOVERY_UTILITY_H
