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

void AddPartialTimeDerivative(const ComponentVarType& rScalarComponent, const DoubleVarType& rTimeDerivativeVariable) override;

void AddPartialTimeDerivative(const ArrayVarType& rVectorVariable, const ArrayVarType& rTimeDerivativeVariable) override;

void CalculateGradient(const DoubleVarType& rScalarVariable, const ArrayVarType& rGradientVariable) override;

void CalculateGradient(const ComponentVarType& rScalarVariable, const ArrayVarType& rGradientVariable) override;

void CalculateGradient(const ArrayVarType& rVectorVariable, const TensorVarType& rGradientVariable) override;

void CalculateDivergence(const ArrayVarType& rVectorVariable, const DoubleVarType& rDivergenceVariable) override;

void CalculateDivergence(const ArrayVarType& rVectorVariable, const ComponentVarType& rDivergenceVariable) override;

void CalculateRotational(const ArrayVarType rVectorVariable, const ArrayVarType& rRotationalVariable) override;

void CalculateMaterialDerivative(const DoubleVarType& rVariable, const DoubleVarType& rMaterialDerivativeVariable) override;

void CalculateMaterialDerivative(const ComponentVarType& rVariable, const DoubleVarType& rMaterialDerivativeVariable) override;

void CalculateMaterialDerivative(const ArrayVarType& rScalarComponent, const ArrayVarType& rMaterialDerivativeVariable) override;

void CalculateLaplacian(const DoubleVarType& rVariable, const DoubleVarType& rLaplacianVariable) override;

void CalculateLaplacian(const ComponentVarType& rScalarComponent, const DoubleVarType& rLaplacianVariable) override;

void CalculateLaplacian(const ArrayVarType& rVectorComponent, const ArrayVarType& rLaplacianVariable) override;

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

    template<class TScalarVariable>
    void AddScalarPartialTimeDerivative(const TScalarVariable& rScalarVariable, const DoubleVarType& rGradientVariable);

    template<class TScalarVariable>
    void CalculateScalarGradient(const TScalarVariable& rScalarVariable, const ArrayVarType& rGradientVariable);

    template<class TScalarVariable>
    void CalculateDivergenceAsScalar(const ArrayVarType& rVariable, const TScalarVariable& rDivergenceVariable);

    template<class TScalarVariable>
    void CalculateScalarMaterialDerivative(const TScalarVariable& rScalarVariable, const DoubleVarType& rMaterialDerivativeVariable);

    template<class TScalarVariable>
    void CalculateScalarLaplacian(const TScalarVariable& rScalarVariable, const DoubleVarType& rLaplacianVariable);

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
