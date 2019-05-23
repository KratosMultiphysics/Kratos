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

#ifndef KRATOS_DERIVATIVE_RECOVERY_UTILITY_H
#define KRATOS_DERIVATIVE_RECOVERY_UTILITY_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "utilities/variable_utils.h"

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

class KRATOS_API(SWIMMING_DEM_APPLICATION) RecoveryVariablesContainer
{
    ///@name Type Definitions
    ///@{

    public:

    typedef array_1d<double, 3> Vector3;
    typedef Variable<Matrix> MatrixVarType;
    typedef VariableComponent<VectorComponentAdaptor<Vector3> > ComponentVarType;

    typedef std::map<std::string, std::string> OperatorsToDerivativeVariablesMapType;
    typedef std::map<std::string, OperatorsToDerivativeVariablesMapType> VariablesToRecoveriesMapType;
    typedef std::map<std::string, bool> RecoveriesStatusMapType;
    typedef std::map<std::string, RecoveriesStatusMapType> VariableRecoveriesStatusMapType;
    static const std::set<std::string> smOperators;

    KRATOS_CLASS_POINTER_DEFINITION(RecoveryVariablesContainer);
    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Kratos parameters.
    RecoveryVariablesContainer()
    {
    }

    void AddRecoveryPair(const std::string OperatorName,
                         const std::string VariableName,
                         const std::string DerivativeVariableName)
    {
        mVariables[VariableName][OperatorName] = DerivativeVariableName;
        for (const auto& op : smOperators){
            mRecoveredStatus[VariableName][op] = false;
        }
    }

    const VariablesToRecoveriesMapType& GetVariables() const
    {
        return mVariables;
    }

    bool CheckIfRecovered(const std::string& variable_name, const std::string& operator_name)
    {
        return mRecoveredStatus[variable_name][operator_name];
    }

    void SetRecoveredStatus(const std::string& variable_name,
                            const std::string& operator_name,
                            const bool status)
    {
        mRecoveredStatus[variable_name][operator_name] = status;
    }

    private:

    VariablesToRecoveriesMapType mVariables;
    VariableRecoveriesStatusMapType mRecoveredStatus;
};

class KRATOS_API(SWIMMING_DEM_APPLICATION) DerivativeRecoveryUtility
{
public:
    ///@name Type Definitions
    ///@{
    typedef array_1d<double, 3> Vector3;
    typedef array_1d<double, 9> Tensor3;
    typedef Variable<double> DoubleVarType;
    typedef VariableComponent<VectorComponentAdaptor<Vector3> > ComponentVarType;
    typedef VariableComponent<VectorComponentAdaptor<Tensor3> > TensorComponentVarType;
    typedef Variable<Vector3> ArrayVarType;
    typedef Variable<Tensor3> TensorVarType;
    using OperatorsToDerivativeVariablesMapType = RecoveryVariablesContainer::OperatorsToDerivativeVariablesMapType;


    /// Pointer definition of DerivativeRecoveryUtility
    KRATOS_CLASS_POINTER_DEFINITION(DerivativeRecoveryUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Kratos parameters.
    DerivativeRecoveryUtility(
        ModelPart& rModelPart,
        Parameters rParameters);

    /// Constructor with Kratos model
    DerivativeRecoveryUtility(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    virtual ~DerivativeRecoveryUtility() = default;

    ///@}
    ///@name Operators
    ///@{

    virtual void Initialize(){};
    virtual void InitializeRecovery(){};
    virtual void FinalizeRecovery(){};
    virtual void Recover();

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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "DerivativeRecoveryUtility" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "DerivativeRecoveryUtility";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:

    bool mStoreFullGradient;
    ModelPart& mrModelPart;
    RecoveryVariablesContainer mVariablesContainer;

    virtual void AddPartialTimeDerivative(const DoubleVarType& rScalarVariable, const DoubleVarType& rTimeDerivativeVariable){}

    virtual void AddPartialTimeDerivative(const ComponentVarType& rScalarComponent, const DoubleVarType& rTimeDerivativeVariable){}

    virtual void AddPartialTimeDerivative(const ArrayVarType& rVectorVariable, const ArrayVarType& rTimeDerivativeVariable){}

    virtual void CalculateGradient(const DoubleVarType& rScalarVariable, const ArrayVarType& rGradientVariable){}

    virtual void CalculateGradient(const ComponentVarType& rScalarComponent, const ArrayVarType& rGradientVariable){}

    virtual void CalculateGradient(const ArrayVarType& rVectorVariable, const TensorVarType& rGradientVariable){}

    virtual void CalculateDivergence(const ArrayVarType& rVectorVariable, const DoubleVarType& rDivergenceVariable){}

    virtual void CalculateDivergence(const ArrayVarType& rVectorVariable, const ComponentVarType& rDivergenceVariable){}

    virtual void CalculateRotational(const ArrayVarType rVectorVariable, const ArrayVarType& rRotationalVariable){}

    virtual void CalculateMaterialDerivative(const DoubleVarType& rVariable, const DoubleVarType& rMaterialDerivativeVariable){}

    virtual void CalculateMaterialDerivative(const ComponentVarType& rScalarComponent, const DoubleVarType& rMaterialDerivativeVariable){}

    virtual void CalculateMaterialDerivative(const ArrayVarType& rVectorVariable, const ArrayVarType& rMaterialDerivativeVariable){}

    virtual void CalculateLaplacian(const DoubleVarType& rScalarVariable, const DoubleVarType& rLaplacianVariable){}

    virtual void CalculateLaplacian(const ComponentVarType& rScalarComponent, const DoubleVarType& rLaplacianVariable){}

    virtual void CalculateLaplacian(const ComponentVarType& rScalarComponent, const ComponentVarType& rLaplacianVariable){}

    virtual void CalculateLaplacian(const ArrayVarType& rVectorVariable, const ArrayVarType& rLaplacianVariable){}

    virtual void CheckDefaultsAndSettings(Parameters rParameters){};

    void SetCurrentVariableRecovery(const std::string& variable_name)
    {
        mCurrentVariable = variable_name;
    }

    const std::string GetCurrentVariableName() const
    {
        return mCurrentVariable;
    }

    const std::string GetDerivativeVariableName(std::string operator_name) const
    {
        const auto& var_name = this->GetCurrentRecovery().at(operator_name);
        return var_name;
    }

    const OperatorsToDerivativeVariablesMapType& GetCurrentRecovery() const
    {
        const auto& recovery = mVariablesContainer.GetVariables().at(mCurrentVariable);
        return recovery;
    }

    bool MustRecover(std::string operator_name)
    {
         auto recovery = this->GetCurrentRecovery();
         return recovery.find(operator_name) != recovery.end();
    }

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{
    std::string mCurrentVariable;

    ///@}p
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    void CheckDefaultVariablesAreInSettings(Parameters rParameters);

    void ReadAllVariablePairs(Parameters rVariablesForRecovery);

    void ReadVariablePairs(std::string OperatorName, Parameters rVariablesForRecovery);

    void CalculateGradient();

    void CalculateDivergence();

    void CalculateRotational();

    void CalculateMaterialDerivative();

    void CalculateLaplacian();


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    template <class TVariable, class TDerivativeVariable>
    bool CalculateGradientIfPossible(const std::string VariableName, const std::string DerivativeVariableName)
    {
        if (KratosComponents<TVariable>::Has(VariableName)){
            const TVariable& variable = KratosComponents<TVariable>::Get(VariableName);
            if (KratosComponents<TDerivativeVariable>::Has(DerivativeVariableName)){
                const TDerivativeVariable& derivative_variable = KratosComponents<TDerivativeVariable>::Get(DerivativeVariableName);
                this->CalculateGradient(variable, derivative_variable);
            }

            else {
                KRATOS_THROW_ERROR(std::invalid_argument, "DerivativeRecoveryUtility: The derivative variable is not registered in Kratos with a type matching that of the variable to derivate: ", DerivativeVariableName);
            }

            return true;
        }

        return false;
    }

    template <class TVariable, class TDerivativeVariable>
    bool CalculateDivergenceIfPossible(const std::string VariableName, const std::string DerivativeVariableName)
    {
        if (KratosComponents<TVariable>::Has(VariableName)){
            const TVariable& variable = KratosComponents<TVariable>::Get(VariableName);
            if (KratosComponents<TDerivativeVariable>::Has(DerivativeVariableName)){
                const TDerivativeVariable& derivative_variable = KratosComponents<TDerivativeVariable>::Get(DerivativeVariableName);
                this->CalculateDivergence(variable, derivative_variable);
            }

            else {
                KRATOS_THROW_ERROR(std::invalid_argument, "DerivativeRecoveryUtility: The derivative variable is not registered in Kratos with a type matching that of the variable to derivate: ", DerivativeVariableName);
            }

            return true;
        }

        return false;
    }

    template <class TVariable, class TDerivativeVariable>
    bool CalculateRotationalIfPossible(const std::string VariableName, const std::string DerivativeVariableName)
    {
        if (KratosComponents<TVariable>::Has(VariableName)){
            const TVariable& variable = KratosComponents<TVariable>::Get(VariableName);
            if (KratosComponents<TDerivativeVariable>::Has(DerivativeVariableName)){
                const TDerivativeVariable& derivative_variable = KratosComponents<TDerivativeVariable>::Get(DerivativeVariableName);
                this->CalculateRotational(variable, derivative_variable);
            }

            else {
                KRATOS_THROW_ERROR(std::invalid_argument, "DerivativeRecoveryUtility: The derivative variable is not registered in Kratos with a type matching that of the variable to derivate: ", DerivativeVariableName);
            }

            return true;
        }

        return false;
    }

    template <class TVariable, class TDerivativeVariable>
    bool CalculateMaterialDerivativeIfPossible(const std::string VariableName, const std::string DerivativeVariableName)
    {
        if (KratosComponents<TVariable>::Has(VariableName)){
            const TVariable& variable = KratosComponents<TVariable>::Get(VariableName);
            if (KratosComponents<TDerivativeVariable>::Has(DerivativeVariableName)){
                const TDerivativeVariable& derivative_variable = KratosComponents<TDerivativeVariable>::Get(DerivativeVariableName);
                this->CalculateMaterialDerivative(variable, derivative_variable);
            }

            else {
                KRATOS_THROW_ERROR(std::invalid_argument, "DerivativeRecoveryUtility: The derivative variable is not registered in Kratos with a type matching that of the variable to derivate: ", DerivativeVariableName);
            }

            return true;
        }

        return false;
    }

    template <class TVariable, class TDerivativeVariable>
    bool CalculateLaplacianIfPossible(const std::string VariableName, const std::string DerivativeVariableName)
    {
        if (KratosComponents<TVariable>::Has(VariableName)){
            const TVariable& variable = KratosComponents<TVariable>::Get(VariableName);
            if (KratosComponents<TDerivativeVariable>::Has(DerivativeVariableName)){
                const TDerivativeVariable& derivative_variable = KratosComponents<TDerivativeVariable>::Get(DerivativeVariableName);
                this->CalculateLaplacian(variable, derivative_variable);
            }

            else {
                KRATOS_THROW_ERROR(std::invalid_argument, "DerivativeRecoveryUtility: The derivative variable is not registered in Kratos with a type matching that of the variable to derivate: ", DerivativeVariableName);
            }

            return true;
        }

        return false;
    }

    /// Default constructor.
    DerivativeRecoveryUtility() = delete;

    /// Assignment operator.
    DerivativeRecoveryUtility& operator=(DerivativeRecoveryUtility const& rOther) = delete;

    /// Copy constructor.
    DerivativeRecoveryUtility(DerivativeRecoveryUtility const& rOther) = delete;

    ///@}

}; // Class DerivativeRecoveryUtility

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_DERIVATIVE_RECOVERY_UTILITY_H
