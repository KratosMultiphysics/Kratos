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

// class KRATOS_API(SWIMMING_DEM_APPLICATION) DifferentiableVariablePairWrapperBase
// {
//     public:

//     KRATOS_CLASS_POINTER_DEFINITION(DifferentiableVariablePairWrapperBase);
//     virtual ~DifferentiableVariablePairWrapperBase(){}
//     DifferentiableVariablePairWrapperBase(const VariableData& rVariable,
//                                           const VariableData& rDifferentiatedVariable)
//     {
//         // std::set<std::string> possible_types = {"scalar", "scalar_component", "vector"};

//         // if (possible_types.find(mVariableType) == possible_types.end()){
//         //     std::ostringstream oss;
//         //     oss << "The provided variable type (" << mVariableType << ") ";

//         //     oss << "is not among the available options (";

//         //     for (auto type : possible_types){
//         //         oss << " " << type << ",";
//         //     }

//         //     oss << ").";
//         //     const auto error_message = oss.str();
//         //     KRATOS_THROW_ERROR(std::invalid_argument, error_message, "");
//         // }
//     }

//     template<class TVariable>
//     const Variable<TVariable>& GetVariableToDerivate() const;

//     template<class TDerivatedVariable>
//     const Variable<TDerivatedVariable>& GetDerivativeVariable() const;
// };

// template <class TVariable, class TDerivatedVariable>
// class DifferentiableVariablePairWrapper : public DifferentiableVariablePairWrapperBase
// {
//     public:

//     KRATOS_CLASS_POINTER_DEFINITION(DifferentiableVariablePairWrapper);

//     DifferentiableVariablePairWrapper(const Variable<TVariable>& rVariable,
//                                       const Variable<TDerivatedVariable>& rDifferentiatedVariable):
//         mVariable(rVariable), mDerivativeVariable(rDifferentiatedVariable){}

//     const Variable<TVariable>& GetVariableToDerivate() const
//     {
//         return mVariable;
//     }

//     const Variable<TDerivatedVariable>& GetDerivativeVariable() const
//     {
//         return mDerivativeVariable;
//     }

//     private:

//     const Variable<TVariable>& mVariable;
//     const Variable<TDerivatedVariable>& mDerivativeVariable;
// };

// template<class TVariable> const Variable<TVariable>& DifferentiableVariablePairWrapperBase::GetVariableToDerivate() const
// {
//     return dynamic_cast<const Variable<TVariable>&>(*this).GetVariableToDerivate();
// }

// template<class TDerivatedVariable> const Variable<TDerivatedVariable>& DifferentiableVariablePairWrapperBase::GetDerivativeVariable() const
// {
//     return dynamic_cast<const Variable<TDerivatedVariable>&>(*this).GetDerivativeVariable();
// }


class KRATOS_API(SWIMMING_DEM_APPLICATION) RecoveryVariablesContainer
{
    ///@name Type Definitions
    ///@{

    public:
    typedef array_1d<double, 3> Vector3;
    typedef Variable<Matrix> MatrixVarType;
    typedef VariableComponent<VectorComponentAdaptor<Vector3> > ComponentVarType;
    typedef std::pair<std::string, std::string> PairOfVariablesType;
    typedef std::vector<PairOfVariablesType> VariablePairsVectorType;
    typedef std::map<std::string, VariablePairsVectorType> VariablesMapType;

    // template <typename T>
    // struct TypeName
    // {
    //     static const char* Get()
    //     {
    //         const std::string name = typeid(T).name();
    //         return name;
    //     }
    // };

    //RecoveryVariablesContainer
    //RecoveryVariablesContainer

    KRATOS_CLASS_POINTER_DEFINITION(RecoveryVariablesContainer);
    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Kratos parameters.
    RecoveryVariablesContainer()
    {
        mVariables["gradient"] = VariablePairsVectorType();
        mVariables["divergence"] = VariablePairsVectorType();
        mVariables["laplacian"] = VariablePairsVectorType();
        mVariables["rotational"] = VariablePairsVectorType();
        mVariables["material_derivative"] = VariablePairsVectorType();
    }

    // template <class TVariable, class TDerivedVariable>
    // void AddRecoveryPair(const std::string OperatorName, const Variable<TVariable>& rOriginVariable, const Variable<TDerivedVariable>& rRecoveredDerivativeVariable)
    // {
    //     const auto variable_type = TypeName<TVariable>::Get();
    //     const auto new_pair = DifferentiableVariablePairWrapper<TVariable, TDerivedVariable>(rOriginVariable, rRecoveredDerivativeVariable);
    //     mVariables[OperatorName].push_back(new_pair);
    //     CheckThatVariableIsNotAlreadyInUse(OperatorName, rOriginVariable, rRecoveredDerivativeVariable);
    // }

    // void AddGradientRecoveryPair(const DoubleVarType& rOriginVariable, const DoubleVarType& rRecoveredDerivativeVariable)
    // {
    //     this->AddRecoveryPair<double, double>("gradient", rOriginVariable, rRecoveredDerivativeVariable);
    // }

    // // void AddGradientRecoveryPair(const ArrayVarType& rOriginVariable, const ArrayVarType& rRecoveredDerivativeVariable)
    // // {
    // //     this->AddRecoveryPair<ArrayVarType, DoubleVarType>("gradient", rOriginVariable, rRecoveredDerivativeVariable);
    // // }

    // void AddDivergenceRecoveryPair(const ArrayVarType& rOriginVariable, const DoubleVarType& rRecoveredDerivativeVariable)
    // {
    //     this->AddRecoveryPair<Vector3, double>("divergence", rOriginVariable, rRecoveredDerivativeVariable);
    // }

    // void AddLaplacianRecoveryPair(const DoubleVarType& rOriginVariable, const DoubleVarType& rRecoveredDerivativeVariable)
    // {
    //     this->AddRecoveryPair<double, double>("laplacian", rOriginVariable, rRecoveredDerivativeVariable);
    // }

    // void AddLaplacianRecoveryPair(const ArrayVarType& rOriginVariable, const ArrayVarType& rRecoveredDerivativeVariable)
    // {
    //     this->AddRecoveryPair<Vector3, Vector3>("laplacian", rOriginVariable, rRecoveredDerivativeVariable);
    // }

    // void AddRotationalRecoveryPair(const ArrayVarType& rOriginVariable, const ArrayVarType& rRecoveredDerivativeVariable)
    // {
    //     this->AddRecoveryPair<Vector3, Vector3>("rotational", rOriginVariable, rRecoveredDerivativeVariable);
    // }

    // void AddMaterialDerivativeRecoveryPair(const DoubleVarType& rOriginVariable, const DoubleVarType& rRecoveredDerivativeVariable)
    // {
    //     this->AddRecoveryPair<double, double>("material_derivative", rOriginVariable, rRecoveredDerivativeVariable);
    // }

    // void AddMaterialDerivativeRecoveryPair(const ArrayVarType& rOriginVariable, const ArrayVarType& rRecoveredDerivativeVariable)
    // {
    //     this->AddRecoveryPair<Vector3, Vector3>("material_derivative", rOriginVariable, rRecoveredDerivativeVariable);
    // }

    void AddRecoveryPair(const std::string OperatorName,
                         const std::string VariableName,
                         const std::string DerivativeVariableName)
    {
        const auto new_pair = PairOfVariablesType(VariableName, DerivativeVariableName);
        mVariables[OperatorName].push_back(new_pair);
        //CheckThatVariableIsNotAlreadyInUse(OperatorName, rOriginVariable, rRecoveredDerivativeVariable);
    }

    const VariablePairsVectorType& GetVariables(const std::string Operator)
    {
        return mVariables[Operator];
    }

    private:

    VariablesMapType mVariables;


    void CheckThatVariableIsNotAlreadyInUse(std::string MyOperator, const VariableData& rOriginVariable, const VariableData& rRecoveredDerivativeVariable)
    {
        KRATOS_TRY
        // bool already_used = false;
        // std::string in_use_operator;
        // std::string variable_using_it;

        // loop over all operators
        // for (auto& operator_entry : mVariables) {
        //      // loop over all pairs for the operator
        //     for (auto pair_wrapper : operator_entry.second){
        //         auto pair = pair_wrapper.GetVariablePair();
        //         if (pair.second.Name() == rRecoveredDerivativeVariable.Name()){
        //             already_used = true;
        //             in_use_operator = pair.second.Name();
        //             variable_using_it = pair.first.Name();
        //             break;
        //         }
        //     }
        // }

        // if (already_used){
        //     std::ostringstream message;
        //     message << "The variable provided for storing the " << MyOperator
        //             << " of the " << rOriginVariable.Name() << " (" << rRecoveredDerivativeVariable.Name() << ") "
        //             << "is already in use to calculate the " << in_use_operator << " of the " << variable_using_it << ".";
        //     auto error_msg = message.str();
        //     KRATOS_THROW_ERROR(std::invalid_argument, error_msg, "");
        // }

        KRATOS_CATCH("")
    }
};

// template <>
// struct RecoveryVariablesContainer::TypeName<RecoveryVariablesContainer::DoubleVarType>
// {
//     static const std::string Get()
//     {
//         const std::string name = "scalar";
//         return name;
//     }
// };

// template <>
// struct RecoveryVariablesContainer::TypeName<RecoveryVariablesContainer::ComponentVarType>
// {
//     static const std::string Get()
//     {
//         const std::string name = "scalar_component";
//         return name;
//     }
// };

// template <>
// struct RecoveryVariablesContainer::TypeName<RecoveryVariablesContainer::ArrayVarType>
// {
//     static const std::string Get()
//     {
//         const std::string name = "vector";
//         return name;
//     }
// };

class KRATOS_API(SWIMMING_DEM_APPLICATION) DerivativeRecoveryUtility
{
public:
    ///@name Type Definitions
    ///@{
    // using DoubleVarType = DoubleVarType;
    // using ArrayVarType = ArrayVarType;
    // using MatrixVarType = RecoveryVariablesContainer::MatrixVarType;
    typedef array_1d<double, 3> Vector3;
    typedef Variable<double> DoubleVarType;
    typedef VariableComponent<VectorComponentAdaptor<Vector3> > ComponentVarType;
    typedef Variable<array_1d<double, 3>> ArrayVarType;
    typedef Variable<Vector> VectorVarType;

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

virtual void AddPartialTimeDerivative(const DoubleVarType& rVariable, const DoubleVarType& rTimeDerivativeVariable){}

// virtual void AddPartialTimeDerivative(const DoubleVarType& rScalarVariable, const DoubleVarType& rTimeDerivativeVariable){}

// virtual void AddPartialTimeDerivative(const ComponentVarType& rScalarComponent, const DoubleVarType& rTimeDerivativeVariable){}

virtual void AddPartialTimeDerivative(const ArrayVarType& rVectorVariable, const ArrayVarType& rTimeDerivativeVariable){}

virtual void CalculateGradient(const DoubleVarType& rScalarVariable, const ArrayVarType& rGradientVariable){}

virtual void CalculateGradient(const ComponentVarType& rScalarComponent, const ArrayVarType& rGradientVariable){}

// // virtual void CalculateGradient(const ArrayVarType& rVectorVariable,
// //                                const ArrayVarType& rComponent0GradientVariable,
// //                                const ArrayVarType& rComponent1GradientVariable,
// //                                const ArrayVarType& rComponent2GradientVariable){}

// virtual void CalculateDivergence(const VariableData& rVectorVariable, const VariableData& rDivergenceVariable){}

virtual void CalculateDivergence(const ArrayVarType& rVectorVariable, const DoubleVarType& rDivergenceVariable){}

virtual void CalculateDivergence(const ArrayVarType& rVectorVariable, const ComponentVarType& rDivergenceVariable){}

virtual void CalculateLaplacian(const DoubleVarType& rScalarVariable, const DoubleVarType& rLaplacianVariable){}

virtual void CalculateLaplacian(const ComponentVarType& rScalarComponent, const DoubleVarType& rLaplacianVariable){}

virtual void CalculateLaplacian(const ArrayVarType& rVectorVariable, const ArrayVarType& rLaplacianVariable){}

virtual void CalculateMaterialDerivative(const DoubleVarType& rVariable, const DoubleVarType& rMaterialDerivativeVariable){}

virtual void CalculateMaterialDerivative(const ComponentVarType& rScalarComponent, const DoubleVarType& rMaterialDerivativeVariable){}

virtual void CalculateMaterialDerivative(const ArrayVarType& rVectorVariable, const ArrayVarType& rMaterialDerivativeVariable){}

virtual void CalculateRotational(const ArrayVarType rVectorVariable, const ArrayVarType& rRotationalVariable){}

virtual void CheckDefaultsAndSettings(Parameters rParameters){};

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    void CheckDefaultVariablesAreInSettings(Parameters rParameters);

    void ReadAllVariablePairs(Parameters rVariablesForRecovery);

    void ReadVariablePairs(std::string OperatorName, Parameters rVariablesForRecovery);

    void CalculateGradient(const std::string VariableName, const std::string DerivativeVariableName);

    void CalculateDivergence(const std::string VariableName, const std::string DerivativeVariableName);

    void CalculateRotational(const std::string VariableName, const std::string DerivativeVariableName);

    void CalculateMaterialDerivative(const std::string VariableName, const std::string DerivativeVariableName);

    void CalculateLaplacian(const std::string VariableName, const std::string DerivativeVariableName);


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    // template <class TVariable, class TDerivativeVariable>
    // bool CalculateDerivativeIfPossible(const std::string OperatorName, const std::string VariableName, const std::string DerivativeVariableName)
    // {
    //     KRATOS_THROW_ERROR(std::invalid_argument, "No perator corresponds to the provided combination of types: ", VariableName + ", " + DerivativeVariableName);
    //     return false;
    // }

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
