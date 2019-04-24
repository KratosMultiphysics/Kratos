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
    typedef Variable<double> ScalarVariableType;
    typedef Variable<Vector3> VectorVariableType;
    typedef Variable<double> DoubleVarType;
    typedef Variable<Vector3> ArrayVarType;
    typedef Variable<Matrix> MatrixVarType;
    typedef VariableComponent<VectorComponentAdaptor<Vector3> > ComponentVariableType;
    typedef std::pair<std::string, std::string> PairOfVariablesType;
    typedef std::vector<PairOfVariablesType> VariablePairsVectorType;
    typedef std::map<std::string, VariablePairsVectorType> VariablesMapType;

    template <typename T>
    struct TypeName
    {
        static const char* Get()
        {
            const std::string name = typeid(T).name();
            return name;
        }
    };

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

    // void AddGradientRecoveryPair(const ScalarVariableType& rOriginVariable, const ScalarVariableType& rRecoveredDerivativeVariable)
    // {
    //     this->AddRecoveryPair<double, double>("gradient", rOriginVariable, rRecoveredDerivativeVariable);
    // }

    // // void AddGradientRecoveryPair(const VectorVariableType& rOriginVariable, const VectorVariableType& rRecoveredDerivativeVariable)
    // // {
    // //     this->AddRecoveryPair<VectorVariableType, ScalarVariableType>("gradient", rOriginVariable, rRecoveredDerivativeVariable);
    // // }

    // void AddDivergenceRecoveryPair(const VectorVariableType& rOriginVariable, const ScalarVariableType& rRecoveredDerivativeVariable)
    // {
    //     this->AddRecoveryPair<Vector3, double>("divergence", rOriginVariable, rRecoveredDerivativeVariable);
    // }

    // void AddLaplacianRecoveryPair(const ScalarVariableType& rOriginVariable, const ScalarVariableType& rRecoveredDerivativeVariable)
    // {
    //     this->AddRecoveryPair<double, double>("laplacian", rOriginVariable, rRecoveredDerivativeVariable);
    // }

    // void AddLaplacianRecoveryPair(const VectorVariableType& rOriginVariable, const VectorVariableType& rRecoveredDerivativeVariable)
    // {
    //     this->AddRecoveryPair<Vector3, Vector3>("laplacian", rOriginVariable, rRecoveredDerivativeVariable);
    // }

    // void AddRotationalRecoveryPair(const VectorVariableType& rOriginVariable, const VectorVariableType& rRecoveredDerivativeVariable)
    // {
    //     this->AddRecoveryPair<Vector3, Vector3>("rotational", rOriginVariable, rRecoveredDerivativeVariable);
    // }

    // void AddMaterialDerivativeRecoveryPair(const ScalarVariableType& rOriginVariable, const ScalarVariableType& rRecoveredDerivativeVariable)
    // {
    //     this->AddRecoveryPair<double, double>("material_derivative", rOriginVariable, rRecoveredDerivativeVariable);
    // }

    // void AddMaterialDerivativeRecoveryPair(const VectorVariableType& rOriginVariable, const VectorVariableType& rRecoveredDerivativeVariable)
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
        bool already_used = false;
        std::string in_use_operator;
        std::string variable_using_it;

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

template <>
struct RecoveryVariablesContainer::TypeName<RecoveryVariablesContainer::ScalarVariableType>
{
    static const std::string Get()
    {
        const std::string name = "scalar";
        return name;
    }
};

template <>
struct RecoveryVariablesContainer::TypeName<RecoveryVariablesContainer::ComponentVariableType>
{
    static const std::string Get()
    {
        const std::string name = "scalar_component";
        return name;
    }
};

template <>
struct RecoveryVariablesContainer::TypeName<RecoveryVariablesContainer::VectorVariableType>
{
    static const std::string Get()
    {
        const std::string name = "vector";
        return name;
    }
};

class KRATOS_API(SWIMMING_DEM_APPLICATION) DerivativeRecoveryUtility
{
public:
    ///@name Type Definitions
    ///@{
    using DoubleVarType = RecoveryVariablesContainer::DoubleVarType;
    using ArrayVarType = RecoveryVariablesContainer::ArrayVarType;
    using MatrixVarType = RecoveryVariablesContainer::MatrixVarType;
    using ScalarVariableType = RecoveryVariablesContainer::ScalarVariableType;
    using ComponentVariableType = RecoveryVariablesContainer::ComponentVariableType;
    using VectorVariableType = RecoveryVariablesContainer::VectorVariableType;

    /// Pointer definition of DerivativeRecoveryUtility
    KRATOS_CLASS_POINTER_DEFINITION(DerivativeRecoveryUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Kratos parameters.
    DerivativeRecoveryUtility(
        ModelPart& rModelPart,
        Parameters rParameters,
        RecoveryVariablesContainer& rVariablesContainer);

    /// Constructor with Kratos model
    DerivativeRecoveryUtility(
        Model& rModel,
        Parameters rParameters,
        RecoveryVariablesContainer& rVariablesContainer);

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
RecoveryVariablesContainer& mrVariablesContainer;

virtual void AddPartialTimeDerivative(const ScalarVariableType& rVariable, const ScalarVariableType& rTimeDerivativeVariable){}

// virtual void AddPartialTimeDerivative(const ScalarVariableType& rScalarVariable, const ScalarVariableType& rTimeDerivativeVariable){}

// virtual void AddPartialTimeDerivative(const ComponentVariableType& rScalarComponent, const ScalarVariableType& rTimeDerivativeVariable){}

virtual void AddPartialTimeDerivative(const VectorVariableType& rVectorVariable, const VectorVariableType& rTimeDerivativeVariable){}

// virtual void CalculateGradient(const VariableData& rVariable, const VariableData& rGradientVariable){}

// // virtual void CalculateGradient(const ScalarVariableType& rScalarVariable, const VectorVariableType& rGradientVariable){}

// // virtual void CalculateGradient(const ComponentVariableType& rScalarComponent, const VectorVariableType& rGradientVariable){}

// // virtual void CalculateGradient(const VectorVariableType& rVectorVariable,
// //                                const VectorVariableType& rComponent0GradientVariable,
// //                                const VectorVariableType& rComponent1GradientVariable,
// //                                const VectorVariableType& rComponent2GradientVariable){}

// virtual void CalculateDivergence(const VariableData& rVectorVariable, const VariableData& rDivergenceVariable){}

// // virtual void CalculateDivergence(const VectorVariableType& rVectorVariable, const ScalarVariableType& rDivergenceVariable){}

// virtual void CalculateLaplacian(const VariableData& rVariable, const VariableData& rLaplacianVariable){}

// // virtual void CalculateLaplacian(const ScalarVariableType& rScalarVariable, const ScalarVariableType& rLaplacianVariable){}

// // virtual void CalculateLaplacian(const ComponentVariableType& rScalarComponent, const ScalarVariableType& rLaplacianVariable){}

// // virtual void CalculateLaplacian(const VectorVariableType& rVectorVariable, const VectorVariableType& rLaplacianVariable){}


virtual void CalculateMaterialDerivative(const ScalarVariableType& rVariable, const ScalarVariableType& rMaterialDerivativeVariable){}

// // virtual void CalculateMaterialDerivative(const ScalarVariableType& rScalarVariable, const ScalarVariableType& rMaterialDerivativeVariable){}

// // virtual void CalculateMaterialDerivative(const ComponentVariableType& rScalarComponent, const ScalarVariableType& rMaterialDerivativeVariable){}

virtual void CalculateMaterialDerivative(const VectorVariableType& rScalarComponent, const VectorVariableType& rMaterialDerivativeVariable){}

// virtual void CalculateRotational(const VariableData rVectorVariable, const VariableData& rRotationalVariable){}

// virtual void CalculateRotational(const VectorVariableType rVectorVariable, const VectorVariableType& rRotationalVariable){}

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

    void CalculateMaterialDerivative(const std::string VariableName, const std::string DerivativeVariableName);


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    template <class TVariable, class TDerivedVariable>
    void AddPartialTimeDerivative(const TVariable& rVariable, const TDerivedVariable& rTimeDerivativeVariable)
    {
        KRATOS_THROW_ERROR(std::invalid_argument, "Wrong combination.", "");
    }

    template <class TVariable, class TDerivedVariable>
    void CalculateMaterialDerivative(const TVariable& rVariable, const TDerivedVariable& rTimeDerivativeVariable)
    {
        KRATOS_THROW_ERROR(std::invalid_argument, "Wrong combination.", "");
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
