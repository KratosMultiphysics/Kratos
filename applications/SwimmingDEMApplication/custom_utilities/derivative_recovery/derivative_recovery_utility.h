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

class KRATOS_API(SWIMMING_DEM_APPLICATION) RecoveryVariablesContainer
{
    ///@name Type Definitions
    ///@{
    typedef std::pair<VariableData, VariableData> PairOfVariablesType;
    typedef std::vector<PairOfVariablesType> VariablesPairsVectorType;
    typedef std::map<std::string, VariablesPairsVectorType> VariablesMapType;

    public:
    KRATOS_CLASS_POINTER_DEFINITION(RecoveryVariablesContainer);
    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Kratos parameters.
    RecoveryVariablesContainer()
    {
        mVariables["gradient"] = VariablesPairsVectorType();
        mVariables["divergence"] = VariablesPairsVectorType();
        mVariables["Laplacian"] = VariablesPairsVectorType();
        mVariables["rotational"] = VariablesPairsVectorType();
        mVariables["material derivative"] = VariablesPairsVectorType();
    }

    void AddRecoveryPair(const std::string OperatorName, const VariableData& rOriginVariable, const VariableData& rRecoveredDerivativeVariable)
    {
        auto new_pair = std::make_pair(rOriginVariable, rRecoveredDerivativeVariable);
        mVariables[OperatorName].push_back(new_pair);
        CheckThatVariableIsNotAlreadyInUse(OperatorName, rOriginVariable, rRecoveredDerivativeVariable);
    }

    void AddGradientRecoveryPair(const VariableData& rOriginVariable, const VariableData& rRecoveredDerivativeVariable)
    {
        this->AddRecoveryPair("gradient", rOriginVariable, rRecoveredDerivativeVariable);
    }

    void AddDivergenceRecoveryPair(const VariableData& rOriginVariable, const VariableData& rRecoveredDerivativeVariable)
    {
        this->AddRecoveryPair("divergence", rOriginVariable, rRecoveredDerivativeVariable);
    }

    void AddLaplacianRecoveryPair(const VariableData& rOriginVariable, const VariableData& rRecoveredDerivativeVariable)
    {
        this->AddRecoveryPair("Laplacian", rOriginVariable, rRecoveredDerivativeVariable);
    }

    void AddRotationalRecoveryPair(const VariableData& rOriginVariable, const VariableData& rRecoveredDerivativeVariable)
    {
        this->AddRecoveryPair("rotational", rOriginVariable, rRecoveredDerivativeVariable);
    }

    void AddMaterialDerivativeRecoveryPair(const VariableData& rOriginVariable, const VariableData& rRecoveredDerivativeVariable)
    {
        this->AddRecoveryPair("material derivative", rOriginVariable, rRecoveredDerivativeVariable);
    }

    const VariablesPairsVectorType& GetVariables(const std::string Operator)
    {
        return mVariables["Operator"];
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
        for (auto& operator_entry : mVariables) {
             // loop over all pairs for the operator
            for (auto pair : operator_entry.second){
                if (pair.second.Name() == rRecoveredDerivativeVariable.Name()){
                    already_used = true;
                    in_use_operator = pair.second.Name();
                    variable_using_it = pair.first.Name();
                    break;
                }
            }
        }

        if (already_used){
            std::ostringstream message;
            message << "The variable provided for storing the " << MyOperator
                    << " of the " << rOriginVariable.Name() << " (" << rRecoveredDerivativeVariable.Name() << ") "
                    << "is already in use to calculate the " << in_use_operator << " of the " << variable_using_it << ".";
            auto error_msg = message.str();
            KRATOS_THROW_ERROR(std::invalid_argument, error_msg, "");
        }

        KRATOS_CATCH("")
    }
};

class KRATOS_API(SWIMMING_DEM_APPLICCheckThatVariableIsNotAlreadyInUseATION) DerivativeRecoveryUtility
{
public:
    ///@name Type Definitions
    ///@{
    typedef Variable<double> ScalarVariableType;
    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > ComponentVariableType;
    typedef Variable<array_1d<double, 3> > VectorVariableType;
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

virtual void AddPartialTimeDerivative(const VariableData& rScalarVariable, const VariableData& rTimeDerivativeVariable){}

// virtual void AddPartialTimeDerivative(const ScalarVariableType& rScalarVariable, const ScalarVariableType& rTimeDerivativeVariable){}

// virtual void AddPartialTimeDerivative(const ComponentVariableType& rScalarComponent, const ScalarVariableType& rTimeDerivativeVariable){}

// virtual void AddPartialTimeDerivative(const VectorVariableType& rVectorVariable, const VectorVariableType& rTimeDerivativeVariable){}

virtual void CalculateGradient(const VariableData& rScalarVariable, const VariableData& rGradientVariable){}

// virtual void CalculateGradient(const ScalarVariableType& rScalarVariable, const VectorVariableType& rGradientVariable){}

// virtual void CalculateGradient(const ComponentVariableType& rScalarComponent, const VectorVariableType& rGradientVariable){}

// virtual void CalculateGradient(const VectorVariableType& rVectorVariable,
//                                const VectorVariableType& rComponent0GradientVariable,
//                                const VectorVariableType& rComponent1GradientVariable,
//                                const VectorVariableType& rComponent2GradientVariable){}

virtual void CalculateDivergence(const VariableData& rVectorVariable, const VariableData& rDivergenceVariable){}

// virtual void CalculateDivergence(const VectorVariableType& rVectorVariable, const ScalarVariableType& rDivergenceVariable){}

virtual void CalculateLaplacian(const VariableData& rScalarVariable, const VariableData& rLaplacianVariable){}

// virtual void CalculateLaplacian(const ScalarVariableType& rScalarVariable, const ScalarVariableType& rLaplacianVariable){}

// virtual void CalculateLaplacian(const ComponentVariableType& rScalarComponent, const ScalarVariableType& rLaplacianVariable){}

// virtual void CalculateLaplacian(const VectorVariableType& rVectorVariable, const VectorVariableType& rLaplacianVariable){}

virtual void CalculateMaterialDerivative(const VariableData& rScalarVariable, const VariableData& rMaterialDerivativeVariable){}

// virtual void CalculateMaterialDerivative(const ScalarVariableType& rScalarVariable, const ScalarVariableType& rMaterialDerivativeVariable){}

// virtual void CalculateMaterialDerivative(const ComponentVariableType& rScalarComponent, const ScalarVariableType& rMaterialDerivativeVariable){}

// virtual void CalculateMaterialDerivative(const VectorVariableType& rScalarComponent, const VectorVariableType& rMaterialDerivativeVariable){}

virtual void CalculateRotational(const VariableData rVectorVariable, const VariableData& rRotationalVariable){}

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
