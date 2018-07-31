// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//

#ifndef ADJOINT_STRUCTURAL_RESPONSE_FUNCTION_H
#define ADJOINT_STRUCTURAL_RESPONSE_FUNCTION_H

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "structural_mechanics_application_variables.h"


namespace Kratos
{
///@addtogroup StructuralMechanicsApplication
///@{

///@name Kratos Classes
///@{

/** \brief AdjointStructuralResponseFunction
*
* This is the response base class for responses in structural mechanics.
* It is designed to be used in adjoint sensitivity analysis.
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) AdjointStructuralResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(AdjointStructuralResponseFunction);

    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> VariableComponentType;

    typedef Variable<array_1d<double, 3>> VariableWithComponentsType;

    typedef ModelPart::NodesContainerType NodesContainerType;

    typedef ModelPart::ElementsContainerType ElementsContainerType;

    typedef ModelPart::ConditionsContainerType ConditionsContainerType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    ///@}
    ///@name Pointer Definitions

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    AdjointStructuralResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings);

    /// Destructor.
    virtual ~AdjointStructuralResponseFunction();

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    ModelPart& GetModelPart();

    ModelPart& GetModelPart() const;

    virtual void Initialize();

    virtual void InitializeSolutionStep(){};

    virtual void FinalizeSolutionStep(){};

    virtual void Check();

    virtual void Clear();

    virtual void CalculateGradient(const Element& rAdjointElem,
                                   const Matrix& rAdjointMatrix,
                                   Vector& rResponseGradient,
                                   ProcessInfo& rProcessInfo);

    virtual void CalculateGradient(const Condition& rAdjointCondition,
                                   const Matrix& rAdjointMatrix,
                                   Vector& rResponseGradient,
                                   ProcessInfo& rProcessInfo);

    virtual void CalculateFirstDerivativesGradient(const Element& rAdjointElem,
                                                   const Matrix& rAdjointMatrix,
                                                   Vector& rResponseGradient,
                                                   ProcessInfo& rProcessInfo);

    virtual void CalculateFirstDerivativesGradient(const Condition& rAdjointCondition,
                                                   const Matrix& rAdjointMatrix,
                                                   Vector& rResponseGradient,
                                                   ProcessInfo& rProcessInfo);

    virtual void CalculateSecondDerivativesGradient(const Element& rAdjointElem,
                                                    const Matrix& rAdjointMatrix,
                                                    Vector& rResponseGradient,
                                                    ProcessInfo& rProcessInfo);

    virtual void CalculateSecondDerivativesGradient(const Condition& rAdjointCondition,
                                                    const Matrix& rAdjointMatrix,
                                                    Vector& rResponseGradient,
                                                    ProcessInfo& rProcessInfo);

    virtual void UpdateSensitivities();

    virtual double CalculateValue(ModelPart& rModelPart);

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    ModelPart& mrModelPart;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{
    template <typename TDataType>
    void UpdateNodalSensitivities(Variable<TDataType> const& rSensitivityVariable, Variable<TDataType> const& rOutputVariable);
    template <typename TDataType>
    void UpdateElementSensitivities(Variable<TDataType> const& rSensitivityVariable, Variable<TDataType> const& rOutputVariable);
    template <typename TDataType>
    void UpdateConditionSensitivities(Variable<TDataType> const& rSensitivityVariable, Variable<TDataType> const& rOutputVariable);

    virtual void CalculateSensitivityGradient(Element& rAdjointElem,
                                              const Variable<double>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo);

    virtual void CalculateSensitivityGradient(Condition& rAdjointCondition,
                                              const Variable<double>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo);

    virtual void CalculateSensitivityGradient(Element& rAdjointElem,
                                              const Variable<array_1d<double,3>>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo);

    virtual void CalculateSensitivityGradient(Condition& rAdjointCondition,
                                              const Variable<array_1d<double,3>>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo);

    void AssembleNodalSensitivityContribution(Variable<double> const& rSensitivityVariable,
                                              Vector const& rSensitivityVector,
                                              Element::GeometryType& rGeom);

    void AssembleNodalSensitivityContribution(Variable<array_1d<double, 3>> const& rSensitivityVariable,
                                              Vector const& rSensitivityVector,
                                              Element::GeometryType& rGeom);

    void AssembleElementSensitivityContribution(Variable<double> const& rSensitivityVariable,
                                                Vector const& rSensitivityVector,
                                                Element& rElem);

    void AssembleElementSensitivityContribution(Variable<array_1d<double, 3>> const& rSensitivityVariable,
                                                Vector const& rSensitivityVector,
                                                Element& rElem);

    void AssembleConditionSensitivityContribution(Variable<double> const& rSensitivityVariable,
                                              Vector const& rSensitivityVector,
                                              Element::GeometryType& rGeom);

    void AssembleConditionSensitivityContribution(Variable<array_1d<double, 3>> const& rSensitivityVariable,
                                              Vector const& rSensitivityVector,
                                              Element::GeometryType& rGeom);

    void ReadDesignVariables(std::vector<std::vector<Variable<double>>>& rScalarDesignVariables,
        std::vector<std::vector<Variable<array_1d<double,3>>>>& rVectorDesignVariables, Parameters DesignVariableSettings);
    ///@}

private:
    ///@name Member Variables
    ///@{

    std::string mSensitivityModelPartName;
    unsigned int mGradientMode;
    double mDelta;

    std::vector<std::vector<Variable<double>>> mNodalSensitivityScalarVariables;
    std::vector<std::vector<Variable<double>>> mElementSensitivityScalarVariables;
    std::vector<std::vector<Variable<double>>> mConditionSensitivityScalarVariables;
    std::vector<std::vector<Variable<array_1d<double,3>>>> mNodalSensitivityVectorVariables;
    std::vector<std::vector<Variable<array_1d<double,3>>>> mElementSensitivityVectorVariables;
    std::vector<std::vector<Variable<array_1d<double,3>>>> mConditionSensitivityVectorVariables;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
};

///@} // Kratos Classes

///@} //Structural Mechanics Application group

} /* namespace Kratos.*/

#endif /* KRATOS_STRUCTURAL_RESPONSE_FUNCTION defined */
