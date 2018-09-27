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

#ifndef ADJOINT_POSTPROCESS_H
#define ADJOINT_POSTPROCESS_H

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "structural_mechanics_application_variables.h"
#include "adjoint_structural_response_function.h"


namespace Kratos
{
///@addtogroup StructuralMechanicsApplication
///@{

///@name Kratos Classes
///@{

/** \brief AdjointPostprocess
*
* This class is responsible to perform the post-processing step in which the
* sensitvities are computed.
* It is designed to be used in adjoint sensitivity analysis.
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) AdjointPostprocess
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(AdjointPostprocess);

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
    AdjointPostprocess(ModelPart& rModelPart, AdjointStructuralResponseFunction& rResponseFunction, Parameters ResponseSettings);

    /// Destructor.
    virtual ~AdjointPostprocess();

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

    virtual void UpdateSensitivities();

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    ModelPart& mrModelPart;
    AdjointStructuralResponseFunction& mrResponseFunction;

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

#endif /* KRATOS_ADJOINT_POSTPROCESS defined */
