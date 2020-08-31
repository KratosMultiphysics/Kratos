//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_RANS_LINE_OUTPUT_PROCESS_H_INCLUDED)
#define KRATOS_RANS_LINE_OUTPUT_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "processes/process.h"

namespace Kratos
{
///@addtogroup RANSApplication
///@{

namespace LineOutputProcessUtilities
{
///@name Type Definitions
///@{

using NodeType = ModelPart::NodeType;
using GeometryType = ModelPart::ElementType::GeometryType;
using SizeType = std::size_t;
using IndicesVector = std::vector<SizeType>;

///@}

///@name Kratos Classes
///@{

/**
 * @brief Class to get variable information
 *
 * @tparam TDataType data type of variable
 */
template <class TDataType>
class VariableDataCollector
{
public:
    /**
     * @brief Updates TDataType values starting from StartIndex
     *
     * @param rValuesVector     Output vector
     * @param rValue            Value to be added
     * @param StartIndex        Start position
     */
    static void inline AddToValuesVector(
        std::vector<double>& rValuesVector,
        const TDataType& rValue,
        const SizeType StartIndex)
    {
        KRATOS_TRY

        KRATOS_DEBUG_ERROR_IF(StartIndex >= rValuesVector.size() + rValue.size())
            << "rValuesVector size is not enough to allocate values. "
               "rValuesVector vector should have atleast "
            << rValue.size()
            << " elements from(including) StartIndex value. [ StartIndex = " << StartIndex
            << ", rValuesVector.size() = " << rValuesVector.size() << " ].\n";

        for (SizeType i = 0; i < rValue.size(); ++i) {
            rValuesVector[StartIndex + i] += rValue[i];
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Updates variable names from given variable
     *
     * @param rNamesVector      Output strings with variable names
     * @param rVariable         Variable
     * @param rValue            Example value to identify dynamic properties of Variable
     * @param StartIndex        Start position
     */
    static void inline AddNamesToVector(
        std::vector<std::string>& rNamesVector,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue,
        const SizeType StartIndex)
    {
        KRATOS_TRY

        KRATOS_DEBUG_ERROR_IF(StartIndex >= rNamesVector.size() + rValue.size())
            << "rNamesVector size is not enough to allocate values. "
               "rNamesVector vector should have atleast "
            << rValue.size()
            << " elements from(including) StartIndex value. [ StartIndex = " << StartIndex
            << ", rNamesVector.size() = " << rNamesVector.size() << " ].\n";

        for (SizeType i = 0; i < rValue.size(); ++i) {
            rNamesVector[StartIndex + i] =
                rVariable.Name() + "_" + std::to_string(i + 1);
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Get the Variable Data Length
     *
     * @param rValue        Value to determin dynamic variable lengths
     * @return SizeType     Length of the rValue (in flattened number of doubles)
     */
    static SizeType inline GetVariableDataLength(
        const TDataType& rValue)
    {
        return rValue.size();
    }
};

///@}
///@name Operations
///@{

/**
 * @brief Modifies names from list of variables
 *
 * The signature of TValueGetterFunction should accept const ModelPart::NodeType&, and const Variable<TDataType>&.
 *
 * @tparam TDataType                    Data type of the variable
 * @tparam TValueGetterFunction         Getter function prototype
 * @param rNamesList                    Output variable names list
 * @param rNode                         Example node to identify dynamic variable properties
 * @param rVariablesList                List of variables
 * @param rVariableValuesStartIndex     Starting positions of each variable
 * @param pValueGetterFunction          Function pointer to get variable value from nodes
 */
template <class TDataType, class TValueGetterFunction>
void inline AddVariablesListNamesToVector(
    std::vector<std::string>& rNamesList,
    const NodeType& rNode,
    const std::vector<const Variable<TDataType>*>& rVariablesList,
    const IndicesVector& rVariableValuesStartIndex,
    TValueGetterFunction* pValueGetterFunction)
{
    for (std::size_t i = 0; i < rVariablesList.size(); ++i) {
        const auto& r_variable = *(rVariablesList[i]);
        VariableDataCollector<TDataType>::AddNamesToVector(
            rNamesList, r_variable, (*pValueGetterFunction)(rNode, r_variable),
            rVariableValuesStartIndex[i]);
    }
}

/**
 * @brief Get the Historical Value
 *
 * Returns rVariable from historical data value container of nodes.
 *
 * @tparam TDataType    Data type
 * @param rNode         Node
 * @param rVariable     Variable to read
 * @return TDataType    Value
 */
template <class TDataType>
TDataType inline GetHistoricalValue(
    const NodeType& rNode,
    const Variable<TDataType>& rVariable)
{
    return rNode.FastGetSolutionStepValue(rVariable);
}

/**
 * @brief Get the Historical Value
 *
 * Returns rVariable from non historical data value container of nodes.
 *
 * @tparam TDataType    Data type
 * @param rNode         Node
 * @param rVariable     Variable to read
 * @return TDataType    Value
 */
template <class TDataType>
TDataType inline GetNonHistoricalValue(
    const NodeType& rNode,
    const Variable<TDataType>& rVariable)
{
    return rNode.GetValue(rVariable);
}

/**
 * @brief Calculates variable start indices
 *
 * @tparam TDataType                Data type
 * @tparam TValueGetterFunction     Function prototype to retrieve nodal values
 * @param rNode                     Example node to identify dynamic properties of variable values
 * @param rVariablesList            List of variables to calculate start positions
 * @param pValueGetterFunction      Function pointer to read values from nodes
 * @param Offset                    Initial offset for flattned data. (This value will be overwritten by new offset)
 * @return IndicesVector    Vector containing all start positions for each variable
 */
template <class TDataType, class TValueGetterFunction>
IndicesVector GetVariableDataStartIndices(
    const NodeType& rNode,
    const std::vector<const Variable<TDataType>*>& rVariablesList,
    TValueGetterFunction* pValueGetterFunction,
    SizeType& Offset)
{
    const auto number_of_variables = rVariablesList.size();
    IndicesVector start_indices(number_of_variables + 1);

    for (SizeType i = 0; i < number_of_variables; ++i) {
        start_indices[i] = Offset;
        Offset += VariableDataCollector<TDataType>::GetVariableDataLength(
            (*pValueGetterFunction)(rNode, *(rVariablesList[i])));
    }

    start_indices[number_of_variables] = Offset;

    return start_indices;
}

/**
 * @brief Adds interpolation contributions for given variables list
 *
 * This method adds interpolation contributions from the given node for list of variables
 * to a flat double vector which is already sized correctly.
 *
 * The signature of TValueGetterFunction should accept const ModelPart::NodeType&, and const Variable<TDataType>&.
 *
 * The size of rValuesList should be able to carry all the values coming from rVariablesList for all TDataType for all nodes.
 * Example: In the case of 3 node model part with DENSITY and VELOCITY variables, this size should be 12.
 *
 * The size of rVariableValuesStartIndex should be one greater than the size of rVariablesList, last position holding the total
 * size of rSamplePointVariableValuesList.
 *
 * @tparam TDataType                        Data type of variables
 * @tparam TValueGetterFunction             Function prototype to retrive nodal values from.
 * @param rValuesList                       Correctly sized double values vector
 * @param rNode                             Node from where nodal values are retrieved from
 * @param ShapeFunctionValue                Shape function value for sampling point
 * @param pValueGetterFunction              Function pointer to retrieve nodal values from
 * @param rVariableValuesStartIndex         List of indices to indicate where corresponding variable should start to write.
 * @param rVariablesList                    List of variables
 * @param StartIndexOffset                  Start array offset index
 */
template <class TDataType, class TValueGetterFunction>
void inline AddInterpolationContributions(
    std::vector<double>& rValuesList,
    const ModelPart::NodeType& rNode,
    const double ShapeFunctionValue,
    TValueGetterFunction* pValueGetterFunction,
    const IndicesVector& rVariableValuesStartIndex,
    const std::vector<const Variable<TDataType>*>& rVariablesList,
    const SizeType StartIndexOffset)
{
    using variable_data_collector =
        LineOutputProcessUtilities::VariableDataCollector<TDataType>;

    for (SizeType i = 0; i < rVariablesList.size(); ++i) {
        variable_data_collector::AddToValuesVector(
            rValuesList, (*pValueGetterFunction)(rNode, *(rVariablesList[i])) * ShapeFunctionValue,
            StartIndexOffset + rVariableValuesStartIndex[i]);
    }
}

/**
 * @brief Interpolates variables for given variable list tuples
 *
 * This method interpolates different types of variables given in the rVariableInfoTuplesList, using getter function
 * and shape function.
 *
 * TVariableInfoTuplesList is a tuple, which should contain 3 tuple_elements.
 *      1. IndicesVector to hold start indices of each variable.
 *      2. std::vector<const Variable<TDataType>*> to hold list of variable pointers.
 *      3. Value getter function pointer returning TDataType and accepts const ModelPart::NodeType&, and const Variable<TDataType>& input args
 *
 * @tparam TValueGetterFunction             Function prototype to retrive nodal values from.
 * @tparam TVariableInfoTuplesList          List of tuple argument prototype
 * @param rGeometry                         Geometry where sample point is inside of
 * @param rSamplingPointShapeFunctions      Shape function values of the sampling point
 * @param LocalSamplePointValuesOffset      Sample point value offset based on local_index of sample point
 * @param rVariableInfoTuplesList           List of tuple_elements
 */
template <class... TVariableInfoTuplesList>
void InterpolateVariables(
    std::vector<double>& rValuesList,
    const GeometryType& rGeometry,
    const Vector& rSamplingPointShapeFunctions,
    const SizeType LocalSamplePointValuesOffset,
    const TVariableInfoTuplesList&... rVariableInfoTuplesList)
{
    KRATOS_TRY

    const SizeType number_of_nodes = rGeometry.PointsNumber();

    KRATOS_DEBUG_ERROR_IF(number_of_nodes != rSamplingPointShapeFunctions.size())
        << "number_of_nodes != rSamplingPointShapeFunctions.size().";

    for (SizeType i_node = 0; i_node < number_of_nodes; ++i_node) {
        const auto& r_node = rGeometry[i_node];
        const double shape_function_value = rSamplingPointShapeFunctions[i_node];

        int dummy[sizeof...(rVariableInfoTuplesList)] = {(
            AddInterpolationContributions(rValuesList, r_node, shape_function_value,
                                          std::get<2>(rVariableInfoTuplesList),
                                          std::get<0>(rVariableInfoTuplesList),
                                          std::get<1>(rVariableInfoTuplesList), LocalSamplePointValuesOffset),
            0)...};

        *dummy = 0;
    }

    KRATOS_CATCH("");
}

///@}

} // namespace LineOutputProcessUtilities

///@name Kratos Classes
///@{

/**
 * @brief Line output process
 *
 * This process outputs double/array_1d<double, 3> variables in historical/non-historical nodal
 * data value containers to files. Values will be sampled along a lint from starting point to end point
 * with given number of sampling points. All the sampled values will be written into one file per time step.
 *
 * Output frequency can be controlled by output control parameters. But, if a string is used as output control
 * variable, then it is assumed every step is output step.
 *
 * Output is given in Comma Seperated Values (CSV) format.
 *
 */
class KRATOS_API(RANS_APPLICATION) RansLineOutputProcess
: public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansLineOutputProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansLineOutputProcess);

    using NodeType = ModelPart::NodeType;
    using SizeType = std::size_t;
    using IndicesVector = std::vector<SizeType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansLineOutputProcess(Model& rModel, Parameters rParameters);

    /// Destructor.
    ~RansLineOutputProcess() override = default;

    /// Assignment operator.
    RansLineOutputProcess& operator=(RansLineOutputProcess const& rOther) = delete;

    /// Copy constructor.
    RansLineOutputProcess(RansLineOutputProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

    void ExecuteInitialize() override;

    void ExecuteInitializeSolutionStep() override;

    void ExecuteFinalizeSolutionStep() override;

    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}

private:
    ///@name Member Variables
    ///@{

    Model& mrModel;

    std::string mModelPartName;
    std::vector<std::string> mVariableNames;
    bool mWriteHeader;
    array_1d<double, 3> mStartPoint;
    array_1d<double, 3> mEndPoint;

    double mOutputStepInterval;
    double mCurrentStepCount = 0.0;
    double mPreviousStepValue;
    std::string mOutputFileName;
    std::string mOutputStepControlVariableName;

    bool mIsHistoricalValue;
    bool mUpdatePointsEachStep;

    SizeType mNumberOfSamplingPoints;
    std::vector<double> mSamplingPoints;
    std::vector<int> mSamplingPointElementIds;
    std::vector<Vector> mSamplingPointElementShapeFunctions;

    IndicesVector mSamplePointLocalIndexList;
    std::vector<IndicesVector> mSamplePointLocalIndexListMaster;

    template <class TDataType>
    using variables_vector_type = std::vector<const Variable<TDataType>*>;

    variables_vector_type<double> mDoubleVariablesList;
    variables_vector_type<array_1d<double, 3>> mArray3VariablesList;
    variables_vector_type<array_1d<double, 4>> mArray4VariablesList;
    variables_vector_type<array_1d<double, 6>> mArray6VariablesList;
    variables_vector_type<array_1d<double, 9>> mArray9VariablesList;
    variables_vector_type<Vector> mVectorVariablesList;
    variables_vector_type<Matrix> mMatrixVariablesList;

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Updates positions of sample points
     *
     * This method is used to recompute and re-find elements where sample points are located.
     *
     */
    void UpdateSamplePoints();

    /**
     * @brief Checks and adds variables to lists
     *
     * This method is used to check whether the given variables are in historical variables list in the case
     * of mIsHistoricalValue = true. Then it adds the proper variable to proper variable list.
     *
     * @tparam TDataType        Data type of variable
     * @param rVariablesList    List of variables, if the variable found it will be added here
     * @param rModelPart        Model part to check the variable exists in historical variable list
     * @param rVariableName     Name of the variable
     * @return true             rVariableName is found under TDataType variables
     * @return false            rVariableName is not found under TDataType variables
     */
    template <class TDataType>
    bool CheckAndAddVariableToList(
        variables_vector_type<TDataType>& rVariablesList,
        const ModelPart& rModelPart,
        const std::string& rVariableName)
    {
        KRATOS_TRY

        if (KratosComponents<Variable<TDataType>>::Has(rVariableName)) {
            const auto& r_variable =
                KratosComponents<Variable<TDataType>>::Get(rVariableName);
            KRATOS_ERROR_IF(mIsHistoricalValue && !rModelPart.HasNodalSolutionStepVariable(r_variable))
                << rVariableName << " is not found in nodal solution step variables list of "
                << rModelPart.Name() << ".";

            rVariablesList.push_back(&r_variable);

            return true;
        }

        return false;

        KRATOS_CATCH("");
    }

    /**
     * @brief Get Variable value from process info
     *
     * This method checks output control variable is in TDataType variables, and if found it
     * returns the value corresponding from ProcessInfo, if not found it returns false and blank string
     *
     * @tparam TOutputDataType                  Output data type
     * @tparam TDataType                        Data type of the output control variable
     * @param rIsValueObtained                  True if value found, false if not
     * @param rOutputValue                      Holding output variable value
     * @param rVariableName                     Output control variable
     * @return bool                             If variable found true otherwise false
     */
    template <class TOutputDataType, class TDataType>
    void GetVariableValueFromProcessInfo(
        bool& rIsValueObtained,
        TOutputDataType& rOutputValue,
        const std::string& rVariableName) const
    {
        KRATOS_TRY

        if (KratosComponents<Variable<TDataType>>::Has(rVariableName)) {
            const auto& r_process_info =
                mrModel.GetModelPart(mModelPartName).GetProcessInfo();
            const auto& r_variable =
                KratosComponents<Variable<TDataType>>::Get(rVariableName);

            if (r_process_info.Has(r_variable)) {
                rOutputValue = static_cast<TOutputDataType>(r_process_info[r_variable]);
                rIsValueObtained = true;
            }
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Checks and returns rVariableName from process info
     *
     * This method checks whether process info has rVariableName in the list of TDataTypes args list,
     * if available, then value will be casted to double and then returned.
     *
     * @tparam TDataTypes           List of data types to be checked whether variable is available
     * @param rVariableName         Name of the variable
     * @return double               Variable value
     */
    template <class... TDataTypes>
    double CheckAndGetOutputVariableValue(
        const std::string& rVariableName) const
    {
        KRATOS_TRY

        bool is_variable_found = false;
        double value = 0.0;

        int dummy[sizeof...(TDataTypes)] = {(
            GetVariableValueFromProcessInfo<double, TDataTypes>(is_variable_found, value, rVariableName),
            0)...};

        *dummy = 0;

        KRATOS_ERROR_IF(!is_variable_found)
            << "Output step control variable name not found in variables list. "
               "[ "
               "output_step_control_variable_name = "
            << rVariableName << " ].\n";

        return value;

        KRATOS_CATCH("");
    }

    /**
     * @brief Checks whether current step is an output step
     *
     * @return true     If it is an output step
     * @return false    If it is not an output step
     */
    bool IsOutputStep();

    /**
     * @brief Write output file
     *
     */
    void WriteOutputFile() const;

    /**
     * @brief Write output file header
     *
     * @param rOutputFileStream File stream in which the header will be written to
     */
    void WriteOutputFileHeader(
        std::ofstream& rOutputFileStream) const;

    /**
     * @brief Get the Output File Name
     *
     * @return std::string Output filename
     */
    std::string GetOutputFileName() const;

    ///@}

}; // Class RansLineOutputProcess

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansLineOutputProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_LINE_OUTPUT_PROCESS_H_INCLUDED defined
