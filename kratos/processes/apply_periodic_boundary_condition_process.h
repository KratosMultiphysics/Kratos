//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//
//

#ifndef APPLY_PERIODIC_CONDITION_PROCESS_H
#define APPLY_PERIODIC_CONDITION_PROCESS_H

// System includes
// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

namespace Kratos
{

class KRATOS_API(KRATOS_CORE) ApplyPeriodicConditionProcess : public Process
{

  public:
    /// Pointer definition of ApplyPeriodicConditionProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyPeriodicConditionProcess);

    typedef Node<3>                                         NodeType;
    typedef Variable<double>                                VariableType;
    typedef NodeType::IndexType                             IndexType;
    typedef ModelPart::NodeIterator                         NodeIteratorType;
    typedef Matrix                                          MatrixType;
    typedef Vector                                          VectorType;
    typedef Geometry<NodeType>                              GeometryType;

    /**
     * @brief Constructor of the process to apply periodic boundary condition
     * @param rMasterModelPart The master model part for the constraints. Constraints are added on here.
     * @param rSlaveModelPart The slave model part for the constraints.
     * @param Settings parameters for the periodic condition to be applied
     */
    ApplyPeriodicConditionProcess(ModelPart &rMasterModelPart, ModelPart &rSlaveModelPart,
                                  Parameters Settings);

    /**
     * @brief Destructor of the process class
     */
    ~ApplyPeriodicConditionProcess();

    /**
     * @brief Function initializes the process
     */
    void ExecuteInitialize() override;

    /**
     * @brief Function initializes the solution step
     */
    void ExecuteInitializeSolutionStep() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

    /**
     * @brief Function to print the information about this current process
     */
    void PrintInfo(std::ostream& rOStream) const override;

    ///@name Type Definitions
    ///@{
    enum class TransformationType {
        TRANSLATION = 1,
        ROTATION = 2
    };
    ///@}

  private:
    MatrixType mTransformationMatrix; // This can be either for rotating or for translating the slave geometry to Master geometry
    MatrixType mTransformationMatrixVariable; // This can be either for rotating or for translating the master variable to slave geometry
    ModelPart &mrMasterModelPart;       // the master modelpart to which the master-slave constraints are added.
    ModelPart &mrSlaveModelPart;
    Parameters mParameters;          // parameters
    double mAngleOfRotation;
    DenseVector<double> mCenterOfRotation;
    DenseVector<double> mAxisOfRotationVector;
    ApplyPeriodicConditionProcess::TransformationType mTransformationType;
    double mDistance;
    DenseVector<double> mDirOfTranslation;
    double mSearchTolerance;
    IndexType mSearchMaxResults;

    /**
     * @brief  Function to remove the common nodes of slave and master modelparts from the slave modelpart
     */
    void RemoveCommonNodesFromSlaveModelPart();

    /**
     * @brief   The function to figure out how the master and slave model parts relate together and add master-slave constraints
     *          to the rModelPart
     */
    template <int TDim>
    void ApplyConstraintsForPeriodicConditions();

    /**
     * @brief   The function adds the linear master-slave constraint to mrMasterModelPart. This function is specifically for applying
     *          periodic conditions for vector variable. This distinction is because, for vector variables the variable should also be
     *          transformed according to the transfromation specified in the settings.
     * @param rSlaveNode The slave node which is to be connected to the rHostGeometry.
     * @param rHostedGeometry the Host geometry which has the rSlaveNode.
     * @param rWeights The weights with which the rSlaveNode is connected to the rHostedGeometry's nodes.
     * @param rVarName The name of the vector variable on which periodic boundary condition can be applied.
     */
    template <int TDim>
    void ConstraintSlaveNodeWithConditionForVectorVariable(NodeType& rSlaveNode, const GeometryType& rHostedGeometry, const VectorType& rWeights, const std::string& rVarName);

    /**
     * @brief   The function adds the linear master-slave constraint to mrMasterModelPart. This function is specifically for applying
     *          periodic conditions for scalar variable. This distinction is because, for scalar variables the variable need NOT be
     *          transformed according to the transfromation specified in the settings.
     * @param rSlaveNode The slave node which is to be connected to the rHostGeometry.
     * @param rHostedGeometry the Host geometry which has the rSlaveNode.
     * @param rWeights The weights with which the rSlaveNode is connected to the rHostedGeometry's nodes.
     * @param rVarName The name of the scalar variable on which periodic boundary condition can be applied.
     */
    template <int TDim>
    void ConstraintSlaveNodeWithConditionForScalarVariable(NodeType& rSlaveNode, const GeometryType& rHostedGeometry, const VectorType& rWeights, const std::string& rVarName);

    /**
     * @brief Calculate the transformation matrix to account for the moving the two periodic condition modelparts together.
     */
    void CalculateTransformationMatrix();

    /*
     * @brief Transform the point(node_cords) using the mTransformationMatrix calculated in CalculateTransformationMatrix function
     * @param rCoordinates The original coordinates which have to be transformed
     * @param rTransformedCoordinates The new coordinates which are transformed with rTransformationMatrix.
     */
    void TransformNode(const array_1d<double, 3 >& rCoordinates, array_1d<double, 3 >& rTransformedCoordinates) const;



}; // Class


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // APPLY_PERIODIC_CONDITION_PROCESS_H  defined
