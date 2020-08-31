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

#ifndef APPLY_SLIDING_INTERFACE_PROCESS_H
#define APPLY_SLIDING_INTERFACE_PROCESS_H

// System includes
// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "utilities/binbased_fast_point_locator_conditions.h"

namespace Kratos
{

template <int TDim>
class KRATOS_API(CHIMERA_APPLICATION) SlidingInterfaceProcess : public Process
{

  public:
    /// Pointer definition of SlidingInterfaceProcess
    KRATOS_CLASS_POINTER_DEFINITION(SlidingInterfaceProcess);

    typedef Node<3>                                         NodeType;
    typedef Variable<double>                                VariableType;
    typedef NodeType::IndexType                             IndexType;
    typedef ModelPart::NodeIterator                         NodeIteratorType;
    typedef Matrix                                          MatrixType;
    typedef Vector                                          VectorType;
    typedef Geometry<NodeType>                              GeometryType;

    /**
     * @brief Constructor of the process to apply periodic boundary condition.
     * @param rMasterModelPart The master model part for the constraints. Constraints are added on here.
     *                          This is assumed to be the stationary modelpart on the interface.
     * @param rSlaveModelPart The slave model part for the constraints. This is the moving modelpart of
     *                          the interface.
     * @param Settings parameters for the periodic condition to be applied
     */
    SlidingInterfaceProcess(ModelPart &rMasterModelPart, ModelPart &rSlaveModelPart,
                                  Parameters Settings);

    /**
     * @brief Destructor of the process class
     */
    ~SlidingInterfaceProcess();

    /**
     * @brief Function initializes the process
     */
    void ExecuteInitialize() override;

    /**
     * @brief Function finalizes the process
     */
    void ExecuteFinalize() override;

    /**
     * @brief Function initializes the solution step
     */
    void ExecuteInitializeSolutionStep() override;

    /**
     * @brief Function finalizes the solution step
     */
    void ExecuteFinalizeSolutionStep() override;

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

    ///@}

  private:
    ModelPart &mrMasterModelPart;    // the master modelpart
    ModelPart &mrSlaveModelPart;     // the slave modelpart
    Parameters mParameters;          // parameters
    double mSearchTolerance;
    IndexType mSearchMaxResults;
    std::string mSearchModelPartName;
    typename BinBasedFastPointLocatorConditions<TDim>::Pointer mpPointLocator;

    /**
     * @brief   The function to figure out how the master and slave model parts relate together and add master-slave constraints
     *          to the root modelpart of the mrMasterModelPart.
     */
    void ApplyConstraintsForSlidingInterface();

    /**
     * @brief   The function adds the linear master-slave constraint to mrMasterModelPart. This function is specifically for applying periodic conditions for vector variable. This distinction is because, for vector variables the variable should also be transformed according to the transfromation specified in the settings.
     * @param rSlaveNode The slave node which is to be connected to the rHostGeometry.
     * @param rHostedGeometry the Host geometry which has the rSlaveNode.
     * @param rWeights The weights with which the rSlaveNode is connected to the rHostedGeometry's nodes.
     * @param rVarName The name of the vector variable on which periodic boundary condition can be applied.
     */
    void ConstraintSlaveNodeWithConditionForVariable(NodeType& rSlaveNode, const GeometryType& rHostedGeometry, const VectorType& rWeights, const std::string& rVarName);

    /**
     * @brief   The function checks if this is a MPI run, in case it is, it returns a reference to gathered_modelpart  else
     *          a reference to the mrMasterModelPart
     */
    void MakeSearchModelpart();

}; // Class


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // APPLY_SLIDING_INTERFACE_PROCESS_H  defined
