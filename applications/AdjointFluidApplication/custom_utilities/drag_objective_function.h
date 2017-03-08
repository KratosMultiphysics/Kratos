//  KratosAdjointFluidApplication
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_DRAG_OBJECTIVE_FUNCTION)
#define KRATOS_DRAG_OBJECTIVE_FUNCTION

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/process_info.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "includes/kratos_parameters.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "custom_utilities/objective_function.h"

namespace Kratos
{
///@addtogroup AdjointFluidApplication
///@{

///@name Kratos Classes
///@{

/// An objective function for drag.
template <unsigned int TDim>
class DragObjectiveFunction : public ObjectiveFunction
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(DragObjectiveFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    DragObjectiveFunction(Parameters& rParameters)
    {
        Parameters DefaultParams(R"(
        {
            "structure_model_part_name": "PLEASE_SPECIFY_MODEL_PART",
            "drag_direction": [1.0, 0.0, 0.0]
        })");

        rParameters.ValidateAndAssignDefaults(DefaultParams);

        mStructureModelPartName = rParameters["structure_model_part_name"].GetString();

        if (rParameters["drag_direction"].IsArray() == false ||
            rParameters["drag_direction"].size() != 3)
        {
            KRATOS_THROW_ERROR(std::runtime_error,
                               "drag_direction vector is not a vector or does "
                               "not have size 3:",
                               rParameters.PrettyPrintJsonString())
        }

        for (unsigned int d = 0; d < TDim; ++d)
            mDragDirection[d] = rParameters["drag_direction"][d].GetDouble();

        if (std::abs(norm_2(mDragDirection) - 1.0) > 1e-3)
        {
            const double magnitude = norm_2(mDragDirection);
            if (magnitude == 0.0)
                KRATOS_THROW_ERROR(std::runtime_error,
                                   "drag_direction is not properly defined.",
                                   "")

            std::cout
                << "WARNING: non unit vector detected in \"drag_direction\": "
                << rParameters.PrettyPrintJsonString() << std::endl;
            std::cout << "normalizing \"drag_direction\"..." << std::endl;

            for (unsigned int d = 0; d < TDim; d++)
                mDragDirection[d] /= magnitude;
        }
    }

    /// Destructor.
    virtual ~DragObjectiveFunction()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    virtual void Initialize(ModelPart& rModelPart)
    {
        if (rModelPart.HasSubModelPart(mStructureModelPartName) == false)
        {
            KRATOS_THROW_ERROR(
                std::runtime_error,
                "invalid parameters \"structure_model_part_name\": ",
                mStructureModelPartName)
        }

// initialize the variables to zero.
#pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodesBegin, NodesEnd);

            for (auto it = NodesBegin; it != NodesEnd; ++it)
                it->Set(STRUCTURE, false);
        }

        ModelPart& rStructureModelPart = rModelPart.GetSubModelPart(mStructureModelPartName);
        for (auto it = rStructureModelPart.NodesBegin();
             it != rStructureModelPart.NodesEnd();
             ++it)
            it->Set(STRUCTURE, true);

        // allocate auxiliary memory
        int NumThreads = OpenMPUtils::GetNumThreads();
        mDragFlagVector.resize(NumThreads);
        mElementIds.resize(NumThreads);

        // initialize drag flag and element id vectors
        Element& rElem = rModelPart.GetElement(0);
        for (unsigned int k = 0; k < mElementIds.size(); ++k)
        {
            mElementIds[k] = 1; // force initialization
            this->GetDragFlagVector(rElem);
        }
    }

    virtual void CalculateAdjointVelocityContribution(const Element& rElem,
                                                      const Matrix& rAdjointMatrix,
                                                      Vector& rRHSContribution,
                                                      ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rRHSContribution.size() != rAdjointMatrix.size1())
            rRHSContribution.resize(rAdjointMatrix.size1(), false);

        Vector& rDragFlagVector = this->GetDragFlagVector(rElem);
        noalias(rRHSContribution) = prod(rAdjointMatrix, rDragFlagVector);

        KRATOS_CATCH("")
    }

    virtual void CalculateAdjointAccelerationContribution(const Element& rElem,
                                                          const Matrix& rAdjointMatrix,
                                                          Vector& rRHSContribution,
                                                          ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rRHSContribution.size() != rAdjointMatrix.size1())
            rRHSContribution.resize(rAdjointMatrix.size1(), false);

        Vector& rDragFlagVector = this->GetDragFlagVector(rElem);
        noalias(rRHSContribution) = prod(rAdjointMatrix, rDragFlagVector);

        KRATOS_CATCH("")
    }

    virtual void CalculateSensitivityContribution(const Element& rElem,
                                                  const Matrix& rDerivativesMatrix,
                                                  Vector& rRHSContribution,
                                                  ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rRHSContribution.size() != rDerivativesMatrix.size1())
            rRHSContribution.resize(rDerivativesMatrix.size1(), false);

        Vector& rDragFlagVector = this->GetDragFlagVector(rElem);
        noalias(rRHSContribution) = prod(rDerivativesMatrix, rDragFlagVector);

        KRATOS_CATCH("")
    }

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}

private:
    ///@name Member Variables
    ///@{

    std::string mStructureModelPartName;
    array_1d<double, TDim> mDragDirection;
    std::vector<Vector> mDragFlagVector;
    std::vector<unsigned int> mElementIds;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    Vector& GetDragFlagVector(const Element& rElement)
    {
        int k = OpenMPUtils::ThisThread();

        // if needed, compute the drag flag vector for this element
        if (rElement.Id() != mElementIds[k])
        {
            const unsigned int NumNodes = rElement.GetGeometry().PointsNumber();
            const unsigned int LocalSize = (TDim + 1) * NumNodes;

            if (mDragFlagVector[k].size() != LocalSize)
                mDragFlagVector[k].resize(LocalSize, false);

            unsigned int LocalIndex = 0;
            for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
            {
                if (rElement.GetGeometry()[iNode].Is(STRUCTURE))
                {
                    for (unsigned int d = 0; d < TDim; d++)
                        mDragFlagVector[k][LocalIndex++] = mDragDirection[d];
                }
                else
                {
                    for (unsigned int d = 0; d < TDim; d++)
                        mDragFlagVector[k][LocalIndex++] = 0.0;
                }

                mDragFlagVector[k][LocalIndex++] = 0.0; // pressure dof
            }

            mElementIds[k] = rElement.Id();
        }

        return mDragFlagVector[k];
    }

    ///@}
};

///@} // Kratos Classes

///@} // Adjoint Fluid Application group

} /* namespace Kratos.*/

#endif /* KRATOS_DRAG_OBJECTIVE_FUNCTION defined */
