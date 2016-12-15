//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//   License:        BSD License
//   Kratos default license: kratos/license.txt
//
//   Project Name:        $AdjointFluidApplication        $
//   Last modified by:    $Author: michael.andre@tum.de   $
//   Date:                $Date:         November  2016   $
//   Revision:            $Revision:                0.0   $


#if !defined(KRATOS_CALCULATE_BOSSAK_DRAG_SENSITIVITY_PROCESS_H_INCLUDED )
#define KRATOS_CALCULATE_BOSSAK_DRAG_SENSITIVITY_PROCESS_H_INCLUDED

// System includes
#include <cmath>
#include <iostream>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{

///@addtogroup AdjointFluidApplication
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

/// Process to calculate drag sensitivity from Bossak adjoint solution.
class CalculateBossakDragSensitivityProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(CalculateBossakDragSensitivityProcess);

    typedef ModelPart::ElementIterator ElementIterator;

    typedef ModelPart::NodeIterator NodeIterator;

    typedef Element::MatrixType MatrixType;

    typedef Element::VectorType VectorType;

    typedef unsigned int IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    CalculateBossakDragSensitivityProcess(ModelPart& rModelPart, Parameters& rParameters)
    : Process(), mrModelPart(rModelPart)
    {
        KRATOS_TRY

        Parameters DefaultParams(R"(
        {
            "model_part_name": "PLEASE_SPECIFY_MODEL_PART",
            "structure_model_part_name": "PLEASE_SPECIFY_MODEL_PART",
            "boundary_model_part_name": "PLEASE_SPECIFY_MODEL_PART",
            "alpha_bossak": -0.3,
            "adjoint_start_time": 0.0,
            "adjoint_end_time": 1.0,
            "drag_direction": [1.0, 0.0, 0.0]
        })");

        rParameters.ValidateAndAssignDefaults(DefaultParams);

        mStructureModelPartName = rParameters["structure_model_part_name"].GetString();
        mBoundaryModelPartName = rParameters["boundary_model_part_name"].GetString();

        mAlphaBossak = rParameters["alpha_bossak"].GetDouble();

        ProcessInfo& rProcessInfo = mrModelPart.GetProcessInfo();
        rProcessInfo[START_TIME] = rParameters["adjoint_start_time"].GetDouble();
        rProcessInfo[END_TIME] = rParameters["adjoint_end_time"].GetDouble();

        if (rProcessInfo[START_TIME] >= rProcessInfo[END_TIME])
        {
            KRATOS_THROW_ERROR(std::runtime_error, "invalid parameters: adjoint_start_time >= adjoint_end_time", rParameters.PrettyPrintJsonString())
        }

        if(rParameters["drag_direction"].IsArray() == true && rParameters["drag_direction"].size() != 3)
        {
            KRATOS_THROW_ERROR(std::runtime_error, "drag_direction vector is not a vector or does not have size 3:",rParameters.PrettyPrintJsonString())
        }

        array_1d<double, 3>& rDragDirection = rProcessInfo[DRAG_DIRECTION];
        rDragDirection[0] = rParameters["drag_direction"][0].GetDouble();
        rDragDirection[1] = rParameters["drag_direction"][1].GetDouble();
        rDragDirection[2] = rParameters["drag_direction"][2].GetDouble();

        double magnitude = 0.0;
        for (IndexType d = 0; d < 3; d++)
        {
            magnitude += rDragDirection[d] * rDragDirection[d];
        }
        magnitude = std::sqrt(magnitude);

        if (std::abs(magnitude - 1.0) > 1e-3)
        {
            std::cout << "WARNING: non unit vector detected in \"drag_direction\": " << rParameters.PrettyPrintJsonString() << std::endl;
            std::cout << "normalizing \"drag_direction\"..." << std::endl;

            for (IndexType d = 0; d < 3; d++)
            {
                rDragDirection[d] /= magnitude;
            }
        }

        rProcessInfo[WINDOW_FUNCTION] = 1.0; // constant window function by default

        KRATOS_CATCH("")
    }

    /// Destructor.
    virtual ~CalculateBossakDragSensitivityProcess() {}

    ///@}
    ///@name Operators
    ///@{

    /// This operator simply calls the Execute method.
    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    virtual void Execute() {}

    virtual void ExecuteInitialize()
    {
        KRATOS_TRY

        if (mrModelPart.HasSubModelPart(mStructureModelPartName) == false)
        {
            KRATOS_THROW_ERROR(std::runtime_error, "invalid parameters \"structure_model_part_name\": ", mStructureModelPartName)
        }

        if (mrModelPart.HasSubModelPart(mBoundaryModelPartName) == false)
        {
            KRATOS_THROW_ERROR(std::runtime_error, "invalid parameters \"boundary_model_part_name\": ", mBoundaryModelPartName)
        }

        // Initialize the variables to zero.
        #pragma omp parallel
        {
          NodeIterator NodesBegin;
          NodeIterator NodesEnd;
          OpenMPUtils::PartitionedIterators(mrModelPart.Nodes(),NodesBegin,
                                            NodesEnd);

          const array_1d<double,3> Zero(3,0.0);

          for (NodeIterator it = NodesBegin; it != NodesEnd; ++it)
          {
              it->FastGetSolutionStepValue(SHAPE_SENSITIVITY) = Zero;
              it->FastGetSolutionStepValue(ADJOINT_VELOCITY) = Zero;
              it->FastGetSolutionStepValue(ADJOINT_PRESSURE) = 0.0;
              it->FastGetSolutionStepValue(ADJOINT_ACCELERATION) = Zero;
              it->Set(STRUCTURE,false);
              it->Set(BOUNDARY,false);
          }
        }

        ModelPart& rStructureModelPart = mrModelPart.GetSubModelPart(mStructureModelPartName);
        for (NodeIterator it = rStructureModelPart.NodesBegin(); it != rStructureModelPart.NodesEnd(); it++)
        {
            it->Set(STRUCTURE,true);
        }

        ModelPart& rBoundaryModelPart = mrModelPart.GetSubModelPart(mBoundaryModelPartName);
        for (NodeIterator it = rBoundaryModelPart.NodesBegin(); it != rBoundaryModelPart.NodesEnd(); it++)
        {
            it->Set(BOUNDARY,true);
        }


        KRATOS_CATCH("")
    }

    virtual void ExecuteBeforeSolutionLoop()
    {
    }

    virtual void ExecuteInitializeSolutionStep()
    {
    }

    virtual void ExecuteFinalizeSolutionStep()
    {
        KRATOS_TRY

        ProcessInfo& rProcessInfo = mrModelPart.GetProcessInfo();
        const int NumThreads = OpenMPUtils::GetNumThreads();
        std::vector< Vector > DragFlagVector(NumThreads);
        std::vector< Vector > FluidAuxVector(NumThreads);
        std::vector< Vector > CoordAuxVector(NumThreads);
        std::vector< Matrix > ShapeDerivativesMatrix(NumThreads);

        double DeltaTime = -rProcessInfo[DELTA_TIME]; // DELTA_TIME < 0
        const IndexType DomainSize = static_cast<IndexType>(rProcessInfo[DOMAIN_SIZE]);

        if (DeltaTime <= 0)
        {
            KRATOS_THROW_ERROR(std::runtime_error, "detected for adjoint solution DELTA_TIME >= 0", "")
        }

        double weight = rProcessInfo[WINDOW_FUNCTION] * DeltaTime / (rProcessInfo[END_TIME] - rProcessInfo[START_TIME]);

        OpenMPUtils::PartitionVector Partition;
        OpenMPUtils::DivideInPartitions(mrModelPart.NumberOfElements(), NumThreads, Partition);
        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            ModelPart::ElementIterator ElementsBegin = mrModelPart.ElementsBegin() + Partition[k];
            ModelPart::ElementIterator ElementsEnd = mrModelPart.ElementsBegin() + Partition[k + 1];

            for (auto itElem = ElementsBegin; itElem != ElementsEnd; itElem++)
            {
                this->GetDragFlagVector(DragFlagVector[k], *(itElem.base()), rProcessInfo);

                if (norm_1(DragFlagVector[k]) == 0.0) // true for most elements
                {
                    continue;
                }

                itElem->Calculate(SHAPE_DERIVATIVE_MATRIX_2,ShapeDerivativesMatrix[k],rProcessInfo);

                itElem->GetFirstDerivativesVector(FluidAuxVector[k],0);

                noalias(FluidAuxVector[k]) += DragFlagVector[k];

                if (CoordAuxVector[k].size() != ShapeDerivativesMatrix[k].size1())
                    CoordAuxVector[k].resize(ShapeDerivativesMatrix[k].size1());

                noalias(CoordAuxVector[k]) = prod(ShapeDerivativesMatrix[k], FluidAuxVector[k]);

                // Carefully write results to nodal variables
                IndexType CoordIndex = 0;
                for (IndexType iNode = 0; iNode < itElem->GetGeometry().PointsNumber(); ++iNode)
                {
                    if (itElem->GetGeometry()[iNode].Is(BOUNDARY))
                    {
                        itElem->GetGeometry()[iNode].SetLock();
                        array_1d<double,3>& rSensitivity =
                                itElem->GetGeometry()[iNode].FastGetSolutionStepValue(SHAPE_SENSITIVITY);
                         for (IndexType d = 0; d < DomainSize; ++d)
                         {
                             rSensitivity[d] += weight * CoordAuxVector[k][CoordIndex++];
                         }
                         itElem->GetGeometry()[iNode].UnSetLock();
                    }
                    else
                    {
                        // Skip this node block.
                        CoordIndex += DomainSize;
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }

    virtual void ExecuteBeforeOutputStep()
    {
    }

    virtual void ExecuteAfterOutputStep()
    {
    }

    virtual void ExecuteFinalize()
    {
    }

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
        return "CalculateBossakDragSensitivityProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << this->Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    double mAlphaBossak;
    std::string mStructureModelPartName;
    std::string mBoundaryModelPartName;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void GetDragFlagVector(Vector& rOutput, Element::Pointer pElement,
            ProcessInfo& rProcessInfo)
        {
            const IndexType DomainSize = static_cast<IndexType>(rProcessInfo[DOMAIN_SIZE]);
            const IndexType NumNodes = pElement->GetGeometry().PointsNumber();
            const IndexType LocalSize = (DomainSize + 1) * NumNodes;

            if (rOutput.size() != LocalSize)
            {
                rOutput.resize(LocalSize, false);
            }

            array_1d<double, 3>& rDragDirection = rProcessInfo[DRAG_DIRECTION];
            IndexType LocalIndex = 0;
            for (IndexType iNode = 0; iNode < NumNodes; ++iNode)
            {
                if (pElement->GetGeometry()[iNode].Is(STRUCTURE))
                {
                    for (IndexType d = 0; d < DomainSize; d++)
                    {
                        rOutput[LocalIndex++] = rDragDirection[d];
                    }
                }
                else
                {
                    for (IndexType d = 0; d < DomainSize; d++)
                    {
                        rOutput[LocalIndex++] = 0.0;
                    }
                }

                rOutput[LocalIndex++] = 0.0; // pressure dof
            }
        }

    ///@}
    ///@name Private Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
}; /* Class CalculateBossakDragSensitivityProcess */

///@}

///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  CalculateBossakDragSensitivityProcess& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const CalculateBossakDragSensitivityProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Adjoint Fluid Application group

}  /* namespace Kratos */

#endif /* KRATOS_CALCULATE_BOSSAK_DRAG_SENSITIVITY_PROCESS_H_INCLUDED defined */
