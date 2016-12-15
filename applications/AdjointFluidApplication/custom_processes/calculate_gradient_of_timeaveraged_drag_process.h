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


#if !defined(KRATOS_CALCULATE_GRADIENT_OF_TIMEAVERAGED_DRAG_PROCESS_H_INCLUDED )
#define KRATOS_CALCULATE_GRADIENT_OF_TIMEAVERAGED_DRAG_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/process_info.h"
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

/// Process to apply time-averaged drag as an objective function.
/**
 * The objective function process calculates \f$ \partial_{\mathbf{u}} (J^n)^T\f$
 * required for the adjoint solution and \f$ \partial_{\mathbf{s}} (J^n)^T\f$
 * required by the sensitivity calculation process, respectively.
 *
 * The time-averaged drag objective function is defined as:
 *
 * \f[
 *  \bar{J} = \frac{1}{t_{end} - t_{start}}\Sigma_{n=1}^N J^n\Delta t
 * \f]
 *
 * with \f$J^n\f$ the drag at step n. \f$t_{start}\f$ and \f$t_{end}\f$ are
 * defined by process info variables ADJOINT_START_TIME, ADJOINT_END_TIME,
 * respectively.
 *
 * @see CalculateBossakSensitivityProcess
 * @see AdjointBossakScheme
 */
class CalculateGradientOfTimeAveragedDragProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(CalculateGradientOfTimeAveragedDragProcess);

    typedef ModelPart::ElementIterator ElementIterator;

    typedef ModelPart::NodeIterator NodeIterator;

    typedef Element::MatrixType MatrixType;

    typedef Element::VectorType VectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    CalculateGradientOfTimeAveragedDragProcess(
            ModelPart::Pointer pComputeModelPart,
            ModelPart::Pointer pDragModelPart,
            ModelPart::Pointer pShapeModelPart,
            Parameters::Pointer pParameters)
    : Process(),
      mpComputeModelPart(pComputeModelPart),
      mpDragModelPart(pDragModelPart),
      mpShapeModelPart(pShapeModelPart),
      mAdjointStartTime(0),
      mAdjointEndTime(0),
      mAlphaBossak(0)
    {
        KRATOS_TRY

        Parameters DefaultParams(R"(
        {
            "compute_model_part_name": "PLEASE_SPECIFY_MODEL_PART",
            "drag_model_part_name": "PLEASE_SPECIFY_MODEL_PART",
            "shape_model_part_name": "PLEASE_SPECIFY_MODEL_PART",
            "adjoint_start_time": 0.0,
            "adjoint_end_time":1.0,
            "alpha_bossak":-0.3
        })");

        // try accessing parameters without defaults so that an error is thrown
        // if they don't exist
        pParameters->GetValue("compute_model_part_name");
        pParameters->GetValue("drag_model_part_name");
        pParameters->GetValue("shape_model_part_name");

        pParameters->ValidateAndAssignDefaults(DefaultParams);

        mAdjointStartTime = pParameters->GetValue("adjoint_start_time").GetDouble();
        mAdjointEndTime = pParameters->GetValue("adjoint_end_time").GetDouble();
        mAlphaBossak = pParameters->GetValue("alpha_bossak").GetDouble();

        if (mAdjointEndTime <= mAdjointStartTime)
        {
            KRATOS_THROW_ERROR(std::logic_error, "Detected adjoint_end_time <= adjoint_end_time in input parameters", "")
        }


        KRATOS_CATCH("")
    }

    /// Destructor.
    virtual ~CalculateGradientOfTimeAveragedDragProcess() {}

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
    }

    virtual void ExecuteBeforeSolutionLoop()
    {
        KRATOS_TRY

        ProcessInfo& rCurrentProcessInfo = mpComputeModelPart->GetProcessInfo();

        for (auto itNode = std::begin(mpComputeModelPart->Nodes());
                itNode != std::end(mpComputeModelPart->Nodes()); itNode++)
        {
            itNode->Set(STRUCTURE,false);
            itNode->Set(BOUNDARY,false);
        }

        // the drag is computed on the STRUCTURE nodes
        for (auto itNode = std::begin(mpDragModelPart->Nodes());
                itNode != std::end(mpDragModelPart->Nodes()); itNode++)
        {
            itNode->Set(STRUCTURE,true);
        }

        // the sensitivity is computed on the BOUNDARY nodes
        for (auto itNode = std::begin(mpShapeModelPart->Nodes());
                itNode != std::end(mpShapeModelPart->Nodes()); itNode++)
        {
            itNode->Set(BOUNDARY,true);
        }

        KRATOS_CATCH("")
    }

    /**
     * Calculate:
     *
     * \f[
     *  \frac{1}{t_{end} - t_{start}} \partial_{\mathbf{u}} (J^n)^T \Delta t
     * \f]
     */
    virtual void ExecuteInitializeSolutionStep()
    {
        KRATOS_TRY

        ProcessInfo& rCurrentProcessInfo = mpComputeModelPart->GetProcessInfo();
        const double DeltaTime = -rCurrentProcessInfo[DELTA_TIME]; // DELTA_TIME < 0

        if (DeltaTime <= 0)
        {
            KRATOS_THROW_ERROR(std::logic_error, "detected DELTA_TIME >= 0", "")
        }

        const array_1d<double,3> Zero(3,0.0);
        const double Weight = DeltaTime / (mAdjointEndTime - mAdjointStartTime);

        // zero objective gradient
        for (auto itNode = std::begin(mpShapeModelPart->Nodes());
                itNode != std::end(mpShapeModelPart->Nodes()); itNode++)
        {
            itNode->SetValue(OBJECTIVE_FUNCTION_GRADIENT,Zero);
        }

        //  \mathbf{r}^n = \frac{1}{\gamma - 1} \partial_{\mathbf{w}}(\mathbf{M w}^n)^T(\dot{\lambda}^n - \dot{\lambda}^{n+1}) + (\partial_{\mathbf{w}}\mathbf{f}^n)^T\lambda^n + (\partial_{\mathbf{w}}J^n)^T

        MatrixType DampingMatrix;
        MatrixType MassMatrix;
        for (auto itElem = std::begin(mpShapeModelPart->Elements());
                itElem != std::end(mpShapeModelPart->Elements()); itNode++)
        {
            //  (\partial_{\mathbf{w}}\mathbf{f}^n)^T
            itElem->CalculateDampingMatrix(DampingMatrix, rCurrentProcessInfo);


        }



        KRATOS_CATCH("")
    }

    /**
     * Calculate:
     *
     * \f[
     *  \frac{1}{t_{end} - t_{start}} \partial_{\mathbf{s}} (J^n)^T \Delta t
     * \f]
     */
    virtual void ExecuteFinalizeSolutionStep()
    {
        KRATOS_TRY

        ProcessInfo& rCurrentProcessInfo = mpComputeModelPart->GetProcessInfo();
        const double DeltaTime = -rCurrentProcessInfo[DELTA_TIME]; // DELTA_TIME < 0

        const double Weight = DeltaTime / (mAdjointEndTime - mAdjointStartTime);



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
        return "CalculateGradientOfTimeAveragedDragProcess";
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

    ModelPart::Pointer mpComputeModelPart;
    ModelPart::Pointer mpDragModelPart;
    ModelPart::Pointer mpShapeModelPart;
    double mAdjointStartTime;
    double mAdjointEndTime;
    double mAlphaBossak;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
}; /* Class CalculateGradientOfTimeAveragedDragProcess */

///@}

///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  CalculateGradientOfTimeAveragedDragProcess& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const CalculateGradientOfTimeAveragedDragProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Adjoint Fluid Application group

}  /* namespace Kratos */

#endif /* KRATOS_CALCULATE_GRADIENT_OF_TIMEAVERAGED_DRAG_PROCESS_H_INCLUDED defined */
