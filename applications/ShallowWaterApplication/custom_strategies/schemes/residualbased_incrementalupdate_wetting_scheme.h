//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_RESIDUALBASED_INCREMENTALUPDATE_WETTING_SCHEME_H_INCLUDED
#define KRATOS_RESIDUALBASED_INCREMENTALUPDATE_WETTING_SCHEME_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"


namespace Kratos
{
///@addtogroup ShallowWaterApplication
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

/// Short class definition.
/** Detail class definition.
*/
template<class TSparseSpace, class TDenseSpace >
class ResidualBasedIncrementalUpdateWettingScheme : public ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ResidualBasedIncrementalUpdateWettingScheme
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedIncrementalUpdateWettingScheme);

    typedef ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace> BaseType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Default constructor.
    ResidualBasedIncrementalUpdateWettingScheme()
        : BaseType()
    {}

    /// @brief Constructor with wetting drying model.
    ResidualBasedIncrementalUpdateWettingScheme(Process::Pointer pWettingModel)
        : BaseType()
        , mpWettingModel(pWettingModel)
    {}

    /// @brief Copy constructor.
    explicit ResidualBasedIncrementalUpdateWettingScheme(ResidualBasedIncrementalUpdateWettingScheme& rOther)
        :BaseType(rOther)
    {
    }

    /// Destructor.
    virtual ~ResidualBasedIncrementalUpdateWettingScheme(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Function called once at the beginning of each solution step.
     * @details The basic operations to be carried in there are the following:
     * - managing variables to be kept constant over the time step (for example time-Scheme constants depending on the actual time step)
     * @param rModelPart The model part of the problem to solve
     * @param rA LHS matrix
     * @param rDx Incremental update of primary variables
     * @param rb RHS Vector
     */
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        if (mpWettingModel != nullptr) {
            mpWettingModel->ExecuteInitializeSolutionStep();
        }
        BaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);
    }

    /**
     * @brief Function called once at the end of a solution step, after convergence is reached if an iterative process is needed
     * @param rModelPart The model part of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        BaseType::FinalizeSolutionStep(rModelPart, rA, rDx, rb);
        if (mpWettingModel != nullptr) {
            mpWettingModel->ExecuteFinalizeSolutionStep();
        }
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
    std::string Info() const override
    {
        return "ResidualBasedIncrementalUpdateWettingScheme";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    Process::Pointer mpWettingModel;

    ///@}
    ///@name Private Operators
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

    ///@}

}; // Class ResidualBasedIncrementalUpdateWettingScheme

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_RESIDUALBASED_INCREMENTALUPDATE_WETTING_SCHEME_H_INCLUDED  defined


