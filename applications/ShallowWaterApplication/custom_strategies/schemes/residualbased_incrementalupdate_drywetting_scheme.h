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

#ifndef KRATOS_RESIDUALBASED_INCREMENTALUPDATE_DRYWETTING_SCHEME_H_INCLUDED
#define KRATOS_RESIDUALBASED_INCREMENTALUPDATE_DRYWETTING_SCHEME_H_INCLUDED


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
class ResidualBasedIncrementalUpdateDryWettingScheme : public ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ResidualBasedIncrementalUpdateDryWettingScheme
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedIncrementalUpdateDryWettingScheme);

    typedef ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace> BaseType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Default constructor.
    ResidualBasedIncrementalUpdateDryWettingScheme()
        : BaseType()
    {}

    /// @brief Constructor with wetting drying model.
    ResidualBasedIncrementalUpdateDryWettingScheme(Process::Pointer pDryWettingModel)
        : BaseType()
        , mpDryWettingModel(pDryWettingModel)
    {}

    /// @brief Copy constructor.
    explicit ResidualBasedIncrementalUpdateDryWettingScheme(ResidualBasedIncrementalUpdateDryWettingScheme& rOther)
        :BaseType(rOther)
    {
    }

    /// Destructor.
    virtual ~ResidualBasedIncrementalUpdateDryWettingScheme(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief It initializes a non-linear iteration
     * @param rModelPart The model of the problem to solve
     * @param rA LHS matrix
     * @param rDx Incremental update of primary variables
     * @param rb RHS Vector
     */
    void InitializeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        if (mpDryWettingModel != 0) {
            mpDryWettingModel->Execute();
        }
        BaseType::InitializeNonLinIteration(rModelPart, rA, rDx, rb);
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
        if (mpDryWettingModel != 0) {
            mpDryWettingModel->ExecuteFinalizeSolutionStep();
        }
        BaseType::FinalizeSolutionStep(rModelPart, rA, rDx, rb);
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
        return "ResidualBasedIncrementalUpdateDryWettingScheme";
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

    Process::Pointer mpDryWettingModel;

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

}; // Class ResidualBasedIncrementalUpdateDryWettingScheme

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_RESIDUALBASED_INCREMENTALUPDATE_DRYWETTING_SCHEME_H_INCLUDED  defined


