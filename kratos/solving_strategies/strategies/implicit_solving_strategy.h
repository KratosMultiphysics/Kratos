//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

#pragma once

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategy.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{

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

/**
 * @brief Implicit solving strategy base class
 * This is the base class from which we will derive all the implicit strategies (line-search, NR, etc...)
 * @tparam TSparseSpace Linear algebra sparse space
 * @tparam TDenseSpace Linear algebra dense space
 * @tparam TLinearSolver Linear solver type
 */
template<class TSparseSpace, class TDenseSpace, class TLinearSolver>
class ImplicitSolvingStrategy : public SolvingStrategy<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    typedef SolvingStrategy<TSparseSpace, TDenseSpace>                              BaseType;

    typedef typename BaseType::TDataType                                           TDataType;

    typedef typename BaseType::TSystemMatrixType                           TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType                           TSystemVectorType;

    typedef typename BaseType::TSystemMatrixPointerType             TSystemMatrixPointerType;

    typedef typename BaseType::TSystemVectorPointerType             TSystemVectorPointerType;

    typedef typename BaseType::LocalSystemMatrixType                   LocalSystemMatrixType;

    typedef typename BaseType::LocalSystemVectorType                   LocalSystemVectorType;

    typedef Scheme<TSparseSpace, TDenseSpace>                                    TSchemeType;

    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> TBuilderAndSolverType;

    typedef ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>      ClassType;

    typedef typename BaseType::TDofType                                             TDofType;

    typedef typename BaseType::DofsArrayType                                   DofsArrayType;

    typedef typename BaseType::NodesArrayType                                 NodesArrayType;

    typedef typename BaseType::ElementsArrayType                           ElementsArrayType;

    typedef typename BaseType::ConditionsArrayType                       ConditionsArrayType;

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(ImplicitSolvingStrategy);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    explicit ImplicitSolvingStrategy() { }

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit ImplicitSolvingStrategy(
        ModelPart& rModelPart,
        Parameters ThisParameters)
        : BaseType(rModelPart, ThisParameters)
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    /**
     * @brief Default constructor.
     * @param rModelPart The model part to be computed
     * @param MoveMeshFlag The flag to set if the mesh is moved or not
     */
    explicit ImplicitSolvingStrategy(
        ModelPart& rModelPart,
        bool MoveMeshFlag = false)
        : BaseType(rModelPart, MoveMeshFlag)
    {
    }

    /** Destructor.
     */
    virtual ~ImplicitSolvingStrategy(){}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Create method
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    typename BaseType::Pointer Create(
        ModelPart& rModelPart,
        Parameters ThisParameters) const override
    {
        return Kratos::make_shared<ClassType>(rModelPart, ThisParameters);
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                         : "implicit_solving_strategy",
            "build_level"                  : 2
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);

        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "implicit_solving_strategy";
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief This sets the build level
     * @param Level The build level
     * @details
     * {
     * 0 -> Build StiffnessMatrix just once
     * 1 -> Build StiffnessMatrix at the beginning of each solution step
     * 2 -> build StiffnessMatrix at each iteration
     * }
     */
    void SetRebuildLevel(int Level) override
    {
        mRebuildLevel = Level;
        mStiffnessMatrixIsBuilt = false;
    }

    /**
     * @brief This returns the build level
     * @details
     * {
     * 0 -> Build StiffnessMatrix just once
     * 1 -> Build StiffnessMatrix at the beginning of each solution step
     * 2 -> build StiffnessMatrix at each iteration
     * }
     * @return The build level
     */
    int GetRebuildLevel() const override
    {
        return mRebuildLevel;
    }

    /**
     * @brief This method sets the flag mStiffnessMatrixIsBuilt
     * @param StiffnessMatrixIsBuilt The flag that tells if the stiffness matrix is built
     */
    void SetStiffnessMatrixIsBuilt(const bool StiffnessMatrixIsBuilt)
    {
        mStiffnessMatrixIsBuilt = StiffnessMatrixIsBuilt;
    }

    /**
     * @brief This method gets the flag mStiffnessMatrixIsBuilt
     * @return mStiffnessMatrixIsBuilt: The flag that tells if the stiffness matrix is built
     */
    bool GetStiffnessMatrixIsBuilt() const
    {
        return mStiffnessMatrixIsBuilt;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ImplicitSolvingStrategy";
    }

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    // Settings for the rebuilding of the stiffness matrix
    int mRebuildLevel;            /// The current rebuild level
    bool mStiffnessMatrixIsBuilt; /// A flag indicating if the stiffness matrix has been built

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        // Add base strategy settings
        BaseType::AssignSettings(ThisParameters);

        // By default the matrices are rebuilt at each iteration
        mRebuildLevel = ThisParameters["build_level"].GetInt();
    }

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{

private:

    ///@}
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


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

    /** Copy constructor.
     */
    ImplicitSolvingStrategy(const ImplicitSolvingStrategy& Other);

    ///@}

}; /* Class ImplicitSolvingStrategy */

///@}

///@name Type Definitions
///@{


///@}

} /* namespace Kratos.*/
