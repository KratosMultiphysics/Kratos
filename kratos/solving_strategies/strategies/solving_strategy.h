//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                    
//

#if !defined(KRATOS_SOLVING_STRATEGY )
#define  KRATOS_SOLVING_STRATEGY

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"

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

/** Solving strategy base class
This is the base class from which we will derive all the strategies (line-search, NR, etc...)

\URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

\URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

\URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

\URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}

\URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

\URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

\URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

\URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}
 */

template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class SolvingStrategy
{
public:
    ///@name Type Definitions
    ///@{
    
//     typedef std::set<Dof::Pointer,ComparePDof>                                    DofSetType;

    typedef typename TSparseSpace::DataType                                        TDataType;
    
    typedef typename TSparseSpace::MatrixType                              TSystemMatrixType;
    
    typedef typename TSparseSpace::VectorType                              TSystemVectorType;

    typedef typename TSparseSpace::MatrixPointerType                TSystemMatrixPointerType;
    
    typedef typename TSparseSpace::VectorPointerType                TSystemVectorPointerType;

    typedef typename TDenseSpace::MatrixType                           LocalSystemMatrixType;
    
    typedef typename TDenseSpace::VectorType                           LocalSystemVectorType;

    typedef Scheme<TSparseSpace, TDenseSpace>                                    TSchemeType;
    
    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> TBuilderAndSolverType;

    typedef typename ModelPart::DofType                                             TDofType;
    
    typedef typename ModelPart::DofsArrayType                                  DofsArrayType;
    
//     typedef Dof<TDataType>                                                          TDofType;
    
//     typedef PointerVectorSet<TDofType, IdentityFunction<TDofType> >            DofsArrayType;
    
//     typedef PointerVectorSet<TDofType, IndexedObject>                          DofsArrayType;
    
    typedef typename DofsArrayType::iterator                                 DofIteratorType;
    
    typedef typename DofsArrayType::const_iterator                   DofConstantIteratorType;
    
    typedef ModelPart::NodesContainerType                                     NodesArrayType;
    
    typedef ModelPart::ElementsContainerType                               ElementsArrayType;
    
    typedef ModelPart::ConditionsContainerType                           ConditionsArrayType;
    
    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(SolvingStrategy);
    
    ///@}
    ///@name Life Cycle
    ///@{

    /** Constructor.
     * @param rModelPart: The model part to be computed
     * @param MoveMeshFlag: The flag to set if the mesh is moved or not
     */

    SolvingStrategy(
        ModelPart& rModelPart, 
        bool MoveMeshFlag = false
    )
        : mrModelPart(rModelPart)
    {
        SetMoveMeshFlag(MoveMeshFlag);
    }
    
    ///@}

    /** Destructor.
     */

    ///@{
    virtual ~SolvingStrategy(){}
    
    ///@}

    /**OPERATIONS ACCESSIBLE FROM THE INPUT:*/
    
    ///@{

    /**
     * Operation to predict the solution ... if it is not called a trivial predictor is used in which the
     * values of the solution step of interest are assumed equal to the old values
     */
    virtual void Predict()
    {
    }

    /**
     * Initialization of member variables and prior operations
     */
    virtual void Initialize()
    {
    }

    /**
     * The problem of interest is solved.
     * This function calls sequentially: Initialize(), InitializeSolutionStep(), Predict(), SolveSolutionStep() and FinalizeSolutionStep().
     * All those functions can otherwise be called separately.
     */
    virtual double Solve()
    {
        Initialize();
        InitializeSolutionStep();
        Predict();
        SolveSolutionStep();
        FinalizeSolutionStep();
        
        return 0.0; 
        
    }

    /**
     * Clears the internal storage
     */
    virtual void Clear()
    {
    }

    /**
     * This should be considered as a "post solution" convergence check which is useful for coupled analysis
     * - the convergence criteria used is the one used inside the "solve" step
     */
    virtual bool IsConverged()
    {
        return true;
    }

    /**
     * This operations should be called before printing the results when non trivial results (e.g. stresses)
     * need to be calculated given the solution of the step
     * This operations should be called only when needed, before printing as it can involve a non negligible cost
     */
    virtual void CalculateOutputData()
    {
    }

    /**
     * Performs all the required operations that should be done (for each step) before solving the solution step.
     * A member variable should be used as a flag to make sure this function is called only once per step.
     */
    virtual void InitializeSolutionStep()
    {
    }

    /**
     * Performs all the required operations that should be done (for each step) after solving the solution step.
     * A member variable should be used as a flag to make sure this function is called only once per step.
     */
    virtual void FinalizeSolutionStep()
    {
    }

    /**
     * Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    virtual bool SolveSolutionStep()
    {
        return true;
    }

    /**
     * This sets the level of echo for the solving strategy
     * @param Level of echo for the solving strategy:
     * 0 -> Mute... no echo at all
     * 1 -> Printing time and basic informations
     * 2 -> Printing linear solver data
     * 3 -> Print of debug informations: Echo of stiffness matrix, Dx, b...
     */
    virtual void SetEchoLevel(int Level)
    {
        mEchoLevel = Level;
    }

    /**
     * This returns the level of echo for the solving strategy
     * @return Level of echo for the solving strategy:
     * 0 -> Mute... no echo at all
     * 1 -> Printing time and basic informations
     * 2 -> Printing linear solver data
     * 3 -> Print of debug informations: Echo of stiffness matrix, Dx, b...
     */
    virtual int GetEchoLevel()
    {
        return mEchoLevel;
    }

    /**
     * This sets the build level
     * @param Level: The build level:
     * 0 -> Build StiffnessMatrix just once
     * 1 -> Build StiffnessMatrix at the beginning of each solution step
     * 2 -> build StiffnessMatrix at each iteration
     */
    virtual void SetRebuildLevel(int Level)
    {
        mRebuildLevel = Level;
        mStiffnessMatrixIsBuilt = false;
    }

    /**
     * This returns the build level
     * @return The build level:
     * 0 -> Build StiffnessMatrix just once
     * 1 -> Build StiffnessMatrix at the beginning of each solution step
     * 2 -> build StiffnessMatrix at each iteration
     */
    virtual int GetRebuildLevel()
    {
        return mRebuildLevel;
    }

    /**
     * This function sets the flag that says if the mesh is moved
     * @param Flag: True if the mesh is moved, false otherwise
     */
    void SetMoveMeshFlag(bool Flag)
    {
        mMoveMeshFlag = Flag;
    }

    /**
     * This function returns the flag that says if the mesh is moved
     * @return True if the mesh is moved, false otherwise
     */
    bool MoveMeshFlag()
    {
        return mMoveMeshFlag;
    }

    /**
     * This function is designed to move the mesh
     * NOTE: Be careful it just consider displacements, derive this method to adapt to your own strategies (ALE, FSI, etc...)
     */
    virtual void MoveMesh()
    {
        KRATOS_TRY

        if (GetModelPart().NodesBegin()->SolutionStepsDataHas(DISPLACEMENT_X) == false)
        {
            KRATOS_ERROR << "It is impossible to move the mesh since the DISPLACEMENT var is not in the Model Part. Either use SetMoveMeshFlag(False) or add DISPLACEMENT to the list of variables" << std::endl;
        }   

        NodesArrayType& NodesArray = GetModelPart().Nodes();
        const int numNodes = static_cast<int>(NodesArray.size());

        #pragma omp parallel for
        for(int i = 0; i < numNodes; i++)  
        {
            auto itNode = NodesArray.begin() + i;

            noalias(itNode->Coordinates()) = itNode->GetInitialPosition().Coordinates();
            noalias(itNode->Coordinates()) += itNode->FastGetSolutionStepValue(DISPLACEMENT);
        }
        
        if (this->GetEchoLevel() != 0 && GetModelPart().GetCommunicator().MyPID() == 0 )
        {
            std::cout<<" MESH MOVED "<<std::endl;
        }

        KRATOS_CATCH("")
    }

    /**
     * Operations to get the pointer to the model
     * @return mrModelPart: The model part member variable
     */
    inline ModelPart& GetModelPart()
    {
        return mrModelPart;
    };

    /**
     * Operations to get the residual norm
     * @return The residual norm
     */
    virtual double GetResidualNorm()
    {
        return 0.0;
    }

    /**
     * Function to perform expensive checks.
     * It is designed to be called ONCE to verify that the input is correct.
     */
    virtual int Check()
    {
        KRATOS_TRY

        // Check if displacement var is needed
        if (mMoveMeshFlag == true)
        {
            for (ModelPart::NodesContainerType::iterator itNode = GetModelPart().NodesBegin();
                 itNode != GetModelPart().NodesEnd(); itNode++)
            {
                if (itNode->SolutionStepsDataHas(DISPLACEMENT) == false)
                {
                    std::cout << "problem on node with Id " << itNode->Id() << std::endl;
                    KRATOS_ERROR << "It is impossible to move the mesh since the DISPLACEMENT var is not in the rModelPart. Either use SetMoveMeshFlag(False) or add DISPLACEMENT to the list of variables" << std::endl;
                }
            }
        }

        for (ModelPart::ElementsContainerType::iterator itElem = GetModelPart().ElementsBegin();
             itElem != GetModelPart().ElementsEnd(); itElem++)
        {
            itElem->Check(GetModelPart().GetProcessInfo());
        }

        for (ModelPart::ConditionsContainerType::iterator itCond = GetModelPart().ConditionsBegin();
             itCond != GetModelPart().ConditionsEnd(); itCond++)
        {
            itCond->Check(GetModelPart().GetProcessInfo());
        }
        
        return 0;
        
        KRATOS_CATCH("")
    }

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    // Level of echo for the solving strategy
    int mEchoLevel;

    // Settings for the rebuilding of the stiffness matrix
    int mRebuildLevel;
    bool mStiffnessMatrixIsBuilt;

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
    ///@name Protected LifeCycle
    ///@{

private:

    ///@}
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;

    bool mMoveMeshFlag;

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
    SolvingStrategy(const SolvingStrategy& Other);


    ///@}

}; /* Class NewSolvingStrategy */

///@}

///@name Type Definitions
///@{


///@}

} /* namespace Kratos.*/

#endif /* KRATOS_SOLVING_STRATEGY  defined */

