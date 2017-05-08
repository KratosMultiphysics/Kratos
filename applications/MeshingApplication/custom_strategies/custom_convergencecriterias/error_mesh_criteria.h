// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____ 
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _ 
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_ERROR_MESH_CRITERIA_H)
#define  KRATOS_ERROR_MESH_CRITERIA_H

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "meshing_application.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

namespace Kratos
{
///@addtogroup MeshingApplication
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

/** @brief Custom convergence criteria for the mortar condition 
 */
template<class TSparseSpace, class TDenseSpace>
class ErrorMeshCriteria : public virtual  ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ErrorMeshCriteria );

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

    typedef TSparseSpace                              SparseSpaceType;

    typedef typename BaseType::TDataType                    TDataType;

    typedef typename BaseType::DofsArrayType            DofsArrayType;

    typedef typename BaseType::TSystemMatrixType    TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType    TSystemVectorType;
    
    typedef ModelPart::ConditionsContainerType    ConditionsArrayType;
    
    typedef ModelPart::NodesContainerType              NodesArrayType;
    
    typedef std::size_t KeyType;

    ///@}
    ///@name Life Cycle
    ///@{
    
    /// Default constructors
    ErrorMeshCriteria(double ErrorTolerance)
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >(),
          mErrorTolerance(ErrorTolerance)
    {
    }

    ///Copy constructor 
    ErrorMeshCriteria( ErrorMeshCriteria const& rOther )
      :BaseType(rOther)
      ,mErrorTolerance(rOther.mErrorTolerance)
    {
    }

    /// Destructor
    virtual ~ErrorMeshCriteria() {}

    ///@}
    ///@name Operators
    ///@{

    /**
     * Compute relative and absolute error.
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     * @return true if convergence is achieved, false otherwise
     */

    bool PostCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    )
    {
        // We initialize the check
        bool ConvergedError = true;
        
        // We initialize the total error
        double TotalErrorPow2 = 0.0;
        double CurrentSolPow2 = 0.0;
        
        // Iterate in the nodes
        NodesArrayType& pNode = rModelPart.Nodes();
        int numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for 
        for(int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            double& NodalError = itNode->GetValue(NODAL_ERROR);
            double MainDoFNodalError = 0.0;           
            double OtherDoFNodalError = 0.0;        
            
            Node<3>::DofsContainerType& NodeDofs = (itNode)->GetDofs();

            std::size_t DoFId;
//             TDataType DofValue;
            
            for (typename Node<3>::DofsContainerType::const_iterator itDof = NodeDofs.begin(); itDof != NodeDofs.end(); itDof++)
            {
                if (itDof->IsFree())
                {
                    DoFId = itDof->EquationId();
//                     DofValue = itDof->GetSolutionStepValue(0);
                    
                    KeyType CurrVar = itDof->GetVariable().Key();
                    if ((CurrVar == DISPLACEMENT_X) || (CurrVar == DISPLACEMENT_Y) || (CurrVar == DISPLACEMENT_Z) || (CurrVar == VELOCITY_X) || (CurrVar == VELOCITY_Y) || (CurrVar == VELOCITY_Z))
                    {
                        MainDoFNodalError += b[DoFId] * b[DoFId];
                    }
                    else
                    {
                        OtherDoFNodalError += b[DoFId] * b[DoFId];
                    }
                }
                else
                {
//                     const double Reaction = itDof->GetSolutionStepReactionValue();
//                     #pragma omp atomic
//                     CurrentSolPow2 += Reaction * Reaction;
                    
                    DoFId = itDof->EquationId();
            
                    #pragma omp atomic
                    CurrentSolPow2 += b[DoFId] * b[DoFId];
                }
            }
        
            const double& NodalH = itNode->FastGetSolutionStepValue(NODAL_H);
            NodalError = NodalH * NodalH * std::sqrt(MainDoFNodalError) + NodalH * std::sqrt(OtherDoFNodalError);
            #pragma omp atomic
            TotalErrorPow2 += NodalError * NodalError;
        }
        
        // Setting the average nodal error
        rModelPart.GetProcessInfo()[AVERAGE_NODAL_ERROR] = mErrorTolerance * std::sqrt((CurrentSolPow2 + TotalErrorPow2)/numNodes);
        
        // Final check
        if (std::sqrt(TotalErrorPow2/(CurrentSolPow2 + TotalErrorPow2)) > mErrorTolerance)
        {
            ConvergedError = false;
        }
    
        if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
        {
            if (ConvergedError == true)
            {
                std::cout << "The error due to the mesh size is under the tolerance prescribed: " << mErrorTolerance << ". No remeshing required" << std::endl;
            }
            else
            {
                std::cout << "The error due to the mesh size is bigger than the tolerance prescribed: " << mErrorTolerance << ". Remeshing required" << std::endl;
            }
        }
        
        return ConvergedError;
    }
    
    /**
     * This function initialize the convergence criteria
     * @param rModelPart: The model part of interest
     */ 
    
    void Initialize(ModelPart& rModelPart)
    {
        BaseType::mConvergenceCriteriaIsInitialized = true;
    }
    
    /**
     * This function initializes the solution step
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */
    
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    )
    {
   
    }
    
    /**
     * This function finalizes the solution step
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */
    
    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    ) 
    {
        
    }

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Acces
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

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
    ///@name Protected LifeCycle
    ///@{
    ///@}

private:
    ///@name Static Member Variables
    ///@{
    
    ///@}
    ///@name Member Variables
    ///@{

    const double mErrorTolerance;
    
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

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}

}; // Class ErrorMeshCriteria 

///@name Explicit Specializations
///@{

}  // namespace Kratos 

#endif /* KRATOS_ERROR_MESH_CRITERIA_H  defined */

