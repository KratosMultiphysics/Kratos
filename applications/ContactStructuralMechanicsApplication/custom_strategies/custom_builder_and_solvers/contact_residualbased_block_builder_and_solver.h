//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
#if !defined(KRATOS_CONTACT_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER )
#define  KRATOS_CONTACT_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER

/* System includes */

/* External includes */

/* Project includes */
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

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

/** Short class definition.

Detail class definition.

Current class provides an implementation for contact builder and solving operations.

the RHS is constituted by the unbalanced loads (residual)

Degrees of freedom are reordered putting the restrained degrees of freedom at
the end of the system ordered in reverse order with respect to the DofSet.

Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
this information.

Calculation of the reactions involves a cost very similiar to the calculation of the total residual

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
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ContactResidualBasedBlockBuilderAndSolver
    : public ResidualBasedBlockBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    ///@name Type Definitions
    ///@{
    
    KRATOS_CLASS_POINTER_DEFINITION(ContactResidualBasedBlockBuilderAndSolver);

    typedef ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /** Constructor.
     */
    ContactResidualBasedBlockBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : ResidualBasedBlockBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >(pNewLinearSystemSolver)
    {
    }

    /** Destructor.
     */
    ~ContactResidualBasedBlockBuilderAndSolver() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    void ApplyDirichletConditions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        ) override
    {
        FixIsolatedNodes(rModelPart);
        
        BaseType::ApplyDirichletConditions(pScheme, rModelPart, A, Dx, b);
        
        FreeIsolatedNodes(rModelPart);
    }

    //**************************************************************************
    //**************************************************************************

    void BuildRHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& b
        ) override
    {
        FixIsolatedNodes(rModelPart);
        
        BaseType::BuildRHS(pScheme, rModelPart, b);
        
        FreeIsolatedNodes(rModelPart);
    }
    
    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

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
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables 
    ///@{

    ///@}
    ///@name Member Variables 
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    /**
     * This method check the ISOLATED nodes and it fixes
     * @param rModelPart The model part to compute
     */
    void FixIsolatedNodes(ModelPart& rModelPart)
    {
        KRATOS_ERROR_IF(!(rModelPart.HasSubModelPart("Contact"))) << "ERROR:: CONTACT MODEL PART NOT CREATED" << std::endl;
        KRATOS_ERROR_IF(!(rModelPart.HasSubModelPart("ComputingContact"))) << "ERROR:: CONTACT COMPUTING MODEL PART NOT CREATED" << std::endl;
        ModelPart& contact_model_part = rModelPart.GetSubModelPart("Contact"); 
        ModelPart& computing_contact_model_part = rModelPart.GetSubModelPart("ComputingContact"); 
        
        // We reset the flag
        auto& nodes_array = contact_model_part.Nodes();
        #pragma omp parallel for 
        for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
            (nodes_array.begin() + i)->Set(VISITED, false);
            (nodes_array.begin() + i)->Set(ISOLATED, false);
        }
        
        // Now we set the flag in the nodes
        auto& conditions_array = computing_contact_model_part.Conditions();
        
        #pragma omp parallel for 
        for(int i = 0; i < static_cast<int>(conditions_array.size()); ++i) {
            auto it_cond = conditions_array.begin() + i;
            auto& geom = it_cond->GetGeometry();
            for (std::size_t i_node = 0; i_node < geom.size(); ++i_node) {
                geom[i_node].SetLock();
                if (geom[i_node].Is(VISITED) == false) {
                    geom[i_node].Set(ISOLATED, it_cond->Is(ISOLATED));
                    geom[i_node].Set(VISITED, true);
                } else {
                    geom[i_node].Set(ISOLATED, geom[i_node].Is(ISOLATED) && it_cond->Is(ISOLATED));
                }
                geom[i_node].UnSetLock();
            }
        }
        
        // We fix the LM
        #pragma omp parallel for 
        for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
            auto it_node = nodes_array.begin() + i;
            if (it_node->Is(ISOLATED) == true) {
                if (it_node->SolutionStepsDataHas(NORMAL_CONTACT_STRESS) == true)
                    it_node->Fix(NORMAL_CONTACT_STRESS);
                else { // We will assume that VECTOR_LAGRANGE_MULTIPLIER
                    it_node->Fix(VECTOR_LAGRANGE_MULTIPLIER_X);
                    it_node->Fix(VECTOR_LAGRANGE_MULTIPLIER_Y);
                    it_node->Fix(VECTOR_LAGRANGE_MULTIPLIER_Z);
                }
            }
        }
    }
    
    /**
     * This method releases the ISOLATED nodes 
     * @param rModelPart The model part to compute
     */
    void FreeIsolatedNodes(ModelPart& rModelPart)
    {
        KRATOS_ERROR_IF(!(rModelPart.HasSubModelPart("Contact"))) << "ERROR:: CONTACT MODEL PART NOT CREATED" << std::endl;
        ModelPart& contact_model_part = rModelPart.GetSubModelPart("Contact");
        
        // We release the LM
        auto& nodes_array = contact_model_part.Nodes();
        #pragma omp parallel for 
        for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
            auto it_node = nodes_array.begin() + i;
            if (it_node->Is(ISOLATED) == true) {
                if (it_node->SolutionStepsDataHas(NORMAL_CONTACT_STRESS) == true)
                    it_node->Free(NORMAL_CONTACT_STRESS);
                else { // We will assume that VECTOR_LAGRANGE_MULTIPLIER
                    it_node->Free(VECTOR_LAGRANGE_MULTIPLIER_X);
                    it_node->Free(VECTOR_LAGRANGE_MULTIPLIER_Y);
                    it_node->Free(VECTOR_LAGRANGE_MULTIPLIER_Z);
                }
            }
        }
        
    }
    
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

}; /* Class ContactResidualBasedBlockBuilderAndSolver */

///@}

///@name Type Definitions */
///@{


///@}

} /* namespace Kratos.*/

#endif /* KRATOS_CONTACT_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER  defined */
