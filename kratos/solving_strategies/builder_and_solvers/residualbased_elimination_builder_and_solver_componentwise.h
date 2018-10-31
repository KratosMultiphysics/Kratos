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
#if !defined(KRATOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVERCOMPONENTWISE )
#define  KRATOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVERCOMPONENTWISE


/* System includes */
#include <set>

#ifdef _OPENMP
#include <omp.h>
#endif


/* External includes */


/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"


namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */

/*@} */


/**@name  Enum's */
/*@{ */


/*@} */
/**@name  Functions */
/*@{ */



/*@} */
/**@name Kratos Classes */
/*@{ */

/** Short class definition.

Detail class definition.

This is a specialization of the standard buliding strategy to the case in which a single variable is to be used in the
building.

the creation of the DofList and the construction of the system matrix is in this case much faster
as the neighborhood relationships are considered to be known


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
         class TDenseSpace ,
         class TLinearSolver,
         class TVariableType
         >
class ResidualBasedEliminationBuilderAndSolverComponentwise
    : public ResidualBasedEliminationBuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver >
{
public:
    /**@name Type Definitions */
    /*@{ */
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedEliminationBuilderAndSolverComponentwise );


    typedef BuilderAndSolver<TSparseSpace,TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;


    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    typedef typename BaseType::ElementsContainerType ElementsContainerType;

    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /** Constructor.
    */
    ResidualBasedEliminationBuilderAndSolverComponentwise(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,TVariableType const& Var)
        : ResidualBasedEliminationBuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver >(pNewLinearSystemSolver)
        , rVar(Var)
    {

        /* 			std::cout << "using the standard builder and solver " << std::endl; */

    }


    /** Destructor.
    */
    ~ResidualBasedEliminationBuilderAndSolverComponentwise() override {}


    /*@} */
    /**@name Operators
    */
    /*@{ */



    //**************************************************************************
    //**************************************************************************
    void Build(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& b) override
    {
        KRATOS_TRY
        if(!pScheme)
            KRATOS_THROW_ERROR(std::runtime_error, "No scheme provided!", "");

        //getting the elements from the model
        ElementsArrayType& pElements = r_model_part.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

        //resetting to zero the vector of reactions
        TSparseSpace::SetToZero( *(BaseType::mpReactionsVector) );

        //create a partition of the element array
        int number_of_threads = OpenMPUtils::GetNumThreads();

#ifdef _OPENMP
        int A_size = A.size1();

        //creating an array of lock variables of the size of the system matrix
        std::vector< omp_lock_t > lock_array(A.size1());

        for(int i = 0; i<A_size; i++)
            omp_init_lock(&lock_array[i]);
#endif

        DenseVector<unsigned int> element_partition;
        CreatePartition(number_of_threads, pElements.size(), element_partition);
        if (this->GetEchoLevel()>0)
        {
            KRATOS_WATCH( number_of_threads );
            KRATOS_WATCH( element_partition );
        }


        double start_prod = OpenMPUtils::GetCurrentTime();

        #pragma omp parallel for firstprivate(number_of_threads) schedule(static,1)
        for(int k=0; k<number_of_threads; k++)
        {
            //contributions to the system
            LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0,0);
            LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

            //vector containing the localization in the system of the different
            //terms
            Element::EquationIdVectorType EquationId;
            ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
            typename ElementsArrayType::ptr_iterator it_begin=pElements.ptr_begin()+element_partition[k];
            typename ElementsArrayType::ptr_iterator it_end=pElements.ptr_begin()+element_partition[k+1];

            unsigned int pos = (r_model_part.Nodes().begin())->GetDofPosition(rVar);


            // assemble all elements
            for (typename ElementsArrayType::ptr_iterator it=it_begin; it!=it_end; ++it)
            {

                //calculate elemental contribution
                (*it)->InitializeNonLinearIteration(CurrentProcessInfo);
                (*it)->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

                Geometry< Node<3> >& geom = (*it)->GetGeometry();
                if(EquationId.size() != geom.size()) EquationId.resize(geom.size(),false);

                for(unsigned int i=0; i<geom.size(); i++)
                    EquationId[i] = geom[i].GetDof(rVar,pos).EquationId();

                //assemble the elemental contribution
#ifdef USE_LOCKS_IN_ASSEMBLY
                this->Assemble(A,b,LHS_Contribution,RHS_Contribution,EquationId,lock_array);
#else
                this->Assemble(A,b,LHS_Contribution,RHS_Contribution,EquationId);
#endif
            }
        }

        DenseVector<unsigned int> condition_partition;
        CreatePartition(number_of_threads, ConditionsArray.size(), condition_partition);

        #pragma omp parallel for firstprivate(number_of_threads) schedule(static,1)
        for(int k=0; k<number_of_threads; k++)
        {
            //contributions to the system
            LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0,0);
            LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

            Condition::EquationIdVectorType EquationId;

            ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

            typename ConditionsArrayType::ptr_iterator it_begin=ConditionsArray.ptr_begin()+condition_partition[k];
            typename ConditionsArrayType::ptr_iterator it_end=ConditionsArray.ptr_begin()+condition_partition[k+1];

            unsigned int pos = (r_model_part.Nodes().begin())->GetDofPosition(rVar);

            // A all elements
            for (typename ConditionsArrayType::ptr_iterator it=it_begin; it!=it_end; ++it)
            {

                //calculate elemental contribution
                (*it)->InitializeNonLinearIteration(CurrentProcessInfo);
                (*it)->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

                Geometry< Node<3> >& geom = (*it)->GetGeometry();
                if(EquationId.size() != geom.size()) EquationId.resize(geom.size(),false);

                for(unsigned int i=0; i<geom.size(); i++)
                {
                    EquationId[i] = geom[i].GetDof(rVar,pos).EquationId();
                }

#ifdef USE_LOCKS_IN_ASSEMBLY
                this->Assemble(A,b,LHS_Contribution,RHS_Contribution,EquationId,lock_array);
#else
                this->Assemble(A,b,LHS_Contribution,RHS_Contribution,EquationId);
#endif
            }
        }
        if (this->GetEchoLevel()>0)
        {
            double stop_prod = OpenMPUtils::GetCurrentTime();
            std::cout << "parallel building time: " << stop_prod - start_prod << std::endl;
        }

#ifdef _OPENMP
        for(int i = 0; i<A_size; i++)
            omp_destroy_lock(&lock_array[i]);
#endif

        KRATOS_CATCH("")

    }




    //**************************************************************************
    //**************************************************************************
    void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part
    ) override
    {
        KRATOS_TRY


        //fills a list of "active" nodes defined as nodes which have neighbours
        // AND no fixed pressure
        mActiveNodes.clear();
        mActiveNodes.reserve(r_model_part.Nodes().size() );
        for (typename NodesArrayType::iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); ++it)
        {
            if( (it->GetValue(NEIGHBOUR_NODES)).size() != 0 )
            {
                mActiveNodes.push_back(*(it.base() ));
            }
        }

        //fills the DofList and give a unique progressive tag to each node
        BaseType::mDofSet.clear();
        BaseType::mDofSet.reserve(mActiveNodes.size() );

        for(WeakPointerVector< Node<3> >::iterator iii = mActiveNodes.begin(); iii!=mActiveNodes.end(); iii++)
        {
            BaseType::mDofSet.push_back( iii->pGetDof(rVar).get() );
        }

        //throws an execption if there are no Degrees of freedom involved in the analysis
        if (BaseType::mDofSet.size()==0)
            KRATOS_THROW_ERROR(std::logic_error, "No degrees of freedom!", "");

        BaseType::mDofSetIsInitialized = true;


    // If reactions are to be calculated, we check if all the dofs have reactions defined
    // This is tobe done only in debug mode

    #ifdef KRATOS_DEBUG

    if(BaseType::GetCalculateReactionsFlag())
    {
        for(auto dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
        {
                KRATOS_ERROR_IF_NOT(dof_iterator->HasReaction()) << "Reaction variable not set for the following : " <<std::endl
                    << "Node : "<<dof_iterator->Id()<< std::endl
                    << "Dof : "<<(*dof_iterator)<<std::endl<<"Not possible to calculate reactions."<<std::endl;
        }
    }
    #endif


        KRATOS_CATCH("")
    }


    //**************************************************************************
    //**************************************************************************
    void ResizeAndInitializeVectors(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixPointerType& pA,
        TSystemVectorPointerType& pDx,
        TSystemVectorPointerType& pb,
        ModelPart& rModelPart
    ) override
    {
        KRATOS_TRY

        if(pA == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(0,0) );
            pA.swap(pNewA);
        }
        if(pDx == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewDx = TSystemVectorPointerType(new TSystemVectorType(0) );
            pDx.swap(pNewDx);
        }
        if(pb == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewb = TSystemVectorPointerType(new TSystemVectorType(0) );
            pb.swap(pNewb);
        }
        if(BaseType::mpReactionsVector == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewReactionsVector = TSystemVectorPointerType(new TSystemVectorType(0) );
            BaseType::mpReactionsVector.swap(pNewReactionsVector);
        }

        TSystemMatrixType& A  = *pA;
        TSystemVectorType& Dx = *pDx;
        TSystemVectorType& b  = *pb;

        //resizing the system vectors and matrix
        if (A.size1() == 0 || BaseType::GetReshapeMatrixFlag() == true) //if the matrix is not initialized
        {
            A.resize(BaseType::mEquationSystemSize,BaseType::mEquationSystemSize,false);
#ifdef _OPENMP
            ParallelConstructGraph(A);
#else
            ConstructGraph(A);
#endif
        }
        else
        {
            if(A.size1() != BaseType::mEquationSystemSize || A.size2() != BaseType::mEquationSystemSize)
            {
                //KRATOS_WATCH("it should not come here!!!!!!!! ... this is SLOW");
                KRATOS_ERROR <<"The equation system size has changed during the simulation. This is not permited."<<std::endl;
                A.resize(BaseType::mEquationSystemSize,BaseType::mEquationSystemSize,true);
#ifdef _OPENMP
                ParallelConstructGraph(A);
#else
                ConstructGraph(A);
#endif
            }
        }
        if(Dx.size() != BaseType::mEquationSystemSize)
            Dx.resize(BaseType::mEquationSystemSize,false);
        if(b.size() != BaseType::mEquationSystemSize)
            b.resize(BaseType::mEquationSystemSize,false);

        //


        //if needed resize the vector for the calculation of reactions
        if(BaseType::mCalculateReactionsFlag == true)
        {
            unsigned int ReactionsVectorSize = BaseType::mDofSet.size();
            if(BaseType::mpReactionsVector->size() != ReactionsVectorSize)
                BaseType::mpReactionsVector->resize(ReactionsVectorSize,false);
        }

        //swapping pointers
// 				pA.swap(pNewA);
// 				pDx.swap(pNewDx);
// 				pb.swap(pNewb);
#ifndef __SUNPRO_CC
        KRATOS_CATCH("")
#endif

    }

    //**************************************************************************
    //**************************************************************************
    void Clear() override
    {
        this->mDofSet = DofsArrayType();

        if(this->mpReactionsVector != NULL)
        {
            TSparseSpace::Clear( (this->mpReactionsVector) );
        }
// 			*(this->mpReactionsVector) = TSystemVectorType();

        if (this->GetEchoLevel()>1)
        {
            KRATOS_WATCH("ResidualBasedEliminationBuilderAndSolver Clear Function called");
        }
    }
    /*@} */
    /**@name Operations */
    /*@{ */


    /*@} */
    /**@name Access */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ResidualBasedEliminationBuilderAndSolverComponentwise";
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

    /*@} */
    /**@name Friends */
    /*@{ */


    /*@} */

protected:
    /**@name Protected static Member Variables */
    /*@{ */


    /*@} */
    /**@name Protected member Variables */
    /*@{ */


    /*@} */
    /**@name Protected Operators*/
    /*@{ */
    //**************************************************************************
    //**************************************************************************
    //**************************************************************************
    //**************************************************************************
    void ConstructGraph(TSystemMatrixType& A)
    {
        KRATOS_TRY

        std::vector< std::vector<int> > index_list(BaseType::mEquationSystemSize);

        int total_size = 0;

        unsigned int pos = (mActiveNodes.begin())->GetDofPosition(rVar);
        //constructing the system matrix row by row
        int index_i;
        for(WeakPointerVector< Node<3> >::iterator in = mActiveNodes.begin();
                in!=mActiveNodes.end(); in++)
        {
            Node<3>::DofType& current_dof = in->GetDof(rVar,pos);
            if( current_dof.IsFixed() == false)
            {
                index_i = (current_dof).EquationId();
                WeakPointerVector< Node<3> >& neighb_nodes = in->GetValue(NEIGHBOUR_NODES);

                std::vector<int>& indices = index_list[index_i];
                indices.reserve(neighb_nodes.size()+1);

                //filling the first neighbours list
                indices.push_back(index_i);
                for( WeakPointerVector< Node<3> >::iterator i =	neighb_nodes.begin();
                        i != neighb_nodes.end(); i++)
                {

                    Node<3>::DofType& neighb_dof = i->GetDof(rVar,pos);
                    if(neighb_dof.IsFixed() == false )
                    {
                        int index_j = (neighb_dof).EquationId();
                        indices.push_back(index_j);
                    }
                }

                //sorting the indices and elminating the duplicates
                std::sort(indices.begin(),indices.end());
                typename std::vector<int>::iterator new_end = std::unique(indices.begin(),indices.end());

                indices.erase(new_end,indices.end());

                total_size += indices.size();
            }
        }

        A.reserve(total_size,false);

        //setting to zero the matrix (and the diagonal matrix)
        for(unsigned int i=0; i<BaseType::mEquationSystemSize; i++)
        {
            std::vector<int>& indices = index_list[i];
            for(unsigned int j=0; j<indices.size(); j++)
            {
                A.push_back(i,indices[j] , 0.00);
            }
        }

        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************
    //**************************************************************************
    //**************************************************************************
#ifdef _OPENMP
    void ParallelConstructGraph(TSystemMatrixType& A)
    {
#ifndef __SUNPRO_CC
        KRATOS_TRY
#endif
        std::vector< std::vector<int> > index_list(BaseType::mEquationSystemSize);

        int number_of_threads = omp_get_max_threads();

        unsigned int pos = (mActiveNodes.begin())->GetDofPosition(rVar);
        //constructing the system matrix row by row

        DenseVector<unsigned int> partition;
        DenseVector<unsigned int> local_sizes(number_of_threads);
        for(int i=0; i<number_of_threads; i++)
            local_sizes[i] = 0;

        CreatePartition(number_of_threads, mActiveNodes.size(), partition);

        #pragma omp parallel for firstprivate(number_of_threads,pos) schedule(static,1)
        for(int k=0; k<number_of_threads; k++)
        {
            WeakPointerVector< Node<3> >::iterator it_begin = mActiveNodes.begin()+partition[k];
            WeakPointerVector< Node<3> >::iterator it_end = mActiveNodes.begin()+partition[k+1];

            for(WeakPointerVector< Node<3> >::iterator in = it_begin;
                    in!=it_end; in++)
            {
                Node<3>::DofType& current_dof = in->GetDof(rVar,pos);
                if( current_dof.IsFixed() == false)
                {
                    int index_i = (current_dof).EquationId();
                    WeakPointerVector< Node<3> >& neighb_nodes = in->GetValue(NEIGHBOUR_NODES);

                    std::vector<int>& indices = index_list[index_i];
                    indices.reserve(neighb_nodes.size()+1);

                    //filling the first neighbours list
                    indices.push_back(index_i);
                    for( WeakPointerVector< Node<3> >::iterator i =	neighb_nodes.begin();
                            i != neighb_nodes.end(); i++)
                    {

                        Node<3>::DofType& neighb_dof = i->GetDof(rVar,pos);
                        if(neighb_dof.IsFixed() == false )
                        {
                            int index_j = (neighb_dof).EquationId();
                            indices.push_back(index_j);
                        }
                    }

                    //sorting the indices and elminating the duplicates
                    std::sort(indices.begin(),indices.end());
                    typename std::vector<int>::iterator new_end = std::unique(indices.begin(),indices.end());
                    indices.erase(new_end,indices.end());

                    local_sizes[k] += indices.size();
                }
            }
        }

        //calculate the total size of the system
        int total_size = 0.0;
        for(int i=0; i<number_of_threads; i++)
            total_size += local_sizes[i];

        A.reserve(total_size,false);

        //setting to zero the matrix (and the diagonal matrix)
        for(unsigned int i=0; i<BaseType::mEquationSystemSize; i++)
        {
            std::vector<int>& indices = index_list[i];
            for(unsigned int j=0; j<indices.size(); j++)
            {
                A.push_back(i,indices[j] , 0.00);
            }
        }
#ifndef __SUNPRO_CC
        KRATOS_CATCH("")
#endif
    }
#endif



    /*@} */
    /**@name Protected Operations*/
    /*@{ */


    /*@} */
    /**@name Protected  Access */
    /*@{ */


    /*@} */
    /**@name Protected Inquiry */
    /*@{ */


    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */



    /*@} */

private:
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */
    TVariableType const & rVar;
    WeakPointerVector<Node<3> > mActiveNodes;

    /*@} */
    /**@name Private Operators*/
    /*@{ */
    //******************************************************************************************
    //******************************************************************************************
    inline void CreatePartition(unsigned int number_of_threads,const int number_of_rows, DenseVector<unsigned int>& partitions)
    {
        partitions.resize(number_of_threads+1);
        int partition_size = number_of_rows / number_of_threads;
        partitions[0] = 0;
        partitions[number_of_threads] = number_of_rows;
        for(unsigned int i = 1; i<number_of_threads; i++)
            partitions[i] = partitions[i-1] + partition_size ;
    }





    /*@} */
    /**@name Private Operations*/
    /*@{ */






    /*@} */
    /**@name Private  Access */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */


    /*@} */

}; /* Class ResidualBasedEliminationBuilderAndSolverComponentwise */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVERCOMPONENTWISE  defined */
