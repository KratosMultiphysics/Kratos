/* *********************************************************
*
*   Last Modified by:    $Author: rrossi $
*   Date:                $Date: 2008-11-10 14:23:32 $
*   Revision:            $Revision: 1.17 $
*
* ***********************************************************/


#if !defined(KRATOS_RESIDUAL_BASED_ELIMINATION_DISCRETE_LAPLACIAN_BUILDER_AND_SOLVER )
#define  KRATOS_RESIDUAL_BASED_ELIMINATION_DISCRETE_LAPLACIAN_BUILDER_AND_SOLVER


/* System includes */
#include <set>


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "utilities/geometry_utilities.h"


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
         int TDim
         >
class ResidualBasedEliminationDiscreteLaplacianBuilderAndSolver
    : public BuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver >
{
public:
    /**@name Type Definitions */
    /*@{ */
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedEliminationDiscreteLaplacianBuilderAndSolver );


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
    ResidualBasedEliminationDiscreteLaplacianBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver, int estimated_number_of_second_neighb = 125
    )
        : BuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver >(pNewLinearSystemSolver)
    {
        mestimated_number_of_second_neighb = estimated_number_of_second_neighb;

    }


    /** Destructor.
    */
    virtual ~ResidualBasedEliminationDiscreteLaplacianBuilderAndSolver() {}


    /*@} */
    /**@name Operators
    */
    /*@{ */


    //**************************************************************************
    //**************************************************************************
    void SystemSolve(
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    )
    {
        KRATOS_TRY

        double norm_b = TSparseSpace::TwoNorm(b);

        if(norm_b != 0.00)
            BaseType::mpLinearSystemSolver->Solve(A,Dx,b);
        else
            TSparseSpace::SetToZero(Dx);

        //prints informations about the current time
        if (BaseType::GetEchoLevel()>1)
        {
            std::cout << *(BaseType::mpLinearSystemSolver) << std::endl;
        }

        KRATOS_CATCH("")

    }

    //**************************************************************************
    //**************************************************************************
    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        KRATOS_TRY

        boost::timer building_time;

        Build(pScheme,r_model_part,A,b);

        if(BaseType::GetEchoLevel()>0)
        {
            std::cout << "Building Time : " << building_time.elapsed() << std::endl;
        }

        if (BaseType::GetEchoLevel()== 3)
        {
            std::cout << "before the solution of the system" << std::endl;
            std::cout << "System Matrix = " << A << std::endl;
            std::cout << "unknowns vector = " << Dx << std::endl;
            std::cout << "RHS vector = " << b << std::endl;
        }

        boost::timer solve_time;

        SystemSolve(A,Dx,b);

        if(BaseType::GetEchoLevel()>0)
        {
            std::cout << "System Solve Time : " << solve_time.elapsed() << std::endl;
        }

        KRATOS_CATCH("")
    }



    //**************************************************************************
    //**************************************************************************
    //this is done in a purely nodal way taking advantage of the neighbour relatinoships
    //which are assumed to be calculated separately
    void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part
    )
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

        WeakPointerVector<Node<3> > aux;
        aux.resize(mActiveNodes.size());

        //getting the dof position
        //unsigned int dof_position = (mActiveNodes.begin())->GetDofPosition(PRESSURE);

        //fills the DofList and give a unique progressive tag to each node
        BaseType::mDofSet.clear();
        BaseType::mDofSet.reserve(mActiveNodes.size() );
        int free_id = 0;
        int fix_id = mActiveNodes.size();
        for(unsigned int i=0; i<mActiveNodes.size(); i++)
        {
            typename Dof< double >::Pointer pDof = mActiveNodes[i].pGetDof(PRESSURE);
            int index;
            if(pDof->IsFixed())
                index = --fix_id;
            else
                index = free_id++;
            BaseType::mDofSet.push_back( pDof.get() );
            pDof->SetEquationId(index);
            aux(index) = mActiveNodes(i);

        }

        //reordering the active nodes according to their fixity
        for(unsigned int i=0; i<mActiveNodes.size(); i++)
        {
            mActiveNodes(i) = aux(i);
        }

        this->mEquationSystemSize = mActiveNodes.size();
        mActiveSize = free_id;


        //throws an execption if there are no Degrees of freedom involved in the analysis
        if (BaseType::mDofSet.size()==0)
            KRATOS_THROW_ERROR(std::logic_error, "No degrees of freedom!", "");

        BaseType::mDofSetIsInitialized = true;

        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************
    void SetUpSystem(
        ModelPart& r_model_part
    )
    {
        KRATOS_TRY

        KRATOS_CATCH("");
    }

    //**************************************************************************
    //**************************************************************************
    void ResizeAndInitializeVectors( typename TSchemeType::Pointer pScheme,
        TSystemMatrixPointerType& pA,
        TSystemVectorPointerType& pDx,
        TSystemVectorPointerType& pb,
        ModelPart& r_model_part
    )
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
        if (A.size1() == 0 || this->GetReshapeMatrixFlag() == true) //if the matrix is not initialized
        {
            A.resize(this->mEquationSystemSize,this->mEquationSystemSize,false);
            ConstructMatrixStructure(A);
        }
        else
        {
            if(A.size1() != this->mEquationSystemSize || A.size2() != 	this->mEquationSystemSize)
            {
                KRATOS_THROW_ERROR(std::logic_error,"it should not resize the matrix!!","")
            }
        }
        if(Dx.size() != this->mEquationSystemSize)
            Dx.resize(this->mEquationSystemSize,false);
        if(b.size() != this->mEquationSystemSize)
            b.resize(this->mEquationSystemSize,false);

        KRATOS_CATCH("")

    }

    //**************************************************************************
    //**************************************************************************
    void Build(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& b)
    {
        KRATOS_TRY
        //typedef typename unsigned int size_type;
        //typedef typename double value_type;
        if (r_model_part.NodesBegin()->SolutionStepsDataHas(ARRHENIUS)==false )
            KRATOS_THROW_ERROR(std::logic_error,"Add  ----ARRHENIUS---- variable!!!!!! ERROR","");


        const Vector& BDFcoeffs = r_model_part.GetProcessInfo()[BDF_COEFFICIENTS];

        //resetting the nodal area and the TAU
        for(WeakPointerVector< Node<3> >::iterator in = mActiveNodes.begin();
                in!=mActiveNodes.end(); in++)
        {
            in->FastGetSolutionStepValue(NODAL_MASS) = 0.00;
        }

        boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim> DN_DX;
        boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim+1> elemental_stabilization;
        array_1d<double,TDim+1> N;
        array_1d<unsigned int ,TDim+1> local_indices;
        array_1d<double,3> vel;
        array_1d<double,TDim> proj_aux;
        array_1d<double,TDim+1> rhs_contribution;

        //getting the dof position
        unsigned int dof_position = (mActiveNodes.begin())->GetDofPosition(PRESSURE);

        //creating gradient matrix
        TSystemMatrixType G;
        ConstructG_Structure(G);

        double aaa = 1.0/(TDim+1.0);

        double dynamic_tau = r_model_part.GetProcessInfo()[DYNAMIC_TAU];

        array_1d<double,3> vg;
        for(ModelPart::ElementsContainerType::iterator i = r_model_part.ElementsBegin();
                i!=r_model_part.ElementsEnd(); i++)
        {

            Geometry< Node<3> >& geom = i->GetGeometry();

            //calculating elemental values
            double Volume;
            //CalculateGeometryData3D(geom,DN_De,J,Jinv,DN_DX,N,Volume);
            GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);

            //determining elemental density and viscosity
            double density = geom[0].FastGetSolutionStepValue(DENSITY);
            double nu = geom[0].FastGetSolutionStepValue(VISCOSITY);
            for(int ii = 1; ii<TDim+1; ii++)
            {
                density += geom[ii].FastGetSolutionStepValue(DENSITY);
                nu += geom[ii].FastGetSolutionStepValue(VISCOSITY);
            }
            density *= aaa;
            nu *= aaa;

            //finiding local indices
            for(int ii = 0; ii<TDim+1; ii++)
            {
                local_indices[ii] = geom[ii].GetDof(PRESSURE,dof_position).EquationId();
            }

            //building matrix G (gradient integrated by parts)
            for(unsigned int ii = 0; ii<TDim+1; ii++)
            {
                unsigned int i1 = local_indices[ii];
                for(unsigned int kkk = 0; kkk<TDim; kkk++)
                {
                    unsigned int row_pos = i1*TDim + kkk;
                    for(unsigned int jj = 0; jj<TDim+1; jj++)
                    {
                        unsigned int i2 = local_indices[jj];
                        double temp = Volume*N[jj];
                        G(row_pos,i2) += temp * DN_DX(ii,kkk);
                    }
                }
            }

            //calculating elemental area, and "gauss point velocity"
            noalias(vg) = ZeroVector(3);; //velocity on the gauss points
            for(int kk = 0; kk<TDim+1; kk++)
            {
                //adding the elemental contribution to the nodal volume
                geom[kk].FastGetSolutionStepValue(NODAL_MASS) += density*Volume*N[kk];

                //calculating the velocity on the gauss point
                noalias(vg) +=	N[kk]* ( geom[kk].FastGetSolutionStepValue(FRACT_VEL) -
                                         geom[kk].FastGetSolutionStepValue(MESH_VELOCITY) );
            }

            //*******************************************************************************************
            //adding stabilization
            //*******************************************************************************************
            double h;

//				double dt_contrib_to_tau = 0.0;
//				if( muse_dt_in_stabilization == true)
//					dt_contrib_to_tau = 1.0/BDFcoeffs[0];

            if(TDim == 2)
                h = sqrt(2.0*Volume);
            else
                h = pow(6.00*Volume,0.33333333333);
            double c1 = 4.00;
            double c2 = 2.00;
            double norm_u = norm_2(vg);
            double tau = 1.00 / (dynamic_tau*BDFcoeffs[0] + c1*nu/(h*h) + c2*norm_u/h );



            //calculating stabilization laplacian LHS
            double laplacian_coeff = (tau)/density;
            noalias(elemental_stabilization) =  laplacian_coeff*prod(DN_DX,trans(DN_DX));

            const array_1d<double,3>& proj_temp = geom[0].FastGetSolutionStepValue(PRESS_PROJ);
            for(int iii = 0; iii<TDim; iii++)
                proj_aux[iii] = N[0]*proj_temp[iii];
            for(int kk = 1; kk<TDim+1; kk++)
            {
                const array_1d<double,3>& proj_temp = geom[kk].FastGetSolutionStepValue(PRESS_PROJ);
                for(int iii = 0; iii<TDim; iii++)
                    proj_aux[iii] += N[kk]*proj_temp[iii];
            }
            proj_aux *= tau;
            noalias(rhs_contribution) = prod(DN_DX , proj_aux);

            double Gaux = 0.00;
            for(int kk = 0; kk<TDim+1; kk++)
            {
                for(int tt = 0; tt<TDim; tt++)
                {
                    const array_1d<double,3>& fv = geom[kk].FastGetSolutionStepValue(FRACT_VEL);
                    Gaux += DN_DX(kk,tt)*fv[tt];
                    Gaux -= geom[kk].FastGetSolutionStepValue(ARRHENIUS);
                }
            }
            noalias(rhs_contribution) -= Gaux * N;



            //adding the stabilization laplacian to the system matrix and the stabilizing contribution to the LHS
            for(int ii = 0; ii<TDim+1; ii++)
            {
                int i1 = local_indices[ii];
                for(int jj = 0; jj<TDim+1; jj++)
                {
                    int i2 = local_indices[jj];
                    A(i1,i2) += Volume * elemental_stabilization(ii,jj);
                }
                b[i1] += Volume * rhs_contribution[ii];
            }
        }

        //auxiliary vector definitions
        TSystemVectorType aux_vect_large(this->mEquationSystemSize*TDim);
        TSystemVectorType aux_vect_small(this->mEquationSystemSize);

        //computing the diagonal matrices to be used in the calculation of DMinvG
        TSystemVectorType Minv_dt(this->mEquationSystemSize*TDim);
        TSystemVectorType Minv_stabilized(this->mEquationSystemSize*TDim);
        ComputeMinvMatrices(BDFcoeffs[0],Minv_dt,Minv_stabilized);

        ////adding pold term
        // RHS += Gtrans*diag*G* pold BUT no other terms added to the LHS
        AddToRHS_DMGdt(Minv_dt,G,b,aux_vect_small,aux_vect_large);



        //RHS - (L + DMinvG)*p
        //note that up to this point A contains only the Laplacian part, not yet the DMinvG
        ModifyForDirichlet(A,G,b,Minv_dt,aux_vect_small,aux_vect_large);
//			ModifyForDirichlet(A,G,b,Minv_stabilized,aux_vect_small,aux_vect_large);

        //modifying the G Matrix to include the Dt effects.
        //G is modified so that Gtrans*G = Gtrans * diag * G
        //which implies that G = sqrt(diag)*G
        modify_G(G,Minv_stabilized);

        ////adding the term Gtrans * dt/area * G to the LHS (this can be implemented as Gtrans * G as G was modified to include the diag)
        add_symmetric_prod(G,A);


        KRATOS_CATCH("");
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
    WeakPointerVector<Node<3> > mActiveNodes;
    unsigned int mActiveSize;
    unsigned int mestimated_number_of_second_neighb;
    bool muse_dt_in_stabilization;

    /*@} */
    /**@name Private Operators*/
    /*@{ */


    /*@} */
    /**@name Private Operations*/
    /*@{ */

    //**************************************************************************
    void ConstructMatrixStructure(
        TSystemMatrixType& A
    )
    {
        KRATOS_TRY

        boost::timer temp;

        std::vector< int > work_array;
        work_array.reserve(1000);

        //guessing the total size of the matrix

        A.reserve(mActiveNodes.size() * mestimated_number_of_second_neighb, false);

        //getting the dof position
        unsigned int dof_position = (mActiveNodes.begin())->GetDofPosition(PRESSURE);

        //building up the matrix graph row by row
        int total_size = 0;
        int row_number = 0;
        for(WeakPointerVector< Node<3> >::iterator in = mActiveNodes.begin();
                in!=mActiveNodes.end(); in++)
        {
            int index_i = in->GetDof(PRESSURE,dof_position).EquationId();

            WeakPointerVector< Node<3> >& neighb_nodes = in->GetValue(NEIGHBOUR_NODES);

            //filling the first neighbours list
            work_array.push_back(index_i);
            for( WeakPointerVector< Node<3> >::iterator i =	neighb_nodes.begin();
                    i != neighb_nodes.end(); i++)
            {
                int index_j = i->GetDof(PRESSURE,dof_position).EquationId();
                work_array.push_back(index_j);
            }

            //adding the second neighbours
            for( WeakPointerVector< Node<3> >::iterator i =	neighb_nodes.begin();
                    i != neighb_nodes.end(); i++)
            {
                WeakPointerVector< Node<3> >& second_neighbours = i->GetValue(NEIGHBOUR_NODES);
                for( WeakPointerVector< Node<3> >::iterator jjj =	second_neighbours.begin();
                        jjj != second_neighbours.end(); jjj++)
                {
                    int second_neighb_index = (jjj->GetDof(PRESSURE,dof_position)).EquationId();
                    work_array.push_back(second_neighb_index);
                }
            }


            //sorting the indices and elminating the duplicates
            std::sort(work_array.begin(),work_array.end());
            typename std::vector<int>::iterator new_end = std::unique(work_array.begin(),work_array.end());
            unsigned int number_of_entries = new_end - work_array.begin();

            //filling up the matrix
            for(unsigned int j=0; j<number_of_entries; j++)
            {
                A.push_back(row_number,work_array[j] , 0.00);
            }
            row_number += 1;

            //clearing the array for the next step
            work_array.erase(work_array.begin(),work_array.end());
            total_size += number_of_entries;
        }

        double measured_avg_second_neighbours = double(total_size) / double(mActiveNodes.size() );
        std::cout << "mesured average number of second neighbours" << measured_avg_second_neighbours << std::endl;
        std::cout << "estimated average number of second neighbours" << mestimated_number_of_second_neighb << std::endl;

        if(measured_avg_second_neighbours > mestimated_number_of_second_neighb)
        {
            std::cout << "**********************************************************************" << std::endl;
            std::cout << "attention the number of second neighbours estimated is too low...INEFFICIENT!!!" << std::endl;
            std::cout << "**********************************************************************" << std::endl;
        }

        //setting the estimate neighbour numbers to a "safe value" for the next steps
        mestimated_number_of_second_neighb = int(measured_avg_second_neighbours * 1.1);

        std::cout << "construction time " << temp.elapsed() << std::endl;
        KRATOS_CATCH("")
    }

    //*************************************************************************************************
    //*************************************************************************************************
    //**************************************************************************
    void ConstructG_Structure(
        TSystemMatrixType& G
    )
    {
        KRATOS_TRY

        G.resize(this->mEquationSystemSize*TDim,this->mEquationSystemSize,false);
        //getting the dof position
        unsigned int dof_position = (mActiveNodes.begin())->GetDofPosition(PRESSURE);

        std::vector<int> work_array;
        work_array.reserve(200);

        //calculate total size of G
        int tot_size = 0;
        for(WeakPointerVector< Node<3> >::iterator in = mActiveNodes.begin();
                in!=mActiveNodes.end(); in++)
        {
            tot_size += (in->GetValue(NEIGHBOUR_NODES)).size();
        }
        tot_size *= TDim;

        G.reserve(tot_size,false);

        int row_number = 0;
        for(WeakPointerVector< Node<3> >::iterator in = mActiveNodes.begin();
                in!=mActiveNodes.end(); in++)
        {

            int index_i = in->GetDof(PRESSURE,dof_position).EquationId();
            work_array.push_back(index_i);

            WeakPointerVector< Node<3> >& neighb_nodes = in->GetValue(NEIGHBOUR_NODES);

            //filling the first neighbours list
            for( WeakPointerVector< Node<3> >::iterator i =	neighb_nodes.begin();
                    i != neighb_nodes.end(); i++)
            {
                int index_j = i->GetDof(PRESSURE,dof_position).EquationId();
                work_array.push_back(index_j);
            }

            //sorting the indices and elminating the duplicates
            std::sort(work_array.begin(),work_array.end());

            for(unsigned int k=0; k<TDim; k++)
            {
                unsigned int base_index = row_number*TDim;
                for(unsigned int j=0; j<work_array.size(); j++)
                {
                    G.push_back(base_index + k,work_array[j], 0.00);
                }
            }
            row_number += 1;


            //clearing the array for the next step
            work_array.erase(work_array.begin(),work_array.end());
        }

        KRATOS_CATCH("")
    }




    //*************************************************************************************************
    //*************************************************************************************************
    //output += trans(input)*input
    void add_symmetric_prod (const TSystemMatrixType& input,
                             TSystemMatrixType& output)
    {
        KRATOS_TRY
        typedef  unsigned int size_type;
        typedef  double value_type;
        double temp;

        for (size_type k = 0; k < input.size1 (); ++ k)
        {
            size_type begin = input.index1_data () [k];
            size_type end = input.index1_data () [k + 1];


            for (size_type i = begin; i < end; ++ i)
            {
                unsigned int index_i = input.index2_data () [i];
                value_type data_i = input.value_data()[i];

                //diagonal term
                output(index_i,index_i) += data_i * data_i;

                //take advantage of symmetry in the calculation of the outdiagonal terms
                for (size_type j = i+1; j < end; ++ j)
                {
                    unsigned int index_j = input.index2_data () [j];
                    temp = data_i * input.value_data()[j];

                    output(index_i,index_j) += temp;
                    output(index_j,index_i) += temp;
                }
            }
        }
        KRATOS_CATCH("");
    }

    //*************************************************************************************************
    //*************************************************************************************************
    void ModifyForDirichlet (TSystemMatrixType& A,
                             const TSystemMatrixType& G,
                             TSystemVectorType& b,
                             const TSystemVectorType& Minv_dt,
                             TSystemVectorType& aux_vect_small,
                             TSystemVectorType& aux_vect_large
                            )
    {
        KRATOS_TRY

        //getting the dof position
        unsigned int dof_position = (mActiveNodes.begin())->GetDofPosition(PRESSURE);

        typedef  unsigned int size_type;
        //typedef  double value_type;

        for(WeakPointerVector< Node<3> >::iterator in = mActiveNodes.begin();
                in!=mActiveNodes.end(); in++)
        {
            int index_i = in->GetDof(PRESSURE,dof_position).EquationId();

            //aux_vect_small[index_i] = in->FastGetSolutionStepValue(PRESSURE,1);
            aux_vect_small[index_i] = in->FastGetSolutionStepValue(PRESSURE);
        }

        //RHS - L * p
        axpy_prod(-aux_vect_small,A,b,false);

        // here we do RHS - DMinv_DtG p) through a number of different steps:
        // aux = -G*p
        axpy_prod(G,-aux_vect_small,aux_vect_large,true);

        // aux *= Minv_dt
        for(unsigned int i = 0; i<b.size(); i++)
        {
            aux_vect_large[i] *= Minv_dt[i];
        }

        // b = b - D * aux
        axpy_prod(aux_vect_large,G, b,false);



        //set to zero all the terms where
        //index2 > mActiveSize && index1<mActiveSize
        size_type begin = A.index1_data () [0];
        size_type end = A.index1_data () [mActiveSize];
        for (size_type i = begin; i < end; ++ i)
        {
            unsigned int index_j = A.index2_data () [i];
            if(index_j >= mActiveSize)
            {
                A.value_data()[i] = 0.00;
            }
        }

        //set to zero all the remaining lines
        begin = A.index1_data () [mActiveSize];
        end = A.index1_data () [this->mEquationSystemSize];
        for (size_type i = begin; i < end; ++ i)
        {
            A.value_data()[i] = 0.00;
        }


        //calculate the average of the diagonal numbers
        double avg_diag = 0.00;
        for(unsigned int i = 0; i < mActiveSize ; i++)
        {
            avg_diag += A(i,i);
        }
        avg_diag /= mActiveSize;

        //put a 1 on the diagonals and 0 on the RHS
        for(unsigned int i = mActiveSize ; i<this->mEquationSystemSize; i++)
        {
            A(i,i) = avg_diag; //1.00;
            b[i] = 0.00;
        }


        KRATOS_CATCH("");
    }

    //*************************************************************************************************
    //*************************************************************************************************
    //Calculate the inverse "mass" matrices, one containing only the dt and the stabilization and the other one also the stabilization for FSI
    void ComputeMinvMatrices(const double& BDF0, TSystemVectorType& Minv_dt, TSystemVectorType& Minv_stab)
    {
        //getting the dof position
        unsigned int dof_position = (mActiveNodes.begin())->GetDofPosition(PRESSURE);
        array_1d<double,3> vg;

        //forming an auxiliary vector with the base
        //Minv_dt = (dt/(density*Area))
        if(TDim == 2)
        {
            for(WeakPointerVector< Node<3> >::iterator in = mActiveNodes.begin();
                    in!=mActiveNodes.end(); in++)
            {
                int index_i = in->GetDof(PRESSURE,dof_position).EquationId();
                int base_index = index_i*TDim;
                double nodal_mass = in->FastGetSolutionStepValue(NODAL_MASS);
                double tmp = (1.00/(BDF0*nodal_mass) );

                if(in->IsFixed(VELOCITY_X))
                {
                    Minv_dt[base_index] = 0.0;
                }
                else
                {
                    Minv_dt[base_index] = tmp ;
                }

                if(in->IsFixed(VELOCITY_Y))
                {
                    Minv_dt[base_index+1] = 0.0;
                }
                else
                {
                    Minv_dt[base_index+1] = tmp;
                }

                //stabilize for FSI as needed
                Minv_stab[base_index] = Minv_dt[base_index] ;
                Minv_stab[base_index+1] = Minv_dt[base_index+1];

                if(in->FastGetSolutionStepValue(IS_INTERFACE) > 0)
                {
                    Minv_stab[base_index] += tmp;
                    Minv_stab[base_index+1] += tmp;
                }
            }
        }
        else if (TDim == 3)
        {
            for(WeakPointerVector< Node<3> >::iterator in = mActiveNodes.begin();
                    in!=mActiveNodes.end(); in++)
            {
                int index_i = in->GetDof(PRESSURE,dof_position).EquationId();
                int base_index = index_i*TDim;
                double nodal_mass = in->FastGetSolutionStepValue(NODAL_MASS);
                double tmp = (1.00/(BDF0*nodal_mass) );

                if(in->IsFixed(VELOCITY_X))	Minv_dt[base_index] = 0.0;
                else	Minv_dt[base_index] = tmp ;

                if(in->IsFixed(VELOCITY_Y))	Minv_dt[base_index+1] = 0.0;
                else	Minv_dt[base_index+1] = tmp;

                if(in->IsFixed(VELOCITY_Z))	Minv_dt[base_index+2] = 0.0;
                else	Minv_dt[base_index+2] = tmp;

                //stabilize for FSI as needed
                Minv_stab[base_index] = Minv_dt[base_index] ;
                Minv_stab[base_index+1] = Minv_dt[base_index+1];
                Minv_stab[base_index+2] = Minv_dt[base_index+2];
                if(in->FastGetSolutionStepValue(IS_INTERFACE) > 0)
                {
                    Minv_stab[base_index] += tmp;
                    Minv_stab[base_index+1] += tmp;
                    Minv_stab[base_index+2] += tmp;
                }
            }
        }

    }

    //*************************************************************************************************
    //*************************************************************************************************
    //G is modified so that Gtrans*G = Gtrans * diag * G
    //    ====>   G = sqrt(diag)*G
    void modify_G(TSystemMatrixType& G,const TSystemVectorType& Minv_stab)
    {
        //modifying G taking in account the diag vector
        //G = diag * G
        typedef unsigned int size_type;
        //typedef double value_type;

        for (size_type k = 0; k < G.size1 (); ++ k)
        {
            size_type begin = G.index1_data () [k];
            size_type end = G.index1_data () [k + 1];
            double diag_term = sqrt(Minv_stab[k]);

            for (size_type i = begin; i < end; ++ i)
            {
                G.value_data()[i] *= diag_term;
            }
        }

    }

    //*************************************************************************************************
    //*************************************************************************************************
    void AddToRHS_DMGdt(	TSystemVectorType& Minv_dt,
                            TSystemMatrixType& G,
                            TSystemVectorType& b,
                            TSystemVectorType& aux_vect_small,
                            TSystemVectorType& aux_vect_large)
    {
        //getting the dof position
        unsigned int dof_position = (mActiveNodes.begin())->GetDofPosition(PRESSURE);

        ////adding pold term
        // RHS += Gtrans*diag*G* pold BUT no other terms added to the LHS
        for(WeakPointerVector< Node<3> >::iterator in = mActiveNodes.begin();
                in!=mActiveNodes.end(); in++)
        {
            int index_i = in->GetDof(PRESSURE,dof_position).EquationId();

            //aux_vect_small[index_i] = in->FastGetSolutionStepValue(PRESSURE,1);
            aux_vect_small[index_i] = in->FastGetSolutionStepValue(PRESSURE_OLD_IT);
        }
        axpy_prod(G,aux_vect_small,aux_vect_large,true);

        //
        for(unsigned int i = 0; i<b.size(); i++)
        {
            aux_vect_large[i] *= Minv_dt[i];
        }

        //
        axpy_prod(aux_vect_large,G, b,false);

    }


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

}; /* Class ResidualBasedEliminationDiscreteLaplacianBuilderAndSolver */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_ELIMINATION_DISCRETE_LAPLACIAN_BUILDER_AND_SOLVER  defined */

