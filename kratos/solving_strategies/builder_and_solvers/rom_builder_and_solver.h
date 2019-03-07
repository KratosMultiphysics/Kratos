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
#if !defined(KRATOS_ROM_BUILDER_AND_SOLVER )
#define  KRATOS_ROM_BUILDER_AND_SOLVER

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"

//default linear solver
//#include "linear_solvers/linear_solver.h"

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
         class TDenseSpace, // = DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ROMBuilderAndSolver : public BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>
{
public:

    /**
     * This struct is used in the component wise calculation only
     * is defined here and is used to declare a member variable in the component wise builder and solver
     * private pointers can only be accessed by means of set and get functions
     * this allows to set and not copy the Element_Variables and Condition_Variables
     * which will be asked and set by another strategy object
     */

    //pointer definition

    KRATOS_CLASS_POINTER_DEFINITION(ROMBuilderAndSolver);

    // The size_t types
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    /// Definition of the classes from the base class
    typedef BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;
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

    /// Additional definitions
    typedef PointerVectorSet<Element, IndexedObject> ElementsContainerType;
    typedef Element::EquationIdVectorType EquationIdVectorType;
    typedef Element::DofsVectorType DofsVectorType;
    typedef boost::numeric::ublas::compressed_matrix<double> CompressedMatrixType;

    /// DoF types definition
    typedef Node<3> NodeType;
    typedef typename NodeType::DofType DofType;
    typedef typename DofType::Pointer DofPointerType;


    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /**
     * @brief Default constructor. (with parameters)
     */
    explicit ROMBuilderAndSolver(typename TLinearSolver::Pointer pNewLinearSystemSolver, Parameters ThisParameters)
    : BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(pNewLinearSystemSolver)
    {
        // Validate default parameters
        Parameters default_parameters = Parameters(R"(
        {
            "nodal_unknowns" : []
            "number_of_rom_dofs" : 10
        })" );

        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        // We set the other member variables
        mpLinearSystemSolver = pNewLinearSystemSolver;

        for(auto it=ThisParameters["nodal_dofs"].begin(); it!=ThisParameters["nodal_dofs"].end(); it++)
        {
            mNodalVariablesNames.push_back(it->GetString());
        }

        mnodal_dofs = mNodalVariablesNames.size();
        mrom_dofs = ThisParameters["number_of_rom_dofs"].GetInt();
    }



    /** Destructor.
     */
    virtual ~ROMBuilderAndSolver()
    {
    }


    /*@} */
    /**@name Operators
     */
    virtual void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part
    ) override
    {
        mDofList.clear();
        for(auto& node : r_model_part.Nodes())
        {
            for(const std::string& var_name : mNodalVariablesNames)
            {
                if( KratosComponents< Variable<double> >::Has( var_name ) ) //case of double variable
                {
                    auto pdof = node.pGetDof( KratosComponents< Variable<double> >::Get( var_name ) );
                    mDofList.push_back(pdof);
                }
                else if( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(var_name) ) //case of component variable
                {
                    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
                    component_type var_component = KratosComponents< component_type >::Get(var_name);
                    auto pdof = node.pGetDof( var_component.GetSourceVariable() );
                    mDofList.push_back( pdof );
                }
            }
        }
    }

    /**
            organises the dofset in order to speed up the building phase
     */
    virtual void SetUpSystem(
        ModelPart& r_model_part
    ) override
    {
        //deliberately not setting the equation id
    }

    Vector ProjectToReducedBasis(const TSystemVectorType& rX, ModelPart::NodesContainerType& rNodes)
    {
        Vector rom_unknowns = ZeroVector(mrom_dofs);
        for(auto& node : rNodes)
        {
            unsigned int node_aux_id = node.GetValue(AUX_ID);
            const Matrix& nodal_rom_basis = node.GetValue(ROM_BASIS);
            for(int i=0; i<mrom_dofs; ++i)
                for(int j=0; j<mnodal_dofs; ++j)
                    rom_unknowns[i] += nodal_rom_basis(i,j)*rX(node_aux_id*mnodal_dofs+j);
        }
        return rom_unknowns;
    }

    void ProjectToFineBasis(const TSystemVectorType& rRomUnkowns, ModelPart::NodesContainerType& rNodes,  TSystemVectorType& rX)
    {
        TSparseSpace::SetToZero(rX);
        for(auto& node : rNodes)
        {
            unsigned int node_aux_id = node.GetValue(AUX_ID);
            const Matrix& nodal_rom_basis = node.GetValue(ROM_BASIS);
            Vector tmp = prod(nodal_rom_basis, rRomUnkowns );

            for(unsigned int i=0; i<tmp.size(); ++i)
                rX[node_aux_id*mnodal_dofs+i] = tmp[i];
        }
    }

    void GetDofValues(const std::vector<DofPointerType>& rDofList, TSystemVectorType& rX)
    {
        unsigned int i=0;
        for(auto& dof : rDofList)
            rX[i++] = dof->GetSolutionStepValue();
    }

    /*@{ */

    /**
            Function to perform the building and solving phase at the same time.
            It is ideally the fastest and safer function to use when it is possible to solve
            just after building
     */
    virtual void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        //define a dense matrix to hold the reduced problem
        Matrix Arom = ZeroMatrix(mrom_dofs,mrom_dofs);
        Vector brom = ZeroVector(mrom_dofs);
        TSystemVectorType x(Dx.size());

        //find the rom basis
        this->GetDofValues(mDofList,x);
        Vector xrom = this->ProjectToReducedBasis(x, rModelPart.Nodes());
        
        //build the system matrix by looping over elements and conditions and ensembling to A
        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        // Getting the elements from the model
        const int nelements = static_cast<int>(rModelPart.Elements().size());

        // Getting the array of the conditions
        const int nconditions = static_cast<int>(rModelPart.Conditions().size());

        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();
        ModelPart::ConditionsContainerType::iterator cond_begin = rModelPart.ConditionsBegin();

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        double start_build = OpenMPUtils::GetCurrentTime();

//       #pragma omp parallel firstprivate(nelements,nconditions, LHS_Contribution, RHS_Contribution, EquationId )
        {
//            # pragma omp for  schedule(guided, 512) nowait
            for (int k = 0; k < nelements; k++)
            {
                ModelPart::ElementsContainerType::iterator it = el_begin + k;

                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool element_is_active = true;
                if ((it)->IsDefined(ACTIVE))
                    element_is_active = (it)->Is(ACTIVE);

                if (element_is_active)
                {
                    //calculate elemental contribution
                    pScheme->CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                    //assemble the elemental contribution - here is where the ROM acts
                    //compute the elemental reduction matrix T
                    const auto& geom = it->GetGeometry();
                    Matrix Telemental(geom.size()*mnodal_dofs, mrom_dofs);
                    for(unsigned int i=0; i<geom.size(); ++i)
                    {
                        const Matrix& rom_nodal_basis = geom[i].GetValue(ROM_BASIS);
                        for(unsigned int k=0; k<rom_nodal_basis.size1(); ++k)
                        {
                            row(Telemental, i*mnodal_dofs+k) = row(rom_nodal_basis,k);
                        }
                    }
                    
                    //compute LHS_ROM = Ttrans*LHS_Contribution*T
                    //and     RHS_ROM = Ttrans*b
                    //and sum such contributions to Arom and brom
                    Matrix aux = prod(LHS_Contribution, Telemental);
                    noalias(Arom) += prod(trans(Telemental), aux);
                    noalias(brom) += prod(trans(Telemental), RHS_Contribution);
 
                    // clean local elemental me overridemory
                    pScheme->CleanMemory(*(it.base()));
                } 
 
            } 
 
            // #pragma omp for  schedule(guided , 512)
            for (int k = 0; k < nconditions;  k++)
            {
                ModelPart::ConditionsContainerType::iterator it = cond_begin + k;

                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool condition_is_active = true;
                if ((it)->IsDefined(ACTIVE))
                    condition_is_active = (it)->Is(ACTIVE);

                if (condition_is_active)
                {
                    //calculate elemental contribution
                    pScheme->Condition_CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                    //assemble the elemental contribution - here is where the ROM acts
                    //compute the elemental reduction matrix T
                    const auto& geom = it->GetGeometry();
                    Matrix Telemental(geom.size()*mnodal_dofs, mrom_dofs);
                    for(unsigned int i=0; i<geom.size(); ++i)
                    {
                        const Matrix& rom_nodal_basis = geom[i].GetValue(ROM_BASIS);
                        for(unsigned int k=0; k<rom_nodal_basis.size1(); ++k)
                        {
                            row(Telemental, i*mnodal_dofs+k) = row(rom_nodal_basis,k);
                        }
                    }
                    
                    //compute LHS_ROM = Ttrans*LHS_Contribution*T
                    //and     RHS_ROM = Ttrans*b
                    //and sum such contributions to Arom and brom
                    Matrix aux = prod(LHS_Contribution, Telemental);
                    noalias(Arom) += prod(trans(Telemental), aux);
                    noalias(brom) += prod(trans(Telemental), RHS_Contribution);

                    // clean local elemental memory
                    pScheme->CleanMemory(*(it.base()));
                }
            }
        }

        const double stop_build = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Build time: " << stop_build - start_build << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Finished parallel building" << std::endl;

        //solve for the rom unkowns dunk = Arom^-1 * brom
        Vector dxrom(xrom.size());
        MathUtils<double>::Solve(Arom,dxrom,brom);

        //update database
        noalias(xrom) += dxrom;
        KRATOS_WATCH(xrom)

        //update fine Dx
        ProjectToFineBasis(dxrom, rModelPart.Nodes(),  Dx);
    }

    void ResizeAndInitializeVectors(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixPointerType& pA,
        TSystemVectorPointerType& pDx,
        TSystemVectorPointerType& pb,
        ModelPart& rModelPart
    ) override
    {
    }

    void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
    }

    void FinalizeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
    }

    virtual void CalculateReactions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
    }

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param r_model_part
     * @return 0 all ok
     */
    virtual int Check(ModelPart& r_model_part) override
    {
        KRATOS_TRY

        return 0;
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

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        return "ROMBuilderAndSolver";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
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

    /** Pointer to the Model.
     */
    typename TLinearSolver::Pointer mpLinearSystemSolver;

    DofsArrayType mDofSet;
    std::vector<DofPointerType> mDofList; 

    bool mReshapeMatrixFlag = false;

    /// flag taking care if the dof set was initialized ot not
    bool mDofSetIsInitialized = false;

    /// flag taking in account if it is needed or not to calculate the reactions
    bool mCalculateReactionsFlag = false;

    /// number of degrees of freedom of the problem to be solve
    unsigned int mEquationSystemSize;
    /*@} */
    /**@name Protected Operators*/
    /*@{ */

    int mEchoLevel = 0;

    TSystemVectorPointerType mpReactionsVector;

    std::vector< std::string > mNodalVariablesNames;
    int mnodal_dofs;
    int mrom_dofs;



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

    /*@} */
    /**@name Private Operators*/
    /*@{ */


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

}; /* Class ROMBuilderAndSolver */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_ROM_BUILDER_AND_SOLVER  defined */

