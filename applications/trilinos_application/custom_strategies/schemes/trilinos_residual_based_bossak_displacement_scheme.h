//  KRATOS  _____     _ _ _                 
//         |_   _| __(_) (_)_ __   ___  ___ 
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__ 
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:             BSD License 
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//      

#if !defined(KRATOS_TRILINOS_RESIDUAL_BASED_BOSSAK_DISPLACEMENT_SCHEME )
#define  KRATOS_TRILINOS_RESIDUAL_BASED_BOSSAK_DISPLACEMENT_SCHEME

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"
#include "Epetra_Import.h"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp"
#include "includes/variables.h"

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

/** @brief TRrilinos Bossak integration scheme (for dynamic problems)
 */
template<class TSparseSpace,  class TDenseSpace >
class TrilinosResidualBasedBossakDisplacementScheme: public ResidualBasedBossakDisplacementScheme< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( TrilinosResidualBasedBossakDisplacementScheme );

    typedef ResidualBasedBossakDisplacementScheme<TSparseSpace,TDenseSpace> BaseType;

    typedef typename BaseType::TDataType                                   TDataType;

    typedef typename BaseType::DofsArrayType                           DofsArrayType;

    typedef typename Element::DofsVectorType                          DofsVectorType;

    typedef typename BaseType::TSystemMatrixType                   TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType                   TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType           LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType           LocalSystemMatrixType;

    typedef ModelPart::NodesContainerType                             NodesArrayType;

    typedef ModelPart::ElementsContainerType                       ElementsArrayType;

    typedef ModelPart::ConditionsContainerType                   ConditionsArrayType;

    typedef typename BaseType::Pointer                               BaseTypePointer;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     * The bossak method
     */
    TrilinosResidualBasedBossakDisplacementScheme(double rAlpham = 0.0)
        :ResidualBasedBossakDisplacementScheme<TSparseSpace,TDenseSpace>(rAlpham),
        mImporterIsInitialized(false)
    {}

    /** Destructor.
     */
    virtual ~TrilinosResidualBasedBossakDisplacementScheme
    () {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * Performing the update of the solution
     * Incremental update within newton iteration. It updates the state variables at the end of the time step: u_{n+1}^{k+1}= u_{n+1}^{k}+ \Delta u
     * @param rModelPart: The model of the problem to solve
     * @param rDofSet: Set of all primary variables
     * @param A: LHS matrix
     * @param Dx: incremental update of primary variables
     * @param b: RHS Vector
     */

    void Update(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b 
        )
    {
        KRATOS_TRY;

        // std::cout << " Update " << std::endl;

        if (DofImporterIsInitialized() == false)
        {
            this->InitializeDofImporter(rDofSet,Dx);
        }

        const int system_size = TSparseSpace::Size1(A);
        
        const unsigned int NumThreads = OpenMPUtils::GetNumThreads();

        // Defining a temporary vector to gather all of the values needed
        Epetra_Vector temp( mpDofImporter->TargetMap() );
        
        // Importing in the new temp vector the values
        const unsigned int ierr = temp.Import(Dx,*mpDofImporter,Insert);
        if(ierr != 0) 
        {
            KRATOS_THROW_ERROR(std::logic_error,"Epetra failure found","");
        }
        
        double* temp_values; //DO NOT make delete of this one!!
        temp.ExtractView( &temp_values );

        Dx.Comm().Barrier();
        
        // Update of displacement (by DOF)
        OpenMPUtils::PartitionVector DofPartition;
        OpenMPUtils::DivideInPartitions(rDofSet.size(), NumThreads, DofPartition);

        const int ndof = static_cast<int>(rDofSet.size());
        typename DofsArrayType::iterator DofBegin = rDofSet.begin();

        #pragma omp parallel for firstprivate(DofBegin)
        for(int i = 0;  i < ndof; i++)
        {
            typename DofsArrayType::iterator itDof = DofBegin + i;
            int global_id = itDof->EquationId();

            if(global_id < system_size)
            {
                if (itDof->IsFree() )
                {
                    const double auxDx = temp[mpDofImporter->TargetMap().LID(global_id)];
                    itDof->GetSolutionStepValue() += auxDx;
                }
            }
        }

        // Updating time derivatives (nodally for efficiency)
        OpenMPUtils::PartitionVector NodePartition;
        OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);

        const int nnodes = static_cast<int>(rModelPart.Nodes().size());
        NodesArrayType::iterator NodeBegin = rModelPart.Nodes().begin();

        #pragma omp parallel for firstprivate(NodeBegin)
        for(int i = 0;  i < nnodes; i++)
        {
            array_1d<double, 3 > DeltaDisplacement;

            NodesArrayType::iterator itNode = NodeBegin + i;

            noalias(DeltaDisplacement) = (itNode)->FastGetSolutionStepValue(DISPLACEMENT) - (itNode)->FastGetSolutionStepValue(DISPLACEMENT, 1);

            array_1d<double, 3 > & CurrentVelocity            = (itNode)->FastGetSolutionStepValue(VELOCITY, 0);
            const array_1d<double, 3 > & PreviousVelocity     = (itNode)->FastGetSolutionStepValue(VELOCITY, 1);

            array_1d<double, 3 > & CurrentAcceleration        = (itNode)->FastGetSolutionStepValue(ACCELERATION, 0);
            const array_1d<double, 3 > & PreviousAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION, 1);

            BaseType::UpdateVelocity     (CurrentVelocity,     DeltaDisplacement, PreviousVelocity, PreviousAcceleration);

            BaseType::UpdateAcceleration (CurrentAcceleration, DeltaDisplacement, PreviousVelocity, PreviousAcceleration);
        }

        KRATOS_CATCH( "" );
    }

    ///@}
    ///@name Operations
    ///@{
    
    void Clear()
    {
        BaseType::Clear();
        
        mpDofImporter.reset();
        mImporterIsInitialized = false;
    }
    
    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    bool DofImporterIsInitialized()
    {
        return mImporterIsInitialized;
    }
    
    ///@}
    ///@name Input and output
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

    /**
     * This initilizes the DofImporter
     */
    virtual void InitializeDofImporter(
        DofsArrayType& rDofSet,
        TSystemVectorType& Dx
        )
    {
        int system_size = TSparseSpace::Size(Dx);
        int number_of_dofs = rDofSet.size();
        std::vector< int > index_array(number_of_dofs);

        // Filling the array with the global ids
        unsigned int counter = 0;
        for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
            int id = i_dof->EquationId();
            if( id < system_size )
            {
                index_array[counter] = id;
                counter += 1;
            }
        }

        std::sort(index_array.begin(),index_array.end());
        std::vector<int>::iterator NewEnd = std::unique(index_array.begin(),index_array.end());
        index_array.resize(NewEnd-index_array.begin());

        int check_size = -1;
        int tot_update_dofs = index_array.size();
        Dx.Comm().SumAll(&tot_update_dofs,&check_size,1);
        if ( (check_size < system_size) &&  (Dx.Comm().MyPID() == 0) )
        {
            std::stringstream Msg;
            Msg << "Dof count is not correct. There are less dofs then expected." << std::endl;
            Msg << "Expected number of active dofs = " << system_size << " dofs found = " << check_size << std::endl;
            KRATOS_THROW_ERROR(std::runtime_error,Msg.str(),"")
        }

        // Defining a map as needed
        Epetra_Map dof_update_map(-1,index_array.size(), &(*(index_array.begin())),0,Dx.Comm() );

        // Defining the importer class
        boost::shared_ptr<Epetra_Import> pDofImporter( new Epetra_Import(dof_update_map,Dx.Map()) );
        mpDofImporter.swap(pDofImporter);

        mImporterIsInitialized = true;
    }
    
    ///@}
    ///@name Protected  Access
    ///@{

    /**
     * Get pointer Epetra_Import instance that can be used to import values from Dx to the owner
     * of each Dof.
     * @note Important: always check that the Importer is initialized before calling using
     * DofImporterIsInitialized or initialize it with InitializeDofImporter.
     * @return Importer
     */
    
    boost::shared_ptr<Epetra_Import> pGetImporter()
    {
        return mpDofImporter;
    }
    
    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@{

private:

    ///@name Static Member Variables
    ///@{
    
    ///@}
    ///@name Member Variables
    ///@{

    bool mImporterIsInitialized;

    boost::shared_ptr<Epetra_Import> mpDofImporter;
    
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
    ///@name Un accessible methods
    ///@{
    ///@}
}; /* Class TrilinosResidualBasedBossakDisplacementScheme */
///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}
}  /* namespace Kratos.*/

#endif /* KRATOS_TRILINOS_RESIDUAL_BASED_BOSSAK_DISPLACEMENT_SCHEME defined */
