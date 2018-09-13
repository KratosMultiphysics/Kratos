//   
//   Project Name:        KratosDamApplication   $
//   Last Modified by:    $Author:Ignasi de Pouplana $
//   Date:                $Date:    February 2017$
//   Revision:            $Revision:         1.0 $
//

#if !defined(KRATOS_TRILINOS_DAM_UP_SCHEME )
#define  KRATOS_TRILINOS_DAM_UP_SCHEME

/* External includes */
#include "Epetra_Import.h"

// Application includes
#include "custom_strategies/schemes/dam_UP_scheme.hpp"
#include "dam_application_variables.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>

class TrilinosDamUPScheme : public DamUPScheme<TSparseSpace,TDenseSpace>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( TrilinosDamUPScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>     BaseType;
    typedef typename BaseType::DofsArrayType     DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    TrilinosDamUPScheme(double beta, double gamma, double rayleigh_m ,double rayleigh_k)
        : DamUPScheme<TSparseSpace,TDenseSpace>(beta, gamma, rayleigh_m ,rayleigh_k),
        mImporterIsInitialized(false) {}
    
    //------------------------------------------------------------------------------------
    
    ///Destructor
    virtual ~TrilinosDamUPScheme() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Update(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        KRATOS_TRY
        
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
        
        this->UpdateVariablesDerivatives(r_model_part);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Clear()
    {
        BaseType::Clear();
        
        mpDofImporter.reset();
        mImporterIsInitialized = false;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    bool DofImporterIsInitialized()
    {
        return mImporterIsInitialized;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:
    
    /// Member Variables
    

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
        Kratos::shared_ptr<Epetra_Import> pDofImporter = Kratos::make_shared<Epetra_Import>(dof_update_map,Dx.Map());
        mpDofImporter.swap(pDofImporter);

        mImporterIsInitialized = true;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Kratos::shared_ptr<Epetra_Import> pGetImporter()
    {
        return mpDofImporter;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    /// Member Variables
    
    bool mImporterIsInitialized;
    Kratos::shared_ptr<Epetra_Import> mpDofImporter;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


}; // Class TrilinosDamUPScheme
}  // namespace Kratos

#endif // KRATOS_TRILINOS_DAM_UP_SCHEME defined