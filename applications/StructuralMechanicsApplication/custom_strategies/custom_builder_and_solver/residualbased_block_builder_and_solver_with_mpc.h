//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala / Navaneeth K Narayanan
//
//

#ifndef KRATOS_SOLVING_STRATEGIES_BUILDER_AND_SOLVERS_RESIDUALBASED_BLOCK_BUILDER_AND_SOLVER_WITH_MPC_H_
#define KRATOS_SOLVING_STRATEGIES_BUILDER_AND_SOLVERS_RESIDUALBASED_BLOCK_BUILDER_AND_SOLVER_WITH_MPC_H_


/* System includes */

#include "utilities/openmp_utils.h"


#include <unordered_set>
/* External includes */
#include "boost/smart_ptr.hpp"

#include "utilities/timer.h"

/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "includes/model_part.h"

// #define USE_GOOGLE_HASH
#ifdef USE_GOOGLE_HASH
#include "sparsehash/dense_hash_set" //included in external libraries
#endif
// #define USE_LOCKS_IN_ASSEMBLY
// #include <iostream>
// #include <fstream>

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

Current class provides an implementation for standard builder and solving operations.

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
class ResidualBasedBlockBuilderAndSolverWithMpc
		: public ResidualBasedBlockBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
	/**@name Type Definitions */
	/*@{ */
	KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBlockBuilderAndSolverWithMpc);


	typedef ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

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
    typedef Dof<double>::Pointer DofPointerType;
    typedef Dof<double> DofType;


	/*@} */
	/**@name Life Cycle
	 */
	/*@{ */

	/** Constructor.
	 */
	ResidualBasedBlockBuilderAndSolverWithMpc(
			typename TLinearSolver::Pointer pNewLinearSystemSolver)
	: ResidualBasedBlockBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >(pNewLinearSystemSolver)
	  {
	  }

	/** Destructor.
	 */
	virtual ~ResidualBasedBlockBuilderAndSolverWithMpc()
	{
	}



    //**************************************************************************
    //**************************************************************************
	// This is modified to include the MPC information.
	//
    void Build(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& b) override
    {
        KRATOS_TRY
        if (!pScheme)
            KRATOS_THROW_ERROR(std::runtime_error, "No scheme provided!", "");

        //getting the elements from the model
        const int nelements = static_cast<int>(r_model_part.Elements().size());

        //getting the array of the conditions
        const int nconditions = static_cast<int>(r_model_part.Conditions().size());

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
        ModelPart::ElementsContainerType::iterator el_begin = r_model_part.ElementsBegin();
        ModelPart::ConditionsContainerType::iterator cond_begin = r_model_part.ConditionsBegin();

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        double start_build = OpenMPUtils::GetCurrentTime();

        #pragma omp parallel firstprivate(nelements,nconditions, LHS_Contribution, RHS_Contribution, EquationId )
        {
            # pragma omp for  schedule(guided, 512) nowait
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
                    // Modifying the local contributions for MPC
                    this->Element_ApplyMultipointConstraints(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                    //assemble the elemental contribution
    #ifdef USE_LOCKS_IN_ASSEMBLY
                    this->Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, BaseType::mlock_array);
    #else
                    this->Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
    #endif
                    // clean local elemental memory
                    pScheme->CleanMemory(*(it.base()));
                }

            }


            //#pragma omp parallel for firstprivate(nconditions, LHS_Contribution, RHS_Contribution, EquationId ) schedule(dynamic, 1024)
            #pragma omp for  schedule(guided, 512)
            for (int k = 0; k < nconditions; k++)
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

                    // Modifying the local contributions for MPC
                    this->Condition_ApplyMultipointConstraints(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                    //assemble the elemental contribution
    #ifdef USE_LOCKS_IN_ASSEMBLY
                    this->Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, BaseType::mlock_array);
    #else
                    this->Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
    #endif

                    // clean local elemental memory
                    pScheme->CleanMemory(*(it.base()));
                }
            }
        }

        // equation_ids.close();
        // KRATOS_THROW_ERROR(std::logic_error,"i want to stop here :-D","")

        double stop_build = OpenMPUtils::GetCurrentTime();
        if (this->GetEchoLevel() >= 1 && r_model_part.GetCommunicator().MyPID() == 0)
            std::cout << "build time: " << stop_build - start_build << std::endl;

        //for (int i = 0; i < A_size; i++)
        //    omp_destroy_lock(&lock_array[i]);
        if (this->GetEchoLevel() > 2 && r_model_part.GetCommunicator().MyPID() == 0)
        {
            KRATOS_WATCH("finished parallel building");
        }

        KRATOS_CATCH("")

    }



protected:


    void ConstructMatrixStructure(
        TSystemMatrixType& A,
        ElementsContainerType& rElements,
        ConditionsArrayType& rConditions,
        ProcessInfo& CurrentProcessInfo) override
    {
        //filling with zero the matrix (creating the structure)
        Timer::Start("MatrixStructure");

        const std::size_t equation_size = BaseType::mDofSet.size();

#ifdef USE_GOOGLE_HASH
        std::vector<google::dense_hash_set<std::size_t> > indices(equation_size);
        const std::size_t empty_key = 2*equation_size + 10;
#else
        std::vector<std::unordered_set<std::size_t> > indices(equation_size);
#endif

        #pragma omp parallel for firstprivate(equation_size)
        for (int iii = 0; iii < static_cast<int>(equation_size); iii++)
        {
#ifdef USE_GOOGLE_HASH
        indices[iii].set_empty_key(empty_key);
#else
        indices[iii].reserve(40);
#endif
        }
        Element::EquationIdVectorType ids(3, 0);

        const int nelements = static_cast<int>(rElements.size());
        #pragma omp parallel for firstprivate(nelements, ids)
        for(int iii=0; iii<nelements; iii++)
        {
            typename ElementsContainerType::iterator i_element = rElements.begin() + iii;
            (i_element)->EquationIdVector(ids, CurrentProcessInfo);

            // Modifying the equation IDs of this element to suit MPCs
            this->Element_ModifyEquationIdsForMPC(*(i_element.base()), ids, CurrentProcessInfo);

            for (std::size_t i = 0; i < ids.size(); i++)
            {
#ifdef _OPENMP
                omp_set_lock(&(BaseType::mlock_array[ids[i]]));
#endif
                auto& row_indices = indices[ids[i]];
                row_indices.insert(ids.begin(), ids.end());

#ifdef _OPENMP
                omp_unset_lock(&(BaseType::mlock_array[ids[i]]));
#endif
            }

        }
        const int nconditions = static_cast<int>(rConditions.size());
        #pragma omp parallel for firstprivate(nconditions, ids)
        for (int iii = 0; iii<nconditions; iii++)
        {
            typename ConditionsArrayType::iterator i_condition = rConditions.begin() + iii;
            (i_condition)->EquationIdVector(ids, CurrentProcessInfo);
            // Modifying the equation IDs of this element to suit MPCs
            this->Condition_ModifyEquationIdsForMPC(*(i_condition.base()), ids, CurrentProcessInfo);

            for (std::size_t i = 0; i < ids.size(); i++)
            {
#ifdef _OPENMP
                omp_set_lock(&(BaseType::mlock_array[ids[i]]));
#endif
                auto& row_indices = indices[ids[i]];
                row_indices.insert(ids.begin(), ids.end());
#ifdef _OPENMP
                omp_unset_lock(&(BaseType::mlock_array[ids[i]]));
#endif
            }
        }
        //count the row sizes
        unsigned int nnz = 0;
        for (unsigned int i = 0; i < indices.size(); i++)
            nnz += indices[i].size();

        A = boost::numeric::ublas::compressed_matrix<double>(indices.size(), indices.size(), nnz);

        double* Avalues = A.value_data().begin();
        std::size_t* Arow_indices = A.index1_data().begin();
        std::size_t* Acol_indices = A.index2_data().begin();

        //filling the index1 vector - DO NOT MAKE PARALLEL THE FOLLOWING LOOP!
        Arow_indices[0] = 0;
        for (int i = 0; i < static_cast<int>(A.size1()); i++)
            Arow_indices[i+1] = Arow_indices[i] + indices[i].size();


        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(A.size1()); i++)
        {
            const unsigned int row_begin = Arow_indices[i];
            const unsigned int row_end = Arow_indices[i+1];
            unsigned int k = row_begin;
            for (auto it = indices[i].begin(); it != indices[i].end(); it++)
            {
                Acol_indices[k] = *it;
                Avalues[k] = 0.0;
                k++;
            }

            indices[i].clear(); //deallocating the memory

            std::sort(&Acol_indices[row_begin], &Acol_indices[row_end]);

        }

        A.set_filled(indices.size()+1, nnz);
        Timer::Stop("MatrixStructure");
    }


protected:
    /*@} */
    /**@name Protected member Variables */
    /*@{ */

    /*
     * This function modifies the provided equation ID vector to accommodate MPC constraints
     */
    void Element_ModifyEquationIdsForMPC(Element::Pointer rCurrentElement,
    		Element::EquationIdVectorType& EquationId,
			ProcessInfo& CurrentProcessInfo){
    	const unsigned int number_of_nodes = rCurrentElement->GetGeometry().PointsNumber();

    	// For each node check if it is a slave or not If it is .. we change the Transformation matrix
    	for ( unsigned int j = 0; j < number_of_nodes; j++ )
    	{
    		if( rCurrentElement->GetGeometry()[j].GetValue(IS_SLAVE) ){ //temporary, will be checked once at the beginning only
    			//slaveFound = true;
    			// Necessary data for iterating and modifying the matrix
    			MpcData &mpcData       = rCurrentElement->GetGeometry()[j].GetValue(SLAVES);
    			std::vector<DofType> slaveDOFs = mpcData.GetSlaveDOFs();
    			int totalNumberOfSlaves  = slaveDOFs.size();

    			for(int slaveI=0; slaveI<totalNumberOfSlaves; slaveI++){ // For each of the Slave DOF
    				std::vector<DofType> masterDOFs = mpcData.GetMasterDOFs(slaveDOFs[slaveI]);

    				int totalNumberOfMasters = masterDOFs.size();
    				for(int i=0; i<totalNumberOfMasters; i++){
    					// We extend the EquationId vector
    					int masterEqID = masterDOFs[i].EquationId();
    					EquationId.push_back(masterEqID);
    				}
    			} // If the node is slave

    		}
    	}

    }

    /*
     * This function modifies the provided equation ID vector to accommodate MPC constraints
     */
    void Condition_ModifyEquationIdsForMPC(Condition::Pointer rCurrentCondition,
    		Element::EquationIdVectorType& EquationId,
			ProcessInfo& CurrentProcessInfo){

    	const unsigned int number_of_nodes = rCurrentCondition->GetGeometry().PointsNumber();

    	// For each node check if it is a slave or not If it is .. we change the Transformation matrix
    	for ( unsigned int j = 0; j < number_of_nodes; j++ )
    	{

    		if( rCurrentCondition->GetGeometry()[j].GetValue(IS_SLAVE) ){ //temporary, will be checked once at the beginning only
    			//slaveFound = true;
    			// Necessary data for iterating and modifying the matrix
    			MpcData &mpcData       = rCurrentCondition->GetGeometry()[j].GetValue(SLAVES);
    			std::vector<DofType> slaveDOFs = mpcData.GetSlaveDOFs();
    			int totalNumberOfSlaves  = slaveDOFs.size();

    			for(int slaveI=0; slaveI<totalNumberOfSlaves; slaveI++){ // For each of the Slave DOF
    				std::vector<DofType> masterDOFs = mpcData.GetMasterDOFs(slaveDOFs[slaveI]);
    				int totalNumberOfMasters = masterDOFs.size();
    				for(int i=0; i<totalNumberOfMasters; i++){
    					// We extend the EquationId vector
    					int masterEqID = masterDOFs[i].EquationId();
    					EquationId.push_back(masterEqID);
    				}
    			} // If the node is slave

    		}
    	}
    }

    /*
     * This function changes/extends the element LHS and RHS to apply MPC
     */

    void Element_ApplyMultipointConstraints(Element::Pointer rCurrentElement,
   		 LocalSystemMatrixType& LHS_Contribution,
			 LocalSystemVectorType& RHS_Contribution,
			 Element::EquationIdVectorType& EquationId,
			 ProcessInfo& CurrentProcessInfo){

    KRATOS_TRY
     
   	 bool slaveFound = false;
   	 const unsigned int number_of_nodes = rCurrentElement->GetGeometry().PointsNumber();
   	 int totalNumberOfSlaves = 0;
   	 int totalNumberOfMasters = 0;
   	 for ( unsigned int j = 0; j < number_of_nodes; j++ ){
   		 if( rCurrentElement->GetGeometry()[j].GetValue(IS_SLAVE) ){ //temporary, will be checked once at the beginning only
   			 slaveFound = true;
   			 break;
   		 }
   	 }

   	 // If no slave is found no need of going on
   	 if(!slaveFound){
   		 return;
   	 }

   	 // Loop Over all nodes
   	 for ( unsigned int j = 0; j < number_of_nodes; ++j ){
   		 // If the node has a slave DOF
   		 if( rCurrentElement->GetGeometry()[j].GetValue(IS_SLAVE) ){
   			 // Necessary data for iterating and modifying the matrix
   			 MpcData &mpcData       = rCurrentElement->GetGeometry()[j].GetValue(SLAVES);
             //mpcData.GetInfo();
   			 totalNumberOfMasters = mpcData.GetTotalNumbeOfMasterDOFs();
   			 totalNumberOfSlaves  = mpcData.GetSlaveDOFs().size();
   			 std::vector<DofType> slaveDOFs = mpcData.GetSlaveDOFs();
   			 std::vector<std::size_t>::iterator it;
   	
   			 // We resize the LHS and RHS contribution with the master sizes
   			 int currentSysSize = LHS_Contribution.size1();
   			 int lhsSize1 = currentSysSize + totalNumberOfMasters;
   			 int lhsSize2 = currentSysSize + totalNumberOfMasters;
   			 LHS_Contribution.resize(lhsSize1,lhsSize2,true); //true for Preserving the data and resizing the matrix
   			 RHS_Contribution.resize(lhsSize1,true);

   			 // Making the extra part of matrx
   			 for(int m=currentSysSize; m<lhsSize1; m++){
   				 for(int n=0; n<lhsSize1; n++){
   					 LHS_Contribution(m,n) = 0.0;
   					 LHS_Contribution(n,m) = 0.0;
   				 }
   				 RHS_Contribution(m) = 0.0;
   			 }

   			 int currentNumberOfMastersProcessed = 0;
   			 // Looping over all the slaves that current element has
   			 for(int slaveI=0; slaveI<totalNumberOfSlaves; ++slaveI){ // For each of the Slave DOF

   				 // Obtaining the local dof number for the slave.
   				 int localSlaveDOF = -1;
   				 int slaveEqId = slaveDOFs[slaveI].EquationId();
   				 it = std::find (EquationId.begin(), EquationId.end(), slaveEqId);
   				 if (it != EquationId.end()){
   					 std::size_t pos = std::distance(EquationId.begin(), it);
   					 localSlaveDOF = pos;
   				 }
   				 //std::cout<<"localSlaveDOF :: "<<localSlaveDOF<<" Global DOF :: "<<slaveDOFs[slaveI]<<std::endl;

   				 if(localSlaveDOF < 0){
   					 KRATOS_THROW_ERROR(std::logic_error,"","");
   				 }

   				 // Master data for current slave
   				 std::vector<DofType> masterDOFs = mpcData.GetMasterDOFs(slaveDOFs[slaveI]);
   				 int numMastersForThisSlave  = masterDOFs.size();
   				 std::vector<double> masterWeights = mpcData.GetMasterWeightsForSlave(slaveDOFs[slaveI]);

   				 // Looping over all the masters for this slave
   				 for(int masterI=0; masterI<numMastersForThisSlave; ++masterI){
                     int totalNumberOfDOFsIncurrentElement =  EquationId.size();
   					 int localMasterDOF = currentNumberOfMastersProcessed + currentSysSize;
                     ++currentNumberOfMastersProcessed;
   					 // Loop over all the actual DOFS
   					 for(int actualDofI=0; actualDofI<totalNumberOfDOFsIncurrentElement; actualDofI++){
   						 // For K(m,u) and K(u,m)
   						 if(actualDofI != localSlaveDOF){
   							 LHS_Contribution(localMasterDOF, actualDofI) += masterWeights[masterI] * LHS_Contribution(localSlaveDOF, actualDofI);
   							 LHS_Contribution(actualDofI, localMasterDOF) += masterWeights[masterI] * LHS_Contribution(actualDofI, localSlaveDOF);

   							 // Making K(s,u) and K(u,s) zero
   							 LHS_Contribution(localSlaveDOF, actualDofI) = 0.0;
   							 LHS_Contribution(actualDofI, localSlaveDOF) = 0.0;
   						 }

   						 // For K(m,s) and K(s,m)
   						 LHS_Contribution(localMasterDOF, localSlaveDOF) = masterWeights[masterI]    * LHS_Contribution(localSlaveDOF, localSlaveDOF);
   						 LHS_Contribution(localSlaveDOF, localMasterDOF) = masterWeights[masterI]    * LHS_Contribution(localSlaveDOF, localSlaveDOF);

   					 }// Loop over all the actual DOFS
   					 LHS_Contribution(localSlaveDOF, localSlaveDOF) = -1*LHS_Contribution(localSlaveDOF, localSlaveDOF);
   					 EquationId.push_back(masterDOFs[masterI].EquationId());
   					 // Changing the RHS side of the equation
   					 RHS_Contribution(localMasterDOF) =  RHS_Contribution(localMasterDOF) + masterWeights[masterI]*RHS_Contribution(localSlaveDOF);   					 

   				 }// Looping over all the masters for this slave

   				 RHS_Contribution(localSlaveDOF) = 0.0;
   			 }// Looping over all the slaves that current node ha

   		 } // If the node has slave DOF

   	 }// Loop Over all nodes
   		KRATOS_CATCH("Applying Multipoint constraints failed ..");
    } // End of function

    void Condition_ApplyMultipointConstraints( Condition::Pointer rCurrentElement,
                 LocalSystemMatrixType& LHS_Contribution,
                 LocalSystemVectorType& RHS_Contribution,
                 Element::EquationIdVectorType& EquationId,
                 ProcessInfo& CurrentProcessInfo){


    KRATOS_TRY
     bool slaveFound = false;
     const unsigned int number_of_nodes = rCurrentElement->GetGeometry().PointsNumber();
     int totalNumberOfSlaves = 0;
     int totalNumberOfMasters = 0;
     for ( unsigned int j = 0; j < number_of_nodes; j++ ){
         if( rCurrentElement->GetGeometry()[j].GetValue(IS_SLAVE) ){ //temporary, will be checked once at the beginning only
             slaveFound = true;
             break;
         }
     }

     // If no slave is found no need of going on
     if(!slaveFound){
         return;
     }

     // Loop Over all nodes
     for ( unsigned int j = 0; j < number_of_nodes; ++j ){
         // If the node has a slave DOF
         if( rCurrentElement->GetGeometry()[j].GetValue(IS_SLAVE) ){
             // Necessary data for iterating and modifying the matrix
             MpcData &mpcData       = rCurrentElement->GetGeometry()[j].GetValue(SLAVES);
             //mpcData.GetInfo();
             totalNumberOfMasters = mpcData.GetTotalNumbeOfMasterDOFs();
             totalNumberOfSlaves  = mpcData.GetSlaveDOFs().size();
             std::vector<DofType> slaveDOFs = mpcData.GetSlaveDOFs();
             std::vector<std::size_t>::iterator it;
    
             // We resize the LHS and RHS contribution with the master sizes
             int currentSysSize = LHS_Contribution.size1();
             int lhsSize1 = currentSysSize + totalNumberOfMasters;
             int lhsSize2 = currentSysSize + totalNumberOfMasters;
             LHS_Contribution.resize(lhsSize1,lhsSize2,true); //true for Preserving the data and resizing the matrix
             RHS_Contribution.resize(lhsSize1,true);

             // Making the extra part of matrx
             for(int m=currentSysSize; m<lhsSize1; m++){
                 for(int n=0; n<lhsSize1; n++){
                     LHS_Contribution(m,n) = 0.0;
                     LHS_Contribution(n,m) = 0.0;
                 }
                 RHS_Contribution(m) = 0.0;
             }

             int currentNumberOfMastersProcessed = 0;
             // Looping over all the slaves that current element has
             for(int slaveI=0; slaveI<totalNumberOfSlaves; ++slaveI){ // For each of the Slave DOF

                 // Obtaining the local dof number for the slave.
                 int localSlaveDOF = -1;
                 int slaveEqId = slaveDOFs[slaveI].EquationId();
                 it = std::find (EquationId.begin(), EquationId.end(), slaveEqId);
                 if (it != EquationId.end()){
                     std::size_t pos = std::distance(EquationId.begin(), it);
                     localSlaveDOF = pos;
                 }
                 //std::cout<<"localSlaveDOF :: "<<localSlaveDOF<<" Global DOF :: "<<slaveDOFs[slaveI]<<std::endl;

                 if(localSlaveDOF < 0){
                     KRATOS_THROW_ERROR(std::logic_error,"","");
                 }

                 // Master data for current slave
                 std::vector<DofType> masterDOFs = mpcData.GetMasterDOFs(slaveDOFs[slaveI]);
                 int numMastersForThisSlave  = masterDOFs.size();
                 std::vector<double> masterWeights = mpcData.GetMasterWeightsForSlave(slaveDOFs[slaveI]);

                 // Looping over all the masters for this slave
                 for(int masterI=0; masterI<numMastersForThisSlave; ++masterI){
                     int totalNumberOfDOFsIncurrentElement =  EquationId.size();
                     int localMasterDOF = currentNumberOfMastersProcessed + currentSysSize;
                     // Loop over all the actual DOFS
                     for(int actualDofI=0; actualDofI<totalNumberOfDOFsIncurrentElement; actualDofI++){
                         // For K(m,u) and K(u,m)
                         if(actualDofI != localSlaveDOF){
                             LHS_Contribution(localMasterDOF, actualDofI) += masterWeights[masterI] * LHS_Contribution(localSlaveDOF, actualDofI);
                             LHS_Contribution(actualDofI, localMasterDOF) += masterWeights[masterI] * LHS_Contribution(actualDofI, localSlaveDOF);

                             // Making K(s,u) and K(u,s) zero
                             LHS_Contribution(localSlaveDOF, actualDofI) = 0.0;
                             LHS_Contribution(actualDofI, localSlaveDOF) = 0.0;
                         }

                         // For K(m,s) and K(s,m)
                         LHS_Contribution(localMasterDOF, localSlaveDOF) = masterWeights[masterI]    * LHS_Contribution(localSlaveDOF, localSlaveDOF);
                         LHS_Contribution(localSlaveDOF, localMasterDOF) = masterWeights[masterI]    * LHS_Contribution(localSlaveDOF, localSlaveDOF);

                     }// Loop over all the actual DOFS
                     LHS_Contribution(localSlaveDOF, localSlaveDOF) = -1*LHS_Contribution(localSlaveDOF, localSlaveDOF);
                     EquationId.push_back(masterDOFs[masterI].EquationId());
                     // Changing the RHS side of the equation
                     RHS_Contribution(localMasterDOF) =  RHS_Contribution(localMasterDOF) + masterWeights[masterI]*RHS_Contribution(localSlaveDOF);
                     ++currentNumberOfMastersProcessed;

                 }// Looping over all the masters for this slave

                 RHS_Contribution(localSlaveDOF) = 0.0;
             }// Looping over all the slaves that current node ha

         } // If the node has slave DOF

     }// Loop Over all nodes
         //std::exit(-1);
        if(0){
         std::cout<<"After modification for MPC :: "<<std::endl;
         std::cout<<"Element ID :: "<<rCurrentElement->Id()<<std::endl;
            /*std::cout<<"LHS :: "<<std::endl;
            std::cout<<LHS_Contribution<<std::endl;
            std::cout<<"RHS :: "<<std::endl;
            std::cout<<RHS_Contribution<<std::endl;*/
            std::cout<<"EquationId :: "<<std::endl;
            for( const auto& it : EquationId )
            {
             std::cout<<it<<std::endl;
            }
        }
        //std::cout<<rCurrentElement->Id()<<" :: Finished appying element MPC .. ! "<<std::endl;
        KRATOS_CATCH("Applying Multipoint constraints failed ..");
    } // End of the function



};

}




#endif /* KRATOS_SOLVING_STRATEGIES_BUILDER_AND_SOLVERS_RESIDUALBASED_BLOCK_BUILDER_AND_SOLVER_WITH_MPC_H_ */
