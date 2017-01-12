/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
/*
==============================================================================
KratosMultiScaleApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-06-06 10:37:00 $
//   Revision:            $Revision: 1.00 $
//
//


#if !defined(KRATOS_MUMPS_LINEAR_SOLVER_H_INCLUDED )
#define  KRATOS_MUMPS_LINEAR_SOLVER_H_INCLUDED

//#define USE_MUMPS_STATIC_LIB

// External includes

// Project includes
#include "utilities/timer.h"
#include "includes/define.h"
#include "linear_solvers/direct_solver.h"
#include <iostream>

extern "C" {
#ifdef USE_MUMPS_STATIC_LIB
#include "smumps_c.h"
#include "dmumps_c.h"
#else
#include "mumpsAggregate.h"
#endif // USE_MUMPS_STATIC_LIB
}

namespace Kratos
{

// ===============================================================================================
//
// MUMPS ENUMS
//
// ===============================================================================================

#define JOB_INIT			-1
#define JOB_END				-2
#define JOB_ANALYSIS		1
#define JOB_FACTORIZATION	2
#define JOB_SOLUTION		3
#define USE_COMM_WORLD	-987654

/*
mumps par%PAR (integer) must be initialized by the user on all processors and is accessed by MUMPS
only during the initialization phase (JOB = –1). It is not altered by MUMPS and its value is
communicated internally to the other phases as required.
Possible values for PAR are:
	0 host is not involved in factorization/solve phases
	1 host is involved in factorization/solve phases
*/
#define HOST_NOT_INVOLVED_IN_FACTORIZATION_AND_SOLVE 0
#define HOST_INVOLVED_IN_FACTORIZATION_AND_SOLVE 1

/*
mumps par%SYM (integer) must be initialized by the user on all processors and is accessed by
MUMPS only during the initialization phase (JOB = –1). It is not altered by MUMPS. Its value
is communicated internally to the other phases as required. Possible values for SYM are:

	0 A is unsymmetric
	
	1 A is suitable for symmetric positive definite since numerical pivoting is not performed and
	pivots are taken directly from the diagonal. In case ScaLAPACK is called, P POTRF is
	used, which assumes positive diagonal pivots (an error -40 is returned in INFOG(1)). In case
	ScaLAPACK is not used (ICNTL(13)>0), this option will also work for more general classes
	of matrices, typically symmetric negative matrices. If the user thinks his matrix is positive
	definite, he/she may want to check that the number of negative pivots (INFOG(12)) is zero
	on exit. Another approach to suppress numerical pivoting which works with ScaLAPACK
	for both positive definite and negative definite matrices consists in setting SYM=2 and
	CNTL(1)=0.0D0 (recommended strategy).
	
	2 A is general symmetric
	
Other values are treated as 0. Note that the value of SYM should be identical on all processors; if
this is not the case, the value on processor 0 is used by the package. For the complex version, the
value SYM=1 is currently treated as SYM=2. We do not have a version for Hermitian matrices in
this release of MUMPS.
*/
#define MATRIX_TYPE_UNSYMMETRIC						0
#define MATRIX_TYPE_SYMMETRIC_POSITIVE_DEFINITE		1
#define MATRIX_TYPE_GENERAL_SYMMETRIC				2


// ===============================================================================================
//
// MUMPS INTERFACE
//
// ===============================================================================================


template< typename TReal >
class MUMPSInterface 
{
};

template<>
class MUMPSInterface<float>
{
public:
	typedef float RealType;
	typedef int IndexType;
	typedef SMUMPS_STRUC_C StructureType;
	static void RunMUMPS(StructureType &s) 
	{ 
#ifdef USE_MUMPS_STATIC_LIB
		smumps_c(&s);
#else
		RunMUMPS_s(&s);
#endif // USE_MUMPS_STATIC_LIB
	}
};

template<>
class MUMPSInterface<double>
{
public:
	typedef double RealType;
	typedef int IndexType;
	typedef DMUMPS_STRUC_C StructureType;
	static void RunMUMPS(StructureType &s) 
	{
		//#pragma omp critical
		//{
#ifdef USE_MUMPS_STATIC_LIB
		dmumps_c(&s);
#else
		RunMUMPS_d(&s);
#endif // USE_MUMPS_STATIC_LIB
		//}
	}
};


// ===============================================================================================
//
// MUMPS SOLVER
//
// ===============================================================================================

#define MUMPS_V_2

template<class TSparseSpaceType, class TDenseSpaceType,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class MUMPSLinearSolver
    : public DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:

    KRATOS_CLASS_POINTER_DEFINITION( MUMPSLinearSolver );

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

	typedef typename TSparseSpaceType::DataType DataType;
	
	typedef typename SparseMatrixType::const_iterator1 SparseIteratorType1;

	typedef typename SparseMatrixType::const_iterator2 SparseIteratorType2;

	typedef MUMPSInterface< DataType > MUMPSInterfaceType;

	typedef typename MUMPSInterfaceType::StructureType MUMPSStructureType;

	typedef typename MUMPSInterfaceType::IndexType MUMPSIndexType;

	typedef typename MUMPSInterfaceType::RealType MUMPSRealType;

public:
	
    MUMPSLinearSolver()
		: m_nnz(0)
		, m_rows(NULL)
		, m_cols(NULL)
		, m_vals(NULL)
		, m_user_max_memory(false)
		, m_max_memory_MB(0)
	{
		m_done = false;
		m_mumps_lib_initialized = false;
	}

	MUMPSLinearSolver(int workspace_memory)
		: m_nnz(0)
		, m_rows(NULL)
		, m_cols(NULL)
		, m_vals(NULL)
		, m_user_max_memory(true)
		, m_max_memory_MB(workspace_memory)
	{

	}

    virtual ~MUMPSLinearSolver()
	{
		this->ClearParams();
	}

private:

	inline std::string GetError(int errorcode)
	{
		switch (errorcode) 
		{
			case -6: case -10:
				return "MUMPSSolver::Error( Matrix is singular )";
			case -13 :
				return "MUMPSSolver::Error( Not enough memory )";
			case -9:
			{
				return "MUMPSSolver::Error( increase ICNTL(14) )";
			}
			default:
			{
				std::stringstream msg;
				msg << "MUMPSSolver::Error( " << errorcode << " )";
				return msg.str();
			}
		}
	}

	inline void ClearParams()
	{
		m_nnz = 0;
		if (m_rows != NULL) delete[] m_rows;
		if (m_cols != NULL) delete[] m_cols;
		if (m_vals != NULL) delete[] m_vals;
		m_rows = NULL;
		m_cols = NULL;
		m_vals = NULL;

		m_params.n = 0;
		m_params.nz = 0;
		m_params.irn = NULL; 
		m_params.jcn = NULL; 
		m_params.a   = NULL; 
	}

	inline void LoadUnsymMatrix(const SparseMatrixType & rA)
	{
		size_t newprobsize = rA.size1();
		size_t newnnz = rA.nnz();

		if(newnnz != m_nnz) {
			this->ClearParams();
			m_nnz = newnnz;
			m_rows = new MUMPSIndexType[m_nnz];
			m_cols = new MUMPSIndexType[m_nnz];
			m_vals = new MUMPSRealType[m_nnz];
		}

		size_t counter(0);

		for( SparseIteratorType1 it1 = rA.begin1(); it1 != rA.end1(); it1++) {
			for( SparseIteratorType2 it2 = it1.begin(); it2 != it1.end(); it2++) {
				m_vals[counter]	= *it2;
				m_rows[counter]	=  it2.index1() + 1;
				m_cols[counter]	=  it2.index2() + 1;
				counter++;
			}
		}

		m_params.n = newprobsize;
		m_params.nz = newnnz;
		m_params.irn = m_rows;
		m_params.jcn = m_cols;
		m_params.a = m_vals;
	}
	
	inline bool InitializationPhase()
	{
		m_params.comm_fortran = USE_COMM_WORLD;
		m_params.par = HOST_INVOLVED_IN_FACTORIZATION_AND_SOLVE;
		m_params.job = JOB_INIT;

		m_params.sym = MATRIX_TYPE_UNSYMMETRIC; // TODO: check matrix type!
		m_params.irn = NULL;
		m_params.jcn = NULL;
		m_params.a   = NULL;

		MUMPSInterfaceType::RunMUMPS(m_params);

		if(m_params.info[0] != 0)
		{
			std::stringstream ss;
			ss << "MUMPSLinearSolver - Error during Initialization Phase:" << std::endl;
			ss << this->GetError(m_params.info[0]);
			m_last_error = ss.str();
			return false;
		}

		return true;
	}

	inline bool AnalysisPhase(const SparseMatrixType & rA)
	{
		this->LoadUnsymMatrix(rA);

		m_params.rhs = 0; // set later in Solve()!

		m_params.icntl[0] = -1; // output errors
		m_params.icntl[1] = -1; // output detailed messages
		m_params.icntl[2] = -1; // output info messages
		m_params.icntl[3] = -1; // verbose level
		m_params.icntl[5] =  7; // automatically find out if to scale or permute 

		m_params.icntl[6] = 7; // automatic pivot order for the factorization
		//m_params.icntl[6] = 5; // force the use of METIS

		m_params.icntl[7] = 7; // simultaneous row and colum iterative scaling (for not symmetric matrices)
		//m_params.icntl[7] = 1; // diagonal scaling (for symmetric matrices)

 		//m_params.icntl[0] = 6; // output errors
 		//m_params.icntl[1] = 6; // output detailed messages
 		//m_params.icntl[2] = 6; // output info messages
 		//m_params.icntl[3] = 0; // verbose level

		m_params.job = JOB_ANALYSIS;

		MUMPSInterfaceType::RunMUMPS(m_params);

		if(m_params.info[0] != 0) // Error Handling
		{
			std::stringstream ss;
			ss << "MUMPSLinearSolver - Error during Analysis Phase:" << std::endl;
			ss << "---Error Code: " << m_params.info[0] << std::endl;
			ss << "---Error Description: " << this->GetError(m_params.info[0]) << std::endl;
			ss << "---This error is not handled! see the manual for the description of this error code" << std::endl;
			m_last_error = ss.str();

			return false;
		}

		return true;
	}

	inline bool FactorizationPhase()
	{
		m_params.icntl[22] = 0; // default to 0. MUMPS chooses the max memory required to factorize
		if(m_user_max_memory)
		{
			int max_mem = m_max_memory_MB;
			int lower_bound = m_params.infog[15];
			if(max_mem < lower_bound)
			{
				max_mem = lower_bound;
				std::stringstream ss;
				ss << "MUMPS Error: the user-defined max memory is less than the minimum required for the factorization: " << lower_bound << " MB " << std::endl;
				std::cout << ss.str();
				return false;
			}
			m_params.icntl[22] = max_mem;
		}

		m_params.job = JOB_FACTORIZATION;

		MUMPSInterfaceType::RunMUMPS(m_params);

		if(m_params.info[0] != 0) // Error Handling
		{
			std::stringstream ss;
			ss << "MUMPSLinearSolver - Error during Factorization Phase:" << std::endl;
			ss << "---Error Code: " << m_params.info[0] << std::endl;
			ss << "---Error Description: " << this->GetError(m_params.info[0]) << std::endl;
			ss << "---This error is not handled! see the manual for the description of this error code" << std::endl;
			m_last_error = ss.str();
			return false;

			//if(m_params.info[0] == -9) // We have an handler for this error code
			//{
			//	// info[0] = -9  -> this means that we need to increase the size of m_params.icntl[13]
			//	// the required extra size can be found in info[1].
			//	// if info[1] is negative than the extra size is -info[1] * 1 million.
			//	int num_missing = m_params.info[1];
			//	if(num_missing < 0){
			//		num_missing = -num_missing * 1000000;
			//		std::stringstream ss;
			//		ss << "MUMPS WARNING: mum missing is " << num_missing << std::endl;
			//		m_last_error = ss.str();
			//		std::cout << ss.str();
			//		return false;
			//	}

			//	// try with this settings
			//	m_params.icntl[13] += num_missing;

			//	MUMPSInterfaceType::RunMUMPS(m_params);

			//	if(m_params.info[0] != 0)
			//	{
			//		if(m_params.info[0] == -9)
			//		{
			//			// This should never happen. But just to check...
			//			std::stringstream ss;
			//			ss << "MUMPSLinearSolver - Error during Factorization Phase:" << std::endl;
			//			ss << "---Error Code: " << m_params.info[0] << std::endl;
			//			ss << "---Error Description: " << this->GetError(m_params.info[0]) << std::endl;
			//			ss << "---It seems the the Error Handler for this error code Failed!" << std::endl;
			//			m_last_error = ss.str();
			//		}
			//		else
			//		{
			//			std::stringstream ss;
			//			ss << "MUMPSLinearSolver - Error during Factorization Phase:" << std::endl;
			//			ss << "---Error Code: " << m_params.info[0] << std::endl;
			//			ss << "---Error Description: " << this->GetError(m_params.info[0]) << std::endl;
			//			ss << "---This error is not handled! see the manual for the description of this error code" << std::endl;
			//			m_last_error = ss.str();
			//		}
			//		return false;
			//	}

			//	return true;
			//}
			//else // No Error Handler...
			//{
			//	std::stringstream ss;
			//	ss << "MUMPSLinearSolver - Error during Factorization Phase:" << std::endl;
			//	ss << "---Error Code: " << m_params.info[0] << std::endl;
			//	ss << "---Error Description: " << this->GetError(m_params.info[0]) << std::endl;
			//	ss << "---This error is not handled! see the manual for the description of this error code" << std::endl;
			//	m_last_error = ss.str();

			//	return false;
			//}
		}

		return true;
	}

	inline bool SolvePhase(const VectorType & rB, VectorType & rX)
	{
		rX = rB;

		m_params.job = JOB_SOLUTION;

		m_params.rhs = &rX.data()[0];
		MUMPSInterfaceType::RunMUMPS(m_params);
			
		if(m_params.info[0] != 0) // Error Handling
		{
			std::stringstream ss;
			ss << "MUMPSLinearSolver - Error during Solve Phase:" << std::endl;
			ss << "---Error Code: " << m_params.info[0] << std::endl;
			ss << "---Error Description: " << this->GetError(m_params.info[0]) << std::endl;
			ss << "---This error is not handled! see the manual for the description of this error code" << std::endl;
			m_last_error = ss.str();

			noalias( rX ) = ZeroVector(rB.size());
			m_params.rhs = NULL;
			return false;
		}

		m_params.rhs = NULL;

		return true;
	}

	inline bool FinalizationPhase()
	{
		m_params.job = JOB_END;
		MUMPSInterfaceType::RunMUMPS(m_params);

		if(m_params.info[0] != 0) // Error Handling
		{
			std::stringstream ss;
			ss << "MUMPSLinearSolver - Error during Finalization Phase:" << std::endl;
			ss << "---Error Code: " << m_params.info[0] << std::endl;
			ss << "---Error Description: " << this->GetError(m_params.info[0]) << std::endl;
			ss << "---This error is not handled! see the manual for the description of this error code" << std::endl;
			m_last_error = ss.str();

			return false;
		}

		return true;
	}

public:
	
	void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
    }

    void InitializeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
		m_done = true;

		if(m_mumps_lib_initialized) {
			m_done = this->FinalizationPhase();
			m_mumps_lib_initialized = false;
			if(!m_done) {
				this->ClearParams();
				return;
			}
		}

		m_params = MUMPSStructureType();
		this->ClearParams();

		if(!this->InitializationPhase()) {
			this->ClearParams();
			m_done = false;
			return;
		}

		if(!this->AnalysisPhase(rA)) {
			this->ClearParams();
			m_done = false;
			return;
		}

		if(!this->FactorizationPhase()) {
			this->ClearParams();
			m_done = false;
			return;
		}

		m_mumps_lib_initialized = true;
    }

    void PerformSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
		m_done = true;
		if(m_mumps_lib_initialized) {
			m_done = this->SolvePhase(rB, rX);
			if(!m_done) {
				this->ClearParams();
			}
		}
		else {
			m_done = false;
		}
    }

    void FinalizeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
		m_done = true;

		if(m_mumps_lib_initialized) {
			m_done = this->FinalizationPhase();
		}

		m_mumps_lib_initialized = false;
		this->ClearParams();
    }

    void Clear()
    {
		this->ClearParams();
    }

    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
		m_last_error = "";

        if(IsNotConsistent(rA, rX, rB)) {
			m_last_error = "The size of the system Matrices and Vectors are not consistent";
			return false;
		}

#ifdef MUMPS_V_2
		this->InitializeSolutionStep(rA, rX, rB);
		if(!m_done) return false;
		this->PerformSolutionStep(rA, rX, rB);
		if(!m_done) return false;
		this->FinalizeSolutionStep(rA, rX, rB);
		if(!m_done) return false;
		return true;
#else
		m_params = MUMPSStructureType();
		this->ClearParams();

		if(!this->InitializationPhase()) {
			this->ClearParams();
			return false;
		}

		if(!this->AnalysisPhase(rA)) {
			this->ClearParams();
			return false;
		}

		if(!this->FactorizationPhase()) {
			this->ClearParams();
			return false;
		}

		if(!this->SolvePhase(rB, rX)) {
			this->ClearParams();
			return false;
		}

		if(!this->FinalizationPhase()) {
			this->ClearParams();
			return false;
		}

		this->ClearParams();

		return true;
#endif // MUMPS_V_2

    }

    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
    {
		m_last_error = "";

        if(IsNotConsistent(rA, rX, rB)) {
			m_last_error = "The size of the system Matrices and Vectors are not consistent";
			return false;
		}

#ifdef MUMPS_V_2
		Vector vB(rB.size1());
		Vector vX(rX.size1());

		this->InitializeSolutionStep(rA, vX, vB);
		if(!m_done) return false;

		for(size_t j = 0; j < rB.size2(); j++) {
			noalias(vB) = column(rB, j);
			noalias(vX) = column(rX, j);
			this->PerformSolutionStep(rA, vX, vB);
			if(!m_done) return false;
			column(rX, j) = vX;
		}

		this->FinalizeSolutionStep(rA, vX, vB);
		if(!m_done) return false;
		return true;
#else
		if(!this->InitializationPhase()) {
			this->ClearParams();
			return false;
		}

		if(!this->AnalysisPhase(rA)) {
			this->ClearParams();
			return false;
		}

		if(!this->FactorizationPhase()) {
			this->ClearParams();
			return false;
		}

		Vector vB(rB.size1());
		Vector vX(rX.size1());
		for(size_t j = 0; j < rB.size2(); j++) {
			noalias(vB) = column(rB, j);
			noalias(vX) = column(rX, j);
			if(!this->SolvePhase(vB, vX)) {
				this->ClearParams();
				return false;
			}
			column(rX, j) = vX;
		}

		return true;
#endif // MUMPS_V_2

    }

    void  PrintInfo(std::ostream& rOStream) const
    {
        if(m_last_error.empty()) {
			std::stringstream ss;
			ss << "MUMPS Linear Solver finished without errors." << std::endl;
			rOStream << ss.str();
		}
		else {
			std::stringstream ss;
			ss << "MUMPS Linear Solver finished with ERROR:" << std::endl;
			ss << m_last_error << std::endl;
			rOStream << ss.str();
		}
    }

    void  PrintData(std::ostream& rOStream) const
    {
    }

private:

    MUMPSLinearSolver& operator=(const MUMPSLinearSolver& Other);
	
private:

	size_t			 m_nnz;
	MUMPSIndexType * m_rows;
	MUMPSIndexType * m_cols;
	MUMPSRealType  * m_vals;
	MUMPSStructureType m_params;
	std::string m_last_error;

	bool m_user_max_memory;
	int m_max_memory_MB;

	bool m_done;
	bool m_mumps_lib_initialized;

}; // Class MUMPSLinearSolver


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream, MUMPSLinearSolver<TSparseSpaceType,
                                  TDenseSpaceType,
                                  TReordererType>& rThis)
{
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MUMPSLinearSolver<TSparseSpaceType,
                                  TDenseSpaceType,
                                  TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


}  // namespace Kratos.

#endif // KRATOS_MUMPS_LINEAR_SOLVER_H_INCLUDED  defined 