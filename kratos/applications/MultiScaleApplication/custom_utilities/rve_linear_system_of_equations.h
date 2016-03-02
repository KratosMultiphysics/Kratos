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
//   Date:                $Date: 2013-11-04 12:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#if !defined(RVE_LINEAR_SYSTEM_OF_EQUATIONS_H_INCLUDED)
#define RVE_LINEAR_SYSTEM_OF_EQUATIONS_H_INCLUDED

#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/define.h"
#include "includes/serializer.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/schemes/scheme.h"
#include "rve_constraint_handler.h"
#include "rve_geometry_descriptor.h"
#include "rve_utilities.h"


//#define RVE_SOE_TIMER_ON

namespace Kratos
{


	namespace RveSOEUtils
	{
		inline void AddUnique(std::vector<std::size_t>& v, const std::size_t& candidate)
		{
			std::vector<std::size_t>::iterator i = v.begin();
			std::vector<std::size_t>::iterator endit = v.end();
			while (i != endit && (*i) != candidate)
				i++;
			if (i == endit)
				v.push_back(candidate);
		}
	}




	template<class TSparseSpace, 
			 class TDenseSpace,
			 class TReorderer = Reorderer<TSparseSpace, TDenseSpace> >
	class RveLinearSystemOfEquations
	{

	public:

		KRATOS_CLASS_POINTER_DEFINITION( RveLinearSystemOfEquations );
		typedef LinearSolver<TSparseSpace, TDenseSpace, TReorderer> LinearSolverType;
		typedef typename LinearSolverType::Pointer LinearSolverPointerType;
		typedef typename TSparseSpace::MatrixType SparseMatrixType;
		typedef typename TSparseSpace::VectorType VectorType;
		typedef typename TDenseSpace::MatrixType  DenseMatrixType;
		typedef Scheme<TSparseSpace,TDenseSpace> SchemeType;
		typedef typename SchemeType::Pointer SchemePointerType;
		typedef ModelPart::DofsArrayType DofsArrayType;
		typedef ModelPart::NodesContainerType NodesArrayType;
		typedef ModelPart::ElementsContainerType ElementsArrayType;
		typedef ModelPart::ConditionsContainerType ConditionsArrayType;
		typedef RveConstraintHandler<TSparseSpace, TDenseSpace> RveConstraintHandlerType;
		typedef typename RveConstraintHandlerType::Pointer RveConstraintHandlerPointerType;
		typedef typename RveConstraintHandlerType::IndexContainerType IndexContainerType;

	public:
		
		RveLinearSystemOfEquations(LinearSolverPointerType pLinearSolver)
			: m_lin_solver(pLinearSolver)
			, m_calculate_reactions(true)
			, m_initialized(false)
			, m_equation_system_size(0)
		{
		}

		virtual ~RveLinearSystemOfEquations()
		{
		}

	public:

		virtual std::string GetInfo()const
		{
			std::stringstream ss;
			m_lin_solver->PrintInfo(ss);
			return ss.str();
		}

	public:

		virtual void Begin(ModelPart& mp, 
						   const RveGeometryDescriptor& geom,
						   SchemePointerType& pScheme, 
						   RveConstraintHandlerPointerType& chandler)
		{
			//if(m_initialized)return; //COMMENTED BY STEFANO
			if(m_initialized)
				this->End();

			chandler->SetupDofSet(mp, geom, pScheme, m_dofset);
			chandler->SetUpSystem(mp, geom, m_dofset, m_equation_system_size, m_transformed_equation_ids, m_equation_id_flag);
			this->ResizeAndInitializeVectors(mp);

			m_initialized = true;
		}

		virtual void BuildRHS(ModelPart& mp, 
							  SchemePointerType& pScheme)
		{
			ElementsArrayType& pElements = mp.Elements();
			ConditionsArrayType& ConditionsArray = mp.Conditions();

			TSparseSpace::SetToZero(m_b);
			if(m_calculate_reactions)
				TSparseSpace::SetToZero(m_r);

			Element::EquationIdVectorType EquationId;
			DenseMatrixType LHS_Contribution = DenseMatrixType(0, 0);
			VectorType RHS_Contribution = VectorType(0);

			ProcessInfo& CurrentProcessInfo = mp.GetProcessInfo();

			for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
			{
				pScheme->Calculate_RHS_Contribution(*it, RHS_Contribution, EquationId, CurrentProcessInfo);
				AssembleRHS(RHS_Contribution, EquationId);
			}

			LHS_Contribution.resize(0, 0, false);
			RHS_Contribution.resize(0, false);

			for (typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it)
			{
				pScheme->Condition_Calculate_RHS_Contribution(*it, RHS_Contribution, EquationId, CurrentProcessInfo);
				AssembleRHS(RHS_Contribution, EquationId);
			}

			if(m_calculate_reactions)
				this->CalculateReactions(mp, pScheme);
		}

		virtual void BuildRHS_Reduced(ModelPart& mp, 
							  SchemePointerType& pScheme,
							  const RveGeometryDescriptor::IndexContainerType& elem_subset)
		{
			ElementsArrayType& pElements = mp.Elements();

			TSparseSpace::SetToZero(m_b);
			if(m_calculate_reactions)
				TSparseSpace::SetToZero(m_r);

			Element::EquationIdVectorType EquationId;
			DenseMatrixType LHS_Contribution = DenseMatrixType(0, 0);
			VectorType RHS_Contribution = VectorType(0);

			ProcessInfo& CurrentProcessInfo = mp.GetProcessInfo();

			for(RveGeometryDescriptor::IndexContainerType::const_iterator it = elem_subset.begin(); it != elem_subset.end(); ++it)
			{
				RveGeometryDescriptor::IndexType elemid = *it;
				ModelPart::ElementType::Pointer ielem = mp.pGetElement(elemid);
				pScheme->Calculate_RHS_Contribution(ielem, RHS_Contribution, EquationId, CurrentProcessInfo);
				AssembleRHS(RHS_Contribution, EquationId);
			}

			if(m_calculate_reactions)
				this->CalculateReactions(mp, pScheme);
		}

		virtual void BuildLHS(ModelPart& mp, 
							  SchemePointerType& pScheme)
		{
			ElementsArrayType& pElements = mp.Elements();
			ConditionsArrayType& ConditionsArray = mp.Conditions();

			TSparseSpace::SetToZero(m_A);

			DenseMatrixType LHS_Contribution = DenseMatrixType(0, 0);
			Element::EquationIdVectorType EquationId;

			ProcessInfo& CurrentProcessInfo = mp.GetProcessInfo();

			for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
			{
				pScheme->Calculate_LHS_Contribution(*it, LHS_Contribution, EquationId, CurrentProcessInfo);
				AssembleLHS(LHS_Contribution, EquationId);
				pScheme->CleanMemory(*it);
			}

			LHS_Contribution.resize(0, 0, false);

			for (typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it)
			{
				pScheme->Condition_Calculate_LHS_Contribution(*it, LHS_Contribution, EquationId, CurrentProcessInfo);
				AssembleLHS(LHS_Contribution, EquationId);
			}

			// factorize
			this->m_lin_solver->InitializeSolutionStep(m_A, m_x, m_b);
		}

		virtual void Build(ModelPart& mp, 
						   SchemePointerType& pScheme)
		{

#ifdef RVE_SOE_TIMER_ON
			RveUtilities::RveTimer timer;
#endif // RVE_SOE_TIMER_ONE

			ElementsArrayType& pElements = mp.Elements();
			ConditionsArrayType& ConditionsArray = mp.Conditions();

			TSparseSpace::SetToZero(m_A);
			TSparseSpace::SetToZero(m_b);
			if(m_calculate_reactions)
				TSparseSpace::SetToZero(m_r);

			Element::EquationIdVectorType EquationId;
			DenseMatrixType LHS_Contribution = DenseMatrixType(0, 0);
			VectorType RHS_Contribution = VectorType(0);

			ProcessInfo& CurrentProcessInfo = mp.GetProcessInfo();

#ifdef RVE_SOE_TIMER_ON
			double time_calc = 0.0;
			double time_assm = 0.0;
#endif // RVE_SOE_TIMER_ONE

			for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
			{
#ifdef RVE_SOE_TIMER_ON
				timer.start();
#endif // RVE_SOE_TIMER_ONE
				pScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
#ifdef RVE_SOE_TIMER_ON
				timer.stop();
				time_calc += timer.value();
				timer.start();
#endif // RVE_SOE_TIMER_ONE
				AssembleLHS(LHS_Contribution, EquationId);
				AssembleRHS(RHS_Contribution, EquationId);
#ifdef RVE_SOE_TIMER_ON
				timer.stop();
				time_assm += timer.value();
#endif // RVE_SOE_TIMER_ONE
				pScheme->CleanMemory(*it);
			}

			LHS_Contribution.resize(0, 0, false);
			RHS_Contribution.resize(0, false);

			for (typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it)
			{
				pScheme->Condition_CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
				AssembleLHS(LHS_Contribution, EquationId);
				AssembleRHS(RHS_Contribution, EquationId);
			}

			if(m_calculate_reactions)
				this->CalculateReactions(mp, pScheme);

			// factorize
#ifdef RVE_SOE_TIMER_ON
			timer.start();
#endif // RVE_SOE_TIMER_ONE
			this->m_lin_solver->InitializeSolutionStep(m_A, m_x, m_b);
#ifdef RVE_SOE_TIMER_ON
			timer.stop();
			double time_fact = timer.value();
			double total_time = time_calc + time_assm + time_fact;
			std::stringstream ss;
			ss << "RVE SOE Timings:\n";
			ss << "Calculate:   " << time_calc << "; % = " << time_calc/total_time << std::endl;
			ss << "Assemble:    " << time_assm << "; % = " << time_assm/total_time << std::endl;
			ss << "Factorize:   " << time_fact << "; % = " << time_fact/total_time << std::endl;
			std::cout << ss.str();
#endif // RVE_SOE_TIMER_ONE
		}

		virtual void Solve()
		{
			//TSparseSpace::SetToZero(m_x);
			if(TSparseSpace::Size(m_x) == 0) return;
			//this->m_lin_solver->PerformSolutionStep(m_A,m_x,m_b); //COMMENTED BY STEFANO
			this->m_lin_solver->Solve(m_A,m_x,m_b); // TODO: dire a ricc se si puo splittare initi/solve/fin
		}

		virtual void End()
		{
			//return; //COMMENTED BY STEFANO
			if(m_initialized)
			{
				this->m_lin_solver->FinalizeSolutionStep(m_A,m_x,m_b);

				m_b = VectorType();
				m_x = VectorType();
				m_r = VectorType();
				m_A = SparseMatrixType();
				
				m_dofset = DofsArrayType();
				m_equation_system_size = 0;

				m_initialized = false;
			}
		}

	protected:

		virtual void ResizeAndInitializeVectors(ModelPart& mp)
		{
			m_A.resize(m_equation_system_size, m_equation_system_size, false);
			ConstructMatrixStructure(mp);

			if (m_x.size() != m_equation_system_size)
				m_x.resize(m_equation_system_size, false);
			if (m_b.size() != m_equation_system_size)
				m_b.resize(m_equation_system_size, false);

			if (m_calculate_reactions)
			{
				size_t rection_vector_size = m_dofset.size() - m_equation_system_size;
				if (m_r.size() != rection_vector_size)
					m_r.resize(rection_vector_size, false);
			}
		}

		virtual void ConstructMatrixStructure(ModelPart& mp)
		{
			ElementsArrayType& rElements = mp.Elements();
			ConditionsArrayType& rConditions = mp.Conditions();
			ProcessInfo& CurrentProcessInfo = mp.GetProcessInfo();

			std::size_t equation_size = m_equation_system_size;
			std::vector<std::vector<std::size_t> > indices(equation_size);

			Element::EquationIdVectorType ids(3, 0);
			for (typename ElementsArrayType::iterator i_element = rElements.begin(); i_element != rElements.end(); i_element++)
			{
				(i_element)->EquationIdVector(ids, CurrentProcessInfo);

				for (std::size_t i = 0; i < ids.size(); i++)
					if (ids[i] < equation_size)
					{
						std::vector<std::size_t>& row_indices = indices[ids[i]];
						for (std::size_t j = 0; j < ids.size(); j++)
							if (ids[j] < equation_size)
							{
								RveSOEUtils::AddUnique(row_indices, ids[j]);
							}
					}
			}

			for (typename ConditionsArrayType::iterator i_condition = rConditions.begin(); i_condition != rConditions.end(); i_condition++)
			{
				(i_condition)->EquationIdVector(ids, CurrentProcessInfo);
				for (std::size_t i = 0; i < ids.size(); i++)
					if (ids[i] < equation_size)
					{
						std::vector<std::size_t>& row_indices = indices[ids[i]];
						for (std::size_t j = 0; j < ids.size(); j++)
							if (ids[j] < equation_size)
							{
								RveSOEUtils::AddUnique(row_indices, ids[j]);
							}
					}
			}

			//allocating the memory needed
			int data_size = 0;
			for (std::size_t i = 0; i < indices.size(); i++)
			{
				data_size += indices[i].size();
			}
			m_A.reserve(data_size, false);

			//filling with zero the matrix (creating the structure)
			for (std::size_t i = 0; i < indices.size(); i++)
			{
				std::vector<std::size_t>& row_indices = indices[i];
				std::sort(row_indices.begin(), row_indices.end());

				for (std::vector<std::size_t>::iterator it = row_indices.begin(); it != row_indices.end(); it++)
				{
					m_A.push_back(i, *it, 0.00);
				}
				row_indices.clear();
			}
		}

		virtual void AssembleLHS(DenseMatrixType& LHS_Contribution,
								 Element::EquationIdVectorType& EquationId)
		{
			unsigned int local_size = LHS_Contribution.size1();

			for (unsigned int i_local = 0; i_local < local_size; i_local++)
			{
				unsigned int i_global = EquationId[i_local];
				// MAZ_01 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				i_global = m_transformed_equation_ids[i_global];
				// MAZ_01 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				if (i_global < m_equation_system_size)
				{
					for (unsigned int j_local = 0; j_local < local_size; j_local++)
					{
						unsigned int j_global = EquationId[j_local];
						// MAZ_01 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
						j_global = m_transformed_equation_ids[j_global];
						// MAZ_01 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
						if (j_global < m_equation_system_size)
						{
							m_A(i_global, j_global) += LHS_Contribution(i_local, j_local);
						}
					}
				}
			}
		}

		virtual void AssembleRHS(VectorType& RHS_Contribution,
								 Element::EquationIdVectorType& EquationId)
		{
			unsigned int local_size = RHS_Contribution.size();

			if (!m_calculate_reactions)
			{
				for (unsigned int i_local = 0; i_local < local_size; i_local++)
				{
					unsigned int i_global = EquationId[i_local];
					// MAZ_01 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					i_global = m_transformed_equation_ids[i_global];
					// MAZ_01 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					if (i_global < m_equation_system_size)
					{
						m_b[i_global] += RHS_Contribution[i_local];
					}
				}
			}
			else 
			{
				for (unsigned int i_local = 0; i_local < local_size; i_local++)
				{
					//unsigned int i_global = EquationId[i_local];
					//// MAZ_01 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					//i_global = m_transformed_equation_ids[i_global];
					//// MAZ_01 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					//if (i_global < m_equation_system_size) 
					//{
					//	m_b[i_global] += RHS_Contribution[i_local];
					//}
					//else
					//{
					//	m_r[i_global - m_equation_system_size] -= RHS_Contribution[i_local];
					//}
					// MAZ_01 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					unsigned int i_global_original = EquationId[i_local];
					unsigned int i_global = m_transformed_equation_ids[i_global_original];
					// MAZ_01 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					if (i_global < m_equation_system_size) 
					{
						m_b[i_global] += RHS_Contribution[i_local];
						if(i_global_original != i_global)
						{
							m_r[i_global_original - m_equation_system_size] -= RHS_Contribution[i_local];
						}
					}
					else
					{
						m_r[i_global - m_equation_system_size] -= RHS_Contribution[i_local];
					}
				}
			}
		}

		virtual void CalculateReactions(ModelPart& mp, 
										SchemePointerType& pScheme)
		{
			for (DofsArrayType::ptr_iterator it2 = m_dofset.ptr_begin(); it2 != m_dofset.ptr_end(); ++it2)
			{
				Dof<double>::Pointer& dofp = *it2;
				if (dofp->IsFixed())
				{
					size_t i = dofp->EquationId() - m_equation_system_size;
					dofp->GetSolutionStepReactionValue() = m_r[i];
				}
			}
		}

	public:

		inline const LinearSolverPointerType& GetLinearSolver()const { return m_lin_solver; }

		inline const bool CalculateReactions()const { return m_calculate_reactions; }
		inline bool&      CalculateReactions()      { return m_calculate_reactions; }

		inline const size_t EquationSystemSize()const { return m_equation_system_size; }

		inline const VectorType&       B()const { return m_b; }
		inline const SparseMatrixType& A()const { return m_A; }
		inline const VectorType&       X()const { return m_x; }
		inline const VectorType&       R()const { return m_r; }

		inline VectorType&             B()      { return m_b; }
		inline SparseMatrixType&       A()      { return m_A; }
		inline VectorType&             X()      { return m_x; }
		inline VectorType&             R()      { return m_r; }

		inline const DofsArrayType& DofSet()const { return m_dofset; }
		inline DofsArrayType&       DofSet()      { return m_dofset; }

		inline const IndexContainerType& TransformedEquationIds()const { return m_transformed_equation_ids; }

	protected:

		LinearSolverPointerType m_lin_solver;
		bool m_calculate_reactions;
		bool m_initialized;
		size_t m_equation_system_size;
		VectorType m_b;
		VectorType m_x;
		VectorType m_r;
		SparseMatrixType m_A;
		DofsArrayType m_dofset;
		IndexContainerType m_transformed_equation_ids;
		IndexContainerType m_equation_id_flag;

	};

} // namespace Kratos



#endif // RVE_LINEAR_SYSTEM_OF_EQUATIONS_H_INCLUDED