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

#if !defined(RVE_CONSTRAINT_HANDLER_PBF_SD_THERMAL_H_INCLUDED)
#define RVE_CONSTRAINT_HANDLER_PBF_SD_THERMAL_H_INCLUDED

#include "rve_constraint_handler.h"

namespace Kratos
{

	/** \brief RveConstraintHandler_PBF_SD_THERMAL
	*
	* Rve Constrain Handler for Periodic Boundary Displacement Fluctuations
	* in Small Displacement formulation
	*/
	template<class TSparseSpace, 
			 class TDenseSpace>
	class RveConstraintHandler_PBF_SD_THERMAL : public RveConstraintHandler<TSparseSpace, TDenseSpace>
	{

	public:

		KRATOS_CLASS_POINTER_DEFINITION( RveConstraintHandler_PBF_SD_THERMAL );

		typedef RveConstraintHandler<TSparseSpace, TDenseSpace> BaseType;

		typedef typename ModelPart::DofsArrayType DofsArrayType;

	public:
		
		RveConstraintHandler_PBF_SD_THERMAL()
			: BaseType()
		{
		}

		virtual ~RveConstraintHandler_PBF_SD_THERMAL()
		{
		}

	public:

		virtual void SetUpSystem(ModelPart& mp, 
								 const RveGeometryDescriptor& geom,
								 DofsArrayType& dofset, 
								 size_t& equation_system_size,
								 IndexContainerType& transformed_equation_ids,
								 IndexContainerType& equation_id_flag)
		{
			int free_id = 0;
			int fix_id = dofset.size();
			equation_system_size = 0;
			size_t fixed_size = 0;

			size_t cnt = 0;
			std::vector<RveGeometryDescriptor::IndexType> slave_node_ids(geom.PeriodicNodesIDs().size());
			for(RveGeometryDescriptor::PeriodicIndexContainerType::const_iterator ms_it = 
				geom.PeriodicNodesIDs().begin(); ms_it != geom.PeriodicNodesIDs().end(); ++ms_it)
			{
				slave_node_ids[cnt++] = ms_it->second;
			}
			std::sort(slave_node_ids.begin(), slave_node_ids.end());

			for (typename DofsArrayType::iterator dof_iterator = dofset.begin(); dof_iterator != dofset.end(); ++dof_iterator)
			{
				if (dof_iterator->IsFixed())
				{
					dof_iterator->SetEquationId(--fix_id);
					fixed_size++;
				}
				else 
				{
					if(!std::binary_search(slave_node_ids.begin(), slave_node_ids.end(), dof_iterator->Id()))
					{
						dof_iterator->SetEquationId(free_id++);
						equation_system_size++;
					}
				}
			}

			for(RveGeometryDescriptor::PeriodicIndexContainerType::const_iterator ms_it = 
				geom.PeriodicNodesIDs().begin(); ms_it != geom.PeriodicNodesIDs().end(); ++ms_it)
			{
				RveGeometryDescriptor::IndexType master_id = ms_it->first;
				RveGeometryDescriptor::IndexType slave_id  = ms_it->second;
				ModelPart::NodeType& master_node = mp.GetNode(master_id);
				ModelPart::NodeType& slave_node  = mp.GetNode(slave_id);
				ModelPart::NodeType::DofsContainerType& master_dofs = master_node.GetDofs();
				for(ModelPart::NodeType::DofsContainerType::iterator dof_iter = master_dofs.begin(); dof_iter != master_dofs.end(); ++dof_iter)
				{
					if(slave_node.HasDofFor(dof_iter->GetVariable()))
					{
						ModelPart::DofType::Pointer slave_dof = slave_node.pGetDof( dof_iter->GetVariable() );
						ModelPart::DofType::Pointer& master_dof = *(dof_iter.base());
						ModelPart::DofType::EquationIdType master_dof_eq_id = dof_iter->EquationId();
						ModelPart::DofType::EquationIdType slave_dof_eq_id  = slave_dof->EquationId();
						slave_dof->SetEquationId( master_dof_eq_id );
					}
				}
			}
		}
		
		virtual void ApplyMacroScaleData(ModelPart& mp, 
										 const RveGeometryDescriptor& geom,
										 const RveMacroscaleTemperatureData& macroScaleData)
		{
			const Element::
			const ProcessInfo& pinfo = mp.GetProcessInfo();
			double T_mean = macroScaleData.Mean_Temp();
			Matrix grad_T = macroScaleData.mGrad_T();

			for (ModelPart::NodeIterator it = mp.NodesBegin(); it != mp.NodesEnd(); ++it)
			{
				it->FastGetSolutionStepValue(TEMPERATURE) = Tm;
				ModelPart::NodeType& index = *it;
				index.Fix(TEMPERATURE);
			}
		}

		virtual void FinalizeSolutionStep(ModelPart& mp, 
										 const RveGeometryDescriptor& geom,
										 const RveMacroscaleTemperatureData& macroScaleData)
		{
			const Matrix& grad_T = macroScaleData.mGrad_T();
			array_1d<double, 3> Tm;
			Tm.clear();
			double T_mean = macroScaleData.Mean_Temp();
			for(ModelPart::NodeIterator it = mp.NodesBegin(); it != mp.NodesEnd(); ++it)
			{
				ModelPart::NodeType& inode = *it;
				if(inode.SolutionStepsDataHas(RVE_FULL_TEMPERATURE))
				{
					const array_1d<double, 3>& Tf = inode.FastGetSolutionStepValue(TEMPERATURE);
					Tm(0) = Tf(0) + T_mean + grad_T;

					it->FastGetSolutionStepValue(RVE_FULL_TEMPERATURE) = Tm;
				}
			}
		}

	};

} // namespace Kratos



#endif // RVE_CONSTRAINT_HANDLER_PBF_SD_THERMAL_H_INCLUDED