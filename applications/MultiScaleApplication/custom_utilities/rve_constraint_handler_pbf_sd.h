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

#if !defined(RVE_CONSTRAINT_HANDLER_PBF_SD_H_INCLUDED)
#define RVE_CONSTRAINT_HANDLER_PBF_SD_H_INCLUDED

#include "rve_constraint_handler.h"
#include "../custom_conditions/rve_weak_periodic_condition_2D2N.h"

namespace Kratos
{

	/** \brief RveConstraintHandler_PBF_SD
	*
	* Rve Constrain Handler for Periodic Boundary Displacement Fluctuations
	* in Small Displacement formulation
	*/
	template<class TSparseSpace, 
			 class TDenseSpace>
	class RveConstraintHandler_PBF_SD : public RveConstraintHandler<TSparseSpace, TDenseSpace>
	{

	public:

		KRATOS_CLASS_POINTER_DEFINITION( RveConstraintHandler_PBF_SD );

		typedef RveConstraintHandler<TSparseSpace, TDenseSpace> BaseType;
		typedef typename ModelPart::DofsArrayType DofsArrayType;
		typedef typename BaseType::IndexContainerType IndexContainerType;
		typedef typename BaseType::SparseMatrixType SparseMatrixType;
		typedef typename BaseType::VectorType VectorType;
		typedef typename BaseType::SchemeType SchemeType;

	public:
		
		RveConstraintHandler_PBF_SD()
			: BaseType()
		{
		}

		virtual ~RveConstraintHandler_PBF_SD()
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
			// MAZ_01 +++++++++++++

			int free_id = 0;
			int fix_id = dofset.size();
			equation_system_size = 0;
			size_t fixed_size = 0;

			// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			// create the transformed equation id vector
			// store the same id for free (non master non slave nodes)
			// store increasing ids for free (master and slave nodes)
			transformed_equation_ids.resize(fix_id);
			// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
					// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					transformed_equation_ids[fix_id] = fix_id; // copy as is - note: after --fix_id
					// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					fixed_size++;
				}
				else 
				{
					if(!std::binary_search(slave_node_ids.begin(), slave_node_ids.end(), dof_iterator->Id()))
					{
						// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
						transformed_equation_ids[free_id] = free_id; // copy as is - note: before fix_id++
						// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
				{ // UX
					ModelPart::DofType::Pointer slave_dof = slave_node.pGetDof( DISPLACEMENT_X );
					ModelPart::DofType::Pointer master_dof = master_node.pGetDof( DISPLACEMENT_X );
					ModelPart::DofType::EquationIdType master_dof_eq_id = master_dof->EquationId();
					slave_dof->SetEquationId( --fix_id );
					transformed_equation_ids[fix_id] = master_dof_eq_id; // map to the master eq id - note: after --fix_id
				}
				{ // UY
					ModelPart::DofType::Pointer slave_dof = slave_node.pGetDof( DISPLACEMENT_Y );
					ModelPart::DofType::Pointer master_dof = master_node.pGetDof( DISPLACEMENT_Y );
					ModelPart::DofType::EquationIdType master_dof_eq_id = master_dof->EquationId();
					slave_dof->SetEquationId( --fix_id );
					transformed_equation_ids[fix_id] = master_dof_eq_id; // map to the master eq id - note: after --fix_id
				}
				if(geom.Dimension() == 3)
				{ // UZ
					ModelPart::DofType::Pointer slave_dof = slave_node.pGetDof( DISPLACEMENT_Z );
					ModelPart::DofType::Pointer master_dof = master_node.pGetDof( DISPLACEMENT_Z );
					ModelPart::DofType::EquationIdType master_dof_eq_id = master_dof->EquationId();
					slave_dof->SetEquationId( --fix_id );
					transformed_equation_ids[fix_id] = master_dof_eq_id; // map to the master eq id - note: after --fix_id
				}
			}
		}
		
		virtual void ApplyMacroScaleData(ModelPart& mp, 
										 const RveGeometryDescriptor& geom,
										 const RveMacroscaleData& macroScaleData)
		{
			const ProcessInfo& pinfo = mp.GetProcessInfo();
			Vector E = -macroScaleData.StrainVector();

			std::vector< Vector > strain_array;
			for(ModelPart::ElementIterator elem_iter = mp.ElementsBegin(); elem_iter != mp.ElementsEnd(); ++elem_iter)
			{
				Element& ielem = *elem_iter;
				Element::GeometryType& igeom = ielem.GetGeometry();
				Element::IntegrationMethod intmethod = ielem.GetIntegrationMethod();
				const Element::GeometryType::IntegrationPointsArrayType& ipts = igeom.IntegrationPoints(intmethod);
				if(strain_array.size() != ipts.size())
					strain_array.resize(ipts.size());
				std::fill(strain_array.begin(), strain_array.end(), E);
				ielem.SetValueOnIntegrationPoints(INITIAL_STRAIN, strain_array, pinfo);
			}

			if(geom.Dimension() == 2)
			{
				for(RveGeometryDescriptor::IndexContainerType::const_iterator it =
					geom.CornerNodesIDs().begin(); it != geom.CornerNodesIDs().end(); ++it)
				{
					RveGeometryDescriptor::IndexType index = *it;
					ModelPart::NodeType& bnd_node = mp.GetNode(index);
					bnd_node.Fix(DISPLACEMENT_X);
					bnd_node.Fix(DISPLACEMENT_Y);
				}
			}
			else
			{
				for(RveGeometryDescriptor::IndexContainerType::const_iterator it =
					geom.CornerNodesIDs().begin(); it != geom.CornerNodesIDs().end(); ++it)
				{
					RveGeometryDescriptor::IndexType index = *it;
					ModelPart::NodeType& bnd_node = mp.GetNode(index);
					bnd_node.Fix(DISPLACEMENT_X);
					bnd_node.Fix(DISPLACEMENT_Y);
					bnd_node.Fix(DISPLACEMENT_Z);
				}
			}
		}

		virtual void FinalizeSolutionStep(ModelPart& mp, 
										 const RveGeometryDescriptor& geom,
										 const RveMacroscaleData& macroScaleData)
		{
			const Vector& E = macroScaleData.StrainVector();
			array_1d<double, 3> Um;
			Um.clear();
			if(geom.Dimension() == 2)
			{
				double x0 = geom.Center()[0];
				double y0 = geom.Center()[1];
				double exx = E(0);
				double eyy = E(1);
				double exy = E(2)/2.0;
				for(ModelPart::NodeIterator it = mp.NodesBegin(); it != mp.NodesEnd(); ++it)
				{
					ModelPart::NodeType& inode = *it;
					if(inode.SolutionStepsDataHas(RVE_FULL_DISPLACEMENT))
					{
						double x = inode.X0() - x0;
						double y = inode.Y0() - y0;
						const array_1d<double, 3>& Uf = inode.FastGetSolutionStepValue(DISPLACEMENT);
						Um(0) = Uf(0) + exx * x + exy * y;
						Um(1) = Uf(1) + exy * x + eyy * y;
						inode.FastGetSolutionStepValue(RVE_FULL_DISPLACEMENT) = Um;
					}
				}
			}
			else
			{
				double x0 = geom.Center()[0];
				double y0 = geom.Center()[1];
				double z0 = geom.Center()[2];
				double exx = E(0);
				double eyy = E(1);
				double ezz = E(2);
				double exy = E(3)/2.0;
				double eyz = E(4)/2.0;
				double exz = E(5)/2.0;
				for(ModelPart::NodeIterator it = mp.NodesBegin(); it != mp.NodesEnd(); ++it)
				{
					ModelPart::NodeType& inode = *it;
					if(inode.SolutionStepsDataHas(RVE_FULL_DISPLACEMENT))
					{
						double x = inode.X0() - x0;
						double y = inode.Y0() - y0;
						double z = inode.Z0() - z0;
						const array_1d<double, 3>& Uf = inode.FastGetSolutionStepValue(DISPLACEMENT);
						Um(0) = Uf(0) + exx * x + exy * y + exz * z;
						Um(1) = Uf(1) + exy * x + eyy * y + eyz * z;
						Um(2) = Uf(2) + exz * x + eyz * y + ezz * z;
						inode.FastGetSolutionStepValue(RVE_FULL_DISPLACEMENT) = Um;
					}
				}
			}
		}

		// MAZ_01 +++++++++++++++++++++++++++++++++++++++++++++++
		// this method can be used to modifiy the standard update
		// made by the scheme. 
		// TODO:
		// check if this is better than creating ad-hoc schemes
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++
		virtual void Update(ModelPart& mp,
							const RveGeometryDescriptor& geom,
							const RveMacroscaleData& macroScaleData,
							const IndexContainerType& trasformed_eq_ids,
							SchemeType& scheme,
							DofsArrayType& dofset,
							SparseMatrixType& A,
							VectorType& Dx,
							VectorType& b,
							size_t equation_system_size)
		{
			for(typename DofsArrayType::iterator i_dof = dofset.begin() ; i_dof != dofset.end() ; ++i_dof) {
				if(i_dof->IsFree()) {
					i_dof->GetSolutionStepValue() += Dx[trasformed_eq_ids[i_dof->EquationId()]];
					
				}
			}
		}

	};

} // namespace Kratos



#endif // RVE_CONSTRAINT_HANDLER_PBF_SD_H_INCLUDED