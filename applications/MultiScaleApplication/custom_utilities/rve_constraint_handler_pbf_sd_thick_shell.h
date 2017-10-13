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

#if !defined(RVE_CONSTRAINT_HANDLER_PBF_SD_THICK_SHELL_H_INCLUDED)
#define RVE_CONSTRAINT_HANDLER_PBF_SD_THICK_SHELL_H_INCLUDED

#include "rve_constraint_handler.h"
#include "../custom_conditions/rve_weak_rotation_constraint_condition_3D.h"

//#define RVE_CH_PBF_SHELL_DONT_FIX_CORNER_ROT
#define RVE_SH_PBF_SHELL_FIX_MICRO_DRILLING_ROT

namespace Kratos
{

	/** \brief RveConstraintHandler_PBF_SD_ThickShell
	*
	* Rve Constrain Handler for Periodic Boundary Displacement Fluctuations
	* in Small Displacement Thick Shell formulation
	*/
	template<class TSparseSpace, 
			 class TDenseSpace>
	class RveConstraintHandler_PBF_SD_ThickShell : public RveConstraintHandler<TSparseSpace, TDenseSpace>
	{

	public:

		KRATOS_CLASS_POINTER_DEFINITION( RveConstraintHandler_PBF_SD_ThickShell );

		typedef RveConstraintHandler<TSparseSpace, TDenseSpace> BaseType;
		typedef typename ModelPart::DofsArrayType DofsArrayType;
		typedef typename BaseType::IndexContainerType IndexContainerType;
		typedef typename BaseType::SparseMatrixType SparseMatrixType;
		typedef typename BaseType::VectorType VectorType;
		typedef typename BaseType::SchemeType SchemeType;

	public:
		
		RveConstraintHandler_PBF_SD_ThickShell()
			: BaseType()
			, m_already_fixed(false)
		{
		}

		virtual ~RveConstraintHandler_PBF_SD_ThickShell()
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

			// create the transformed equation id vector
			// store the same id for free (non master non slave nodes)
			// store increasing ids for free (master and slave nodes)
			transformed_equation_ids.resize(fix_id);

			size_t cnt = 0;
#ifdef RVE_CH_PBF_SHELL_DONT_FIX_CORNER_ROT
			std::vector<RveGeometryDescriptor::IndexType> slave_node_ids(geom.PeriodicNodesIDs().size()+3);
#else
			std::vector<RveGeometryDescriptor::IndexType> slave_node_ids(geom.PeriodicNodesIDs().size());
#endif // RVE_CH_PBF_SHELL_DONT_FIX_CORNER_ROT
			for(RveGeometryDescriptor::PeriodicIndexContainerType::const_iterator ms_it = 
				geom.PeriodicNodesIDs().begin(); ms_it != geom.PeriodicNodesIDs().end(); ++ms_it)
			{
				slave_node_ids[cnt++] = ms_it->second;
			}
#ifdef RVE_CH_PBF_SHELL_DONT_FIX_CORNER_ROT
			slave_node_ids[cnt++] = geom.CornerNodesIDs()[1];
			slave_node_ids[cnt++] = geom.CornerNodesIDs()[2];
			slave_node_ids[cnt++] = geom.CornerNodesIDs()[3];
#endif // RVE_CH_PBF_SHELL_DONT_FIX_CORNER_ROT
			std::sort(slave_node_ids.begin(), slave_node_ids.end());

			for (typename DofsArrayType::iterator dof_iterator = dofset.begin(); dof_iterator != dofset.end(); ++dof_iterator)
			{
				if (dof_iterator->IsFixed())
				{
					dof_iterator->SetEquationId(--fix_id);
					transformed_equation_ids[fix_id] = fix_id; // copy as is - note: after --fix_id
					fixed_size++;
				}
				else 
				{
#ifdef RVE_CH_PBF_SHELL_DONT_FIX_CORNER_ROT
					if(dof_iterator->GetVariable() == DISPLACEMENT_X || dof_iterator->GetVariable() == DISPLACEMENT_Y ||
					   dof_iterator->GetVariable() == DISPLACEMENT_Z ||
					   dof_iterator->GetVariable() == ROTATION_X || dof_iterator->GetVariable() == ROTATION_Y)
					{
						if(!std::binary_search(slave_node_ids.begin(), slave_node_ids.end(), dof_iterator->Id()))
						{
							transformed_equation_ids[free_id] = free_id; // copy as is - note: before fix_id++
							dof_iterator->SetEquationId(free_id++);
							equation_system_size++;
						}
					}
					else
					{
							if(!std::binary_search(slave_node_ids.begin(), slave_node_ids.end(), dof_iterator->Id()))
							{
								transformed_equation_ids[free_id] = free_id; // copy as is - note: before fix_id++
								dof_iterator->SetEquationId(free_id++);
								equation_system_size++;
							}
							else
							{
								if( dof_iterator->Id() == geom.CornerNodesIDs()[1] || 
									dof_iterator->Id() == geom.CornerNodesIDs()[2] ||
									dof_iterator->Id() == geom.CornerNodesIDs()[3] )
								{
									transformed_equation_ids[free_id] = free_id; // copy as is - note: before fix_id++
									dof_iterator->SetEquationId(free_id++);
									equation_system_size++;
								}
							}
					}
#else
					if(!std::binary_search(slave_node_ids.begin(), slave_node_ids.end(), dof_iterator->Id()))
					{
						transformed_equation_ids[free_id] = free_id; // copy as is - note: before fix_id++
						dof_iterator->SetEquationId(free_id++);
						equation_system_size++;
					}
#endif // RVE_CH_PBF_SHELL_DONT_FIX_CORNER_ROT
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
				{ // UZ
					ModelPart::DofType::Pointer slave_dof = slave_node.pGetDof( DISPLACEMENT_Z );
					ModelPart::DofType::Pointer master_dof = master_node.pGetDof( DISPLACEMENT_Z );
					ModelPart::DofType::EquationIdType master_dof_eq_id = master_dof->EquationId();
					slave_dof->SetEquationId( --fix_id );
					transformed_equation_ids[fix_id] = master_dof_eq_id; // map to the master eq id - note: after --fix_id
				}	
				{ // RX
					ModelPart::DofType::Pointer slave_dof = slave_node.pGetDof( ROTATION_X );
					ModelPart::DofType::Pointer master_dof = master_node.pGetDof( ROTATION_X );
					ModelPart::DofType::EquationIdType master_dof_eq_id = master_dof->EquationId();
					slave_dof->SetEquationId( --fix_id );
					transformed_equation_ids[fix_id] = master_dof_eq_id; // map to the master eq id - note: after --fix_id
				}
				{ // RY
					ModelPart::DofType::Pointer slave_dof = slave_node.pGetDof( ROTATION_Y );
					ModelPart::DofType::Pointer master_dof = master_node.pGetDof( ROTATION_Y );
					ModelPart::DofType::EquationIdType master_dof_eq_id = master_dof->EquationId();
					slave_dof->SetEquationId( --fix_id );
					transformed_equation_ids[fix_id] = master_dof_eq_id; // map to the master eq id - note: after --fix_id
				}
#ifndef RVE_SH_PBF_SHELL_FIX_MICRO_DRILLING_ROT
				{ // RZ
					ModelPart::DofType::Pointer slave_dof = slave_node.pGetDof( ROTATION_Z );
					ModelPart::DofType::Pointer master_dof = master_node.pGetDof( ROTATION_Z );
					ModelPart::DofType::EquationIdType master_dof_eq_id = master_dof->EquationId();
					slave_dof->SetEquationId( --fix_id );
					transformed_equation_ids[fix_id] = master_dof_eq_id; // map to the master eq id - note: after --fix_id
				}
#endif // !RVE_SH_PBF_SHELL_FIX_MICRO_DRILLING_ROT
			}
#ifdef RVE_CH_PBF_SHELL_DONT_FIX_CORNER_ROT
			for(unsigned int i=1; i<4; i++)
			{
				RveGeometryDescriptor::IndexType master_id = geom.CornerNodesIDs()[0];
				RveGeometryDescriptor::IndexType slave_id  = geom.CornerNodesIDs()[i];
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
				{ // UZ
					ModelPart::DofType::Pointer slave_dof = slave_node.pGetDof( DISPLACEMENT_Z );
					ModelPart::DofType::Pointer master_dof = master_node.pGetDof( DISPLACEMENT_Z );
					ModelPart::DofType::EquationIdType master_dof_eq_id = master_dof->EquationId();
					slave_dof->SetEquationId( --fix_id );
					transformed_equation_ids[fix_id] = master_dof_eq_id; // map to the master eq id - note: after --fix_id
				}
				{ // RX
					ModelPart::DofType::Pointer slave_dof = slave_node.pGetDof( ROTATION_X );
					ModelPart::DofType::Pointer master_dof = master_node.pGetDof( ROTATION_X );
					ModelPart::DofType::EquationIdType master_dof_eq_id = master_dof->EquationId();
					slave_dof->SetEquationId( --fix_id );
					transformed_equation_ids[fix_id] = master_dof_eq_id; // map to the master eq id - note: after --fix_id
				}
				{ // RY
					ModelPart::DofType::Pointer slave_dof = slave_node.pGetDof( ROTATION_Y );
					ModelPart::DofType::Pointer master_dof = master_node.pGetDof( ROTATION_Y );
					ModelPart::DofType::EquationIdType master_dof_eq_id = master_dof->EquationId();
					slave_dof->SetEquationId( --fix_id );
					transformed_equation_ids[fix_id] = master_dof_eq_id; // map to the master eq id - note: after --fix_id
				}
			}
#endif // RVE_CH_PBF_SHELL_DONT_FIX_CORNER_ROT
		}
		
		virtual void AddConditions(ModelPart& mp, 
								   const RveGeometryDescriptor& geom)
		{
			Properties::Pointer cnd_prop(new Properties(0));
			size_t cnd_id = 0;

			// create a constraint condition for each surface element
			ModelPart::NodeIterator node_it = mp.NodesEnd();
			--node_it;
			ModelPart::NodeType::Pointer master_node = mp.pGetNode(node_it->Id());
			for(ModelPart::ElementIterator it = mp.ElementsBegin(); it != mp.ElementsEnd(); ++it)
			{
				Element& elem = *it;
				Element::GeometryType& elem_geom = elem.GetGeometry();
				Element::NodesArrayType new_nodes;
				for(unsigned int i = 0; i < elem_geom.PointsNumber(); i++)
					new_nodes.push_back(elem_geom.pGetPoint(i));
				Element::GeometryType::Pointer cnd_geom = elem_geom.Create(new_nodes);
				// NOTE: use the same ID of the element to have a 1-to-1 match with IDs
				RveWeakRotationCondition3D::Pointer cnd(new RveWeakRotationCondition3D(elem.GetId(), cnd_geom, cnd_prop));
				cnd->SetLagrangianNode(master_node);
				mp.AddCondition(cnd);
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
					strain_array.resize(ipts.size(), E);
				ielem.SetValueOnIntegrationPoints(INITIAL_STRAIN, strain_array, pinfo);
			}

			if(m_already_fixed) return;

#ifdef RVE_CH_PBF_SHELL_DONT_FIX_CORNER_ROT
			/*for(RveGeometryDescriptor::IndexContainerType::const_iterator it =
				geom.CornerNodesIDs().begin(); it != geom.CornerNodesIDs().end(); ++it)
			{
				RveGeometryDescriptor::IndexType index = *it;
				ModelPart::NodeType& bnd_node = mp.GetNode(index);
				bnd_node.Fix(DISPLACEMENT_X);
				bnd_node.Fix(DISPLACEMENT_Y);
				bnd_node.Fix(DISPLACEMENT_Z);
			}*/
			std::vector<unsigned int> fixed_nodes;
			// case x1
			fixed_nodes.push_back(8);// 23 (B) 11(A)
			// case x2
			/*fixed_nodes.push_back(151);
			fixed_nodes.push_back(79);
			fixed_nodes.push_back(43);
			fixed_nodes.push_back(133);*/
			// case x3
			/*fixed_nodes.push_back(347);
			fixed_nodes.push_back(247);
			fixed_nodes.push_back(155);
			fixed_nodes.push_back(321);
			fixed_nodes.push_back(201);
			fixed_nodes.push_back(79);
			fixed_nodes.push_back(300);
			fixed_nodes.push_back(171);
			fixed_nodes.push_back(43);*/

			for(std::vector<unsigned int>::const_iterator it =
				fixed_nodes.begin(); it != fixed_nodes.end(); ++it)
			{
				RveGeometryDescriptor::IndexType index = *it;
				ModelPart::NodeType& bnd_node = mp.GetNode(index);
				bnd_node.Fix(DISPLACEMENT_X);
				bnd_node.Fix(DISPLACEMENT_Y);
				bnd_node.Fix(DISPLACEMENT_Z);
				bnd_node.Fix(ROTATION_X);
				bnd_node.Fix(ROTATION_Y);
			}
#else
			for(RveGeometryDescriptor::IndexContainerType::const_iterator it =
				geom.CornerNodesIDs().begin(); it != geom.CornerNodesIDs().end(); ++it)
			{
				RveGeometryDescriptor::IndexType index = *it;
				ModelPart::NodeType& bnd_node = mp.GetNode(index);
				bnd_node.Fix(DISPLACEMENT_X);
				bnd_node.Fix(DISPLACEMENT_Y);
				bnd_node.Fix(DISPLACEMENT_Z);
				bnd_node.Fix(ROTATION_X);
				bnd_node.Fix(ROTATION_Y);
#ifndef RVE_SH_PBF_SHELL_FIX_MICRO_DRILLING_ROT
				bnd_node.Fix(ROTATION_Z);
#endif // !RVE_SH_PBF_SHELL_FIX_MICRO_DRILLING_ROT
			}
#endif // RVE_CH_PBF_SHELL_DONT_FIX_CORNER_ROT
			
#ifdef RVE_SH_PBF_SHELL_FIX_MICRO_DRILLING_ROT
			for(ModelPart::NodeIterator it = mp.NodesBegin(); it != mp.NodesEnd(); ++it)
			{
				RveGeometryDescriptor::IndexType index = *it;
				ModelPart::NodeType& inode = *it;
				inode.Fix(ROTATION_Z);
			}  
#endif // !RVE_SH_PBF_SHELL_FIX_MICRO_DRILLING_ROT

			m_already_fixed = true;
		}

		virtual void FinalizeSolutionStep(ModelPart& mp, 
										 const RveGeometryDescriptor& geom,
										 const RveMacroscaleData& macroScaleData)
		{
			const Vector& E = macroScaleData.StrainVector();
			array_1d<double, 3> Um;
			array_1d<double, 3> Rm;
			Um.clear();
			Rm.clear();
			
			double x0 = geom.Center()[0];
			double y0 = geom.Center()[1];
			double exx = E(0);
			double eyy = E(1);
			double exy = E(2)/2.0;
			double kxx = E(3);
			double kyy = E(4);
			double kxy = E(5)/2.0;
			double gyz = E(6);
			double gxz = E(7);
			for(ModelPart::NodeIterator it = mp.NodesBegin(); it != mp.NodesEnd(); ++it)
			{
				ModelPart::NodeType& inode = *it;
				if(inode.SolutionStepsDataHas(RVE_FULL_DISPLACEMENT) && inode.SolutionStepsDataHas(RVE_FULL_ROTATION))
				{
					double x = inode.X0() - x0;
					double y = inode.Y0() - y0;
					const array_1d<double, 3>& Uf = inode.FastGetSolutionStepValue(DISPLACEMENT);
					const array_1d<double, 3>& Rf = inode.FastGetSolutionStepValue(ROTATION);
					
					// macro = E*X
					Um(0) = Uf(0) + exx*x + exy*y;
					Um(1) = Uf(1) + exy*x + eyy*y;
					
					// macro = -[K*X]^T*X/2 + G*X = [G-([K*X]^T)/2]*X
					Um(2) = Uf(2) - (x*(kxx*x + kxy*y))/2.0 - (y*(kxy*x + kyy*y))/2.0 + gxz*x + gyz*y;
					
					// R = K*X (without sign convention for rotations)
					// R = P'K*X (with P = [0 1; -1 0], such that P*R = [ry;-rx], K*X = [ry;-rx] -> P'*K*X = [rx; ry])
					Rm(0) = Rf(0) - kxy*x - kyy*y;
					Rm(1) = Rf(1) + kxx*x + kxy*y;
					
					inode.FastGetSolutionStepValue(RVE_FULL_DISPLACEMENT) = Um;
					inode.FastGetSolutionStepValue(RVE_FULL_ROTATION) = Rm;
				}
			}
		}

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

	protected:

		bool m_already_fixed;
	};

} // namespace Kratos



#endif // RVE_CONSTRAINT_HANDLER_PBF_SD_THICK_SHELL_H_INCLUDED