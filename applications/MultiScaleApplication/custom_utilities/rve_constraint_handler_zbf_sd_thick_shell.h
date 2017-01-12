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

#if !defined(RVE_CONSTRAINT_HANDLER_ZBF_SD_THICK_SHELL_H_INCLUDED)
#define RVE_CONSTRAINT_HANDLER_ZBF_SD_THICK_SHELL_H_INCLUDED

#include "rve_constraint_handler.h"
#include "../custom_conditions/rve_weak_rotation_constraint_condition_3D.h"

namespace Kratos
{

	/** \brief RveConstraintHandler_ZBF_SD_ThickShell
	*
	* Rve Constrain Handler for Zero Boundary Displacement Fluctuations
	* in Small Displacement Thick Shell formulation
	*/
	template<class TSparseSpace, 
			 class TDenseSpace>
	class RveConstraintHandler_ZBF_SD_ThickShell : public RveConstraintHandler<TSparseSpace, TDenseSpace>
	{

	public:

		KRATOS_CLASS_POINTER_DEFINITION( RveConstraintHandler_ZBF_SD_ThickShell );

		typedef RveConstraintHandler<TSparseSpace, TDenseSpace> BaseType;
		typedef typename ModelPart::DofsArrayType DofsArrayType;
		typedef typename BaseType::IndexContainerType IndexContainerType;
		typedef typename BaseType::SparseMatrixType SparseMatrixType;
		typedef typename BaseType::VectorType VectorType;
		typedef typename BaseType::SchemeType SchemeType;

	public:
		
		RveConstraintHandler_ZBF_SD_ThickShell()
			: BaseType()
			, m_already_fixed(false)
		{
		}

		virtual ~RveConstraintHandler_ZBF_SD_ThickShell()
		{
		}

	public:
		
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

			for(RveGeometryDescriptor::IndexContainerType::const_iterator it =
				geom.BoundaryNodesIDs().begin(); it != geom.BoundaryNodesIDs().end(); ++it)
			{
				RveGeometryDescriptor::IndexType index = *it;
				ModelPart::NodeType& bnd_node = mp.GetNode(index);
				bnd_node.Fix(DISPLACEMENT_X);
				bnd_node.Fix(DISPLACEMENT_Y);
				bnd_node.Fix(DISPLACEMENT_Z);
				bnd_node.Fix(ROTATION_X);
				bnd_node.Fix(ROTATION_Y);
				bnd_node.Fix(ROTATION_Z);
			}

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

	protected:

		bool m_already_fixed;
	};

} // namespace Kratos



#endif // RVE_CONSTRAINT_HANDLER_ZBF_SD_THICK_SHELL_H_INCLUDED