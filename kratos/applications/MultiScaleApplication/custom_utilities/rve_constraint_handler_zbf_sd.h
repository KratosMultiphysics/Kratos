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

#if !defined(RVE_CONSTRAINT_HANDLER_ZBF_SD_H_INCLUDED)
#define RVE_CONSTRAINT_HANDLER_ZBF_SD_H_INCLUDED

#include "rve_constraint_handler.h"

namespace Kratos
{

	/** \brief RveConstraintHandler_ZBF_SD
	*
	* Rve Constrain Handler for Zero Boundary Displacement Fluctuations
	* in Small Displacement formulation
	*/
	template<class TSparseSpace, 
			 class TDenseSpace>
	class RveConstraintHandler_ZBF_SD : public RveConstraintHandler<TSparseSpace, TDenseSpace>
	{

	public:

		KRATOS_CLASS_POINTER_DEFINITION( RveConstraintHandler_ZBF_SD );

		typedef RveConstraintHandler<TSparseSpace, TDenseSpace> BaseType;

	public:
		
		RveConstraintHandler_ZBF_SD()
			: BaseType()
		{
		}

		virtual ~RveConstraintHandler_ZBF_SD()
		{
		}

	public:

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
					geom.BoundaryNodesIDs().begin(); it != geom.BoundaryNodesIDs().end(); ++it)
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
					geom.BoundaryNodesIDs().begin(); it != geom.BoundaryNodesIDs().end(); ++it)
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

	};

} // namespace Kratos



#endif // RVE_CONSTRAINT_HANDLER_ZBF_SD_H_INCLUDED