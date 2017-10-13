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

#if !defined(RVE_HOMOGENIZER_THICK_SHELL_H_INCLUDED)
#define RVE_HOMOGENIZER_THICK_SHELL_H_INCLUDED

#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/define.h"
#include "includes/serializer.h"
#include "rve_linear_system_of_equations.h"
#include "rve_utilities.h"
#include "rve_homogenizer.h"
#include "rve_config.h"

namespace Kratos
{

	template<class TSparseSpace, 
			 class TDenseSpace,
			 class TReorderer = Reorderer<TSparseSpace, TDenseSpace> >
	class RveHomogenizerThickShell : public RveHomogenizer<TSparseSpace, TDenseSpace>
	{

	public:

		KRATOS_CLASS_POINTER_DEFINITION( RveHomogenizerThickShell );
		typedef RveLinearSystemOfEquations<TSparseSpace, TDenseSpace, TReorderer> RveLinearSystemOfEquationsType;
		typedef typename RveLinearSystemOfEquationsType::Pointer RveLinearSystemOfEquationsPointerType;
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
		typedef RveHomogenizer<TSparseSpace, TDenseSpace> BaseType;

	public:
		
		RveHomogenizerThickShell()
			: BaseType()
		{
		}

		virtual ~RveHomogenizerThickShell()
		{
		}
		
		virtual void HomogenizeStressTensor(ModelPart& mp, 
											const RveGeometryDescriptor& geomDescriptor,
											RveLinearSystemOfEquationsPointerType& soe,
											RveConstraintHandlerPointerType& constraintHandler,
											RveMacroscaleData& macroScaleData,
											Vector& S)
		{
			ProcessInfo& processInfo = mp.GetProcessInfo();

			if(S.size() != 8) S.resize(8, false);
			noalias(S) = ZeroVector(8);

#if RVE_HOMOGENIZER_OPTIMIZATION == 0
			double totalVolume(0.0);

			std::vector< Matrix > generalizedForceTensors;
			std::vector< Matrix > generalizedMomentTensors;
		
			//for(ModelPart::ElementIterator it = mp.ElementsBegin(); it != mp.ElementsEnd(); ++it)
			//{
			//	Element& ielem = *it;
			//	Element::GeometryType& igeom = ielem.GetGeometry();
			//	Element::IntegrationMethod intmethod = ielem.GetIntegrationMethod();
			//	const Element::GeometryType::IntegrationPointsArrayType& ipts = igeom.IntegrationPoints(intmethod);
			//
			//	ielem.GetValueOnIntegrationPoints(SHELL_FORCE_GLOBAL, generalizedForceTensors, processInfo);
			//	if(generalizedForceTensors.size() != ipts.size()) continue;
			//	ielem.GetValueOnIntegrationPoints(SHELL_MOMENT_GLOBAL, generalizedMomentTensors, processInfo);
			//	if(generalizedMomentTensors.size() != ipts.size()) continue;
			//
			//	for(size_t point_id = 0; point_id < ipts.size(); point_id++)
			//	{
			//		double dV = igeom.DeterminantOfJacobian(point_id, intmethod) * ipts[point_id].Weight();
			//		Matrix& igpGeneralizedForceTensor = generalizedForceTensors[point_id];
			//		Matrix& igpGeneralizedMomentTensor = generalizedMomentTensors[point_id];
			//		
			//		if(igpGeneralizedForceTensor.size1() >= 3 && igpGeneralizedForceTensor.size2() >= 3 &&
			//		   igpGeneralizedMomentTensor.size1() >= 3 && igpGeneralizedMomentTensor.size2() >= 3)
			//		{
			//			S(0) += igpGeneralizedForceTensor(0,0) * dV;
			//			S(1) += igpGeneralizedForceTensor(1,1) * dV;
			//			S(2) += igpGeneralizedForceTensor(0,1) * dV;
			//			
			//			S(3) += igpGeneralizedMomentTensor(0,0) * dV;
			//			S(4) += igpGeneralizedMomentTensor(1,1) * dV;
			//			S(5) += igpGeneralizedMomentTensor(0,1) * dV;
			//			
			//			S(6) += igpGeneralizedForceTensor(1,2) * dV;
			//			S(7) += igpGeneralizedForceTensor(0,2) * dV;
			//			
			//			totalVolume += dV;
			//		}
			//	}
			//}

			if(totalVolume == 0.0)
				noalias(S) = ZeroVector(3);
			else
				S /= totalVolume;
			
#else
			double totalVolume = geomDescriptor.DomainSize();
			double X_center = geomDescriptor.Center()[0];
			double Y_center = geomDescriptor.Center()[1];
			for(RveGeometryDescriptor::IndexContainerType::const_iterator it =
				geomDescriptor.BoundaryNodesIDs().begin(); it != geomDescriptor.BoundaryNodesIDs().end(); ++it)
			{
				RveGeometryDescriptor::IndexType index = *it;
				ModelPart::NodeType& bnd_node = mp.GetNode(index);
				double x = bnd_node.X0() - X_center;
				double y = bnd_node.Y0() - Y_center;
				ModelPart::NodeType::DofType& dof_ux = bnd_node.GetDof(DISPLACEMENT_X);
				ModelPart::NodeType::DofType& dof_uy = bnd_node.GetDof(DISPLACEMENT_Y);
				ModelPart::NodeType::DofType& dof_uz = bnd_node.GetDof(DISPLACEMENT_Z);
				ModelPart::NodeType::DofType& dof_rx = bnd_node.GetDof(ROTATION_X);
				ModelPart::NodeType::DofType& dof_ry = bnd_node.GetDof(ROTATION_Y);
				double fx = dof_ux.GetSolutionStepReactionValue();
				double fy = dof_uy.GetSolutionStepReactionValue();
				double fz = dof_uz.GetSolutionStepReactionValue();
				double mx = dof_rx.GetSolutionStepReactionValue();
				double my = dof_ry.GetSolutionStepReactionValue();
				S(0) +=   fx*x;
				S(1) +=   fy*y;
				S(2) +=   (fy*x)/2.0 + (fx*y)/2.0;
				S(3) +=   my*x - (fz*x*x)/2.0;
				S(4) += - mx*y - (fz*y*y)/2.0;
				S(5) +=   (my*y)/2.0 - (mx*x)/2.0 - (fz*x*y)/2.0;
				S(6) +=   fz*y;
				S(7) +=   fz*x;
			}
			if(totalVolume == 0.0)
				noalias(S) = ZeroVector(8);
			else
				S /= totalVolume;
#endif
		}

	};

} // namespace Kratos



#endif // RVE_HOMOGENIZER_THICK_SHELL_H_INCLUDED