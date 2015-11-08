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
//   Last Modified by:    $Author: Stefano Zaghi $
//   Date:                $Date: 2015-03-04 12:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#if !defined(RVE_HOMOGENIZER_THERMAL_H_INCLUDED)
#define RVE_HOMOGENIZER_THERMAL_H_INCLUDED

#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/define.h"
#include "includes/serializer.h"
#include "rve_linear_system_of_equations.h"

namespace Kratos
{

	template<class TSparseSpace, 
			 class TDenseSpace,
			 class TReorderer = Reorderer<TSparseSpace, TDenseSpace> >
	class RveHomogenizerThermal
	{

	public:

		KRATOS_CLASS_POINTER_DEFINITION( RveHomogenizerThermal );
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

	public:
		
		RveHomogenizerThermal()
		{
		}

		virtual ~RveHomogenizerThermal()
		{
		}
		
		virtual void HomogenizeHeatFluxVector(ModelPart& mp, 
											  const RveGeometryDescriptor& geomDescriptor,
											  RveLinearSystemOfEquationsPointerType& soe,
											  RveConstraintHandlerPointerType& constraintHandler,
											  Vector& q)
		{
			if(geomDescriptor.Dimension() == 2)
				this->HomogenizeHeatFluxVector_2D(mp, geomDescriptor, soe, constraintHandler, q);
			else
				this->HomogenizeHeatFluxVector_3D(mp, geomDescriptor, soe, constraintHandler, q);
		}
		
		virtual void HomogenizeTengentConductivityMatrix(ModelPart& mp, 
														 const RveGeometryDescriptor& geomDescriptor,
														 RveLinearSystemOfEquationsPointerType& soe,
														 RveConstraintHandlerPointerType& constraintHandler,
														 RveMacroscaleTemperatureData& macroScaleTempData,
														 SchemePointerType& pScheme,
														 const Vector& q,
														 Matrix& K,
														 bool move_mesh)
		{
			size_t hflux_size = geomDescriptor.Dimension() == 2 ? 3 : 6;

			Vector saved_hflux_vector(hflux_size);
			noalias(saved_hflux_vector) = macroScaleTempData.HFluxVector();

			Vector pert_hflux_vector(hflux_size);

			if (K.size1() != hflux_size || K.size2() != hflux_size)
				K.resize(hflux_size, hflux_size, false);
			noalias(K) = ZeroMatrix(hflux_size, hflux_size);

			// compute the perturbation parameters
			double norm_hflux = norm_1(saved_hflux_vector);
			double perturbation = 1.0E-6*norm_hflux;

			for (size_t j = 0; j < hflux_size; j++)
			{
				// apply perturbed strain vector
				noalias(macroScaleTempData.HFluxVector()) = saved_hflux_vector;
				macroScaleTempData.HFluxVector()(j) += perturbation;
				constraintHandler->ApplyMacroScaleData(mp, geomDescriptor, macroScaleTempData);

				soe->BuildRHS(mp, pScheme);
				soe->Solve();
				pScheme->Update(mp, soe->DofSet(), soe->A(), soe->X(), soe->B());
				if(move_mesh) this->MoveMesh(mp);

				this->HomogenizeHeatFluxVector(mp, geomDescriptor, soe, constraintHandler, pert_hflux_vector);
				
				for (size_t i = 0; i < hflux_size; i++)
					K(i, j) = (pert_hflux_vector(i) - q(i)) / perturbation;
			}
			// reset stored strain vector
			noalias(macroScaleTempData.HFluxVector()) = saved_hflux_vector;
		}
		
		virtual void HomogenizeVariable(ModelPart& mp, 
										const RveGeometryDescriptor& geomDescriptor,
										const Variable<double>& rThisVariable, 
										double& rValue)
		{
			rValue = 0.0;
			std::vector< double > gp_values;
			ProcessInfo& processInfo = mp.GetProcessInfo();
			double totalVolume(0.0);
			for(ModelPart::ElementIterator it = mp.ElementsBegin(); it != mp.ElementsEnd(); ++it)
			{
				Element& ielem = *it;
				Element::GeometryType& igeom = ielem.GetGeometry();
				Element::IntegrationMethod intmethod = ielem.GetIntegrationMethod();
				const Element::GeometryType::IntegrationPointsArrayType& ipts = igeom.IntegrationPoints(intmethod);

				ielem.GetValueOnIntegrationPoints(rThisVariable, gp_values, processInfo);
				if(gp_values.size() != ipts.size()) continue;

				for(size_t point_id = 0; point_id < ipts.size(); point_id++)
				{
					double dV = igeom.DeterminantOfJacobian(point_id, intmethod) * ipts[point_id].Weight();
					rValue += gp_values[point_id] * dV;
					totalVolume += dV;
				}
			}
			if(totalVolume > 0.0)
				rValue /= totalVolume;
			else
				rValue = 0.0;
		}

		virtual void HomogenizeVariable(ModelPart& mp, 
										const RveGeometryDescriptor& geomDescriptor,
										const Variable<Vector>& rThisVariable, 
										Vector& rValue)
		{
			rValue = Vector();
			std::vector< Vector > gp_values;
			ProcessInfo& processInfo = mp.GetProcessInfo();
			double totalVolume(0.0);
			for(ModelPart::ElementIterator it = mp.ElementsBegin(); it != mp.ElementsEnd(); ++it)
			{
				Element& ielem = *it;
				Element::GeometryType& igeom = ielem.GetGeometry();
				Element::IntegrationMethod intmethod = ielem.GetIntegrationMethod();
				const Element::GeometryType::IntegrationPointsArrayType& ipts = igeom.IntegrationPoints(intmethod);

				ielem.GetValueOnIntegrationPoints(rThisVariable, gp_values, processInfo);
				if(gp_values.size() != ipts.size()) continue;

				for(size_t point_id = 0; point_id < ipts.size(); point_id++)
				{
					double dV = igeom.DeterminantOfJacobian(point_id, intmethod) * ipts[point_id].Weight();

					Vector& i_gp_value = gp_values[point_id];
					if(rValue.size() == 0) 
					{
						rValue.resize(i_gp_value.size(),false);
						rValue += i_gp_value * dV;
						totalVolume += dV;
					}
					else
					{
						if(rValue.size() == i_gp_value.size())
						{
							rValue += i_gp_value * dV;
							totalVolume += dV;
						}
					}
				}
			}
			if(totalVolume > 0.0)
				rValue /= totalVolume;
			else
				rValue.clear();
		}

		virtual void HomogenizeVariable(ModelPart& mp, 
										const RveGeometryDescriptor& geomDescriptor,
										const Variable<Matrix>& rThisVariable, 
										Matrix& rValue)
		{
			rValue = Matrix();
			std::vector< Matrix > gp_values;
			ProcessInfo& processInfo = mp.GetProcessInfo();
			double totalVolume(0.0);
			for(ModelPart::ElementIterator it = mp.ElementsBegin(); it != mp.ElementsEnd(); ++it)
			{
				Element& ielem = *it;
				Element::GeometryType& igeom = ielem.GetGeometry();
				Element::IntegrationMethod intmethod = ielem.GetIntegrationMethod();
				const Element::GeometryType::IntegrationPointsArrayType& ipts = igeom.IntegrationPoints(intmethod);

				ielem.GetValueOnIntegrationPoints(rThisVariable, gp_values, processInfo);
				if(gp_values.size() != ipts.size()) continue;

				for(size_t point_id = 0; point_id < ipts.size(); point_id++)
				{
					double dV = igeom.DeterminantOfJacobian(point_id, intmethod) * ipts[point_id].Weight();

					Matrix& i_gp_value = gp_values[point_id];
					if(rValue.size1() == 0 && rValue.size2() == 0) 
					{
						rValue.resize(i_gp_value.size1(), i_gp_value.size2(),false);
						rValue += i_gp_value * dV;
						totalVolume += dV;
					}
					else
					{
						if(rValue.size1() == i_gp_value.size1() && rValue.size2() == i_gp_value.size2())
						{
							rValue += i_gp_value * dV;
							totalVolume += dV;
						}
					}
				}
			}
			if(totalVolume > 0.0)
				rValue /= totalVolume;
			else
				rValue.clear();
		}

		virtual void HomogenizeVariable(ModelPart& mp, 
										const RveGeometryDescriptor& geomDescriptor,
										const Variable< array_1d<double,3 > >& rThisVariable, 
										array_1d<double,3 >& rValue)
		{
			rValue.clear();
			std::vector< array_1d<double,3 > > gp_values;
			ProcessInfo& processInfo = mp.GetProcessInfo();
			double totalVolume(0.0);
			for(ModelPart::ElementIterator it = mp.ElementsBegin(); it != mp.ElementsEnd(); ++it)
			{
				Element& ielem = *it;
				Element::GeometryType& igeom = ielem.GetGeometry();
				Element::IntegrationMethod intmethod = ielem.GetIntegrationMethod();
				const Element::GeometryType::IntegrationPointsArrayType& ipts = igeom.IntegrationPoints(intmethod);

				ielem.GetValueOnIntegrationPoints(rThisVariable, gp_values, processInfo);
				if(gp_values.size() != ipts.size()) continue;

				for(size_t point_id = 0; point_id < ipts.size(); point_id++)
				{
					double dV = igeom.DeterminantOfJacobian(point_id, intmethod) * ipts[point_id].Weight();
					rValue += gp_values[point_id] * dV;
					totalVolume += dV;
				}
			}
			if(totalVolume > 0.0)
				rValue /= totalVolume;
			else
				rValue.clear();
		}

		virtual void HomogenizeVariable(ModelPart& mp, 
										const RveGeometryDescriptor& geomDescriptor,
										const Variable< array_1d<double,6 > >& rThisVariable, 
										array_1d<double,6 >& rValue)
		{
			rValue.clear();
			std::vector< array_1d<double,6 > > gp_values;
			ProcessInfo& processInfo = mp.GetProcessInfo();
			double totalVolume(0.0);
			for(ModelPart::ElementIterator it = mp.ElementsBegin(); it != mp.ElementsEnd(); ++it)
			{
				Element& ielem = *it;
				Element::GeometryType& igeom = ielem.GetGeometry();
				Element::IntegrationMethod intmethod = ielem.GetIntegrationMethod();
				const Element::GeometryType::IntegrationPointsArrayType& ipts = igeom.IntegrationPoints(intmethod);

				ielem.GetValueOnIntegrationPoints(rThisVariable, gp_values, processInfo);
				if(gp_values.size() != ipts.size()) continue;

				for(size_t point_id = 0; point_id < ipts.size(); point_id++)
				{
					double dV = igeom.DeterminantOfJacobian(point_id, intmethod) * ipts[point_id].Weight();
					rValue += gp_values[point_id] * dV;
					totalVolume += dV;
				}
			}
			if(totalVolume > 0.0)
				rValue /= totalVolume;
			else
				rValue.clear();
		}

	protected:

		virtual void HomogenizeHeatFluxVector_2D(ModelPart& mp,
											   const RveGeometryDescriptor& geomDescriptor,
											   RveLinearSystemOfEquationsPointerType& soe,
											   RveConstraintHandlerPointerType& constraintHandler,
											   Vector& q)
		{
			ProcessInfo& processInfo = mp.GetProcessInfo();

			if(S.size() != 3) q.resize(3, false);
			noalias(q) = ZeroVector(3);

			double totalVolume(0.0);

			// std::vector< Matrix > hfluxVector;
		
			for(ModelPart::ElementIterator it = mp.ElementsBegin(); it != mp.ElementsEnd(); ++it)
			{
				Element& ielem = *it;
				Element::GeometryType& igeom = ielem.GetGeometry();
				Element::IntegrationMethod intmethod = ielem.GetIntegrationMethod();
				const Element::GeometryType::IntegrationPointsArrayType& ipts = igeom.IntegrationPoints(intmethod);

				hfluxVector += igeom[i].FastGetSolutionStepValue(rSourceVar);
				ielem.GetValueOnIntegrationPoints(PK2_STRESS_TENSOR, hfluxVector, processInfo);
				if (hfluxVector.size() != ipts.size()) continue;

				for(size_t point_id = 0; point_id < ipts.size(); point_id++)
				{
					double dV = igeom.DeterminantOfJacobian(point_id, intmethod) * ipts[point_id].Weight();
					//Matrix& igpHFluxVector = hfluxVector[point_id];
					
					q(0) += hfluxVector(0) * dV;
					q(1) += hfluxVector(1) * dV;
					q(2) += hfluxVector(2) * dV;

					totalVolume += dV;
				}
			}

			if(totalVolume == 0.0)
				noalias(q) = ZeroVector(3);
			else
				q /= totalVolume;

			/*Vector s2(3,0.0);
			array_1d<double,2> iX;
			array_1d<double,2> iF;
			Matrix ms2(2,2,0.0);
			for(RveGeometryDescriptor::IndexContainerType::const_iterator it =
				geomDescriptor.CornerNodesIDs().begin(); it != geomDescriptor.CornerNodesIDs().end(); ++it)
			{
				ModelPart::NodeType& inode = mp.GetNode(*it);
				iX[0] = inode.X0();
				iX[1] = inode.Y0();
				iF[0] = inode.GetDof(DISPLACEMENT_X).GetSolutionStepReactionValue();
				iF[1] = inode.GetDof(DISPLACEMENT_Y).GetSolutionStepReactionValue();
				ms2 += outer_prod(iF, iX);
			}
			ms2 /= totalVolume;
			s2(0) = ms2(0,0);
			s2(1) = ms2(1,1);
			s2(2) = (ms2(0,1)+ms2(1,0))/2.0;
			KRATOS_WATCH(S);
			KRATOS_WATCH(ms2);*/
		}

		virtual void HomogenizeHeatFluxVector_3D(ModelPart& mp,
											   const RveGeometryDescriptor& geomDescriptor,
											   RveLinearSystemOfEquationsPointerType& soe,
											   RveConstraintHandlerPointerType& constraintHandler,
											   Vector& q)
		{
			ProcessInfo& processInfo = mp.GetProcessInfo();

			if(q.size() != 6) q.resize(6, false);
			noalias(q) = ZeroVector(6);

			double totalVolume(0.0);

			std::vector< Matrix > hfluxVector;
		
			for(ModelPart::ElementIterator it = mp.ElementsBegin(); it != mp.ElementsEnd(); ++it)
			{
				Element& ielem = *it;
				Element::GeometryType& igeom = ielem.GetGeometry();
				Element::IntegrationMethod intmethod = ielem.GetIntegrationMethod();
				const Element::GeometryType::IntegrationPointsArrayType& ipts = igeom.IntegrationPoints(intmethod);

				ielem.GetValueOnIntegrationPoints(PK2_STRESS_TENSOR, hfluxVector, processInfo);
				if (hfluxVector.size() != ipts.size()) continue;

				for(size_t point_id = 0; point_id < ipts.size(); point_id++)
				{
					double dV = igeom.DeterminantOfJacobian(point_id, intmethod) * ipts[point_id].Weight();
					Matrix& igpHFluxVector = hfluxVector[point_id];
					
					q(0) += hfluxVector(0) * dV;
					q(1) += hfluxVector(1) * dV;
					q(2) += hfluxVector(2) * dV;
					q(3) += hfluxVector(3) * dV;
					q(4) += hfluxVector(4) * dV;
					q(5) += hfluxVector(5) * dV;

					totalVolume += dV;
				}
			}

			if(totalVolume == 0.0)
				noalias(q) = ZeroVector(6);
			else
				q /= totalVolume;
		}

		/**
		* Updates the nodal coordinates if necessary
		*/
		virtual void MoveMesh(ModelPart& mp)
		{
			// Not necessary ??
			//for (ModelPart::NodeIterator i = mp.NodesBegin(); i != mp.NodesEnd(); ++i)
			//{
			//	ModelPart::NodeType& node = *i;
			//	noalias(node.GetInitialPosition()) = node.GetInitialPosition() + node.FastGetSolutionStepValue(DISPLACEMENT);
			//}
		}

	};

} // namespace Kratos



#endif // RVE_HOMOGENIZER_THERMAL_H_INCLUDED