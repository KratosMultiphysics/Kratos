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

#if !defined(RVE_HOMOGENIZER_H_INCLUDED)
#define RVE_HOMOGENIZER_H_INCLUDED

#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/define.h"
#include "includes/serializer.h"
#include "rve_linear_system_of_equations.h"
#include "rve_utilities.h"
#include "rve_config.h"

namespace Kratos
{

	template<class TSparseSpace, 
			 class TDenseSpace,
			 class TReorderer = Reorderer<TSparseSpace, TDenseSpace> >
	class RveHomogenizer
	{

	public:

		KRATOS_CLASS_POINTER_DEFINITION( RveHomogenizer );
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
		
		RveHomogenizer()
		{
		}

		virtual ~RveHomogenizer()
		{
		}
		
		virtual void HomogenizeStressTensor(ModelPart& mp, 
											const RveGeometryDescriptor& geomDescriptor,
											RveLinearSystemOfEquationsPointerType& soe,
											RveConstraintHandlerPointerType& constraintHandler,
											RveMacroscaleData& macroScaleData,
											Vector& S)
		{
			if(geomDescriptor.Dimension() == 2)
				this->HomogenizeStressTensor_2D(mp,geomDescriptor,soe,constraintHandler,macroScaleData,S);
			else
				this->HomogenizeStressTensor_3D(mp,geomDescriptor,soe,constraintHandler,macroScaleData,S);
		}
		
		virtual void HomogenizeTengentConstitutiveTensor(ModelPart& mp, 
														 const RveGeometryDescriptor& geomDescriptor,
														 RveLinearSystemOfEquationsPointerType& soe,
														 RveConstraintHandlerPointerType& constraintHandler,
														 RveMacroscaleData& macroScaleData,
														 SchemePointerType& pScheme,
														 const Vector& S,
														 Matrix& C,
														 const Vector& U,
														 bool move_mesh)
		{
			size_t strain_size = macroScaleData.StrainVector().size();
			
			Vector saved_strain_vector(strain_size);
			noalias(saved_strain_vector) = macroScaleData.StrainVector();

			Vector pert_stress_vector(strain_size);

			if(C.size1() != strain_size || C.size2() != strain_size)
				C.resize(strain_size, strain_size, false);
			noalias(C) = ZeroMatrix(strain_size, strain_size);
			
			// compute the perturbation parameters
			double perturbation = 1.0E-6*norm_1(macroScaleData.StrainVectorOld());
			if(perturbation < 1.0e-7)
				perturbation = 1.0e-7;
			perturbation = 1.0e-6;

			for(size_t j = 0; j < strain_size; j++)
			{
				if(j > 0) RveUtilities::RestoreSolutionVector(mp, U);

				// apply perturbed strain vector
				noalias(macroScaleData.StrainVector()) = saved_strain_vector;
				macroScaleData.StrainVector()(j) += perturbation;
				constraintHandler->ApplyMacroScaleData(mp, geomDescriptor, macroScaleData);
				soe->BuildRHS(mp, pScheme);
				soe->Solve();

				constraintHandler->Update(mp,geomDescriptor,macroScaleData, 
				                        soe->TransformedEquationIds(), *pScheme,
				                        soe->DofSet(), soe->A(), soe->X(), soe->B(),
										soe->EquationSystemSize());
				if(move_mesh) this->MoveMesh(mp);

#if RVE_HOMOGENIZER_OPTIMIZATION == 1
				soe->BuildRHS(mp, pScheme);
#elif RVE_HOMOGENIZER_OPTIMIZATION == 2
				soe->BuildRHS_Reduced_WithReactions(mp, pScheme, geomDescriptor.BoundaryElementsIDs());
#else
				// this is just to plot reactions
				soe->BuildRHS_Reduced_WithReactions(mp, pScheme, geomDescriptor.BoundaryElementsIDs());
#endif

				this->HomogenizeStressTensor(mp, geomDescriptor, soe, constraintHandler, macroScaleData, pert_stress_vector);
				for(size_t i = 0; i < strain_size; i++)
					C(i, j) = (pert_stress_vector(i) - S(i))/perturbation;
			}
			// reset stored strain vector
			noalias(macroScaleData.StrainVector()) = saved_strain_vector;
		}
		
		virtual void HomogenizeVariable(ModelPart& mp, 
										const RveGeometryDescriptor& geomDescriptor,
										const Variable<int>& rThisVariable, 
										int& rValue)
		{
			rValue = 0;
			//double d_rValue = 0.0;
			//std::vector< int > gp_values;
			//ProcessInfo& processInfo = mp.GetProcessInfo();
			//double totalVolume(0.0);
			//for(ModelPart::ElementIterator it = mp.ElementsBegin(); it != mp.ElementsEnd(); ++it)
			//{
			//	Element& ielem = *it;
			//	Element::GeometryType& igeom = ielem.GetGeometry();
			//	Element::IntegrationMethod intmethod = ielem.GetIntegrationMethod();
			//	const Element::GeometryType::IntegrationPointsArrayType& ipts = igeom.IntegrationPoints(intmethod);

			//	ielem.GetValueOnIntegrationPoints(rThisVariable, gp_values, processInfo);
			//	if(gp_values.size() != ipts.size()) continue;

			//	for(size_t point_id = 0; point_id < ipts.size(); point_id++)
			//	{
			//		double dV = igeom.DeterminantOfJacobian(point_id, intmethod) * ipts[point_id].Weight();
			//		d_rValue += gp_values[point_id] * dV;
			//		totalVolume += dV;
			//	}
			//}
			//if(totalVolume > 0.0)
			//	d_rValue /= totalVolume;
			//else
			//	d_rValue = 0.0;
			//rValue = (int)rValue;
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
						rValue.clear();
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

			//Vector all_detJ; // <<<<<<<<<<<<<<<

			for(ModelPart::ElementIterator it = mp.ElementsBegin(); it != mp.ElementsEnd(); ++it)
			{
				Element& ielem = *it;
				Element::GeometryType& igeom = ielem.GetGeometry();
				Element::IntegrationMethod intmethod = ielem.GetIntegrationMethod();
				const Element::GeometryType::IntegrationPointsArrayType& ipts = igeom.IntegrationPoints(intmethod);

				ielem.GetValueOnIntegrationPoints(rThisVariable, gp_values, processInfo);
				if(gp_values.size() != ipts.size()) continue;

				//all_detJ = igeom.DeterminantOfJacobian(all_detJ, intmethod); // <<<<<<<<<<<<<<<<<

				for(size_t point_id = 0; point_id < ipts.size(); point_id++)
				{
					double dV = igeom.DeterminantOfJacobian(point_id, intmethod) * ipts[point_id].Weight();
					//double dV = all_detJ[point_id] * ipts[point_id].Weight();

					Matrix& i_gp_value = gp_values[point_id];
					if(rValue.size1() == 0 && rValue.size2() == 0) 
					{
						rValue.resize(i_gp_value.size1(), i_gp_value.size2(),false);
						rValue.clear();
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

		virtual void HomogenizeStressTensor_2D(ModelPart& mp, 
											   const RveGeometryDescriptor& geomDescriptor,
											   RveLinearSystemOfEquationsPointerType& soe,
											   RveConstraintHandlerPointerType& constraintHandler,
											   RveMacroscaleData& macroScaleData,
											   Vector& S)
		{
			ProcessInfo& processInfo = mp.GetProcessInfo();

			if(S.size() != 3) S.resize(3, false);
			noalias(S) = ZeroVector(3);

#if RVE_HOMOGENIZER_OPTIMIZATION == 0
			double totalVolume(0.0);

			std::vector< Matrix > stressTensors;
		
			for(ModelPart::ElementIterator it = mp.ElementsBegin(); it != mp.ElementsEnd(); ++it)
			{
				Element& ielem = *it;
				Element::GeometryType& igeom = ielem.GetGeometry();
				Element::IntegrationMethod intmethod = ielem.GetIntegrationMethod();
				const Element::GeometryType::IntegrationPointsArrayType& ipts = igeom.IntegrationPoints(intmethod);

				ielem.GetValueOnIntegrationPoints(PK2_STRESS_TENSOR, stressTensors, processInfo);
				if(stressTensors.size() != ipts.size()) continue;

				for(size_t point_id = 0; point_id < ipts.size(); point_id++)
				{
					double dV = igeom.DeterminantOfJacobian(point_id, intmethod) * ipts[point_id].Weight();
					Matrix& igpStressTensor = stressTensors[point_id];
					
					if(igpStressTensor.size1() >= 2 && igpStressTensor.size2() >= 2)
					{
						S(0) += igpStressTensor(0,0) * dV;
						S(1) += igpStressTensor(1,1) * dV;
						S(2) += igpStressTensor(0,1) * dV;

						totalVolume += dV;
					}
				}
			}

			if(totalVolume == 0.0)
				noalias(S) = ZeroVector(3);
			else
				S /= totalVolume;
				
			/*{
				double X_center = geomDescriptor.Center()[0];
				double Y_center = geomDescriptor.Center()[1];
				array_1d<double, 2> X;
				array_1d<double, 2> f;
				Matrix Sig(2,2,0.0);
				Vector S2(S.size(),0.0);
				for(RveGeometryDescriptor::IndexContainerType::const_iterator it =
					geomDescriptor.BoundaryNodesIDs().begin(); it != geomDescriptor.BoundaryNodesIDs().end(); ++it)
				{
					RveGeometryDescriptor::IndexType index = *it;
					ModelPart::NodeType& bnd_node = mp.GetNode(index);
					X[0] = bnd_node.X0() - X_center;
					X[1] = bnd_node.Y0() - Y_center;
					ModelPart::NodeType::DofType& dof_x = bnd_node.GetDof(DISPLACEMENT_X);
					ModelPart::NodeType::DofType& dof_y = bnd_node.GetDof(DISPLACEMENT_Y);
					f[0] = dof_x.GetSolutionStepReactionValue();
					f[1] = dof_y.GetSolutionStepReactionValue();
					Sig += outer_prod(f,X);
				}
				S2(0) = Sig(0,0);
				S2(1) = Sig(1,1);
				S2(2) = 0.5*(Sig(0,1)+Sig(1,0));
				if(totalVolume == 0.0)
					noalias(S2) = ZeroVector(3);
				else
					S2 /= totalVolume;
				std::stringstream ss;
				ss << "-------------------------------\n";
				ss << MathHelpers::VectorToString(S,4,std::scientific);
				ss << MathHelpers::VectorToString(S2,4, std::scientific);
				ss << "-------------------------------\n";
				std::cout << ss.str();
			}*/

#else
			double totalVolume = geomDescriptor.DomainSize();
			double X_center = geomDescriptor.Center()[0];
			double Y_center = geomDescriptor.Center()[1];
			array_1d<double, 2> X;
			array_1d<double, 2> f;
			Matrix Sig(2,2,0.0);
			for(RveGeometryDescriptor::IndexContainerType::const_iterator it =
				geomDescriptor.BoundaryNodesIDs().begin(); it != geomDescriptor.BoundaryNodesIDs().end(); ++it)
			{
				RveGeometryDescriptor::IndexType index = *it;
				ModelPart::NodeType& bnd_node = mp.GetNode(index);
				X[0] = bnd_node.X0() - X_center;
				X[1] = bnd_node.Y0() - Y_center;
				ModelPart::NodeType::DofType& dof_x = bnd_node.GetDof(DISPLACEMENT_X);
				ModelPart::NodeType::DofType& dof_y = bnd_node.GetDof(DISPLACEMENT_Y);
				f[0] = dof_x.GetSolutionStepReactionValue();
				f[1] = dof_y.GetSolutionStepReactionValue();
				Sig += outer_prod(f,X);
			}
			S(0) = Sig(0,0);
			S(1) = Sig(1,1);
			S(2) = 0.5*(Sig(0,1)+Sig(1,0));
			if(totalVolume == 0.0)
				noalias(S) = ZeroVector(3);
			else
				S /= totalVolume;
#endif
		}

		virtual void HomogenizeStressTensor_3D(ModelPart& mp, 
											   const RveGeometryDescriptor& geomDescriptor,
											   RveLinearSystemOfEquationsPointerType& soe,
											   RveConstraintHandlerPointerType& constraintHandler,
											   RveMacroscaleData& macroScaleData,
											   Vector& S)
		{
			ProcessInfo& processInfo = mp.GetProcessInfo();

			if(S.size() != 6) S.resize(6, false);
			noalias(S) = ZeroVector(6);

#if RVE_HOMOGENIZER_OPTIMIZATION == 0
			double totalVolume(0.0);

			std::vector< Matrix > stressTensors;

			for(ModelPart::ElementIterator it = mp.ElementsBegin(); it != mp.ElementsEnd(); ++it)
			{
				Element& ielem = *it;
				Element::GeometryType& igeom = ielem.GetGeometry();
				Element::IntegrationMethod intmethod = ielem.GetIntegrationMethod();
				const Element::GeometryType::IntegrationPointsArrayType& ipts = igeom.IntegrationPoints(intmethod);

				ielem.GetValueOnIntegrationPoints(PK2_STRESS_TENSOR, stressTensors, processInfo);
				if(stressTensors.size() != ipts.size()) continue;

				for(size_t point_id = 0; point_id < ipts.size(); point_id++)
				{
					double dV = igeom.DeterminantOfJacobian(point_id, intmethod) * ipts[point_id].Weight();
					Matrix& igpStressTensor = stressTensors[point_id];

					S(0) += igpStressTensor(0,0) * dV;
					S(1) += igpStressTensor(1,1) * dV;
					S(2) += igpStressTensor(2,2) * dV;
					S(3) += igpStressTensor(0,1) * dV;
					S(4) += igpStressTensor(1,2) * dV;
					S(5) += igpStressTensor(0,2) * dV;

					totalVolume += dV;
				}
			}

			if(totalVolume == 0.0)
				noalias(S) = ZeroVector(6);
			else
				S /= totalVolume;

#else

			double totalVolume = geomDescriptor.DomainSize();
			double X_center = geomDescriptor.Center()[0];
			double Y_center = geomDescriptor.Center()[1];
			double Z_center = geomDescriptor.Center()[2];
			array_1d<double, 3> X;
			array_1d<double, 3> f;
			Matrix Sig(3,3,0.0);
			for(RveGeometryDescriptor::IndexContainerType::const_iterator it =
				geomDescriptor.BoundaryNodesIDs().begin(); it != geomDescriptor.BoundaryNodesIDs().end(); ++it)
			{
				RveGeometryDescriptor::IndexType index = *it;
				ModelPart::NodeType& bnd_node = mp.GetNode(index);
				X[0] = bnd_node.X0() - X_center;
				X[1] = bnd_node.Y0() - Y_center;
				X[2] = bnd_node.Z0() - Z_center;
				ModelPart::NodeType::DofType& dof_x = bnd_node.GetDof(DISPLACEMENT_X);
				ModelPart::NodeType::DofType& dof_y = bnd_node.GetDof(DISPLACEMENT_Y);
				ModelPart::NodeType::DofType& dof_z = bnd_node.GetDof(DISPLACEMENT_Z);
				f[0] = dof_x.GetSolutionStepReactionValue();
				f[1] = dof_y.GetSolutionStepReactionValue();
				f[2] = dof_z.GetSolutionStepReactionValue();
				Sig += outer_prod(f,X);
			}
			S(0) = Sig(0,0);
			S(1) = Sig(1,1);
			S(2) = Sig(2,2);
			S(3) = 0.5*(Sig(0,1)+Sig(1,0));
			S(4) = 0.5*(Sig(1,2)+Sig(2,1));
			S(5) = 0.5*(Sig(0,2)+Sig(2,0));

			if(totalVolume == 0.0)
				noalias(S) = ZeroVector(6);
			else
				S /= totalVolume;

#endif

		}

		/**
		* Updates the nodal coordinates if necessary
		*/
		virtual void MoveMesh(ModelPart& mp)
		{
			for (ModelPart::NodeIterator i = mp.NodesBegin(); i != mp.NodesEnd(); ++i)
			{
				ModelPart::NodeType& node = *i;
				noalias(node.GetInitialPosition()) = node.GetInitialPosition() + node.FastGetSolutionStepValue(DISPLACEMENT);
			}
		}

	};

} // namespace Kratos



#endif // RVE_HOMOGENIZER_H_INCLUDED