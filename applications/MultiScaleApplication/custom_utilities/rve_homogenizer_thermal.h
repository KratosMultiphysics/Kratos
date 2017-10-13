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
#include "rve_utilities.h"

//#define RVE_OPTIMIZATION 0 // sigma = 1/V(int(sigma_i))
//#define RVE_OPTIMIZATION 1 // sigma = 1/V(sum(boundary_F X))
#define  RVE_OPTIMIZATION 2 // sigma = 1/V(sum(boundary_F X))(with build rhs reduced)

namespace Kratos
{

	template<class TSparseSpace, 
			 class TDenseSpace,
			 class TReorderer = Reorderer<TSparseSpace, TDenseSpace> >
	class RveHomogenizerThermal : public RveHomogenizer<TSparseSpace, TDenseSpace, TReorderer>
	{

	public:

		KRATOS_CLASS_POINTER_DEFINITION( RveHomogenizerThermal );

		typedef RveHomogenizer < TSparseSpace, TDenseSpace, TReorderer > BaseType;


		typedef typename BaseType::RveLinearSystemOfEquationsType RveLinearSystemOfEquationsType;
		typedef typename BaseType::RveLinearSystemOfEquationsPointerType RveLinearSystemOfEquationsPointerType;
		typedef typename BaseType::SparseMatrixType SparseMatrixType;
		typedef typename BaseType::VectorType VectorType;
		typedef typename BaseType::DenseMatrixType  DenseMatrixType;
		typedef typename BaseType::SchemeType SchemeType;
		typedef typename BaseType::SchemePointerType SchemePointerType;
		typedef typename BaseType::DofsArrayType DofsArrayType;
		typedef typename BaseType::NodesArrayType NodesArrayType;
		typedef typename BaseType::ElementsArrayType ElementsArrayType;
		typedef typename BaseType::ConditionsArrayType ConditionsArrayType;
		typedef typename BaseType::RveConstraintHandlerType RveConstraintHandlerType;
		typedef typename BaseType::RveConstraintHandlerPointerType RveConstraintHandlerPointerType;

	public:
		
		RveHomogenizerThermal()
			: BaseType()
		{
		}

		virtual ~RveHomogenizerThermal()
		{
		}
		
		virtual void HomogenizeStressTensor(ModelPart& mp,
											const RveGeometryDescriptor& geomDescriptor,
											RveLinearSystemOfEquationsPointerType& soe,
											RveConstraintHandlerPointerType& constraintHandler,
											RveMacroscaleData& macroScaleData,
											Vector& S)
		{
			//ProcessInfo& processInfo = mp.GetProcessInfo();

			size_t ndim = geomDescriptor.Dimension(); // 2D or 3D
			if (S.size() != ndim) S.resize(ndim, false);
			noalias(S) = ZeroVector(ndim);

#if RVE_OPTIMIZATION == 0
			double totalVolume(0.0);

			std::vector<Vector> stressTensors;

			for (ModelPart::ElementIterator it = mp.ElementsBegin(); it != mp.ElementsEnd(); ++it)
			{
				Element& ielem = *it;
				Element::GeometryType& igeom = ielem.GetGeometry();
				Element::IntegrationMethod intmethod = ielem.GetIntegrationMethod();
				const Element::GeometryType::IntegrationPointsArrayType& ipts = igeom.IntegrationPoints(intmethod);

				ielem.GetValueOnIntegrationPoints(FLUX_RVE, stressTensors, processInfo);
				if (stressTensors.size() != ipts.size()) continue;

				for (size_t point_id = 0; point_id < ipts.size(); point_id++)
				{
					double dV = igeom.DeterminantOfJacobian(point_id, intmethod) * ipts[point_id].Weight();
					Vector& igpStressTensor = stressTensors[point_id];

					for (size_t i = 0; i < ndim; i++)
						S(i) += igpStressTensor(i) * dV;

					totalVolume += dV;
				}
			}

			if (totalVolume == 0.0)
				noalias(S) = ZeroVector(ndim);
			else
				S /= totalVolume;

#else // TODO: -> MODIFY FOR THE TEMPERATURE REACTION
			double totalVolume = geomDescriptor.DomainSize();
			//ModelPart::NodeType& ref_node = mp.GetNode(geomDescriptor.ReferenceNodeID()); // NOW ALWAYS (0,0)
			array_1d<double, 2> X;
			array_1d<double, 2> f;
			Matrix Sig(2, 2, 0.0);
			for (RveGeometryDescriptor::IndexContainerType::const_iterator it =
				geomDescriptor.BoundaryNodesIDs().begin(); it != geomDescriptor.BoundaryNodesIDs().end(); ++it)
			{
				RveGeometryDescriptor::IndexType index = *it;
				ModelPart::NodeType& bnd_node = mp.GetNode(index);
				X[0] = bnd_node.X0() /*- ref_node.X0()*/;
				X[1] = bnd_node.Y0() /*- ref_node.Y0()*/;
				ModelPart::NodeType::DofType& dof = bnd_node.GetDof(TEMPERATURE);
				f[0] = dof.GetSolutionStepReactionValue();
				f[1] = dof.GetSolutionStepReactionValue();
				Sig += outer_prod(f, X);
			}
			S(0) = Sig(0, 0);
			S(1) = Sig(1, 1);

			if (totalVolume == 0.0)
				noalias(S) = ZeroVector(ndim);
			else
				S /= totalVolume;
#endif
		}

		void SaveSolutionVector(Vector& U, ModelPart& mp)
		{
			size_t counter = 0;
			for (ModelPart::NodeIterator node_iter = mp.NodesBegin();
				node_iter != mp.NodesEnd();
				++node_iter)
			{
				ModelPart::NodeType& iNode = *node_iter;
				for (ModelPart::NodeType::DofsContainerType::iterator dof_iter = iNode.GetDofs().begin();
					dof_iter != iNode.GetDofs().end();
					++dof_iter)
				{
					ModelPart::DofType& iDof = *dof_iter;
					U(counter++) = iDof.GetSolutionStepValue();
				}
			}
		}

		void RestoreSolutionVector(const Vector& U, ModelPart& mp)
		{
			size_t counter = 0;
			for (ModelPart::NodeIterator node_iter = mp.NodesBegin();
				node_iter != mp.NodesEnd();
				++node_iter)
			{
				ModelPart::NodeType& iNode = *node_iter;
				for (ModelPart::NodeType::DofsContainerType::iterator dof_iter = iNode.GetDofs().begin();
					dof_iter != iNode.GetDofs().end();
					++dof_iter)
				{
					ModelPart::DofType& iDof = *dof_iter;
					iDof.GetSolutionStepValue() = U(counter++);
				}
			}
		}

		size_t CalculateTotalNumberOfDofs(ModelPart& mp)
		{
			size_t n(0);
			for (ModelPart::NodeIterator node_iter = mp.NodesBegin();
				node_iter != mp.NodesEnd();
				++node_iter)
			{
				ModelPart::NodeType& iNode = *node_iter;
				n += iNode.GetDofs().size();
			}
			return n;
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
			size_t strain_size = geomDescriptor.Dimension();

			Vector saved_strain_vector(strain_size);
			noalias(saved_strain_vector) = macroScaleData.StrainVector();
			Vector pert_stress_vector(strain_size);

			if (C.size1() != strain_size || C.size2() != strain_size)
				C.resize(strain_size, strain_size, false);
			noalias(C) = ZeroMatrix(strain_size, strain_size);

			// compute the perturbation parameters
			double perturbation = 1.0E-6*norm_1(macroScaleData.StrainVectorOld());
			if(perturbation < 1.0e-7)
				perturbation = 1.0e-7;

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
										soe->R(), soe->EquationSystemSize());
				if (move_mesh) this->MoveMesh(mp);

#if RVE_OPTIMIZATION == 1
				soe->BuildRHS(mp, pScheme);
#elif RVE_OPTIMIZATION == 2
				soe->BuildRHS_Reduced(mp, pScheme, geomDescriptor.BoundaryElementsIDs());
#endif
				constraintHandler->PostUpdate(mp,geomDescriptor,macroScaleData, 
				                        soe->TransformedEquationIds(), *pScheme,
				                        soe->DofSet(), soe->A(), soe->X(), soe->B(),
										soe->R(), soe->EquationSystemSize());
				this->HomogenizeStressTensor(mp, geomDescriptor, soe, constraintHandler, macroScaleData, pert_stress_vector);
				
				for (size_t i = 0; i < strain_size; i++)
					C(i, j) = (pert_stress_vector(i) - S(i)) / perturbation;
			}
			// reset stored strain vector
			noalias(macroScaleData.StrainVector()) = saved_strain_vector;
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

		/**
		* Updates the nodal coordinates if necessary
		*/
		virtual void MoveMesh(ModelPart& mp)
		{
		}

	};

} // namespace Kratos



#endif // RVE_HOMOGENIZER_THERMAL_H_INCLUDED