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

#if !defined(RVE_CONSTRAINT_HANDLER_H_INCLUDED)
#define RVE_CONSTRAINT_HANDLER_H_INCLUDED

#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/define.h"
#include "includes/serializer.h"
#include "solving_strategies/schemes/scheme.h"
#include "rve_geometry_descriptor.h"

namespace Kratos
{

	template<class TSparseSpace, 
			 class TDenseSpace>
	class RveConstraintHandler
	{

	public:

		KRATOS_CLASS_POINTER_DEFINITION( RveConstraintHandler );
		typedef typename TSparseSpace::MatrixType SparseMatrixType;
		typedef typename TSparseSpace::VectorType VectorType;
		typedef typename TDenseSpace::MatrixType  DenseMatrixType;
		typedef Scheme<TSparseSpace,TDenseSpace> SchemeType;
		typedef typename SchemeType::Pointer SchemePointerType;
		typedef ModelPart::DofsArrayType DofsArrayType;
		typedef ModelPart::NodesContainerType NodesArrayType;
		typedef ModelPart::ElementsContainerType ElementsArrayType;
		typedef ModelPart::ConditionsContainerType ConditionsArrayType;
		typedef Dof<double>::IndexType IndexType;
		typedef std::vector< IndexType > IndexContainerType;
		typedef std::vector< double > DoubleContainerType;

	public:
		
		RveConstraintHandler()
		{
		}

		virtual ~RveConstraintHandler()
		{
		}

	public:

		virtual void SetupDofSet(ModelPart& mp, 
								 const RveGeometryDescriptor& geom,
								 SchemePointerType& pScheme, 
								 DofsArrayType& dofset)
		{
			ElementsArrayType& pElements = mp.Elements();
			Element::DofsVectorType ElementalDofList;
			ProcessInfo& CurrentProcessInfo = mp.GetProcessInfo();

			DofsArrayType temp;

			for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
			{
				pScheme->GetElementalDofList(*it, ElementalDofList, CurrentProcessInfo);
				for (typename Element::DofsVectorType::iterator i = ElementalDofList.begin(); i != ElementalDofList.end(); ++i)
					temp.push_back( i->get() );
			}

			ConditionsArrayType& pConditions = mp.Conditions();
			for (typename ConditionsArrayType::ptr_iterator it = pConditions.ptr_begin(); it != pConditions.ptr_end(); ++it)
			{
				pScheme->GetConditionDofList(*it, ElementalDofList, CurrentProcessInfo);
				for (typename Element::DofsVectorType::iterator i = ElementalDofList.begin(); i != ElementalDofList.end(); ++i)
					temp.push_back( i->get() );
			}

			temp.Unique();
			dofset = temp;
		}

		virtual void SetUpSystem(ModelPart& mp, 
								 const RveGeometryDescriptor& geom,
								 DofsArrayType& dofset, 
								 size_t& equation_system_size,
								 IndexContainerType& transformed_equation_ids,
								 IndexContainerType& equation_id_flag)
		{
			int free_id = 0;
			int fix_id = dofset.size();

			transformed_equation_ids.resize(fix_id);
			equation_id_flag.resize(fix_id);

			for (typename DofsArrayType::iterator dof_iterator = dofset.begin(); dof_iterator != dofset.end(); ++dof_iterator)
			{
				if (dof_iterator->IsFixed())
				{
					size_t index = --fix_id;
					dof_iterator->SetEquationId(index);
					equation_id_flag[index] = 0;
					transformed_equation_ids[index] = index;
				}
				else
				{
					size_t index = free_id++;
					dof_iterator->SetEquationId(index);
					equation_id_flag[index] = 0;
					transformed_equation_ids[index] = index;
				}
			}

			equation_system_size = fix_id;
		}

		virtual void AddConditions(ModelPart& mp, 
								   const RveGeometryDescriptor& geom)
		{
			/**
			this method can be used to add conditions to the rve modelpart, for example
			if one wants to apply constraints by means of lagrange multipliers.
			*/
			/*for each elem in rve.Edges()
				rve.addcondition(minimalcondition2D(...))*/
		}

		virtual void ApplyMacroScaleData(ModelPart& mp, 
										 const RveGeometryDescriptor& geom,
										 const RveMacroscaleData& macroScaleData)
		{
			// todo:
			// for each element
			//    element->SetValueOnGaussP(INITIAL_STRAIN, -macroScaleData.StrainVector())
		}
		
		virtual void FinalizeSolutionStep(ModelPart& mp, 
										 const RveGeometryDescriptor& geom,
										 const RveMacroscaleData& macroScaleData)
		{
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
			// standard: call the scheme and do nothing else!
			scheme.Update(mp, dofset, A, Dx, b);
		}

	};

} // namespace Kratos



#endif // RVE_CONSTRAINT_HANDLER_H_INCLUDED