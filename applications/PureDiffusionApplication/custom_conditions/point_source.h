//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Masó Sotomayor
//

#if !defined(KRATOS_POINT_SOURCE_CONDITION_H_INCLUDED)
#define KRATOS_POINT_SOURCE_CONDITION_H_INCLUDED

// System includes 

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h" 

namespace Kratos
{

	class PointSource : public Condition
	{
	public:
		// Counted pointer of PointSource
		KRATOS_CLASS_POINTER_DEFINITION(PointSource);

		// Default constructor
		PointSource(IndexType NewId, GeometryType::Pointer pGeometry);
		PointSource(IndexType NewId, GeometryType::Pointer pGeometry, 
		            PropertiesType::Pointer pProperties);

		// Destructor
		virtual ~PointSource();

		Condition::Pointer Create(IndexType NewId, 
		                          NodesArrayType const& ThisNodes, 
		                          PropertiesType::Pointer pProperties) const;

		void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, 
		                          VectorType& rRightHandSideVector, 
		                          ProcessInfo& rCurrentProcessInfo);

		void CalculateRightHandSide(VectorType& rRightHandSideVector, 
		                            ProcessInfo& rCurrentProcessInfo);

		void EquationIdVector(EquationIdVectorType& rResult, 
		                      ProcessInfo& rCurrentProcessInfo);

		void GetDofList(DofsVectorType& rConditionalDofList, 
		                ProcessInfo& rCurrentProcessInfo);

	protected:


	private:

		friend class Serializer;

		// A private default constructor necessary for serialization
		PointSource(){};

		virtual void save(Serializer& rSerializer) const
		{
			KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
		}

		virtual void load(Serializer& rSerializer)
		{
			KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
		}

	};  // class PointSource

}  // namespace Kratos

#endif  // KRATOS_POINT_SOURCE_CONDITION_H_INCLUDED
