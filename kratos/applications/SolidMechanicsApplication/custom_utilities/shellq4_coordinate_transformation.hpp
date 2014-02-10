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
//   Date:                $Date: 2013-10-03 12:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#if !defined(KRATOS_SHELLQ4_COORDINATE_TRANSFORMATION_H_INCLUDED )
#define  KRATOS_SHELLQ4_COORDINATE_TRANSFORMATION_H_INCLUDED

#include "shellq4_local_coordinate_system.hpp"

namespace Kratos
{

	/** \brief ShellQ4_CoordinateTransformation
	*
	* This class represents a basic (linear) coordinate transformation that can be used
	* by any element whose geometry is a QUAD 4 in 3D space, with 6 D.O.F.s per node.
	* It's main aim is to:
	* 1) Create the local coordinate system
	* 2) Transform the incoming global displacements in local coordinate system
	* 3) Transform the outgoing matrices and vectors in global coordinate system
	*/
	class ShellQ4_CoordinateTransformation
	{

	public:

		KRATOS_CLASS_POINTER_DEFINITION( ShellQ4_CoordinateTransformation );

		typedef Element::GeometryType GeometryType;

		typedef Vector VectorType;

		typedef Matrix MatrixType;

	public:

		ShellQ4_CoordinateTransformation(const GeometryType::Pointer & pGeometry)
			: mpGeometry(pGeometry)
		{
		}

		virtual ~ShellQ4_CoordinateTransformation()
		{
		}

	private:

		ShellQ4_CoordinateTransformation(const ShellQ4_CoordinateTransformation & other);

		ShellQ4_CoordinateTransformation & operator = (const ShellQ4_CoordinateTransformation & other);

	public:

		virtual ShellQ4_CoordinateTransformation::Pointer Create(GeometryType::Pointer pGeometry)const
		{
			return ShellQ4_CoordinateTransformation::Pointer( new ShellQ4_CoordinateTransformation( pGeometry ) );
		}

		virtual void Initialize()
		{
		}

		virtual void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
		{
		}
		
		virtual void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
		{
		}
		
		virtual void InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
		{
		}
		
		virtual void FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
		{
		}

		virtual ShellQ4_LocalCoordinateSystem CreateReferenceCoordinateSystem()const
		{
			const GeometryType & geom = GetGeometry();
			return ShellQ4_LocalCoordinateSystem(geom[0].GetInitialPosition(),
											     geom[1].GetInitialPosition(),
											     geom[2].GetInitialPosition(),
											     geom[3].GetInitialPosition());
		}

		virtual ShellQ4_LocalCoordinateSystem CreateLocalCoordinateSystem()const
		{
			return CreateReferenceCoordinateSystem();
		}

		virtual Vector CalculateLocalDisplacements(const ShellQ4_LocalCoordinateSystem & LCS, 
												   const VectorType & globalDisplacements)
		{
			MatrixType R(24, 24);
			LCS.ComputeTotalRotationMatrix( R );
			if(LCS.IsWarped()) {
				MatrixType W(24, 24);
				LCS.ComputeTotalWarpageMatrix( W );
				R = prod( W, R );
			}
			return prod( R, globalDisplacements );
		}

		virtual void FinalizeCalculations(const ShellQ4_LocalCoordinateSystem & LCS,
										  const VectorType & globalDisplacements,
										  const VectorType & localDisplacements,
										  MatrixType & rLeftHandSideMatrix,
										  VectorType & rRightHandSideVector,
										  const bool RHSrequired,
										  const bool LHSrequired)
		{
			MatrixType R(24, 24);
			LCS.ComputeTotalRotationMatrix( R );
			if(LCS.IsWarped()) {
				MatrixType W(24, 24);
				LCS.ComputeTotalWarpageMatrix( W );
				R = prod( W, R );
			}

			if(LHSrequired) {
				MatrixType temp(24, 24);
				noalias( temp ) = prod( trans( R ), rLeftHandSideMatrix );
				noalias( rLeftHandSideMatrix ) = prod( temp, R );
			}

			if(RHSrequired) {
				rRightHandSideVector = prod( trans( R ), rRightHandSideVector );
			}
		}

		virtual MatrixType GetNodalDeformationalRotationTensor(const ShellQ4_LocalCoordinateSystem & LCS,
			                                                   const Vector& globalDisplacements,
															   size_t nodeid)
		{
			return IdentityMatrix(3,3);
		}

		virtual MatrixType GetNodalDeformationalRotationTensor(const ShellQ4_LocalCoordinateSystem & LCS,
			                                                   const Vector& globalDisplacements,
															   const Vector& N)
		{
			return IdentityMatrix(3,3);
		}

	public:

		inline const GeometryType & GetGeometry()const { return *mpGeometry; }

	protected:

		ShellQ4_CoordinateTransformation()
			: mpGeometry(NULL)
		{
		}

	private:

		GeometryType::Pointer mpGeometry;

	private:

		friend class Serializer;

		virtual void save(Serializer& rSerializer) const
		{
			rSerializer.save("pGeom", mpGeometry);
		}

		virtual void load(Serializer& rSerializer)
		{
			rSerializer.load("pGeom", mpGeometry);
		}

	};

}


#endif // KRATOS_SHELLQ4_COORDINATE_TRANSFORMATION_H_INCLUDED