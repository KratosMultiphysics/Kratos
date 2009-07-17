/*
==============================================================================
KratosStructuralApplication 
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
/* *********************************************************   
*          
*   Last Modified by:    $Author: janosch $
*   Date:                $Date: 2008-01-25 08:37:31 $
*   Revision:            $Revision: 1.10 $
*
* ***********************************************************/

#if !defined(KRATOS_ORTHOTROPIC_2D_H_INCLUDED )
#define  KRATOS_ORTHOTROPIC_2D_H_INCLUDED

// System includes 

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "fluency_criteria/fluency_criteria.h"


namespace Kratos
{
	class Plane_Stress_Damage_Orthotropic_2D : public ConstitutiveLaw<Node<3> >
	{
		public:
		/**
		* Type Definitions
		*/
		typedef ConstitutiveLaw<Node<3> > BaseType;
                typedef FluencyCriteria FluencyCriteriaType; 
		/**
		* Counted pointer of Orthotropic3D
		*/
		typedef boost::shared_ptr<Plane_Stress_Damage_Orthotropic_2D> Pointer;

		/**
		* Life Cycle 
		*/
		/**
		* Default constructor.
		*/
                Plane_Stress_Damage_Orthotropic_2D();
		Plane_Stress_Damage_Orthotropic_2D(FluencyCriteriaType& FluencyCriteria);

		virtual boost::shared_ptr<ConstitutiveLaw<Node<3> > > Clone() const
		{
		boost::shared_ptr<ConstitutiveLaw<Node<3> > > p_clone(new Plane_Stress_Damage_Orthotropic_2D());
		return p_clone;
		}

		/**
		* Destructor.
		*/
		virtual ~Plane_Stress_Damage_Orthotropic_2D();

		/**
		* Operators 
		*/
		/**
		* Operations
		*/
         
		bool Has( const Variable<double>& rThisVariable );
		bool Has( const Variable<Vector>& rThisVariable );
		bool Has( const Variable<Matrix>& rThisVariable );


		double GetValue( const Variable<double>& rThisVariable );
		Vector GetValue( const Variable<Vector>& rThisVariable );
		Matrix GetValue( const Variable<Matrix>& rThisVariable );

		void SetValue( const Variable<double>& rThisVariable, const double rValue, 
		const ProcessInfo& rCurrentProcessInfo );

		void SetValue( const Variable<array_1d<double, 3> >& rThisVariable, 
		const array_1d<double, 3>& rValue, const ProcessInfo& rCurrentProcessInfo );

		void SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
		const ProcessInfo& rCurrentProcessInfo );

		void SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
		const ProcessInfo& rCurrentProcessInfo );


		void InitializeMaterial( const Properties& props,
		const GeometryType& geom,
		const Vector& ShapeFunctionsValues );


		void CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult);

		void CalculateStress( const Vector& StrainVector, Vector& StressVector);
		
		void InitializeSolutionStep( const Properties& props,
		const GeometryType& geom, //this is just to give the array of nodes
		const Vector& ShapeFunctionsValues ,
		const ProcessInfo& CurrentProcessInfo);

		void UpdateMaterial( const Vector& StrainVector,
		const Properties& props,
		const GeometryType& geom, //this is just to give the array of nodes
		const Vector& ShapeFunctionsValues ,
		const ProcessInfo& CurrentProcessInfo);

		void FinalizeSolutionStep( const Properties& props,
		const GeometryType& geom, //this is just to give the array of nodes
		const Vector& ShapeFunctionsValues ,
		const ProcessInfo& CurrentProcessInfo);

		void CalculateCauchyStresses( Vector& Cauchy_StressVector,
		const Matrix& F,
		const Vector& PK2_StressVector,
		const Vector& GreenLagrangeStrainVector);

		void CalculateStressAndTangentMatrix(Vector& StressVector,
                    const Vector& StrainVector,
                    Matrix& algorithmicTangent);


                protected:

		void CalculateElasticMatrix(Matrix& C, const Vector& E, const Vector& NU, const double&  Orthotropic_Angle);
  
		double mIsotropic_Elastic_Limit;
                FluencyCriteriaType mFluencyCriteria;
		
		Vector mInSituStress;
		Matrix mCtangent;
		Vector mCurrentStress;
                Vector mOrtotropic_Elastic_Limit;

		void CalculateTransformationMatrix(Matrix& T, const double& Orthotropic_Angle);


		}; // Class Plane_Stress_Damage_Orthotropic_2D 
}  // namespace Kratos.
#endif 



