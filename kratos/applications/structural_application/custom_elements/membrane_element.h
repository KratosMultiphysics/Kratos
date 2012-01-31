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
*   Last Modified by:    $Author: rrossi $
*   Date:                $Date: 2008-10-13 07:00:53 $
*   Revision:            $Revision: 1.4 $
*
* ***********************************************************/


#if !defined(KRATOS_MEMBRANE_BEPPE_ELEMENT_H_INCLUDED )
#define  KRATOS_MEMBRANE_BEPPE_ELEMENT_H_INCLUDED



// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"


namespace Kratos
{
class MembraneElement
	: public Element
    {
	public:

		// Counted pointer of MembraneElement
		KRATOS_CLASS_POINTER_DEFINITION(MembraneElement);  

		// Constructor using an array of nodes 
		MembraneElement(IndexType NewId, GeometryType::Pointer pGeometry);

		// Constructor using an array of nodes with properties
		MembraneElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

		// Destructor
		virtual ~MembraneElement();

		
		// Name Operations

		Element::Pointer Create(
			IndexType NewId,
			NodesArrayType const& ThisNodes,
			PropertiesType::Pointer pProperties) const;

		void EquationIdVector(
			EquationIdVectorType& rResult,
			ProcessInfo& rCurrentProcessInfo);

		void GetDofList(
			DofsVectorType& ElementalDofList,
			ProcessInfo& rCurrentProcessInfo);

		void Initialize ();

		void CalculateRightHandSide(
			VectorType& rRightHandSideVector,
			ProcessInfo& rCurrentProcessInfo);

		void CalculateLocalSystem(
			MatrixType& rLeftHandSideMatrix,
			VectorType& rRightHandSideVector,
			ProcessInfo& rCurrentProcessInfo);

		void CalculateOnIntegrationPoints(
			const Variable<Matrix>& rVariable,
			std::vector<Matrix>& Output,
			const ProcessInfo& rCurrentProcessInfo);

		void MassMatrix(
			MatrixType& rMassMatrix,
			ProcessInfo& rCurrentProcessInfo);

		void DampMatrix(
			MatrixType& rDampMatrix,
			ProcessInfo& rCurrentProcessInfo);

		void FinalizeSolutionStep(
			ProcessInfo& rCurrentProcessInfo);

		void GetValuesVector(
			Vector& values,
			int Step = 0);

		void GetFirstDerivativesVector(
			Vector& values,
			int Step = 0);

		void GetSecondDerivativesVector(
			Vector& values,
			int Step = 0);
	
		void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
				std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);
	//	///////////serve qui? non l'ho visto nel membrane element$$$$$$$$$$$$$$$$$$$$$ 
//		ConstitutiveLaws::Pointer GetConstitutiveLaw(int GaussPointNumber)
//			{return mConstitutiveLawVector[GaussPointNumber];}



	protected:


	private:
		///@name Static Member Variables 

// 		static Matrix                                             msB;
// 		static boost::numeric::ublas::bounded_matrix<double,3,3>  msQ;
// 		static Matrix											  msD;
// 		static Vector				                  msStrainVector;
// 		static Vector				                  msStressVector;
// 		static boost::numeric::ublas::bounded_matrix<double,2,2>  msC;
// 		static Matrix											  msDN_DX;

		/// privat variables

		std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;
		Geometry< Point<3,double> >::Pointer  mpReferenceGeometry; 

		Vector mDetJ0;
		
		double mTotalDomainInitialSize;
		double mdensity;
		double mThickness0; //thickness in the reference configuration


		Vector mThickness;									//container of thickness
		std::vector< array_1d<double,3> > mStrainsVector;	//container of Strain
		std::vector< array_1d<double,6> > mStressesVector;	//container of Stress
		std::vector< array_1d<double,6> > mCauchyStressesVector;	//container of Stress
//		std::vector<Vector> mCauchyStressesVector;			//container of CauchyVector


		std::vector< array_1d<double,3> >  mV1;
		std::vector< array_1d<double,3> >  mV2;
		std::vector< Matrix >              mG_Vector;



		// privat name Operations

		//double GetElementalPressure(
		//	const ProcessInfo& rCurrentProcessInfo);

		void CalculateAll(
			MatrixType& rLeftHandSideMatrix,
			VectorType& rRightHandSideVector,
			const ProcessInfo& rCurrentProcessInfo,
			bool CalculateStiffnessMatrixFlag,
			bool CalculateResidualVectorFlag);

		void CalculateAndAddKm(
			Matrix& K,
			Matrix& msB,
			Matrix& msD,
			double weight);

		void CalculateAndAddKg(
			Matrix& K,
			boost::numeric::ublas::bounded_matrix<double,3,3>& msQ,
			const Matrix& DN_De,
			Vector& msStressVector,
			double weight);

		void CalculateAndSubKp(
			Matrix& K,
			array_1d<double,3>& ge,
			array_1d<double,3>& gn,
			const Matrix& DN_De,
			const Vector& N,
			double pressure,
			double weight);

		void MakeCrossMatrix(
			boost::numeric::ublas::bounded_matrix<double,3,3>& M,
			array_1d<double,3>& U);

		void CrossProduct(
			array_1d<double,3>& cross,
			array_1d<double,3>& a,
			array_1d<double,3>& b);

		void SubtractMatrix(
			MatrixType& Destination,
			boost::numeric::ublas::bounded_matrix<double,3,3>& InputMatrix, 
			int InitialRow,
			int InitialCol);

		void ExpandReducedMatrix(
			Matrix& Destination,
			Matrix& ReducedMatrix);

		void CalculateQ(
			boost::numeric::ublas::bounded_matrix<double,3,3>& msQ,
			Matrix& msG);

		void CalculateB(
			Matrix& msB,
			boost::numeric::ublas::bounded_matrix<double,3,3>& msQ,
			const Matrix& DN_De,
			array_1d<double,3>& ge,
			array_1d<double,3>& gn);

		void CalculateJ(
			boost::numeric::ublas::bounded_matrix<double,2,2>& j,
			array_1d<double,3>& ge,
			array_1d<double,3>& gn,
			array_1d<double,3>& v3);

		void CalculateStrain(
			Vector& StrainVector,
			boost::numeric::ublas::bounded_matrix<double,2,2>& C);

		void CalculateAndAdd_BodyForce(
			const Vector& N,
			const ProcessInfo& rCurrentProcessInfo,
			array_1d<double,3>& BodyForce,
			VectorType& rRightHandSideVector,
			double weight);

		void CalculateAndAdd_PressureForce(
			VectorType& residualvector,
			const Vector& N,
			const array_1d<double,3>& v3,
			double pressure,
			double weight,
			const ProcessInfo& rCurrentProcessInfo);

		// this function transforms the local stress (with 3 components)
		// to the global one (with 6 components)
		void Calculate_GlobalStressVector(
			array_1d<double,6>& GlobalVector,
			Vector& LocalStressVector,
			array_1d<double,3>& v1,
			array_1d<double,3>& v2);

		//auxiliary function needed in the calculation of output stresses
		inline array_1d<double,6> VoigtTensorComponents(
			array_1d<double,3>& a,
			array_1d<double,3>& b);

	    ///@} 
	    ///@name Serialization
	    ///@{ 
    
	    friend class Serializer; 

	    // A private default constructor necessary for serialization 
	    MembraneElement(){}

	    virtual void save(Serializer& rSerializer) const
	    {
	       KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
	    }

	    virtual void load(Serializer& rSerializer)
	    {
	       KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer,  Element );
	    }

	    ///@} 

	};	// class MembraneElement.

}	// namespace Kratos.

#endif // KRATOS_MEMBRANE_BEPPE_ELEMENT_H_INCLUDED  defined 
