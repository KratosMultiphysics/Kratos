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


#if !defined(KRATOS_EBST_H_INCLUDED )
#define  KRATOS_EBST_H_INCLUDED



// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"


namespace Kratos {

    class Ebst
    : public Element {
    public:

        // Counted pointer of Ebst
        KRATOS_CLASS_POINTER_DEFINITION(Ebst);


        // Constructor void
        Ebst();

        // Constructor using an array of nodes
        Ebst(IndexType NewId, GeometryType::Pointer pGeometry);

        // Constructor using an array of nodes with properties
        Ebst(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

        // Destructor
        virtual ~Ebst();


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

        void Initialize();

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

        //		void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
        //				std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);



    protected:

      void CalculateEquivalentStresses(
          const array_1d<double, 3 > & membrane_strain,
          const array_1d<double, 3 > & bending_strain,
          boost::numeric::ublas::bounded_matrix<double, 3, 3 >& Dmat_m,
          boost::numeric::ublas::bounded_matrix<double, 3, 3 >& Dmat_f,
          array_1d<double, 3 > & membrane_stress,
          array_1d<double, 3 > & bending_stress,
          double& h_on_h0, //  ratio between current thickness and original thickness h/h0
          const ProcessInfo& rCurrentProcessInfo
          );
      
        //cartesian derivatives (reference configuration)
	array_1d<double,3> mvxe;
	array_1d<double,3> mvye;
        boost::numeric::ublas::bounded_matrix<double, 2,3 > mdcgM; //central element
        boost::numeric::ublas::bounded_matrix<double, 2,6 > mdcg1;
        boost::numeric::ublas::bounded_matrix<double, 2,6 > mdcg2;
        boost::numeric::ublas::bounded_matrix<double, 2,6 > mdcg3;

        //area in the reference configuration
        double Area0;
        array_1d<double,3> mK0;

        //vector of constitutive laws across the thickness
        std::vector<ConstitutiveLaw<Node<3> >::Pointer> mConstitutiveLawVector;


    private:
        ///@name Static Member Variables

        
        //double GetElementalPressure(
        //	const ProcessInfo& rCurrentProcessInfo);

        bool HasNeighbour(unsigned int index, const Node < 3 > & neighb);

        unsigned int NumberOfActiveNeighbours(WeakPointerVector< Node < 3 > >& neighbs);

        void CalculateAll(
                MatrixType& rLeftHandSideMatrix,
                VectorType& rRightHandSideVector,
                const ProcessInfo& rCurrentProcessInfo,
                bool CalculateStiffnessMatrixFlag,
                bool CalculateResidualVectorFlag);
		
	void CalculateCartesianDerOnGauss(
	  const double eta1, 
	  const double eta2, 
	  const boost::numeric::ublas::bounded_matrix<double, 6, 3 >& ms_coord,
	  const array_1d<double,3>& vxe,
	  const array_1d<double,3>& vye,
	  boost::numeric::ublas::bounded_matrix<double, 2,2 >& Jinv,
	  boost::numeric::ublas::bounded_matrix<double, 2,6 >& dcg
	  );
	  
	  void CalculatePhiG(
		boost::numeric::ublas::bounded_matrix<double, 2,3 >& phiG,
		const boost::numeric::ublas::bounded_matrix<double, 2, 6 >& dcgG,
		const boost::numeric::ublas::bounded_matrix<double, 6, 3 >& ms_coord
		);

	  void CalculatePhiM(
		boost::numeric::ublas::bounded_matrix<double, 2,3 >& phiM,
		const boost::numeric::ublas::bounded_matrix<double, 2, 3 >& dcgM,
		const boost::numeric::ublas::bounded_matrix<double, 6, 3 >& ms_coord
		);
		
	void CalculateAndAdd_MembraneStrain(
		array_1d<double,3>& strain,
		const boost::numeric::ublas::bounded_matrix<double, 2,3 >& phiG
		);

	void CalculateAndAdd_MembraneB(
		boost::numeric::ublas::bounded_matrix<double, 3, 18 >& B_m,
		const boost::numeric::ublas::bounded_matrix<double, 2, 6 >& dcgG,
		const boost::numeric::ublas::bounded_matrix<double, 2, 3 >& phiG
		);	   

	void Calculate_h_ab(
		array_1d<double, 3 >& h_ab,
		const unsigned int alpha, 
		const unsigned int beta,
		const boost::numeric::ublas::bounded_matrix<double, 2,3 >& phiG1,
		const boost::numeric::ublas::bounded_matrix<double, 2,3 >& phiG2,
		const boost::numeric::ublas::bounded_matrix<double, 2,3 >& phiG3,
		const boost::numeric::ublas::bounded_matrix<double, 2,3 >& dcgM
		);

        void CalculateAndAdd_Membrane_Kg(
		boost::numeric::ublas::bounded_matrix<double, 18, 18 >& K,
		const boost::numeric::ublas::bounded_matrix<double, 2, 6 >& dcgG1,
                const boost::numeric::ublas::bounded_matrix<double, 2, 6 >& dcgG2,
                const boost::numeric::ublas::bounded_matrix<double, 2, 6 >& dcgG3,
		const array_1d<double,3>& membrane_stress
		);

    }; // class KRATOS_EBST_H_INCLUDED.

} // namespace Kratos.

#endif // KRATOS_MEMBRANE_BEPPE_ELEMENT_H_INCLUDED  defined 
