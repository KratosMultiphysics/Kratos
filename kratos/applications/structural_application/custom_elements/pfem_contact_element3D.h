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


#if !defined(KRATOS_PFEM_CONTACT_ELEMENT3D_H_INCLUDED )
#define  KRATOS_PFEM_CONTACT_ELEMENT3D_H_INCLUDED



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


namespace Kratos {

    class PfemContactElement3D
    : public Element {
    public:

      
      
      
        // Counted pointer of Ebst
        KRATOS_CLASS_POINTER_DEFINITION(PfemContactElement3D);

        // Constructor void
	 PfemContactElement3D(){}

        // Constructor using an array of nodes
        PfemContactElement3D(IndexType NewId, GeometryType::Pointer pGeometry);

        // Constructor using an array of nodes with properties
        PfemContactElement3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

        // Destructor
        virtual ~PfemContactElement3D();


        // Name Operations
        Element::Pointer Create(
                IndexType NewId,
                NodesArrayType const& ThisNodes,
                PropertiesType::Pointer pProperties) const;

        virtual void EquationIdVector(
                EquationIdVectorType& rResult,
                ProcessInfo& rCurrentProcessInfo);

        virtual void GetDofList(
                DofsVectorType& ElementalDofList,
                ProcessInfo& rCurrentProcessInfo);

        virtual void Initialize();

        virtual void CalculateRightHandSide(
                VectorType& rRightHandSideVector,
                ProcessInfo& rCurrentProcessInfo);

        virtual void CalculateLocalSystem(
                MatrixType& rLeftHandSideMatrix,
                VectorType& rRightHandSideVector,
                ProcessInfo& rCurrentProcessInfo);

       /* virtual void CalculateOnIntegrationPoints(
                const Variable<Matrix>& rVariable,
                std::vector<Matrix>& Output,
                const ProcessInfo& rCurrentProcessInfo);*/

        virtual void MassMatrix(
                MatrixType& rMassMatrix,
                ProcessInfo& rCurrentProcessInfo);

        virtual void DampMatrix(
                MatrixType& rDampMatrix,
                ProcessInfo& rCurrentProcessInfo);
        
        virtual void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

        virtual void FinalizeSolutionStep(
                ProcessInfo& rCurrentProcessInfo);

        virtual void GetValuesVector(
                Vector& values,
                int Step = 0);

        virtual void GetFirstDerivativesVector(
                Vector& values,
                int Step = 0);

        virtual void GetSecondDerivativesVector(
                Vector& values,
                int Step = 0);
        virtual void Calculate( const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo);


	            virtual void InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo);
	/*virtual void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);*/

        //		void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
        //				std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);
		virtual std::string Info() const
		{
			return "PfemContactElement3D #" ;
		}

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << Info() << Id();
	}


    protected:
		
            virtual void CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
                              ProcessInfo& rCurrentProcessInfo,
                              bool CalculateStiffnessMatrixFlag,
                              bool CalculateResidualVectorFlag);


	//auxiliary function needed in the calculation of output stresses
	inline array_1d<double,6> VoigtTensorComponents(
		array_1d<double,3>& a,
		array_1d<double,3>& b);


	  bool mcontact_is_active;
	  double mpenetration;

          bool Check_image_inside_face(const array_1d<double,3> n, int reference_face, const Geometry< Node<3> >& geom,bool image_inside);
          bool same_side(const array_1d<double,3> p0, const array_1d<double,3> p1,const array_1d<double,3> a,const array_1d<double,3> b);
	  void CalculateOldIterationContactHeight(double& H_zero, const int reference_face);
	  void CalculateMinDistanceAndNormal(double& h,array_1d<double,3>& n,const boost::numeric::ublas::bounded_matrix<double, 4, 3 > ordered_points);
          void DetectContact(Geometry< Node<3> >& geom, boost::numeric::ublas::bounded_matrix<double, 4, 3 > DN_DX, const unsigned int single_node_index, array_1d<double,3>& n, double& h);
          void CheckIsContactMaster(int& flag);
	  void FlagVariableCheckForNonSuitableElements(double& accepted);
     private:
       
       
	    ///@}
	    ///@name Serialization
	    ///@{	
	    friend class Serializer; 
	    

	    virtual void save(Serializer& rSerializer) const
	    {
	    rSerializer.save("Name","PfemContactElement3D");
	    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
	    }

	    virtual void load(Serializer& rSerializer)
	    {
	    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer,  Element );
	    }
       
			
    }; // class KRATOS_EBST_H_INCLUDED.

} // namespace Kratos.

#endif // KRATOS_MEMBRANE_BEPPE_ELEMENT_H_INCLUDED  defined 
