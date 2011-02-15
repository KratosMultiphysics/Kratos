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


#if !defined(KRATOS_PFEM_CONTACT_ELEMENT3D_VEL_H_INCLUDED )
#define  KRATOS_PFEM_CONTACT_ELEMENT3D_VEL_H_INCLUDED



// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "custom_elements/pfem_contact_element3D.h"


namespace Kratos {

    class PfemContactElement3DVel
    : public PfemContactElement3D {
    public:

        // Counted pointer of Ebst
        KRATOS_CLASS_POINTER_DEFINITION(PfemContactElement3DVel);

        // Constructor void
   //     PfemContactElement3D();

        // Constructor using an array of nodes
        PfemContactElement3DVel(IndexType NewId, GeometryType::Pointer pGeometry);

        // Constructor using an array of nodes with properties
        PfemContactElement3DVel(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

        // Destructor
        virtual ~PfemContactElement3DVel();


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

       void CalculateLocalVelocityContribution(MatrixType& rDampMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo);

// 	void Calculate( const Variable<double>& rVariable, 
// 			      double& Output, 
// 			      const ProcessInfo& rCurrentProcessInfo);

	/*virtual void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);*/

        //		void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
        //				std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);
		virtual std::string Info() const
		{
			return "PfemContactElement3DVel #" ;
		}

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << Info() << Id();
	}


    protected:
		


    private:

        
    }; // class KRATOS_EBST_H_INCLUDED.

} // namespace Kratos.

#endif // KRATOS_MEMBRANE_BEPPE_ELEMENT_H_INCLUDED  defined 
