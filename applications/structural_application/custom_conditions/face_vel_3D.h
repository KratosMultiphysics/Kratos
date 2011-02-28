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
*   Date:                $Date: 2007-03-13 15:01:51 $
*   Revision:            $Revision: 1.3 $
*
* ***********************************************************/


#if !defined(KRATOS_FACE_VEL_3D_H_INCLUDED )
#define  KRATOS_FACE_VEL_3D_H_INCLUDED



// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "custom_conditions/face3D.h"


namespace Kratos
{
class FaceVel3D
	: public Face3D
    {
	public:

		// Counted pointer of Face3D
		KRATOS_CLASS_POINTER_DEFINITION(FaceVel3D);  


		// Constructor void
		FaceVel3D();

		// Constructor using an array of nodes 
		FaceVel3D(IndexType NewId, GeometryType::Pointer pGeometry);

		// Constructor using an array of nodes with properties
		FaceVel3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

		// Destructor
		virtual ~FaceVel3D();

		
		// Name Operations

		Condition::Pointer Create(
			IndexType NewId,
			NodesArrayType const& ThisNodes,
			PropertiesType::Pointer pProperties) const;

		void EquationIdVector(
			EquationIdVectorType& rResult,
			ProcessInfo& rCurrentProcessInfo);

		void GetDofList(
			DofsVectorType& ElementalDofList,
			ProcessInfo& rCurrentProcessInfo);

	protected:


	private:
		///@name Static Member Variables 

	///@}
	///@name Serialization
	///@{	
	friend class Serializer;

	virtual void save(Serializer& rSerializer)
	{
	rSerializer.save("Name","FaceVel3D");
	KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Face3D );
	}

	virtual void load(Serializer& rSerializer)
	{
	KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer,  Face3D );
	}
		
		

	};	// class Face3D.

}	// namespace Kratos.

#endif // KRATOS_FACE3D_H_INCLUDED  defined 
