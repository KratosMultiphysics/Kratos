//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.3 $
//
// 



// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d.h"
#include "geometries/nurbs_2d.h"
#include "nurbs_testcase_application.h"
#include "includes/variables.h"
#include "includes/model_part.h"


namespace Kratos
{
//Initiallize Variables
    KRATOS_CREATE_VARIABLE(int, NURBS_ID);
//    KRATOS_CREATE_VARIABLE(double, DISTANCE);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(NURBS_COORDINATES);

    KratosNurbsTestcaseApplication::KratosNurbsTestcaseApplication():



        mNurbsPatchGeometry3D(PointerVector< Node<3> >(),Vector(), Vector(), Vector(),int(),int(),int(),int()),
        mNurbsPoisson2D( 0, Element::GeometryType::Pointer( new NurbsPatchGeometry2D<Node<3> >( ) ) )

    {
//        Node<3>::Pointer pNode = Node<3>::Pointer( new Node<3>(0) );
//        mNurbsPatchGeometry = NurbsPatchGeometry< Node<3> >(pNode,double(), Vector(), Vector(),int(),int(),int(),int());
    }
 	
 	void KratosNurbsTestcaseApplication::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();
 		std::cout << "Initializing KratosNurbsTestcaseApplication... " << std::endl;
        KRATOS_REGISTER_ELEMENT("NurbsPoisson2D", mNurbsPoisson2D);  //and here is our element


        KRATOS_REGISTER_VARIABLE(NURBS_ID);
//        KRATOS_REGISTER_VARIABLE(DISTANCE);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(NURBS_COORDINATES);

 	}

    

}  // namespace Kratos.


