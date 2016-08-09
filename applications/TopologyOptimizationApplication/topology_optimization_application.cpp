// ==============================================================================
/*
 KratosTopologyOptimizationApplication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 (Released on march 05, 2007).

 Copyright (c) 2016: Daniel Baumgaertner
                     daniel.baumgaertner@tum.de
                     Chair of Structural Analysis
                     Technische Universitaet Muenchen
                     Arcisstrasse 21 80333 Munich, Germany

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
*/
//==============================================================================
//
//   Project Name:        KratosTopology                        $
//   Last modified by:	  $Author:   daniel.baumgaertner@tum.de $
// 						  $Co-Author: Octaviano Malfavón Farías $
//   Date:                $Date:                    August 2016 $
//   Revision:            $Revision:                        0.0 $
//
// ==============================================================================

// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/line_2d.h"
#include "includes/variables.h"

#include "topology_optimization_application.h"


// Geometries that must be added when more elements are added into the application (SOLID MECHANICS APPLICATION)
//#include "geometries/tetrahedra_3d_10.h"
//#include "geometries/hexahedra_3d_20.h"
//#include "geometries/hexahedra_3d_27.h"
//#include "geometries/prism_3d_6.h"
//#include "geometries/prism_3d_15.h"


namespace Kratos
{
    //Create Variables with Python connection
    KRATOS_CREATE_VARIABLE( double, E_MIN )
    KRATOS_CREATE_VARIABLE( double, E_0 )
    KRATOS_CREATE_VARIABLE( double, PENAL )
    KRATOS_CREATE_VARIABLE( double, X_PHYS )
    KRATOS_CREATE_VARIABLE( double, X_PHYS_OLD )
    KRATOS_CREATE_VARIABLE( double, DCDX )
    KRATOS_CREATE_VARIABLE( double, DVDX )
    KRATOS_CREATE_VARIABLE( double, SOLID_VOID )
    KRATOS_CREATE_VARIABLE( double, LOCAL_STRAIN_ENERGY )


    KratosTopologyOptimizationApplication::KratosTopologyOptimizationApplication():
		mSmallDisplacementSIMPElement3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ), // dummy element for surface representation
        mSmallDisplacementSIMPElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
        mSmallDisplacementSIMPElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) )

//        Extra elements that can be added in the future
//        mSmallDisplacementSIMPElement3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) ),
//        mSmallDisplacementSIMPElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10 ) ) ) ),
//        mSmallDisplacementSIMPElement3D15N( 0, Element::GeometryType::Pointer( new Prism3D15 <Node<3> >( Element::GeometryType::PointsArrayType( 15 ) ) ) ),
//        mSmallDisplacementSIMPElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20 ) ) ) ),
//        mSmallDisplacementSIMPElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType( 27 ) ) ) )

 	{}
 	
 	void KratosTopologyOptimizationApplication::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();

        std::cout << "     KRATOS|_   _/_ \\| _ \\ _ \\| |  / _ \\/ __\\ \\ / /         " << std::endl;
        std::cout << "             | | (_) |  _/(_) | |_| (_) |(_ |\\ V /               " << std::endl;
        std::cout << "             |_|\\___/|_| \\___/|____\\___/ ___| |_|OPTIMIZATION  " << std::endl;
        std::cout << "Initializing KratosTopologyOptimizationApplication...    " << std::endl;

        //Register small displacement elements
        KRATOS_REGISTER_ELEMENT( "SmallDisplacementSIMPElement3D3N", mSmallDisplacementSIMPElement3D3N ) // dummy element for surface representation
        KRATOS_REGISTER_ELEMENT( "SmallDisplacementSIMPElement3D4N", mSmallDisplacementSIMPElement3D4N )
        KRATOS_REGISTER_ELEMENT( "SmallDisplacementSIMPElement3D8N", mSmallDisplacementSIMPElement3D8N )

//        Extra elements that can be added in the future
//        KRATOS_REGISTER_ELEMENT( "SmallDisplacementSIMPElement3D6N", mSmallDisplacementSIMPElement3D6N )
//        KRATOS_REGISTER_ELEMENT( "SmallDisplacementSIMPElement3D10N", mSmallDisplacementSIMPElement3D10N )
//        KRATOS_REGISTER_ELEMENT( "SmallDisplacementSIMPElement3D15N", mSmallDisplacementSIMPElement3D15N )
//        KRATOS_REGISTER_ELEMENT( "SmallDisplacementSIMPElement3D20N", mSmallDisplacementSIMPElement3D20N )
//        KRATOS_REGISTER_ELEMENT( "SmallDisplacementSIMPElement3D27N", mSmallDisplacementSIMPElement3D27N )

        //Register Variables with Python connection
        KRATOS_REGISTER_VARIABLE( E_MIN )
        KRATOS_REGISTER_VARIABLE( E_0 )
        KRATOS_REGISTER_VARIABLE( PENAL )
        KRATOS_REGISTER_VARIABLE( X_PHYS )
        KRATOS_REGISTER_VARIABLE( X_PHYS_OLD )
        KRATOS_REGISTER_VARIABLE( DCDX )
        KRATOS_REGISTER_VARIABLE( DVDX )
        KRATOS_REGISTER_VARIABLE( SOLID_VOID )
        KRATOS_REGISTER_VARIABLE( LOCAL_STRAIN_ENERGY )

 
 	}

}  // namespace Kratos.


