// ==============================================================================
/*
 KratosALEApllication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 (Released on march 05, 2007).

 Copyright (c) 2016: Pooyan Dadvand, Riccardo Rossi, Andreas Winterstein
                     pooyan@cimne.upc.edu
                     rrossi@cimne.upc.edu
                     a.winterstein@tum.de
- CIMNE (International Center for Numerical Methods in Engineering),
  Gran Capita' s/n, 08034 Barcelona, Spain
- Chair of Structural Analysis, Technical University of Munich
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

/* ****************************************************************************
 *  Projectname:         $KratosALEApplication
 *  Last Modified by:    $Author: A.Winterstein@tum.de $
 *  Date:                $Date: June 2016 $
 *  Revision:            $Revision: 1.5 $
 * ***************************************************************************/



// System includes


// External includes


// Project includes
//#include "includes/define.h"
#include "ale_application.h"

#include "geometries/triangle_2d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/prism_3d_15.h"
#include "geometries/prism_3d_6.h"


namespace Kratos
{
//

KratosALEApplication::KratosALEApplication():
    mLaplacianMeshMovingElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mLaplacianMeshMovingElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mLaplacianMeshMovingElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mLaplacianMeshMovingElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) ),
    mStructuralMeshMovingElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mStructuralMeshMovingElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mStructuralMeshMovingElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mStructuralMeshMovingElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) ),
    mStructuralMeshMovingElement3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) ),
    mStructuralMeshMovingElement3D15N( 0, Element::GeometryType::Pointer( new Prism3D15 <Node<3> >( Element::GeometryType::PointsArrayType( 15 ) ) ) )
{}


void KratosALEApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "KRATOS    ___   __   ____             " << std::endl;
    std::cout << "         / _ | / /  / __/             " << std::endl;
    std::cout << "        / __ |/ /__/ _/               " << std::endl;
    std::cout << "       /_/ |_/____/___/  application  " << std::endl;
    std::cout << "Initializing KratosALEApplication...  " << std::endl;

    // register elements
    KRATOS_REGISTER_ELEMENT("LaplacianMeshMovingElement2D3N", mLaplacianMeshMovingElement2D3N);
    KRATOS_REGISTER_ELEMENT("LaplacianMeshMovingElemtent3D4N", mLaplacianMeshMovingElement3D4N);
    KRATOS_REGISTER_ELEMENT("LaplacianMeshMovingElement2D4N", mLaplacianMeshMovingElement2D4N);
    KRATOS_REGISTER_ELEMENT("LaplacianMeshMovingElemtent3D8N", mLaplacianMeshMovingElement3D8N);
    KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElement2D3N", mStructuralMeshMovingElement2D3N);
    KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElement2D4N", mStructuralMeshMovingElement2D4N);
    KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElement3D4N", mStructuralMeshMovingElement3D4N);
    KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElement3D8N", mStructuralMeshMovingElement3D8N);
    KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElement3D6N", mStructuralMeshMovingElement3D6N);
    KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElement3D15N", mStructuralMeshMovingElement3D15N);


    // register variables
    //MESH_VELOCITY currently put to the core since used in other applications
    //KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MESH_VELOCITY;
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MESH_DISPLACEMENT);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MESH_REACTION);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MESH_RHS);

}

}  // namespace Kratos.
