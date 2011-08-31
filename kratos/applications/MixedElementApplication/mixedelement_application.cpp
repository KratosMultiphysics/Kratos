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
#include "mixedelement_application.h"
#include "includes/variables.h"


namespace Kratos
{
    //Example

    KRATOS_CREATE_VARIABLE(double, SX)
    KRATOS_CREATE_VARIABLE(double, SY)
    KRATOS_CREATE_VARIABLE(double, SZ)
    KRATOS_CREATE_VARIABLE(double, SXY)
    KRATOS_CREATE_VARIABLE(double, SXZ)
    KRATOS_CREATE_VARIABLE(double, SYZ)
    //	KRATOS_CREATE_VARIABLE(double, IS_INTERFACE);
    //	KRATOS_CREATE_VARIABLE(double, NODAL_AREA);
    //

    KratosMixedElementApplication::KratosMixedElementApplication():
     		mSigmaUElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
                mIrriducibleElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>()))))
    // 		mElem3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>()))))
    {
    }
    //

    void KratosMixedElementApplication::Register()
    {
        // calling base class register to register Kratos components
        KratosApplication::Register();
        std::cout << "Initializing KratosMixedElementApplication... " << std::endl;

        KRATOS_REGISTER_VARIABLE( SX)
        KRATOS_REGISTER_VARIABLE( SY)
        KRATOS_REGISTER_VARIABLE( SZ)
        KRATOS_REGISTER_VARIABLE( SXY)
        KRATOS_REGISTER_VARIABLE( SXZ)
        KRATOS_REGISTER_VARIABLE( SYZ)

        KRATOS_REGISTER_ELEMENT("SigmaUElement2D", mSigmaUElement2D);
        KRATOS_REGISTER_ELEMENT("IrriducibleElement2D", mIrriducibleElement2D);
//        KRATOS_REGISTER_ELEMENT("Elemt3D", mElem3D);

    }

} // namespace Kratos.


