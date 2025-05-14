// KRATOS
// _____   __               __  __      _ _   _
//|  __ \ / _|             |  \/  |    | | | (_)
//| |__) | |_ ___ _ __ ___ | \  / | ___| | |_ _ _ __   __ _
//|  ___/|  _/ _ \ '_ ` _ \| |\/| |/ _ \ | __| | '_ \ / _` |
//| |    | ||  __/ | | | | | |  | |  __/ | |_| | | | | (_| |
//|_|    |_| \___|_| |_| |_|_|  |_|\___|_|\__|_|_| |_|\__, |
//                                                     __/ |
//                                                    |___/ APPLICATION
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Julio Marti
//


// System includes

// External includes
//

// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_27.h"
#include "pfem_melting_application.h"
#include "includes/variables.h"

namespace Kratos {



KratosPfemMeltingApplication::KratosPfemMeltingApplication()
    : KratosApplication("PfemMeltingApplication"),

	mLagrangianFluidVMS2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
	mLagrangianFluidVMS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
	//mHYPOELASTICSOLID2D  (0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
	//mHYPOELASTICSOLID3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4))))
	mHYPO2D  (0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
	mHYPO3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
	mQFLUID2D  (0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
	mQFLUID3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),

    mEulerianConvDiffLumped2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
    mEulerianConvDiffLumped3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4))))
      {}

void KratosPfemMeltingApplication::Register() {
    KRATOS_INFO("") <<
    " KRATOS                                                     " << std::endl <<
    " _____   __               __  __      _ _   _               " << std::endl <<
    "|  __ \\ / _|             |  \\/  |    | | | (_)              " << std::endl <<
    "| |__) | |_ ___ _ __ ___ | \\  / | ___| | |_ _ _ __   __ _   " << std::endl <<
    "|  ___/|  _/ _ \\ '_ ` _ \\| |\\/| |/ _ \\ | __| | '_ \\ / _` |  " << std::endl <<
    "| |    | ||  __/ | | | | | |  | |  __/ | |_| | | | | (_| |  " << std::endl <<
    "|_|    |_| \\___|_| |_| |_|_|  |_|\\___|_|\\__|_|_| |_|\\__, |  " << std::endl <<
    "                                                     __/ |  " << std::endl <<
    "                                                    |___/ APPLICATION " << std::endl;

    // Registering variables
    KRATOS_REGISTER_VARIABLE(ACTIVATION_ENERGY)
    KRATOS_REGISTER_VARIABLE(ARRHENIUS_COEFFICIENT)
    KRATOS_REGISTER_VARIABLE(RADIOUS)
    KRATOS_REGISTER_VARIABLE(HEAT_OF_VAPORIZATION)
    KRATOS_REGISTER_VARIABLE(ARRHENIUS_VALUE)
    KRATOS_REGISTER_VARIABLE(NODAL_VOLUME)
    KRATOS_REGISTER_VARIABLE(IS_SOLID)

    KRATOS_REGISTER_VARIABLE(DELTA_SIGMA_XX)
    KRATOS_REGISTER_VARIABLE(DELTA_SIGMA_XY)
    KRATOS_REGISTER_VARIABLE(DELTA_SIGMA_XZ)

    KRATOS_REGISTER_VARIABLE(DELTA_SIGMA_YX)
    KRATOS_REGISTER_VARIABLE(DELTA_SIGMA_YY)
    KRATOS_REGISTER_VARIABLE(DELTA_SIGMA_YZ)

    KRATOS_REGISTER_VARIABLE(DELTA_SIGMA_ZX)
    KRATOS_REGISTER_VARIABLE(DELTA_SIGMA_ZY)
    KRATOS_REGISTER_VARIABLE(DELTA_SIGMA_ZZ)


    KRATOS_REGISTER_VARIABLE(HISTORICAL_SIGMA_XX)
    KRATOS_REGISTER_VARIABLE(HISTORICAL_SIGMA_XY)
    KRATOS_REGISTER_VARIABLE(HISTORICAL_SIGMA_XZ)

    KRATOS_REGISTER_VARIABLE(HISTORICAL_SIGMA_YX)
    KRATOS_REGISTER_VARIABLE(HISTORICAL_SIGMA_YY)
    KRATOS_REGISTER_VARIABLE(HISTORICAL_SIGMA_YZ)

    KRATOS_REGISTER_VARIABLE(HISTORICAL_SIGMA_ZX)
    KRATOS_REGISTER_VARIABLE(HISTORICAL_SIGMA_ZY)
    KRATOS_REGISTER_VARIABLE(HISTORICAL_SIGMA_ZZ)

    KRATOS_REGISTER_VARIABLE(PRESSUREAUX)
    
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(INITIAL_POSITION)
    
    
    //element
        //element
    KRATOS_REGISTER_VARIABLE( TOTAL_CAUCHY_STRESS )
  

    KRATOS_REGISTER_ELEMENT("LagrangianFluidVMS2D", mLagrangianFluidVMS2D);
    KRATOS_REGISTER_ELEMENT("LagrangianFluidVMS3D", mLagrangianFluidVMS3D);
    //KRATOS_REGISTER_ELEMENT("HYPOELASTICSOLID2D", mHYPOELASTICSOLID2D);
    //KRATOS_REGISTER_ELEMENT("HYPOELASTICSOLID3D", mHYPOELASTICSOLID3D);
    KRATOS_REGISTER_ELEMENT("HYPO2D", mHYPO2D);
    KRATOS_REGISTER_ELEMENT("HYPO3D", mHYPO3D);

    KRATOS_REGISTER_ELEMENT("QFLUID2D", mQFLUID2D);
    KRATOS_REGISTER_ELEMENT("QFLUID3D", mQFLUID3D);

    KRATOS_REGISTER_ELEMENT("EulerianConvDiffLumped2D", mEulerianConvDiffLumped2D);
    KRATOS_REGISTER_ELEMENT("EulerianConvDiffLumped3D", mEulerianConvDiffLumped3D);

}

}  // namespace Kratos.
