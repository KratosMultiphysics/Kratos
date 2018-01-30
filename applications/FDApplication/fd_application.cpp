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
#include "fd_application.h"
#include "includes/variables.h"

namespace Kratos {

	KratosFDApplication::KratosFDApplication() {}

	void KratosFDApplication::Register() {
		KratosApplication::Register();
		std::cout << "Initializing KratosFDApplication... " << std::endl;
	}

}  // namespace Kratos.
