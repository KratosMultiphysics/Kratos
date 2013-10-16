//
//   Project Name:        Kratos
//   Last Modified by:    $Author: G.Casas $
//   Date:                $Date: 2011-6-13 08:56:42 $
//   Revision:            $Revision: 1.5 $
//
//
//README::::look to the key word "VERSION" if you want to find all the points where you have to change something so that you can pass from a kdtree to a bin data search structure;

#if !defined(KRATOS_CUSTOM_FUNCTIONS)
#define KRATOS_CUSTOM_FUNCTIONS

// /* External includes */

// System includes

// Project includes
#include "includes/model_part.h"
#include "utilities/timer.h"

//Database includes
#include "custom_utilities/discrete_particle_configure.h"
#include "discrete_particle_configure.h"
//SALVA_BEGINNING
// Project includes
#include "includes/define.h"
#include "../../DEM_application/custom_elements/discrete_element.h"
#include "../../DEM_application/custom_elements/spheric_swimming_particle.h"
//
#include "includes/define.h"
//#include "custom_elements/spheric_swimming_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
//#include "../../DEM_application/DEM_application.h"

#include "../../DEM_application/custom_elements/spheric_particle.h"
#include "../swimming_DEM_application.h"
#include "../../../kratos/utilities/geometry_utilities.h"

//
//SALVA_ENDING
//const double prox_tol = 0.00000000001;

namespace Kratos
{
//SALVA_STARTING

    void CalculateGeometryData2D(
            Element::GeometryType& geom,
            boost::numeric::ublas::bounded_matrix<double, 4, 3 > & DN_DX,
            array_1d<double, 4 > & N,
            double& Area) {
        double x10 = geom[1].X() - geom[0].X();
        double y10 = geom[1].Y() - geom[0].Y();

        double x20 = geom[2].X() - geom[0].X();
        double y20 = geom[2].Y() - geom[0].Y();

        //Jacobian is calculated:
        //  |dx/dxi  dx/deta|	|x1-x0   x2-x0|
        //J=|				|=	|			  |
        //  |dy/dxi  dy/deta|	|y1-y0   y2-y0|


        double detJ = x10 * y20 - y10 * x20;

        DN_DX(0, 0) = -y20 + y10;
        DN_DX(0, 1) = x20 - x10;
        DN_DX(1, 0) = y20;
        DN_DX(1, 1) = -x20;
        DN_DX(2, 0) = -y10;
        DN_DX(2, 1) = x10;

        DN_DX /= detJ;
        N[0] = 0.333333333333333;
        N[1] = 0.333333333333333;
        N[2] = 0.333333333333333;

        Area = 0.5 * detJ;
    }

class CustomFunctionsCalculator
{
public:

        static const std::size_t space_dim                  = 3; ///WARNING: generalize to 2d.
        typedef DiscreteParticleConfigure<space_dim>        Configure;
        typedef Configure::ContainerType                    ParticlePointerVector;
        typedef ParticlePointerVector::iterator             ParticlePointerIterator;
        typedef Configure::IteratorType                     ParticleIterator;


    KRATOS_CLASS_POINTER_DEFINITION(CustomFunctionsCalculator);

  
    CustomFunctionsCalculator() {};
    /// Calculator

    virtual ~CustomFunctionsCalculator() {};

    /// Default calculator
    
    void CalculatePressureGradient(ModelPart& r_model_part) {

            for (ModelPart::NodesContainerType::iterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++) {
                inode->FastGetSolutionStepValue(AUX_DOUBLE_VAR) = 0.0;
                noalias(inode->FastGetSolutionStepValue(PRESSURE_GRADIENT)) = ZeroVector(3);
            }
            const std::size_t TDim = 3;
            array_1d <double, TDim + 1 > elemental_pressures;
            array_1d <double, TDim> grad;
            array_1d <double, TDim + 1 > N; //Shape functions vector//
            boost::numeric::ublas::bounded_matrix<double, TDim + 1, TDim> DN_DX;

            for (ModelPart::ElementIterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ielem++) {
                //compute shape function derivatives
                Geometry< Node < 3 > >& geom = ielem->GetGeometry();
                double Volume;

                if (geom.size() == 4)
                    GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
                else 
                    CalculateGeometryData2D(geom, DN_DX, N, Volume);
                    
                //get the pressure gradients;
                
                for (unsigned int i = 0; i < geom.size(); i++)
                    elemental_pressures[i] = geom[i].FastGetSolutionStepValue(PRESSURE);

                noalias(grad) = prod(trans(DN_DX), elemental_pressures);
                double nodal_area = Volume / static_cast<double>(geom.size());
                grad *= nodal_area;

                for (unsigned int i = 0; i < geom.size(); i++)
                    geom[i].FastGetSolutionStepValue(PRESSURE_GRADIENT) += grad;

                for (unsigned int i = 0; i < geom.size(); i++)
                    geom[i].FastGetSolutionStepValue(AUX_DOUBLE_VAR) += nodal_area;

            }

            for (ModelPart::NodesContainerType::iterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++) {
                inode->FastGetSolutionStepValue(PRESSURE_GRADIENT) /= inode->FastGetSolutionStepValue(AUX_DOUBLE_VAR);
            }
        }

    }; // Class CustomFunctionsCalculator


} // namespace Kratos.

#endif // KRATOS_CREATE_AND_DESTROY  defined


