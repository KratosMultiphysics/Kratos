//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_POROMECHANICS_MATH_UTILITIES )
#define  KRATOS_POROMECHANICS_MATH_UTILITIES

/* System includes */
#include <cmath>
#include "utilities/math_utils.h"

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"

namespace Kratos
{

class PoromechanicsMathUtilities
{
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

public:

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static inline double CalculateVonMises(const Vector& StressVector)
    {
        KRATOS_TRY

        Matrix LocalStressTensor  = MathUtils<double>::StressVectorToTensor(StressVector); //reduced dimension stress tensor

        Matrix StressTensor = ZeroMatrix(3); //3D stress tensor
        for(unsigned int i=0; i<LocalStressTensor.size1(); i++)
        {
            for(unsigned int j=0; j<LocalStressTensor.size2(); j++)
            {
                StressTensor(i,j) = LocalStressTensor(i,j);
            }
        }

        double SigmaEquivalent = 0.5 * ( (StressTensor(0,0)-StressTensor(1,1))*(StressTensor(0,0)-StressTensor(1,1)) +
                                         (StressTensor(1,1)-StressTensor(2,2))*(StressTensor(1,1)-StressTensor(2,2)) +
                                         (StressTensor(2,2)-StressTensor(0,0))*(StressTensor(2,2)-StressTensor(0,0)) +
                                         6.0*( StressTensor(0,1)*StressTensor(0,1)+StressTensor(1,2)*StressTensor(1,2)+StressTensor(2,0)*StressTensor(2,0) ) );

        if( SigmaEquivalent < 0 ) SigmaEquivalent = 0.0;

        SigmaEquivalent = sqrt(SigmaEquivalent);

        return SigmaEquivalent;

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static inline void CalculateRotationMatrix(Element::GeometryType& geom,boost::numeric::ublas::bounded_matrix<double,2,2>& RotationMatrix)
    {   
        //Define mid-plane points for quadrilateral_interface_2d_4
        array_1d<double, 3> pmid0 = 0.5 * (geom.GetPoint( 0 ) + geom.GetPoint( 3 ));
        array_1d<double, 3> pmid1 = 0.5 * (geom.GetPoint( 1 ) + geom.GetPoint( 2 ));
        
        //Unitary vector in local x direction
		array_1d<double, 3> Vx( pmid1 - pmid0 );
        double inv_norm = 1.0/norm_2(Vx);
        Vx = inv_norm*Vx;
        
        //Rotation Matrix
        RotationMatrix(0,0) = Vx[0];
        RotationMatrix(0,1) = Vx[1];
        
        RotationMatrix(1,0) = -Vx[1];
        RotationMatrix(1,1) = Vx[0];
    }
    
    //------------------------------------------------------------------------------------
    
    static inline void CalculateRotationMatrix(Element::GeometryType& geom,boost::numeric::ublas::bounded_matrix<double,3,3>& RotationMatrix)
    {
        array_1d<double, 3> pmid0;
        array_1d<double, 3> pmid1;
        array_1d<double, 3> pmid2;
        
        //Define mid-plane points for prism_interface_3d_6
        if(geom.PointsNumber()==6)
        {
            pmid0 = 0.5 * (geom.GetPoint( 0 ) + geom.GetPoint( 3 ));
            pmid1 = 0.5 * (geom.GetPoint( 1 ) + geom.GetPoint( 4 ));
            pmid2 = 0.5 * (geom.GetPoint( 2 ) + geom.GetPoint( 5 ));
        }
        //Define mid-plane points for hexahedra_interface_3d_8
        else
        {
            pmid0 = 0.5 * (geom.GetPoint( 0 ) + geom.GetPoint( 4 ));
            pmid1 = 0.5 * (geom.GetPoint( 1 ) + geom.GetPoint( 5 ));
            pmid2 = 0.5 * (geom.GetPoint( 2 ) + geom.GetPoint( 6 ));
            //array_1d<double, 3> pmid3 = 0.5 * (geom.GetPoint( 3 ) + geom.GetPoint( 7 ));
        }
        
        //Unitary vector in local x direction
        array_1d<double, 3> Vx( pmid1 - pmid0 );
        double inv_norm_x = 1.0/norm_2(Vx);
        Vx = inv_norm_x*Vx;
        
        //Unitary vector in local z direction
        array_1d<double, 3> Vaux( pmid2 - pmid0 );
        array_1d<double, 3> Vz;
        MathUtils<double>::CrossProduct( Vz, Vx, Vaux);
        double inv_norm_z = 1.0/norm_2(Vz);
        Vz = inv_norm_z*Vz;
        
        //Unitary vector in local y direction
        array_1d<double, 3> Vy;
        MathUtils<double>::CrossProduct( Vy, Vz, Vx);
        
        //Rotation Matrix
        RotationMatrix(0,0) = Vx[0];
        RotationMatrix(0,1) = Vx[1];
        RotationMatrix(0,2) = Vx[2];
        
        RotationMatrix(1,0) = Vy[0];
        RotationMatrix(1,1) = Vy[1];
        RotationMatrix(1,2) = Vy[2];
        
        RotationMatrix(2,0) = Vz[0];
        RotationMatrix(2,1) = Vz[1];
        RotationMatrix(2,2) = Vz[2];
    }

}; /* Class PoromechanicsMathUtilities*/
} /* namespace Kratos.*/

#endif /* KRATOS_POROMECHANICS_MATH_UTILITIES defined */
