//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           Fabruary 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_ELEMENT_UTILITIES )
#define  KRATOS_ELEMENT_UTILITIES

/* System includes */
#include <cmath>
#include "utilities/math_utils.h"

// Project includes
#include "includes/element.h"

namespace Kratos
{

class ElementUtilities
{
    
public:
    
    typedef std::size_t SizeType;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static inline void InterpolateVariableWithComponents(array_1d<double,2>& rVector,const Matrix& Ncontainer, 
                                                        const array_1d<double,8>& VariableWithComponents,const SizeType& GPoint)
    {
        //Quadrilateral_2d_4
        noalias(rVector) = ZeroVector(2);
        
        SizeType index = 0;
        for(SizeType i=0; i<4; i++)
        {
            rVector[0] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
            rVector[1] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
        }
    }
    
    //----------------------------------------------------------------------------------------
    
    static inline void InterpolateVariableWithComponents(array_1d<double,3>& rVector,const Matrix& Ncontainer, 
                                                        const array_1d<double,18>& VariableWithComponents,const SizeType& GPoint)
    {
        //Prism_3d_6
        noalias(rVector) = ZeroVector(3);
        
        SizeType index = 0;
        for(SizeType i=0; i<6; i++)
        {
            rVector[0] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
            rVector[1] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
            rVector[2] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
        }
    }
    
    //----------------------------------------------------------------------------------------

    static inline void InterpolateVariableWithComponents(array_1d<double,3>& rVector,const Matrix& Ncontainer, 
                                                        const array_1d<double,24>& VariableWithComponents,const SizeType& GPoint)
    {
        //Hexahedral_3d_8
        noalias(rVector) = ZeroVector(3);
        
        SizeType index = 0;
        for(SizeType i=0; i<8; i++)
        {
            rVector[0] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
            rVector[1] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
            rVector[2] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static inline void FillArray1dOutput(array_1d<double,3> rOutputValue, const array_1d<double,2>& ComputedValue)
    {
        rOutputValue[0] = ComputedValue[0];
        rOutputValue[1] = ComputedValue[1];
        rOutputValue[2] = 0.0;
    }
    
    //----------------------------------------------------------------------------------------
    
    static inline void FillArray1dOutput(array_1d<double,3> rOutputValue, const array_1d<double,3>& ComputedValue)
    {
        rOutputValue[0] = ComputedValue[0];
        rOutputValue[1] = ComputedValue[1];
        rOutputValue[2] = ComputedValue[2];
    }
    
    //----------------------------------------------------------------------------------------

    static inline void FillArray1dOutput(array_1d<double,6> rOutputValue, const boost::numeric::ublas::bounded_matrix<double,2,2>& ComputedValue)
    {
        rOutputValue[0] = ComputedValue(0,0);
        rOutputValue[1] = ComputedValue(1,1);
        rOutputValue[2] = 0.0;
        rOutputValue[3] = ComputedValue(0,1);
        rOutputValue[4] = 0.0;
        rOutputValue[5] = 0.0;
    }
    
    //----------------------------------------------------------------------------------------

    static inline void FillArray1dOutput(array_1d<double,6> rOutputValue, const boost::numeric::ublas::bounded_matrix<double,3,3>& ComputedValue)
    {
        rOutputValue[0] = ComputedValue(0,0);
        rOutputValue[1] = ComputedValue(1,1);
        rOutputValue[2] = ComputedValue(2,2);
        rOutputValue[3] = ComputedValue(0,1);
        rOutputValue[4] = ComputedValue(1,2);
        rOutputValue[5] = ComputedValue(0,2);
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static inline void GetDisplacementsVector(array_1d<double,8>& rDisplacementVector, const Element::GeometryType& Geom)
    {
        //Quadrilateral_2d_4
        SizeType index = 0;
        for(SizeType i=0; i<4; i++)
        {
            rDisplacementVector[index++] = Geom[i].FastGetSolutionStepValue(DISPLACEMENT_X);
            rDisplacementVector[index++] = Geom[i].FastGetSolutionStepValue(DISPLACEMENT_Y);
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void GetDisplacementsVector(array_1d<double,18>& rDisplacementVector, const Element::GeometryType& Geom)
    {
        //Prism_3d_6
        SizeType index = 0;
        for(SizeType i=0; i<6; i++)
        {
            rDisplacementVector[index++] = Geom[i].FastGetSolutionStepValue(DISPLACEMENT_X);
            rDisplacementVector[index++] = Geom[i].FastGetSolutionStepValue(DISPLACEMENT_Y);
            rDisplacementVector[index++] = Geom[i].FastGetSolutionStepValue(DISPLACEMENT_Z);
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void GetDisplacementsVector(array_1d<double,24>& rDisplacementVector, const Element::GeometryType& Geom)
    {
        //Hexahedral_3d_8
        SizeType index = 0;
        for(SizeType i=0; i<8; i++)
        {
            rDisplacementVector[index++] = Geom[i].FastGetSolutionStepValue(DISPLACEMENT_X);
            rDisplacementVector[index++] = Geom[i].FastGetSolutionStepValue(DISPLACEMENT_Y);
            rDisplacementVector[index++] = Geom[i].FastGetSolutionStepValue(DISPLACEMENT_Z);
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void GetVolumeAccelerationVector(array_1d<double,8>& rVolumeAccelerationVector, const Element::GeometryType& Geom)
    {
        //Quadrilateral_2d_4
        SizeType index = 0;
        for(SizeType i=0; i<4; i++)
        {
            rVolumeAccelerationVector[index++] = Geom[i].FastGetSolutionStepValue(VOLUME_ACCELERATION_X);
            rVolumeAccelerationVector[index++] = Geom[i].FastGetSolutionStepValue(VOLUME_ACCELERATION_Y);
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void GetVolumeAccelerationVector(array_1d<double,18>& rVolumeAccelerationVector, const Element::GeometryType& Geom)
    {
        //Prism_3d_6
        SizeType index = 0;
        for(SizeType i=0; i<6; i++)
        {
            rVolumeAccelerationVector[index++] = Geom[i].FastGetSolutionStepValue(VOLUME_ACCELERATION_X);
            rVolumeAccelerationVector[index++] = Geom[i].FastGetSolutionStepValue(VOLUME_ACCELERATION_Y);
            rVolumeAccelerationVector[index++] = Geom[i].FastGetSolutionStepValue(VOLUME_ACCELERATION_Z);
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void GetVolumeAccelerationVector(array_1d<double,24>& rVolumeAccelerationVector, const Element::GeometryType& Geom)
    {
        //Hexahedral_3d_8
        SizeType index = 0;
        for(SizeType i=0; i<8; i++)
        {
            rVolumeAccelerationVector[index++] = Geom[i].FastGetSolutionStepValue(VOLUME_ACCELERATION_X);
            rVolumeAccelerationVector[index++] = Geom[i].FastGetSolutionStepValue(VOLUME_ACCELERATION_Y);
            rVolumeAccelerationVector[index++] = Geom[i].FastGetSolutionStepValue(VOLUME_ACCELERATION_Z);
        }
    }

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

}; /* Class ElementUtilities*/
} /* namespace Kratos.*/

#endif /* KRATOS_ELEMENT_UTILITIES defined */
