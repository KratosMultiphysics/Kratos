//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//                   Pablo Becker
//

#if !defined(KRATOS_SHALLOW_WATER_PARTICLE_H_INCLUDED )
#define  KRATOS_SHALLOW_WATER_PARTICLE_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>



// External includes


// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "includes/serializer.h"
#include "geometries/point.h"
#include "includes/model_part.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// PFEM Particle class.

/*
@see Geometry
@see Node
@see IntegrationPoint
*/
//template<std::size_t TDimension, class TDataType = double> //always size 3!
class ShallowParticle : public Point
{
public: 
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Point
    typedef ShallowParticle* Pointer;
    
    typedef double TDataType; 

    typedef array_1d<double,3> BaseType;

    typedef Point Type;

    typedef BaseType CoordinatesArrayType;

    typedef std::size_t SizeType;

    typedef std::size_t IndexType;

    ///@}
    ///@name Constants
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ShallowParticle(TDataType const& NewX, TDataType const& NewY, TDataType const& NewZ) : Point(NewX, NewY, NewZ)
    {
        this->ERASE_FLAG = true; //initializing as useless particle
        this->SCALAR1 = 0.0;
        this->VECTOR1 = ZeroVector(3);
    }

    ShallowParticle() : Point(0.0, 0.0, 0.0)
    {
        this->ERASE_FLAG = true;
        this->SCALAR1 = 0.0;
        this->VECTOR1 = ZeroVector(3);
    }

    ~ShallowParticle()
    {
    }

    // Returning references
    float& GetScalar1()
    {
        return this->SCALAR1;
    }

    float& GetVector1(const unsigned int i)
    {
        return this->VECTOR1[i];
    }

    array_1d<float,3>& GetVector1()
    {
        return this->VECTOR1;
    }

    bool& GetEraseFlag()
    {
        return this->ERASE_FLAG;
    }

private:
    float SCALAR1;
    array_1d<float,3> VECTOR1;
    bool ERASE_FLAG;

};

} //namespace Kratos

#endif //KRATOS_SHALLOW_WATER_PARTICLE_H_INCLUDED defined
