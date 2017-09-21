
/*
==============================================================================
KratosShallowWaterApplication
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last Modified by:    Miguel Masó Sotomayor
//   Date:                September 19th 2017
//   Revision:            2.1
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
class ShallowParticle : public Point<3>
{
public: 
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Point
    typedef ShallowParticle* Pointer;
    
    typedef double TDataType; 

    typedef array_1d<double,3> BaseType;

    typedef Point<3, double> Type;

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
    ShallowParticle(TDataType const& NewX, TDataType const& NewY, TDataType const& NewZ) : Point<3>(NewX, NewY, NewZ)
    {
        this->ERASE_FLAG=true; //initializing as useless particle
        this->SCALAR1=0.0;
        this->VECTOR1=ZeroVector(3);
        this->VELOCITY=ZeroVector(3);
    }

    ShallowParticle() : Point<3>(0.0, 0.0, 0.0)
    {
        this->ERASE_FLAG=true;
        this->SCALAR1=0.0;
        this->VECTOR1=ZeroVector(3);
        this->VELOCITY=ZeroVector(3);
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

    float& GetVelocity(const unsigned int i)
    {
        return this->VELOCITY[i];
    }

    array_1d<float,3>& GetVelocity()
    {
        return this->VELOCITY;
    }

    bool& GetEraseFlag()
    {
        return this->ERASE_FLAG;
    }

private:
    float SCALAR1;
    array_1d<float,3> VECTOR1;
    array_1d<float,3> VELOCITY;
    bool ERASE_FLAG;

};

} //namespace Kratos

#endif //KRATOS_SHALLOW_WATER_PARTICLE_H_INCLUDED defined
