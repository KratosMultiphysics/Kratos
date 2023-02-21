/*
==============================================================================
Kratos
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
//   Last Modified by:    $Author: pablo $
//   Date:                $Date: 2013-03-05 11:34:12 $
//   Revision:            $Revision: 1.5 $
//
//


#if !defined(KRATOS_PFEM_PARTICLE_FLUID_H_INCLUDED )
#define  KRATOS_PFEM_PARTICLE_FLUID_H_INCLUDED



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
class PFEM_Particle_Fluid : public Point
{
public: 
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Point
    //KRATOS_CLASS_POINTER_DEFINITION(PFEM_Particle_Fluid);
    typedef PFEM_Particle_Fluid* Pointer;
    
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
    PFEM_Particle_Fluid(TDataType const& NewX, TDataType const& NewY, TDataType const& NewZ) : Point(NewX, NewY, NewZ)
    {
		this->ERASE_FLAG=true; //initializing as useless particle
		this->VELOCITY=ZeroVector(3);
		this->DISTANCE=0.0;
    }
    
    PFEM_Particle_Fluid() : Point(0.0, 0.0, 0.0)
    {
		this->ERASE_FLAG=true;
		this->VELOCITY=ZeroVector(3);
		this->DISTANCE=0.0;
    }
    
    ~PFEM_Particle_Fluid()
		{}
    
    //returning references
	array_1d<float,3>& GetVelocity()
	{
		return this->VELOCITY;
	}
	
	float& GetVelocity(const unsigned int i)
	{
		return this->VELOCITY[i];
	}
	
	float& GetDistance()
	{
		return this->DISTANCE;
	}
	
	//double& GetTemperature()
	//{
	//	return this->TEMPERATURE;
	//}
	
	//double& GetOxygen()
	//{
	//	return this->OXYGEN;
	//}
	/*
	double& GetPressure()
	{
		return this->PRESSURE;
	}
	double& GetGradientDiscontinuity()
	{
		return this->GRADIENT_DISCONTINUITY;
	}
	*/
	/*
	unsigned int& GetElementId()
	{
		return this->ELEMENT_ID;
	}
	*/
	bool& GetEraseFlag()
	{
		return this->ERASE_FLAG;
	}
	
	

    
    
    //Copy constructor. Initialize this point with the coordinates
    //of given particle.
    
    //PFEM_Particle(PFEM_Particle const& rOtherParticle)
     //   : Point(rOtherParticle) {}
      
       
private: 
	array_1d<float,3> VELOCITY;
	float DISTANCE;
	//double TEMPERATURE;
	//double OXYGEN;
	//Element::Pointer ELEMENT_WEAKPOINTER;
	//unsigned int ELEMENT_ID;
	//double GRADIENT_DISCONTINUITY;
	bool ERASE_FLAG;
	


};

//typedef std::vector<PFEM_Particle_Fluid*> FluidParticlePointerVector;


//inline std::ostream& operator << (std::ostream& rOStream,
//                                 const std::vector<PFEM_Particle_Fluid*>& rThis)
//{
//}


} //namespace Kratos

#endif //KRATOS_PFEM_PARTICLE_FLUID_H_INCLUDED defined
