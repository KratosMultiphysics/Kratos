// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___ 
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _| 
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  Pablo Becker
//

#if !defined(KRATOS_CONVECTION_PARTICLE_H_INCLUDED )
#define  KRATOS_CONVECTION_PARTICLE_H_INCLUDED

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
class Convection_Particle : public Point
{
public: 
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Point
    typedef Convection_Particle* Pointer;
    
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
    Convection_Particle(TDataType const& NewX, TDataType const& NewY, TDataType const& NewZ) : Point(NewX, NewY, NewZ)
    {
		this->ERASE_FLAG=true; //initializing as useless particle
		this->SCALAR1=0.0;
    }
    
    Convection_Particle() : Point(0.0, 0.0, 0.0)
    {
		this->ERASE_FLAG=true;
		this->SCALAR1=0.0;
    }
    
    ~Convection_Particle()
		{}
    
	float& GetScalar1()
	{
		return this->SCALAR1;
	}
	
	bool& GetEraseFlag()
	{
		return this->ERASE_FLAG;
	}
       
private: 
	float SCALAR1;
	//double TEMPERATURE;
	//double OXYGEN;
	//Element::Pointer ELEMENT_WEAKPOINTER;
	//unsigned int ELEMENT_ID;
	//double GRADIENT_DISCONTINUITY;
	bool ERASE_FLAG;
	


};

} //namespace Kratos

#endif //KRATOS_CONVECTION_PARTICLE_H_INCLUDED defined
