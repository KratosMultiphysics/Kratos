//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:         BSD License 
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

#if !defined(KRATOS_ELEMENT_SIZE_CALCULATOR_H )
#define  KRATOS_ELEMENT_SIZE_CALCULATOR_H

// External includes 

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_utilities/fluid_element_data.h"


namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

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

/// Auxiliary and specialized functions for elements derived from FluidElement.
template< std::size_t TDim, std::size_t TNumNodes >
class ElementSizeCalculator {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ElementSizeCalculator
    KRATOS_CLASS_POINTER_DEFINITION(ElementSizeCalculator);

    using ShapeDerivativesType = typename FluidElementData<TDim,TNumNodes,false>::ShapeDerivativesType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Deleted default constructor
    ElementSizeCalculator() = delete;
    
    /// Deleted copy constructor.
    ElementSizeCalculator(ElementSizeCalculator const& rOther) = delete;

    /// Destructor.
    ~ElementSizeCalculator();

    ///@}
    ///@name Operators
    ///@{

    /// Deleted assignment operator.
    ElementSizeCalculator& operator=(ElementSizeCalculator const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{
    
    static double MinimumElementSize(const Geometry<Node<3> >& rGeometry);

    ///@}

};  // Class ElementSizeCalculator

///@}

///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_ELEMENT_SIZE_CALCULATOR_H  defined 


