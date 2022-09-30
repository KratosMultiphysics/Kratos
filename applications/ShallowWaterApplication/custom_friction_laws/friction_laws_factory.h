//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_FRICTION_LAWS_FACTORY_H_INCLUDED
#define KRATOS_FRICTION_LAWS_FACTORY_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "friction_law.h"


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

/** 
 * @class FrictionLawsFactory
 * @ingroup ShallowWaterApplication
 * @brief The base class for the bottom and surface friction laws
 * @details This class does nothing, define derived friction laws in order to make use of it
 * @author Miguel Maso Sotomayor
 */
class FrictionLawsFactory
{
public:
    ///@name Type Definitions
    ///@{
    
    typedef Node <3> NodeType;

    typedef Geometry<NodeType> GeometryType;

    ///@}
    ///@name Pointer Definition
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(FrictionLawsFactory);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    FrictionLawsFactory() {}

    /**
     * @brief Destructor
     */
    ~FrictionLawsFactory() {}

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Create a bottom friction law
     */
    FrictionLaw::Pointer CreateBottomFrictionLaw(
        const GeometryType& rGeometry,
        const Properties& rProperty,
        const ProcessInfo& rProcessInfo);

    /**
     * @brief Create a surface friction law
     */
    FrictionLaw::Pointer CreateSurfaceFrictionLaw(
        const GeometryType& rGeometry,
        const Properties& rProperty,
        const ProcessInfo& rProcessInfo);

    ///@}
    ///@name Input and output
    ///@{

    /**
     * @brief Turn back information as a string.
     */
    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "FrictionLawsFactory";
        return buffer.str();
    }

    /**
     * @brief Print information about this object.
     */
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /**
     * @brief Print object's data.
     */
    void PrintData(std::ostream& rOStream) const {}

    ///@}
    ///@name Friends
    ///@{


    ///@}

}; // Class FrictionLawsFactory

///@}

}  // namespace Kratos.

#endif // KRATOS_FRICTION_LAWS_FACTORY_H_INCLUDED  defined
