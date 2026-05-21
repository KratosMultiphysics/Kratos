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

#ifndef KRATOS_FRICTION_LAW_H_INCLUDED
#define KRATOS_FRICTION_LAW_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/node.h"
#include "includes/properties.h"
#include "includes/process_info.h"
#include "geometries/geometry.h"


namespace Kratos
{
///@addtogroup ShallowWaterApplication
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

/** 
 * @class FrictionLaw
 * @ingroup ShallowWaterApplication
 * @brief The base class for the bottom and surface friction laws
 * @details This class does nothing, define derived friction laws in order to make use of it
 * @author Miguel Maso Sotomayor
 */
class FrictionLaw
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(FrictionLaw);
    
    typedef Node NodeType;

    typedef Geometry<NodeType> GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FrictionLaw() {}

    /// Destructor.
    virtual ~FrictionLaw() {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Check on the completeness of the input
     */
    virtual int Check()
    {
        return 0;
    }

    /**
     * @brief Initialize the friction law variables
     */
    virtual void Initialize(
        const GeometryType& rGeometry,
        const Properties& rProperty,
        const ProcessInfo& rProcessInfo) {}

    /**
     * @brief Calculate the LHS coefficient for the given data
     * @param rHeight The layer depth
     * @param rVector The layer velocity
     * @return The LHS coefficient
     */
    virtual double CalculateLHS(const double& rHeight, const array_1d<double,3>& rVelocity)
    {
        return 0.0;
    }

    /**
     * @brief Calculate the RHS coefficient for the given data
     * @param rHeight The layer depth
     * @param rVector The layer velocity
     * @return The components of the RHS coefficients
     */
    virtual array_1d<double,3> CalculateRHS(const double& rHeight, const array_1d<double,3>& rVelocity)
    {
        return ZeroVector(3);
    }

    /**
     * @brief Calculate the LHS coefficient for the given data
     * @param rInnerLayer The current layer velocity
     * @return The LHS coefficient
     */
    virtual double CalculateLHS(const array_1d<double,3>& rInnerLayer)
    {
        return 0.0;
    }

    /**
     * @brief Calculate the RHS coefficient for the given data
     * @param rInnerLayer The current layer velocity
     * @return The components of the RHS coefficients
     */
    virtual array_1d<double,3> CalculateRHS(const array_1d<double,3>& rInnerLayer)
    {
        return ZeroVector(3);
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /**
     * @brief Turn back information as a string.
     */
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "FrictionLaw";
        return buffer.str();
    }

    /**
     * @brief Print information about this object.
     */
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /**
     * @brief Print object's data.
     */
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}
    ///@name Friends
    ///@{


    ///@}

private:

    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    FrictionLaw& operator=(FrictionLaw const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    FrictionLaw(FrictionLaw const& rOther) {}


    ///@}

}; // Class FrictionLaw

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                FrictionLaw& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const FrictionLaw& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FRICTION_LAW_H_INCLUDED  defined
