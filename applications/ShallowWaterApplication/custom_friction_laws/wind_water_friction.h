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

#ifndef KRATOS_WIND_WATER_FRICTION_H_INCLUDED
#define KRATOS_WIND_WATER_FRICTION_H_INCLUDED


// System includes


// External includes


// Project includes
#include "friction_law.h"
#include "shallow_water_application_variables.h"


namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name Kratos Classes
///@{

/** 
 * @class WindWaterFriction
 * @ingroup ShallowWaterApplication
 * @brief The base class for the bottom and surface friction laws
 * @details This class computes the bottom friction according to the Manning law
 * @author Miguel Maso Sotomayor
 */
class WindWaterFriction : public FrictionLaw
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(WindWaterFriction);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default Constructor
     */
    WindWaterFriction() {}

    /**
     * @brief Constructor with data
     */
    WindWaterFriction(
        const GeometryType& rGeometry,
        const Properties& rProperty,
        const ProcessInfo& rProcessInfo);

    /**
     * @brief Destructor
     */
    virtual ~WindWaterFriction() {}

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Initialize the friction law variables
     */
    void Initialize(
        const GeometryType& rGeometry,
        const Properties& rProperty,
        const ProcessInfo& rProcessInfo) override;

    /**
     * @brief Calculate the LHS coefficient for the given data
     * @param rVector The layer velocity
     * @return The LHS coefficient
     */
    double CalculateLHS(const array_1d<double,3>& rVelocity) override;

    /**
     * @brief Calculate the RHS coefficient for the given data
     * @param rVector The layer velocity
     * @return The components of the RHS coefficients
     */
    array_1d<double,3> CalculateRHS(const array_1d<double,3>& rVelocity) override;

    ///@}
    ///@name Input and output
    ///@{

    /**
     * @brief Turn back information as a string.
     */
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "WindWaterFriction";
        return buffer.str();
    }

    ///@}

private:

    ///@name Member variables
    ///@{

    double mAirDensity;
    double mWaterDensity;
    array_1d<double,3> mWind;

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    WindWaterFriction& operator=(WindWaterFriction const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    WindWaterFriction(WindWaterFriction const& rOther) {}

    ///@}

}; // Class WindWaterFriction

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const WindWaterFriction& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_WIND_WATER_FRICTION_H_INCLUDED  defined
