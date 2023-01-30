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

#ifndef KRATOS_MANNING_LAW_H_INCLUDED
#define KRATOS_MANNING_LAW_H_INCLUDED


// System includes


// External includes


// Project includes
#include "friction_law.h"


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
 * @class ManningLaw
 * @ingroup ShallowWaterApplication
 * @brief The base class for the bottom and surface friction laws
 * @details This class computes the bottom friction according to the Manning law
 * @author Miguel Maso Sotomayor
 */
class ManningLaw : public FrictionLaw
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ManningLaw);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default Constructor
     */
    ManningLaw() {}

    /**
     * @brief Constructor with data
     */
    ManningLaw(
        const GeometryType& rGeometry,
        const Properties& rProperty,
        const ProcessInfo& rProcessInfo);

    /**
     * @brief Initialize the friction law variables
     */
    virtual ~ManningLaw() {}

    ///@}
    ///@name Operators
    ///@{


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
     * @param rHeight The layer depth
     * @param rVector The layer velocity or momentum
     * @return The LHS coefficient
     */
    double CalculateLHS(const double& rHeight, const array_1d<double,3>& rVelocity) override;

    /**
     * @brief Calculate the RHS coefficient for the given data
     * @param rHeight The layer depth
     * @param rVector The layer velocity or momentum
     * @return The components of the RHS coefficients
     */
    array_1d<double,3> CalculateRHS(const double& rHeight, const array_1d<double,3>& rVelocity) override;

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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "ManningLaw";
        return buffer.str();
    }

    ///@}

protected:

    ///@name Member variables
    ///@{

    double mManning2;

    double mEpsilon;

    ///@}

private:

    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ManningLaw& operator=(ManningLaw const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    ManningLaw(ManningLaw const& rOther) {}


    ///@}

}; // Class ManningLaw

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const ManningLaw& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MANNING_LAW_H_INCLUDED  defined
