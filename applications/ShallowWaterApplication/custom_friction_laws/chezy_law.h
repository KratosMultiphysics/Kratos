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

#ifndef KRATOS_CHEZY_LAW_H_INCLUDED
#define KRATOS_CHEZY_LAW_H_INCLUDED


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
 * @class ChezyLaw
 * @ingroup ShallowWaterApplication
 * @brief The base class for the bottom and surface friction laws
 * @details This class computes the bottom friction according to the Chezy law
 * @author Miguel Maso Sotomayor
 */
class ChezyLaw : public FrictionLaw
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ChezyLaw);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default Constructor
     */
    ChezyLaw() {}

    /**
     * @brief Constructor with data
     */
    ChezyLaw(
        const GeometryType& rGeometry,
        const Properties& rProperty,
        const ProcessInfo& rProcessInfo);

    /**
     * @brief Destructor
     */
    virtual ~ChezyLaw() {}

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
        buffer << "ChezyLaw";
        return buffer.str();
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}

private:

    ///@name Member variables
    ///@{

    double mCoefficient;

    double mEpsilon;

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ChezyLaw& operator=(ChezyLaw const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    ChezyLaw(ChezyLaw const& rOther) {}


    ///@}

}; // Class ChezyLaw

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const ChezyLaw& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CHEZY_LAW_H_INCLUDED  defined
