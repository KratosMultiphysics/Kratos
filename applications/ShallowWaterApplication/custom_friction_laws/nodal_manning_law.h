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

#ifndef KRATOS_NODAL_MANNING_LAW_H_INCLUDED
#define KRATOS_NODAL_MANNING_LAW_H_INCLUDED


// System includes


// External includes


// Project includes
#include "manning_law.h"


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
 * @class NodalManningLaw
 * @ingroup ShallowWaterApplication
 * @brief This class computes the bottom friction according to the Manning law
 * @author Miguel Maso Sotomayor
 */
class NodalManningLaw : public ManningLaw
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(NodalManningLaw);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    NodalManningLaw() {}

    /**
     * @brief Constructor with data
     */
    NodalManningLaw(
        const GeometryType& rGeometry,
        const Properties& rProperty,
        const ProcessInfo& rProcessInfo);

    /**
     * @brief Destructor
     */
    virtual ~NodalManningLaw() {}

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
        buffer << "NodalManningLaw";
        return buffer.str();
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}

private:

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    NodalManningLaw& operator=(NodalManningLaw const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    NodalManningLaw(NodalManningLaw const& rOther) {}


    ///@}

}; // Class NodalManningLaw

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const NodalManningLaw& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_NODAL_MANNING_LAW_H_INCLUDED  defined
