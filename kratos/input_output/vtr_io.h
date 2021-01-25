//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#if !defined(KRATOS_VTR_IO_H_INCLUDED )
#define  KRATOS_VTR_IO_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/io.h"
#include "input_output/xml_io.h"
#include "processes/internals/cartesian_mesh_colors.h"


namespace Kratos
{
///@addtogroup ApplicationNameApplication
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

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(KRATOS_CORE) VtrIO : public IO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of VtrIO
    KRATOS_CLASS_POINTER_DEFINITION(VtrIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    VtrIO() = delete;


    /// Constructor with filename.
    VtrIO(std::string const& Filename);

    /// Constructor with stream.
    VtrIO(std::iostream* pInputStream);


    /// Destructor.
    virtual ~VtrIO(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    void Read(std::vector<Internals::CartesianMeshColors>& MultiBlockMesh);

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override;


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    XmlIO mXmlIO;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    VtrIO& operator=(VtrIO const& rOther);

    /// Copy constructor.
    VtrIO(VtrIO const& rOther);


    ///@}

}; // Class VtrIO

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                VtrIO& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const VtrIO& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_VTR_IO_H_INCLUDED  defined
