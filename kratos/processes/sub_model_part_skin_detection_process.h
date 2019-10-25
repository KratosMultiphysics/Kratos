//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#if !defined(KRATOS_SUB_MODEL_PART_SKIN_DETECTION_PROCESS_H_INCLUDED)
#define KRATOS_SUB_MODEL_PART_SKIN_DETECTION_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "skin_detection_process.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
template<SizeType TDim>
class SubModelPartSkinDetectionProcess: public SkinDetectionProcess<TDim>
{
public:
///@name Type Definitions
///@{

/// Pointer definition of SubModelPartSkinDetectionProcess
KRATOS_CLASS_POINTER_DEFINITION(SubModelPartSkinDetectionProcess);

///@}
///@name Life Cycle
///@{

/// Constructor
SubModelPartSkinDetectionProcess(ModelPart& rModelPart, Parameters Settings);

/// Deleted default constructor.
SubModelPartSkinDetectionProcess() = delete;

/// Deleted copy constructor.
SubModelPartSkinDetectionProcess(SubModelPartSkinDetectionProcess const &rOther) = delete;

/// Destructor.
~SubModelPartSkinDetectionProcess() override = default;

///@}
///@name Operators
///@{

/// Deleted sssignment operator.
SubModelPartSkinDetectionProcess &operator=(SubModelPartSkinDetectionProcess const &rOther) = delete;

///@}
///@name Operations
///@{

//void Execute() override;

///@}
///@name Input and output
///@{

std::string Info() const override
{
    return "SkinDetectionProcess";
}

/// Print information about this object.
void PrintInfo(std::ostream& rOStream) const override
{
    rOStream << "SkinDetectionProcess";
}

/// Print object's data.
void PrintData(std::ostream& rOStream) const override
{
}

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

Parameters GetDefaultSettings() const override;

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

///@}
///@name Private Operations
///@{

///@}

}; // Class SubModelPartSkinDetectionProcess

///@}
///@name Input and output
///@{

/// input stream function
template<SizeType TDim>
inline std::istream &operator>>(std::istream &rIStream,
                                SubModelPartSkinDetectionProcess<TDim> &rThis);

/// output stream function
template<SizeType TDim>
inline std::ostream &operator<<(std::ostream &rOStream,
                                const SubModelPartSkinDetectionProcess<TDim> &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@}

} // namespace Kratos.

#endif // KRATOS_SUB_MODEL_PART_SKIN_DETECTION_PROCESS_H_INCLUDED  defined
