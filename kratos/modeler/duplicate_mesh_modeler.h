//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                    
//

#if !defined(KRATOS_COPY_MODELER_H_INCLUDED )
#define  KRATOS_COPY_MODELER_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "modeler/modeler.h"

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

/// Short class definition.

/** Detail class definition.
 */
class KRATOS_API(KRATOS_CORE) DuplicateMeshModeler : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DuplicateMeshModeler
    KRATOS_CLASS_POINTER_DEFINITION(DuplicateMeshModeler);

    typedef Modeler BaseType;

    typedef Point PointType;

    typedef Node NodeType;

    typedef Geometry<NodeType> GeometryType;

    typedef PointerVector<NodeType> NodesVectorType;

    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// constructor.

    DuplicateMeshModeler(ModelPart& rSourceModelPart) :
        mrModelPart(rSourceModelPart)
    {
    }

    /// Destructor.

    virtual ~DuplicateMeshModeler()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void GenerateMesh(ModelPart& rThisModelPart, Element const& rReferenceElement, Condition const& rReferenceCondition) override;

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

    virtual std::string Info() const override
    {
        return "DuplicateMeshModeler";
    }

    /// Print information about this object.

    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.

    virtual void PrintData(std::ostream& rOStream) const override
    {
    }


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

    void GenerateSubModelParts(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelpart);


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

    ModelPart& mrModelPart;

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


    ///@}

}; // Class DuplicateMeshModeler

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream & operator >>(std::istream& rIStream,
                                  DuplicateMeshModeler& rThis);

/// output stream function

inline std::ostream & operator <<(std::ostream& rOStream,
                                  const DuplicateMeshModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


} // namespace Kratos.

#endif // KRATOS_COPY_MODELER_H_INCLUDED  defined


