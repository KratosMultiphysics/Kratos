/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:34 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_GEOMETRICAL_OBJECT_H_INCLUDED )
#define  KRATOS_GEOMETRICAL_OBJECT_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>


// External includes


// Project includes
#include "includes/define.h"


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
class GeometricalObject : public IndexedObject
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GeometricalObject
    KRATOS_CLASS_POINTER_DEFINITION(GeometricalObject);
    
    typedef IndexedObject BaseType;
    
    typedef Node < 3 > NodeType;

    typedef Geometry<NodeType> GeometryType;

    typedef std::size_t IndexType;

    typedef std::size_t result_type;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GeometricalObject(IndexType NewId = 0) : BaseType(NewId),
        mpGeometry()
    {}
    
    /// Default constructor.
    GeometricalObject(IndexType NewId, GeometryType::Pointer pGeometry) : BaseType(NewId),
        mpGeometry(pGeometry)
    {}

    /// Destructor.
    virtual ~GeometricalObject() {}

    /// Copy constructor.
    GeometricalObject(GeometricalObject const& rOther) : BaseType(rOther.Id()),
        mpGeometry(rOther.mpGeometry) 
    {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    GeometricalObject& operator=(GeometricalObject const& rOther)
    {
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{
      
    GeometryType::Pointer pGetGeometry()
    {
        return mpGeometry;
    }

    const GeometryType::Pointer pGetGeometry() const
    {
        return mpGeometry;
    }

    GeometryType& GetGeometry()
    {
        return *mpGeometry;
    }

    GeometryType const& GetGeometry() const
    {
        return *mpGeometry;
    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "geometrical object # "
               << Id();
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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
      
    /**
     * pointer to the condition geometry
     */
    GeometryType::Pointer mpGeometry;
    
    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{
      

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, IndexedObject );
        rSerializer.save("Geometry",mpGeometry);
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, IndexedObject );
        rSerializer.load("Geometry",mpGeometry);
    }


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

}; // Class GeometricalObject

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  GeometricalObject& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const GeometricalObject& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_GEOMETRICAL_OBJECT_H_INCLUDED  defined 


