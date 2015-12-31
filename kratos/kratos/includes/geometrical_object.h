// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.




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


