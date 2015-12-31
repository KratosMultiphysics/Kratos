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



#if !defined(KRATOS_NEIGHBOURS_H_INCLUDED )
#define  KRATOS_NEIGHBOURS_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>


// External includes


// Project includes
#include "includes/define.h"
#include "containers/weak_pointer_vector.h"


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
template<class TNodeType, class TElementType>
class Neighbours
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Neighbours
    KRATOS_CLASS_POINTER_DEFINITION(Neighbours);

    typedef std::size_t IndexType;

    typedef boost::weak_ptr<TNodeType> NodeWeakPointer;

    typedef boost::weak_ptr<TElementType> ElementWeakPointer;

    /** An array of pointers to elements. */
    typedef WeakPointerVector<TElementType> NeighbourElementsArrayType;

    /** An array of pointers to nodes. */
    typedef WeakPointerVector<TNodeType> NeighbourNodesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Neighbours() {}

    /// Copy constructor.
    Neighbours(const Neighbours& rOther)
        : mIndex(rOther.mIndex)
        , mpNeighbourElements(rOther.mpNeighbourElements)
        , mpNeighbourNodes(rOther.mpNeighbourNodes)
    {
    }

    Neighbours(IndexType NewIndex,
               typename NeighbourElementsArrayType::Pointer pNewNeighbourElements,
               typename NeighbourNodesArrayType::Pointer pNewNeighbourNodes)
        : mIndex(NewIndex)
        , mpNeighbourElements(pNewNeighbourElements)
        , mpNeighbourNodes(pNewNeighbourNodes)
    {
    }

    Neighbours(IndexType NewIndex)
        :  mIndex(NewIndex)
    {
    }

    /// Destructor.
    virtual ~Neighbours() {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Neighbours& operator=(const Neighbours& rOther)
    {
        mIndex = rOther.mIndex;
        mpNeighbourElements = rOther.mpNeighbourElements;
        mpNeighbourNodes = rOther.mpNeighbourNodes;

        return *this;
    }


    ///@}
    ///@name Operations
    ///@{

    Neighbours Clone()
    {
        typename NeighbourElementsArrayType::Pointer p_neighbour_elements(new NeighbourElementsArrayType(*mpNeighbourElements));
        typename NeighbourNodesArrayType::Pointer p_neighbour_nodes(new NeighbourNodesArrayType(*mpNeighbourNodes));

        return Neighbours(mIndex, p_neighbour_elements, p_neighbour_nodes);
    }


    ///@}
    ///@name Access
    ///@{

    IndexType Index()
    {
        return mIndex;
    }


    const typename NeighbourElementsArrayType::Pointer pNeighbourElements() const
    {
        return mpNeighbourElements;
    }

    typename NeighbourElementsArrayType::Pointer pNeighbourElements()
    {
        return mpNeighbourElements;
    }

    NeighbourElementsArrayType const& NeighbourElements() const
    {
        return *mpNeighbourElements;
    }

    NeighbourElementsArrayType& NeighbourElements()
    {
        return *mpNeighbourElements;
    }

    const typename NeighbourNodesArrayType::Pointer pNeighbourNodes() const
    {
        return mpNeighbourNodes;
    }

    typename NeighbourNodesArrayType::Pointer pNeighbourNodes()
    {
        return mpNeighbourNodes;
    }

    NeighbourNodesArrayType const& NeighbourNodes() const
    {
        return *mpNeighbourNodes;
    }

    NeighbourNodesArrayType& NeighbourNodes()
    {
        return *mpNeighbourNodes;
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
        return "Neighbours";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "Number of neighbour nodes    : " << mpNeighbourNodes->size() << std::endl;
        rOStream << "Number of neighbour elements : " << mpNeighbourElements->size() << std::endl;
// 	  rOStream << "Neighbour nodes : " <<;
// 	  for(
// 	  rOStream << "Neighbour elements : " << mpNeighbourElements->size() << std::endl;

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

    IndexType mIndex;

    typename NeighbourElementsArrayType::Pointer mpNeighbourElements;

    typename NeighbourNodesArrayType::Pointer mpNeighbourNodes;


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

}; // Class Neighbours

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TNodeType, class TElementType>
inline std::istream& operator >> (std::istream& rIStream,
                                  Neighbours<TNodeType, TElementType>& rThis);

/// output stream function
template<class TNodeType, class TElementType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Neighbours<TNodeType, TElementType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_NEIGHBOURS_H_INCLUDED  defined 


