//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_INTERFACE_NODE_INCLUDED_H_INCLUDED )
#define  KRATOS_INTERFACE_NODE_INCLUDED_H_INCLUDED

// System includes

// External includes

// Project includes
#include "interface_object.h"


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

/// Node on the Interface for Searching
/** This class Is the "wrapper" for nodes on the interface. It selects the best result by the closest distance to the 
* point of which neighbor have to be found
* Look into the class description of the MapperCommunicator to see how this Object is used in the application
*/
class InterfaceNode : public InterfaceObject
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InterfaceNode
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceNode);

    ///@}
    ///@name Life Cycle
    ///@{

    // A default constructor necessary for serialization 
    InterfaceNode() : InterfaceObject()
    {
    }
    
    InterfaceNode(Node<3>& rNode, const int EchoLevel) : mpNode(&rNode)
    {
        SetCoordinates();
        mEchoLevel = EchoLevel;
    }

    /// Destructor.
    virtual ~InterfaceNode() { }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Node<3>* pGetBase()
    {
        return mpNode;
    }

    bool EvaluateResult(const array_1d<double, 3>& GlobalCooords,
                        double& rMinDistance, const double Distance,
                        std::vector<double>& rShapeFunctionValues) override   // I am an object in the bins
    {
        bool is_closer = false;

        if (Distance < rMinDistance)
        {
            rMinDistance = Distance;
            is_closer = true;
        }

        return is_closer;
    }

    // Functions used for Debugging
    void PrintNeighbors(const int CommRank) override
    {
        array_1d<double, 3> neighbor_coordinates = mpNode->GetValue(NEIGHBOR_COORDINATES);
        double neighbor_comm_rank = mpNode->GetValue(NEIGHBOR_RANK);

        PrintMatchInfo("InterfaceNode", CommRank,
                       neighbor_comm_rank, neighbor_coordinates);
    }

    void WriteRankAndCoordinatesToVariable(const int CommRank) override
    {
        // This function writes the coordinates and the rank of the
        // InterfaceObject to the variables "NEIGHBOR_COORDINATES"
        // and "NEIGHBOR_RANK", for debugging
        array_1d<double, 3> neighbor_coordinates;
        // TODO exchange with "Coordinates()"
        neighbor_coordinates[0] = this->X();
        neighbor_coordinates[1] = this->Y();
        neighbor_coordinates[2] = this->Z();
        mpNode->SetValue(NEIGHBOR_COORDINATES, neighbor_coordinates);
        mpNode->SetValue(NEIGHBOR_RANK, CommRank);
    }

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
        std::stringstream buffer;
        buffer << "InterfaceNode" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "InterfaceNode";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {}


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

    Node<3>* mpNode;
        
    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;
    
    virtual void save(Serializer& rSerializer) const 
    {
        KRATOS_ERROR << "This object is not supposed to be used with serialization!" << std::endl;        
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, InterfaceObject);
    }
    virtual void load(Serializer& rSerializer) 
    {
        KRATOS_ERROR << "This object is not supposed to be used with serialization!" << std::endl;
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, InterfaceObject);
    }

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void SetCoordinates() override
    {
        this->Coordinates() = mpNode->Coordinates();
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

    /// Assignment operator.
    InterfaceNode& operator=(InterfaceNode const& rOther);

    //   /// Copy constructor.
    //   InterfaceNode(InterfaceNode const& rOther){}


    ///@}

}; // Class InterfaceNode

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  InterfaceNode& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const InterfaceNode& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_NODE_INCLUDED_H_INCLUDED  defined
