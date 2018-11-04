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

#if !defined(KRATOS_INTERFACE_OBJECT_INCLUDED_H_INCLUDED )
#define  KRATOS_INTERFACE_OBJECT_INCLUDED_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "mapper_utilities.h"
#include "../mapping_application_variables.h"


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

/// Base Class for Searching Objects
/** This class provides a set of functions that is used for identifying the nearest object.
* It is needed such that the bin-search can be used with both nodes and elements/conditions
* The bin search is implemented to work with this kind of object
* It implements the function "EvaluateResult", which is used by the local search to determine which
* of the objects in teh vicinity of a point is the best search result. This function has to be
* implemented by all subclasses. It cannot be made pure virtual, because for remote searching,
* objects of this type have to be created
* Look into the class description of the MapperCommunicator to see how this Object is used in the application
*/
class InterfaceObject : public Point
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InterfaceObject
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceObject);

    ///@}
    ///@name  Enum's
    ///@{

    enum PairingStatus
    {
        NoNeighbor = 0,
        Approximation = 1,
        NeighborFound = 2,
    };

    ///@}
    ///@name Life Cycle
    ///@{

    InterfaceObject(double X, double Y, double Z) : Point(X, Y, Z)   // constuct from coordinates
    {
        SetInitialValuesToMembers();
    }

    /// Destructor.
    virtual ~InterfaceObject() { }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Reset()
    {
        SetInitialValuesToMembers();
        SetCoordinates();
    }

    bool IsInBoundingBox(double* pBoundingBox[])
    {
        return MapperUtilities::PointIsInsideBoundingBox(*pBoundingBox, this->Coordinates());;
    }

    void ProcessSearchResult(const double Distance, const int PairingStatus, const int Rank)
    {
        if (mPairingStatus < PairingStatus || (mPairingStatus == PairingStatus
                                               && Distance < mMinDistanceNeighbor))
        {
            mPairingStatus = PairingStatus;
            mMinDistanceNeighbor = Distance;
            mNeighborRank = Rank;
        }
    }

    int GetPairingStatus()
    {
        return mPairingStatus;
    }

    bool NeighborOrApproximationFound()
    {
        if (mPairingStatus >= PairingStatus::Approximation)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    bool HasNeighborOrApproximationInPartition(const int PartitionIndex)
    {
        bool return_value = false;
        if (mPairingStatus >= PairingStatus::Approximation)
        {
            if (mNeighborRank == PartitionIndex)
                return_value = true;
        }
        return return_value;
    }

    int GetNeighborRank()
    {
        return mNeighborRank;
    }

    void SetIsBeingSent()
    {
        mIsBeingSent = true;
    }

    bool GetIsBeingSent()
    {
        return mIsBeingSent;
    }

    virtual Node<3>* pGetBaseNode()
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
        return nullptr;
    }

    virtual Geometry<Node<3>>* pGetBaseGeometry()
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
        return nullptr;
    }

    virtual bool EvaluateResult(const array_1d<double, 3>& rGlobalCoords,
                                double& rMinDistance, const double Distance,
                                std::vector<double>& rShapeFunctionValues)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
        return false;
    }

    virtual bool ComputeApproximation(const array_1d<double, 3>& rGlobalCoords, double& rMinDistance,
                                      std::vector<double>& rShapeFunctionValues)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
        return false;
    }

    // Functions used for Debugging
    virtual void PrintNeighbors(const int CommRank)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    virtual void WriteRankAndCoordinatesToVariable(const int CommRank)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
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
        buffer << "InterfaceObject" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "InterfaceObject";
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

    int mEchoLevel = 0;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    // This constructor is called by its derived classes
    InterfaceObject() : Point(0.0f, 0.0f, 0.0f)
    {
        SetInitialValuesToMembers();
    }

    void SetInitialValuesToMembers()
    {
        mMinDistanceNeighbor = std::numeric_limits<double>::max();
        mPairingStatus = PairingStatus::NoNeighbor;
        mNeighborRank = 0;
        mIsBeingSent = false;
    }

    virtual void SetCoordinates()
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    void PrintMatchInfo(const std::string& rInterfaceObjectType,
                        const int CommRank, const int NeighborCommRank,
                        array_1d<double, 3>& rNeighborCoordinates)
    {

        std::cout << rInterfaceObjectType << " ["
                  << this->X() << " "
                  << this->Y() << " "
                  << this->Z() << "], "
                  << "Rank " << CommRank
                  << " || Neighbor ["
                  << rNeighborCoordinates[0] << " "
                  << rNeighborCoordinates[1] << " "
                  << rNeighborCoordinates[2] << "], Rank "
                  << NeighborCommRank << std::endl;
    }


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

    double mMinDistanceNeighbor;
    int mPairingStatus; // 0 no Neighbor found; 1 approximation (i.e. nearest Node found); 2 match found
    int mNeighborRank;
    bool mIsBeingSent;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_ERROR << "This object is not supposed to be used with serialization!" << std::endl;
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Point);
    }
    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_ERROR << "This object is not supposed to be used with serialization!" << std::endl;
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Point);
    }

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
    InterfaceObject& operator=(InterfaceObject const& rOther);

    //   /// Copy constructor.
    //   InterfaceObject(InterfaceObject const& rOther){}


    ///@}

}; // Class InterfaceObject

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  InterfaceObject& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const InterfaceObject& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_OBJECT_INCLUDED_H_INCLUDED  defined
