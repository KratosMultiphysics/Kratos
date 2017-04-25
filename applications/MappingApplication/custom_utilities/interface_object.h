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

/// Short class definition.
/** Detail class definition.
*/
class InterfaceObject : public Point<3>
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

    // TODO hide constructors that are only called from the derived classes
    InterfaceObject() : Point<3>(0.0f, 0.0f, 0.0f)    // Default Constructor /TODO maybe delete zeros
    {
        SetInitialValuesToMembers();
    }

    InterfaceObject(double X, double Y, double Z) : Point<3>(X, Y, Z)   // constuct from coordinates
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
        // xmax, xmin,  ymax, ymin,  zmax, zmin
        bool is_inside = false;

        if (this->X() < *pBoundingBox[0] && this->X() > *pBoundingBox[1])   // check x-direction
        {
            if (this->Y() < *pBoundingBox[2] && this->Y() > *pBoundingBox[3])   // check y-direction
            {
                if (this->Z() < *pBoundingBox[4] && this->Z() > *pBoundingBox[5])   // check z-direction
                {
                    is_inside = true;
                }
            }
        }
        return is_inside;
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

    // These functions have to be duplicated because virtual templates are not possible in C++
    // Scalars
    virtual double GetObjectValue(const Variable<double>& rVariable,
                                  const Kratos::Flags& rOptions)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    virtual void SetObjectValue(const Variable<double>& rVariable,
                                const double& Value,
                                const Kratos::Flags& rOptions,
                                const double Factor)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    virtual double GetObjectValueInterpolated(const Variable<double>& rVariable,
            const std::vector<double>& rShapeFunctionValues)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    // Vectors
    virtual array_1d<double, 3> GetObjectValue(const Variable< array_1d<double, 3> >& rVariable,
            const Kratos::Flags& rOptions)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    virtual void SetObjectValue(const Variable< array_1d<double, 3> >& rVariable,
                                const array_1d<double, 3>& rValue,
                                const Kratos::Flags& rOptions,
                                const double Factor)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    virtual array_1d<double, 3> GetObjectValueInterpolated(const Variable< array_1d<double, 3> >& rVariable,
            const std::vector<double>& rShapeFunctionValues)
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "InterfaceObject" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "InterfaceObject";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


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
