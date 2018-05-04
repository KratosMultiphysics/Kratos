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

#if !defined(KRATOS_INTERFACE_MESHLESS_POINT_INCLUDED_H_INCLUDED )
#define  KRATOS_INTERFACE_MESHLESS_POINT_INCLUDED_H_INCLUDED

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
class InterfaceMeshlessPoint : public InterfaceObject
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InterfaceMeshlessPoint
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceMeshlessPoint);

    ///@}
    ///@name Life Cycle
    ///@{

    // A default constructor necessary for serialization
    InterfaceMeshlessPoint() : InterfaceObject()
    {
    }

    InterfaceMeshlessPoint(Element& rElement, const int EchoLevel) : mpElement(&rElement)
    {
        SetCoordinates();
        mEchoLevel = EchoLevel;
    }

    /// Destructor.
    virtual ~InterfaceMeshlessPoint() { }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element* pGetBaseElement() override
    {
        return mpElement;
    }

    bool EvaluateResult(const InterfaceObject::Pointer rObject,
                        double& rMinDistance, const double Distance,
                        std::vector<double>& rShapeFunctionValues) override   // I am an object in the bins
    {
        (mpElement->GetValue(CONTACT_ELEMENTS)).push_back(rObject->pGetBaseElement());

        return false;
    }

    bool ComputeApproximation(const array_1d<double, 3>& rGlobalCoords, double& rMinDistance,
                                      std::vector<double>& rShapeFunctionValues) override
    {
        // do nothing
        return false;
    }

    // Functions used for Debugging
    void PrintNeighbors(const int CommRank) override
    {
        array_1d<double, 3> neighbor_coordinates = mpElement->GetValue(NEIGHBOR_COORDINATES);
        double neighbor_comm_rank = mpElement->GetValue(NEIGHBOR_RANK);

        PrintMatchInfo("InterfaceMeshlessPoint", CommRank,
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
        mpElement->SetValue(NEIGHBOR_COORDINATES, neighbor_coordinates);
        mpElement->SetValue(NEIGHBOR_RANK, CommRank);
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
        buffer << "InterfaceMeshlessPoint" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "InterfaceMeshlessPoint";
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

    Element* mpElement;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_ERROR << "This object is not supposed to be used with serialization!" << std::endl;
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, InterfaceObject);
    }
    virtual void load(Serializer& rSerializer) override
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
        // std::cout << *mpElement << std::endl;
        // std::cout << SHAPE_FUNCTION_VALUES << std::endl;
        const auto& r_geom = mpElement->GetGeometry();

        // std::cout << r_geom << std::endl;

        KRATOS_ERROR_IF_NOT(mpElement->Has(INTEGRATION_WEIGHT)) << "No INTEGRATION_WEIGHT available!"
            << std::endl;
        KRATOS_ERROR_IF_NOT(mpElement->Has(SHAPE_FUNCTION_VALUES)) << "No SHAPE_FUNCTION_VALUES available!"
            << std::endl;

        const auto& r_shape_function_values = mpElement->GetValue(SHAPE_FUNCTION_VALUES);

        array_1d<double,3> coords;

        for (int i=0;i<3;++i)
            coords[i] = 0.0;

        for (std::size_t i=0; i<r_geom.PointsNumber(); ++i)
        {
            coords[0] += r_geom[i].X() * r_shape_function_values[i];
            coords[1] += r_geom[i].Y() * r_shape_function_values[i];
            coords[2] += r_geom[i].Z() * r_shape_function_values[i];
        }

        noalias(this->Coordinates()) = coords;
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
    InterfaceMeshlessPoint& operator=(InterfaceMeshlessPoint const& rOther);

    //   /// Copy constructor.
    //   InterfaceMeshlessPoint(InterfaceMeshlessPoint const& rOther){}


    ///@}

}; // Class InterfaceMeshlessPoint

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  InterfaceMeshlessPoint& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const InterfaceMeshlessPoint& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_MESHLESS_POINT_INCLUDED_H_INCLUDED  defined
