//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_CALCULATE_DISTANCE_TO_BOUNDARY_PROCESS_H_INCLUDED
#define KRATOS_CALCULATE_DISTANCE_TO_BOUNDARY_PROCESS_H_INCLUDED


// System includes


// External includes


// Project includes
#include "processes/process.h"
#include "geometries/geometry.h"
#include "includes/kratos_parameters.h"


namespace Kratos
{
///@addtogroup ShallowWaterApplication
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

/// Forward declaration
class ModelPart;

/**
 * @class CalculateDistanceToBoundaryProcess
 * @ingroup ShallowWaterApplication
 * @brief Calculate the minimum distance from all the nodes to a boundary condition in 2D
 * @details The boundary conditions are assumed to be contained in a line
 * @author Miguel Maso Sotomayor
 */
class KRATOS_API(SHALLOW_WATER_APPLICATION) CalculateDistanceToBoundaryProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CalculateDistanceToBoundaryProcess
    KRATOS_CLASS_POINTER_DEFINITION(CalculateDistanceToBoundaryProcess);

    typedef Geometry<Point> GeometryType;

    typedef Node<3> NodeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CalculateDistanceToBoundaryProcess(
        ModelPart& rComputingModelPart,
        ModelPart& rBoundaryModelPart,
        Parameters ThisParameters = Parameters()) :
            Process(),
            mrModelPart(rComputingModelPart)
    {
        ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
        mRSquaredThreshold = ThisParameters["r_squared_threshold"].GetDouble();
        FindApproximatingGeometry(mpBoundary, rBoundaryModelPart);
    }

    /// Destructor.
    virtual ~CalculateDistanceToBoundaryProcess() {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    int Check() override;

    void ExecuteBeforeSolutionLoop() override;

    const Parameters GetDefaultParameters() const override;

    double GetRSquared();

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
        buffer << "CalculateDistanceToBoundaryProcess";
        return buffer.str();
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

private:
    ///@name Static Member Variables
    ///@{

    ModelPart& mrModelPart;
    GeometryType::Pointer mpBoundary;
    double mRSquaredThreshold;
    double mRSquared;

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void FindApproximatingGeometry(GeometryType::Pointer& pEntity, const ModelPart& rBoundaryPart);

    double RSquared(const GeometryType& rLine, const ModelPart& rBoundaryPart);

    double SquaredDistance(const Point& rPointA, const Point& rPointB);

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
    CalculateDistanceToBoundaryProcess& operator=(CalculateDistanceToBoundaryProcess const& rOther) = delete;

    /// Copy constructor.
    CalculateDistanceToBoundaryProcess(CalculateDistanceToBoundaryProcess const& rOther) = delete;


    ///@}

}; // Class CalculateDistanceToBoundaryProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                CalculateDistanceToBoundaryProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const CalculateDistanceToBoundaryProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CALCULATE_DISTANCE_TO_BOUNDARY_PROCESS_H_INCLUDED  defined
