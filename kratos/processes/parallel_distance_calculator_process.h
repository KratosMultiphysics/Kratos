//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

#if !defined(KRATOS_PARALLEL_DISTANCE_CALCULATOR_H_INCLUDED )
#define  KRATOS_PARALLEL_DISTANCE_CALCULATOR_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "process.h"

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
template <unsigned int TDim>
class KRATOS_API(KRATOS_CORE) ParallelDistanceCalculatorProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;

    /// Pointer definition of ParallelDistanceCalculatorProcess
    KRATOS_CLASS_POINTER_DEFINITION(ParallelDistanceCalculatorProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    ParallelDistanceCalculatorProcess(
        ModelPart &rModelPart,
        Parameters Settings);

    ParallelDistanceCalculatorProcess(
        Model &rModel,
        Parameters Settings);

    ParallelDistanceCalculatorProcess(
        ModelPart& rModelPart,
        const Variable<double>& rDistanceVar,
        const Variable<double>& rAreaVar,
        const unsigned int MaxLevels,
        const double MaxDistance,
        const bool CalculateExactDistancesToPlane = false);

    /// Copy constructor.
    ParallelDistanceCalculatorProcess(ParallelDistanceCalculatorProcess<TDim> const& rOther) = delete;

    /// Destructor.
    virtual ~ParallelDistanceCalculatorProcess() {};

    ///Function to calculate a signed distance function suitable for calculations using the Level Set Method
    ///the function assumes given a "signed distance" distributions and recomputes the distances
    ///respecting as accurately as possible the position of the zero of the original distributions
    ///@param rModelPart is the ModelPart on which we will operate
    ///@param rDistanceVar is the Variable that we will use in calculating the distance
    ///@param rAreaVar is the Variable that we will use for L2 projections
    ///@param MaxLevels is the number of maximum "layers" of element that will be used in the calculation of the distances
    ///@param MaxDistance distances will not be computed after reaching this limit
    // void CalculateDistances(
    //     ModelPart& rModelPart,
    //     const Variable<double>& rDistanceVar,
    //     const Variable<double>& rAreaVar,
    //     const unsigned int MaxLevels,
    //     const double MaxDistance,
    //     Flags Options = CALCULATE_EXACT_DISTANCES_TO_PLANE.AsFalse());

    ///Function to calculate a signed distance function suitable for calculations using the Level Set Method
	///The difference of this function with previous one is the fact that it wont recalculate the exact distance
	///in divided elements in order to preserve the current distance.
    ///the function assumes given a "signed distance" distributions and recomputes the distances
    ///respecting as accurately as possible the position of the zero of the original distributions
    ///@param rModelPart is the ModelPart on which we will operate
    ///@param rDistanceVar is the Variable that we will use in calculating the distance
    ///@param rAreaVar is the Variable that we will use for L2 projections
    ///@param MaxLevels is the number of maximum "layers" of element that will be used in the calculation of the distances
    ///@param MaxDistance distances will not be computed after reaching this limit
    // void CalculateInterfacePreservingDistances(
    //     ModelPart& rModelPart,
    //     const Variable<double>& rDistanceVar,
    //     const Variable<double>& rAreaVar,
    //     const unsigned int MaxLevels,
    //     const double MaxDistance);

    /// A simplified version of CalculateDistances to be used when the rDistanceVar == 0 surface is described by a set of nodes
    /**
     * @param rModelPart is the ModelPart on which we will operate
     * @param rDistanceVar is the Variable that we will use in calculating the distance
     * @param rAreaVar is the Variable that we will use for L2 projections
     * @param max_levels is the number of maximum "layers" of element that will be used in the calculation of the distances
     * @param max_distance distances will not be computed after reaching this limit
     * @see ParallelDistanceCalculatorProcess::CalculateDistances
     */
    // void CalculateDistancesLagrangianSurface(
    //     ModelPart& rModelPart,
    //     const Variable<double>& rDistanceVar,
    //     const Variable<double>& rAreaVar,
    //     const unsigned int MaxLevels,
    //     const double MaxDistance);

    // double FindMaximumEdgeSize(ModelPart& rModelPart);

    const Parameters GetDefaultParameters() const override;

    int Check() override;

    void Execute() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


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
        buffer << "ParallelDistanceCalculatorProcess" << TDim << "D";
        return buffer.str();
    };

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ParallelDistanceCalculatorProcess" << TDim << "D";
    };

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {};

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

        ModelPart& mrModelPart;
        const Variable<double>* mpDistanceVar;
        const Variable<double>* mpAreaVar;
        unsigned int mMaxLevels;
        double mMaxDistance;
        bool mCalculateExactDistancesToPlane;


    ///@}
    ///@name Protected Operators
    ///@{

    bool IsDivided(const array_1d<double,TDim+1>& rDistance);

    bool IsActive(const array_1d<double,TDim+1>& rVisited);

    void AddDistanceToNodes(
        Geometry<NodeType>& rGeometry,
        const BoundedMatrix<double,TDim+1,TDim>& rDN_DX,
        const double& Volume);

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


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void ResetVariables();

    void CalculateExactDistancesOnDividedElements();

	void ExtendDistancesByLayer();

    void AssignDistanceSign();

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
}; // Class ParallelDistanceCalculatorProcess

///@}
///@name Input and output
///@{

/// input stream function
template<unsigned int TDim>
inline std::istream& operator >> (
    std::istream& rIStream,
    ParallelDistanceCalculatorProcess<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template<unsigned int TDim>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const ParallelDistanceCalculatorProcess<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_PARALLEL_DISTANCE_CALCULATOR_H_INCLUDED  defined
