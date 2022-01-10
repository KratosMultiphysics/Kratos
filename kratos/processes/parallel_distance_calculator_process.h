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

    const Parameters GetDefaultParameters() const override;

    int Check() override;

    void Execute() override;

    double FindMaximumEdgeSize();

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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "ParallelDistanceCalculatorProcess" << TDim << "D";
        return buffer.str();
    };

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ParallelDistanceCalculatorProcess" << TDim << "D";
    };

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {};

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
