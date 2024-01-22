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


#ifndef KRATOS_ESTIMATE_TIME_STEP_UTILITY_H_INCLUDED
#define KRATOS_ESTIMATE_TIME_STEP_UTILITY_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/node.h"
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

/// Forward declaration of ModelPart
class ModelPart;

/// Utility to estimate the time step in terms of the courant number.
/** The velocity can be the sum of the convective velocity and the wave speed
*/
class KRATOS_API(SHALLOW_WATER_APPLICATION) EstimateTimeStepUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EstimateTimeStepUtility
    KRATOS_CLASS_POINTER_DEFINITION(EstimateTimeStepUtility);

    typedef Node NodeType;

    typedef Geometry<NodeType> GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    EstimateTimeStepUtility(ModelPart& rThisModelPart, Parameters ThisParameters);

    /// Destructor.
    ~EstimateTimeStepUtility(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    double Execute() const;

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
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const;


    ///@}
    ///@name Friends
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    const ModelPart& mrModelPart;
    bool mEstimateDt;
    bool mAdaptiveDt;
    double mConstantDt;
    double mCourant;
    double mMinDt;
    double mMaxDt;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    double EstimateTimeStep() const;

    double ElementCharacteristicTime(const GeometryType& rElement, double Gravity) const;

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
    EstimateTimeStepUtility& operator=(EstimateTimeStepUtility const& rOther);

    /// Copy constructor.
    EstimateTimeStepUtility(EstimateTimeStepUtility const& rOther);


    ///@}

}; // Class EstimateTimeStepUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                EstimateTimeStepUtility& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const EstimateTimeStepUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_ESTIMATE_TIME_STEP_UTILITY_H_INCLUDED  defined
