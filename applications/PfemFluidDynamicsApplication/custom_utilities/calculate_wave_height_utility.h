//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics PfemFluidDynamics Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//                   Miguel Maso Sotomayor
//

#ifndef KRATOS_CALCULATE_WAVE_HEIGHT_UTILITY_H_INCLUDED
#define KRATOS_CALCULATE_WAVE_HEIGHT_UTILITY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/** 
* @ingroup PfemFluidDynamicsApplication
* @author Miguel Maso Sotomayor
* @brief This function computes the wave height at a given point
* @details The direction is taken from the gravity variable in the ProcessInfo
*/
class KRATOS_API(PFEM_FLUID_DYNAMICS_APPLICATION) CalculateWaveHeightUtility
{
public:
    ///@name Pointer Definition
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(CalculateWaveHeightUtility);

    ///@}
    ///@name Type Definitions
    ///@{

    typedef Node NodeType;

    typedef Geometry<NodeType> GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor with ModelPart and Parameters.
     */
    CalculateWaveHeightUtility(ModelPart &rThisModelPart, Parameters ThisParameters);

    /**
     * @brief Destructor.
     */
    ~CalculateWaveHeightUtility() {}

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Calculate the wave height.
     */
    double Calculate(const array_1d<double,3>& rCoordinates) const;

    ///@}
    ///@name Input and output
    ///@{

    ///@brief Turn back information as a string.
    std::string Info() const
    {
      return "CalculateWaveHeightUtility";
    }

    ///@brief Print information about this object.
    void PrintInfo(std::ostream &rOStream) const
    {
      rOStream << "CalculateWaveHeightUtility";
    }

    ///@brief Print information about this object.
    void PrintData(std::ostream &rOStream) const
    {
    }

    ///@}

protected:

    ///@name Protected LifeCycle
    ///@{

    // /// Copy constructor.
    // CalculateWaveHeightUtility(CalculateWaveHeightUtility const &rOther);

    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}

private:

    ///@name Member Variables
    ///@{

    ModelPart &mrModelPart;
    array_1d<double,3> mDirection;
    array_1d<double,3> mCoordinates;
    double mMeanWaterLevel;
    bool mUseLocalElementSize;
    bool mUseNearestNode;
    double mRelativeRadius;
    double mAbsoluteRadius;
    double mMeanElementSize;

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Calculate the averaged wave height from the nearest nodes with a given search radius.
     */
    double CalculateAverage(const array_1d<double,3>& rCoordinates, double SearchRadius) const;


    /**
     * @brief Calculate the averaged wave height from the nearest nodes.
     */
    double CalculateAverage(const array_1d<double,3>& rCoordinates) const;

    /**
     * @brief Calculate the wave height from the nearest node.
     */
    double CalculateNearest(const array_1d<double,3>& rCoordinates) const;

    ///@}
    ///@name Private  Access
    ///@{

    // /// Assignment operator.
    // CalculateWaveHeightUtility &operator=(CalculateWaveHeightUtility const &rOther);

    ///@}
    ///@name Serialization
    ///@{


    ///@}

}; // Class CalculateWaveHeightUtility

///@}
///@name Input and output
///@{

///@brief output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const CalculateWaveHeightUtility &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

} // namespace Kratos.

#endif // KRATOS_CALCULATE_WAVE_HEIGHT_UTILITY_H_INCLUDED  defined
