//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//                   Andreas Winterstein
//


#if !defined(KRATOS_CALCULATE_MESH_VELOCITY_UTILITY_H_INCLUDED )
#define  KRATOS_CALCULATE_MESH_VELOCITY_UTILITY_H_INCLUDED


// System includes
#include<map>
#include<tuple>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"


namespace Kratos
{
///@addtogroup MeshMovingApplication
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
 */
class CalculateMeshVelocityUtility
{
public:
    ///@name  Enum's
    ///@{

    enum IntegrationMethod{
        bdf1,
        bdf2,
        generalized_alpha,
        newmark,
        bossak
    };

    ///@}

    ///@name Type Definitions
    ///@{

    /// Pointer definition of CalculateMeshVelocityUtility
    KRATOS_CLASS_POINTER_DEFINITION(CalculateMeshVelocityUtility);

    typedef std::size_t SizeType;

    typedef std::tuple<IntegrationMethod, SizeType> TupleType;

    typedef std::map<std::string, TupleType> MethodsMapType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CalculateMeshVelocityUtility(ModelPart& rModelPart,
                                 Parameters Settings);

    /// Destructor.
    virtual ~CalculateMeshVelocityUtility() = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    CalculateMeshVelocityUtility& operator=(CalculateMeshVelocityUtility const& rOther) = delete;

    /// Copy constructor.
    CalculateMeshVelocityUtility(CalculateMeshVelocityUtility const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void CalculateMeshVelocities();

    static SizeType GetMinimumBufferSize(const std::string& rIntegrationMethod);

    ///@}
    ///@name Access
    ///@{

    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;

    static MethodsMapType msAvailableMethods;

    IntegrationMethod mIntegrationMethod;

    double mbeta;
    double mgamma;

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateMeshVelocitiesBDF(double DeltaTime);

    void CalculateMeshVelocitiesGeneralizedAlpha(double DeltaTime);

    static const TupleType& GetMethod(const std::string& rIntegrationMethod);

    ///@}
    ///@name Private  Access
    ///@{

    ///@}

}; // Class CalculateMeshVelocityUtility

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CALCULATE_MESH_VELOCITY_UTILITY_H_INCLUDED  defined
