//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:     Ruben Zorrilla
//

#if !defined(KRATOS_ACCELERATION_LIMITATION_UTILITIES_H_INCLUDED )
#define  KRATOS_ACCELERATION_LIMITATION_UTILITIES_H_INCLUDED

// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{
  ///@addtogroup FluidDynamicsApplication
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

  class KRATOS_API(FLUID_DYNAMICS_APPLICATION) AccelerationLimitationUtilities
  {
  public:

    ///@name Type Definitions
    ///@{

    typedef Geometry<Node<3>>                                 GeometryType;
    typedef IntegrationPoint<3>                       IntegrationPointType;
    typedef std::vector<IntegrationPointType>   IntegrationPointsArrayType;

    /// Pointer definition of AccelerationLimitationUtilities
    KRATOS_CLASS_POINTER_DEFINITION(AccelerationLimitationUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    /**
     * @brief Construct a new Acceleration Limitation Utilities object
     *
     * @param ModelPart model part to be controlled
     * @param multipleOfG specification of maximal acceleration in multiples of the gravitational acceleration
     */
    AccelerationLimitationUtilities( ModelPart &ModelPart, double multipleOfG ) : mrModelPart(ModelPart) {

        this->mMaximalAccelaration = multipleOfG;
    };

    /// Destructor.
    ~AccelerationLimitationUtilities() {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Execution based on already specified parameters
     * The velocity at every node is reduced based on a physically possible acceleration
     */
    void Execute();


    /**
     * @brief Set the Limit As Multiple Of Gravitional Acceleration object
     *
     * @param newMaxAcc new maximal acceleration
     */
    void SetLimitAsMultipleOfGravitionalAcceleration( double& newMaxAcc );


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

    ModelPart &mrModelPart;
    // The utilities model part saved as a reference

    double mMaximalAccelaration = 3.0;
    // given in multipples of gravitational acceleration

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
    AccelerationLimitationUtilities& operator=(AccelerationLimitationUtilities const& rOther);

    /// Copy constructor.
    AccelerationLimitationUtilities(AccelerationLimitationUtilities const& rOther);

    ///@}

}; // Class AccelerationLimitationUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const AccelerationLimitationUtilities& rThis);

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_ACCELERATION_LIMITATION_UTILITIES_H_INCLUDED  defined
