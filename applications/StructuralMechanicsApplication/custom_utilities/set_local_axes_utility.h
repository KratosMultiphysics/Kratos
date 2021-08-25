// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                         license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

#if !defined(KRATOS_SET_LOCAL_AXES_UTILITY_H_INCLUDED )
#define  KRATOS_SET_LOCAL_AXES_UTILITY_H_INCLUDED

// System includes
#include "includes/model_part.h"

// Include kratos definitions

// Project includes

// Configures

// External includes

namespace Kratos {
///@addtogroup StructuralMechanicsApplication
///@{

///@name Kratos Classes
///@{

/**
 * @class SetLocalAxesUtility
 * @ingroup StructuralMechanicsApplication
 * @brief setting the local axes of the elements in a modelpart
 * according to a certain criterion
 * @details This class provides a method for computing and setting local
 * axes to elements
 * @author Alejandro Cornejo
 */

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SetLocalAxesUtility
{
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of NodeSearchUtility
      KRATOS_CLASS_POINTER_DEFINITION(SetLocalAxesUtility);


      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      SetLocalAxesUtility() {}

      /// Destructor.
      ~SetLocalAxesUtility(){
      }

      ///@}
      ///@name Operations
      ///@{

      /**
       * @brief Perform Node Search
       * @details Searches for nodes within a given 'Radius' of the current node
       * @param rStructureNodes Nodes Container
       * @param Id NodeId of current node.
       * @param Radius Search radius.
       * @param rResults Results container.
       **/
      void SetLocalAxisCartesianSystem(
          ModelPart &rModelPart,
          Parameters ThisParameters);


      void CheckAndNormalizeVector(
        BoundedVector<double, 3>& rVector)
      {
        const double norm = MathUtils<double>::Norm3(rVector);
        if (norm > std::numeric_limits<double>::epsilon()) {
          rVector /= norm;
        } else {
          KRATOS_ERROR << "The norm of one LOCAL_AXIS is null" << std::endl;
        }
      }

      ///@}
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
      {
          std::stringstream buffer;
          buffer << "SetLocalAxesUtility" ;

          return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const  {rOStream << "SetLocalAxesUtility";}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const  {}

      ///@}

    private:
      ///@name Member Variables
      ///@{

      ///@}
      ///@name Un accessible methods
      ///@{

      /// Assignment operator.
      SetLocalAxesUtility& operator=(SetLocalAxesUtility const& rOther)
      {
          return *this;
      }

      /// Copy constructor.
      SetLocalAxesUtility(SetLocalAxesUtility const& rOther)
      {
          *this = rOther;
      }

      ///@}

    }; // Class SetLocalAxesUtility

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif