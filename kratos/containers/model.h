//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Pooyan Dadvand
//


#if !defined(KRATOS_MODEL_H_INCLUDED )
#define  KRATOS_MODEL_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <unordered_map>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

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
  class KRATOS_API(KRATOS_CORE)  Model
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of Model
      KRATOS_CLASS_POINTER_DEFINITION(Model);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      Model(){};

      /// Destructor.
      virtual ~Model(){};


      ///@}
      ///@name Operators
      ///@{
      void AddModelPart( ModelPart::Pointer pmodel_part); //TODO: change this conveniently
      
      ModelPart& GetModelPart(std::string name);

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
      virtual std::string Info() const;

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const;

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const;


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
      std::unordered_map< std::string, ModelPart* > mflat_map; //TODO: deprecate this
      std::unordered_map< std::string, ModelPart* > mroot_map; 


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
//       Model& operator=(Model const& rOther);

      /// Copy constructor.
//       Model(Model const& rOther);


      ///@}

    }; // Class Model

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    Model& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const Model& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MODEL_H_INCLUDED  defined
