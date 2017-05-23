//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   
//

#if !defined(KRATOS_REORDER_AND_OPTIMIZE_MODELPART_PROCESS_H_INCLUDED )
#define  KRATOS_REORDER_AND_OPTIMIZE_MODELPART_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"


namespace Kratos
{
  ///@addtogroup KratosCore
  ///@{

  ///@name Kratos Classes
  ///@{

  /// Short class definition.
  /** Detail class definition.
  */
  class ReorderAndOptimizeModelPartProcess : public Process
    {
    public:
		using GeometryType = Geometry<Point<3> >;
      ///@name Type Definitions
      ///@{

      /// Pointer definition of ReorderAndOptimizeModelPartProcess
      KRATOS_CLASS_POINTER_DEFINITION(ReorderAndOptimizeModelPartProcess);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor is deleted.
      ReorderAndOptimizeModelPartProcess() = delete;

	  /// Constructor to be used. Takes the geometry to be meshed and ModelPart to be filled
	  ReorderAndOptimizeModelPartProcess( ModelPart& rOutputModelPart, Parameters settings);

	  /// The object is not copyable.
	  ReorderAndOptimizeModelPartProcess(ReorderAndOptimizeModelPartProcess const& rOther) = delete;

      /// Destructor.
      virtual ~ReorderAndOptimizeModelPartProcess() ;

      ///@}
      ///@name Operators
      ///@{

	  /// It is not assignable.
	  ReorderAndOptimizeModelPartProcess& operator=(ReorderAndOptimizeModelPartProcess const& rOther) = delete;

      ///@}
      ///@name Operations
      ///@{

	  void Execute() override;


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
      virtual std::string Info() const override;

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const override;

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const override;


      ///@}
      ///@name Friends
      ///@{


      ///@}

      protected:
      ///@name Static Member Variables
      ///@{


      ///@}
      ///@name Member Variables
      ///@{
		  ModelPart& mrModelPart;


      ///@}
      ///@name Private Operations
      ///@{
      void ActualizeSubModelPart(ModelPart& subpart);
      void OptimizeOrdering();
      
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

    }; // Class ReorderAndOptimizeModelPartProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    ReorderAndOptimizeModelPartProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const ReorderAndOptimizeModelPartProcess& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_REORDER_AND_OPTIMIZE_MODELPART_PROCESS_H_INCLUDED  defined
