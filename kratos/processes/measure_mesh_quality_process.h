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

#if !defined(KRATOS_MEASURE_MESH_QUALITY_PROCESS_H_INCLUDED )
#define  KRATOS_MEASURE_MESH_QUALITY_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"



namespace Kratos
{
  ///@addtogroup Kratos
  ///@{
  ///@name Kratos Classes
  ///@{

  /// Short class definition.
  /** Detail class definition.
  */
  class MeasureMeshQualityProcess : public Process
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of MeasureMeshQualityProcess
      KRATOS_CLASS_POINTER_DEFINITION(MeasureMeshQualityProcess);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Constructor takes the modelpart to apply smoothing to its mesh 0.
      MeasureMeshQualityProcess(ModelPart& rModelPart, std::size_t Dimension);

      /// Destructor.
      ~MeasureMeshQualityProcess() override;


      ///@}
      ///@name Operators
      ///@{


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
      std::string Info() const override;

      /// Print information about this object.
      void PrintInfo(std::ostream& rOStream) const override;

      /// Print object's data.
      void PrintData(std::ostream& rOStream) const override;


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

		ModelPart& mrModelPart;
		std::size_t mDimension;

		std::size_t mNumberOfInvertedElements;
		std::size_t mNumberOfSlivers;
		double mMinArea;
		double mMaxArea;
		double mMinH;
		double mMinAngle;


      ///@}
      ///@name Member Variables
      ///@{


      ///@}
      ///@name Private Operators
      ///@{


      ///@}
      ///@name Private Operations
      ///@{

		void ResetMeasures();


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
      MeasureMeshQualityProcess& operator=(MeasureMeshQualityProcess const& rOther);

      /// Copy constructor.
      MeasureMeshQualityProcess(MeasureMeshQualityProcess const& rOther);


      ///@}

    }; // Class MeasureMeshQualityProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    MeasureMeshQualityProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const MeasureMeshQualityProcess& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MEASURE_MESH_QUALITY_PROCESS_H_INCLUDED  defined
