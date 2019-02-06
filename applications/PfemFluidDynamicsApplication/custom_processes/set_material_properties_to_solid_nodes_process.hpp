//
//   Project Name:        KratosPfemFluidApplication $
//   Created by:          $Author:           AFranci $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:        October 2018 $
//   Revision:            $Revision:             0.0 $
//
//

#if !defined(KRATOS_SET_MATERIAL_PROPERTIES_TO_SOLID_NODES_PROCESS_H_INCLUDED )
#define  KRATOS_SET_MATERIAL_PROPERTIES_TO_SOLID_NODES_PROCESS_H_INCLUDED


// System includes

// External includes

// Project includes


#include "spatial_containers/spatial_containers.h"

#include "custom_processes/set_material_properties_to_solid_nodes_process.hpp"
#include "custom_utilities/mesher_utilities.hpp"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "utilities/math_utils.h"
#include "custom_processes/mesher_process.hpp"

///VARIABLES used:
//Data:     
//StepData: 
//Flags:    (checked) 
//          (set)     
//          (modified)  
//          (reset)   


namespace Kratos
{

  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{
  typedef  ModelPart::NodesContainerType                      NodesContainerType;
  typedef  ModelPart::ElementsContainerType                ElementsContainerType;
  typedef  ModelPart::MeshType::GeometryType::PointsArrayType    PointsArrayType;
  typedef std::size_t SizeType;

 
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
  class SetMaterialPropertiesToSolidNodesProcess
    : public MesherProcess
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SetMaterialPropertiesToSolidNodesProcess
    KRATOS_CLASS_POINTER_DEFINITION( SetMaterialPropertiesToSolidNodesProcess );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SetMaterialPropertiesToSolidNodesProcess(ModelPart& rModelPart)
      : mrModelPart(rModelPart)
    {

    }

    /// Destructor.
    virtual ~SetMaterialPropertiesToSolidNodesProcess()
    {
    }

    void operator()()
    {
      Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    void Execute() override
    {
      KRATOS_TRY

	double density = 0;
      double young_modulus = 0;
      double poisson_ratio = 0;
  
#pragma omp parallel
      {
	  
	ModelPart::ElementIterator ElemBegin;
	ModelPart::ElementIterator ElemEnd;
	OpenMPUtils::PartitionedIterators(mrModelPart.Elements(),ElemBegin,ElemEnd);
	for ( ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem )
	  {
	    ModelPart::PropertiesType &elemProperties=itElem->GetProperties();

	    density = elemProperties[DENSITY];
	    young_modulus = elemProperties[YOUNG_MODULUS];
	    poisson_ratio = elemProperties[POISSON_RATIO];
	    
	    Geometry<Node <3> >& rGeom = itElem->GetGeometry();
	    const SizeType NumNodes = rGeom.PointsNumber();
	    for (SizeType i = 0; i < NumNodes; ++i)
	      {
		rGeom[i].FastGetSolutionStepValue(YOUNG_MODULUS)=young_modulus;
		rGeom[i].FastGetSolutionStepValue(DENSITY)=density;
		rGeom[i].FastGetSolutionStepValue(POISSON_RATIO)=poisson_ratio;
	      }

	    
	  }


      }

      
      KRATOS_CATCH(" ")    
	};

    ///@}
    ///@name Operators
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
    std::string Info() const  override
    {
      return "SetMaterialPropertiesToSolidNodesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "SetMaterialPropertiesToSolidNodesProcess";
    }

    


    

  protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

    
    ///@}
    ///@name Protected  Access
    ///@{
    ModelPart& mrModelPart;


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
    SetMaterialPropertiesToSolidNodesProcess& operator=(SetMaterialPropertiesToSolidNodesProcess const& rOther);

    /// Copy constructor.
    //SetMaterialPropertiesToSolidNodesProcess(SetMaterialPropertiesToSolidNodesProcess const& rOther);


    ///@}

  }; // Class SetMaterialPropertiesToSolidNodesProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    SetMaterialPropertiesToSolidNodesProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const SetMaterialPropertiesToSolidNodesProcess& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}


}  // namespace Kratos.

#endif // KRATOS_SET_MATERIAL_PROPERTIES_TO_SOLID_NODES_PROCESS_H_INCLUDED  defined 

