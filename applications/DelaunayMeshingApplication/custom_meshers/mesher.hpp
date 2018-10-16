//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#if !defined(KRATOS_MESHER_H_INCLUDED )
#define  KRATOS_MESHER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "containers/variables_list_data_value_container.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/openmp_utils.h"
#include "custom_utilities/mesher_utilities.hpp"
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
class KRATOS_API(DELAUNAY_MESHING_APPLICATION) Mesher
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of Mesher
    KRATOS_CLASS_POINTER_DEFINITION( Mesher );

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    typedef MesherUtilities::MeshingInfoParameters  InfoParametersType;
    typedef MesherUtilities::MeshingParameters   MeshingParametersType;
    typedef MesherUtilities::RefiningParameters   RefineParametersType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Mesher() {}

    /// Copy constructor.
    Mesher(Mesher const& rOther)
      :mpMeshingVariables(rOther.mpMeshingVariables)
      ,mPreMeshingProcesses(rOther.mPreMeshingProcesses)
      ,mPostMeshingProcesses(rOther.mPostMeshingProcesses)
      ,mpMesherUtilities(rOther.mpMesherUtilities)
      ,mpDataTransferUtilities(rOther.mpDataTransferUtilities)
      ,mEchoLevel(rOther.mEchoLevel)
    {};

    /// Destructor.
    virtual ~Mesher() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    //************* STARTING METHODS

    /**
     * Called to initialize the mesher
     * Must be called before any calculation is done
     */
    void Initialize();


    /**
     * level of echo for the mesh mesher
     */
    virtual void SetEchoLevel(int Level)
    {
        mEchoLevel = Level;
    }

    int GetEchoLevel()
    {
        return mEchoLevel;
    }

    /**
     * Remesh information is given to the mesher
     */
    void SetMeshingParameters( MeshingParametersType::Pointer& rMeshingParameters );

    /**
     * Pre and Post meshing processes are given to the mesher
     */
    void SetPreMeshingProcess( MesherProcess::Pointer pPreMeshingProcess );

    void SetPostMeshingProcess( MesherProcess::Pointer pPostMeshingProcess );


    void SetPreMeshingProcessVector( std::vector<MesherProcess::Pointer>& rPreMeshingProcessVector );

    void SetPostMeshingProcessVector( std::vector<MesherProcess::Pointer>& rPostMeshingProcessVector );


    /**
     * Utilities for the mesher are given to the mesher
     */
    void SetMesherUtilities( MesherUtilities::Pointer rMesherUtilities );

    /**
     * Transfer utilities are given to the mesher
     */
    void SetDataTransferUtilities( MeshDataTransferUtilities::Pointer rDataTransferUtilities );


    //*******************************************************************************************
    //*******************************************************************************************

    /**
     * Mesher :: Initilize
     */
    virtual void InitializeMesher(ModelPart& rModelPart);


    /**
     * Mesh Generation :: Remesh all ModelPart
     */
    virtual void ExecuteMeshing(ModelPart& rModelPart);

    /**
     * Mesher :: Finalize
     */
    virtual void FinalizeMesher(ModelPart& rModelPart);

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
    virtual std::string Info() const
    {
        return "Mesher";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


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

    MeshingParametersType::Pointer          mpMeshingVariables;

    std::vector<MesherProcess::Pointer>   mPreMeshingProcesses;
    std::vector<MesherProcess::Pointer>  mPostMeshingProcesses;

    MesherUtilities::Pointer                 mpMesherUtilities;

    MeshDataTransferUtilities::Pointer mpDataTransferUtilities;

    int mEchoLevel;

    ///@}
    ///@name Protected Operators
    ///@{

    /// Assignment operator.
    Mesher& operator=(Mesher const& rOther);

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * Mesher :: Start Echo
     */
    virtual void StartEcho(ModelPart& rSubModelPart,
			   std::string GenerationMessage);


    /**
     * Mesher :: End Echo
     */
    virtual void EndEcho(ModelPart& rSubModelPart,
			 std::string GenerationMessage);



    /**
     * Mesher :: Set Nodes to mesh
     */
    virtual void SetNodes(ModelPart& rModelPart,
			  MeshingParametersType& rMeshingVariables);

    /**
     * Mesher :: Set Elements to mesh
     */
    virtual void SetElements(ModelPart& rModelPart,
			     MeshingParametersType& rMeshingVariables);

    /**
     * Mesher :: Set Elements to mesh
     */
    virtual void SetNeighbours(ModelPart& rModelPart,
			       MeshingParametersType& rMeshingVariables);


    /**
     * Mesher :: Process to be done at the begining of the Generation
     */
    virtual void ExecutePreMeshingProcesses();


    /**
     * Mesher :: Process to be done at the end of the Generation
     */
    virtual void ExecutePostMeshingProcesses();



    /**
     * Mesher :: Delaunay Tessellation
     */
    virtual void Generate (ModelPart& rModelPart,
			   MeshingParametersType& rMeshingVariables){};


    /**
     * Mesher :: Set Element Neighbours
     */
    virtual void SetElementNeighbours(ModelPart& rModelPart,
				      MeshingParametersType& rMeshingVariables);

    /**
     * Mesher :: Recover Boundary Position
     */
    virtual void RecoverBoundaryPosition(ModelPart& rModelPart,
					 MeshingParametersType& rMeshingVariables);


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

}; // Class Mesher

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  Mesher& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Mesher& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_MESHER_H_INCLUDED  defined


