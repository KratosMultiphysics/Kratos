//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_MESH_MODELER_H_INCLUDED )
#define  KRATOS_MESH_MODELER_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <stdlib.h>

#include <boost/timer.hpp>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/model_part.h"
#include "containers/variables_list_data_value_container.h"
#include "spatial_containers/spatial_containers.h"

#include "processes/find_nodal_h_process.h"

#include "custom_processes/elemental_neighbours_search_process.hpp"
#include "custom_processes/nodal_neighbours_search_process.hpp"
#include "custom_processes/build_mesh_boundary_process.hpp"

#include "custom_utilities/boundary_normals_calculation_utilities.hpp"
#include "custom_utilities/mesh_error_calculation_utilities.hpp"
#include "custom_utilities/change_tip_elements_utilities.hpp"
#include "custom_utilities/modeler_utilities.hpp"


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
class MeshModeler
{
public:


    ///@name Type Definitions
    ///@{

    /// Pointer definition of MeshModeler
    KRATOS_CLASS_POINTER_DEFINITION( MeshModeler );

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;
  
    typedef ModelerUtilities::InfoParameters         InfoParametersType;
    typedef ModelerUtilities::MeshingParameters   MeshingParametersType;
    typedef ModelerUtilities::RefiningParameters   RefineParametersType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MeshModeler() {}

    /// Copy constructor.
    MeshModeler(MeshModeler const& rOther)
      :mpMeshingVariables(rOther.mpMeshingVariables)
      ,mPreMeshingProcesses(rOther.mPreMeshingProcesses)
      ,mPostMeshingProcesses(rOther.mPostMeshingProcesses)
      ,mpModelerUtilities(rOther.mpModelerUtilities)
      ,mpDataTransferUtilities(rOther.mpDataTransferUtilities)
      ,mEchoLevel(rOther.mEchoLevel)
    {};

    /// Destructor.
    virtual ~MeshModeler() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    //************* STARTING METHODS

    /**
     * Called to initialize the modeler
     * Must be called before any calculation is done
     */
    void Initialize();
  
    
    /**
     * level of echo for the mesh modeler
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
     * Remesh information is given to the modeler
     */
    void SetMeshingParameters( MeshingParametersType::Pointer& rMeshingParameters, ModelPart::IndexType MeshId );

    /**
     * Pre and Post meshing processes are given to the modeler
     */
    void SetPreMeshingProcess( Process::Pointer pPreMeshingProcess );

    void SetPostMeshingProcess( Process::Pointer pPostMeshingProcess );

    void SetPreMeshingProcessVector( std::vector<Process::Pointer>& rPreMeshingProcessVector );

    void SetPostMeshingProcessVector( std::vector<Process::Pointer>& rPostMeshingProcessVector );

    /**
     * Utilities for the mesher are given to the modeler
     */
    void SetModelerUtilities( ModelerUtilities::Pointer rModelerUtilities );
  
    /**
     * Transfer utilities are given to the modeler
     */
    void SetDataTransferUtilities( MeshDataTransferUtilities::Pointer rDataTransferUtilities );


    //*******************************************************************************************
    //*******************************************************************************************


    /**
     * Mesh Modeler :: Initilize
     */
    virtual void InitializeMeshModeler(ModelPart& rModelPart, ModelPart::IndexType MeshId=0);
  
    /**
     * Mesh Generation :: Remesh all ModelPart
     */
    virtual void GenerateMesh(ModelPart& rModelPart, ModelPart::IndexType MeshId=0);
       
    /**
     * Mesh Modeler :: Finalize
     */
    virtual void FinalizeMeshModeler(ModelPart& rModelPart, ModelPart::IndexType MeshId=0);
 
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
        return "MeshModeler";
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
  
    std::vector<Process::Pointer>         mPreMeshingProcesses;

    std::vector<Process::Pointer>        mPostMeshingProcesses;
    
    ModelerUtilities::Pointer               mpModelerUtilities;

    MeshDataTransferUtilities::Pointer mpDataTransferUtilities;  

    int mEchoLevel;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    /// Assignment operator.
    MeshModeler& operator=(MeshModeler const& rOther);

    /**
     * Mesh Modeler :: Start Echo
     */
    virtual void StartEcho(ModelPart& rModelPart,
			   std::string GenerationMessage, 
			   ModelPart::IndexType MeshId=0);


    /**
     * Mesh Modeler :: End Echo
     */
    virtual void EndEcho(ModelPart& rModelPart,
			 std::string GenerationMessage, 
			 ModelPart::IndexType MeshId=0);



    /**
     * Mesh Modeler :: Process to be done at the begining of the Generation
     */
    virtual void InitializeMeshGeneration(ModelPart& rModelPart,
					  MeshingParametersType& rMeshingVariables,
					  ModelPart::IndexType MeshId=0);

    /**
     * Mesh Modeler :: Process to be done at the end of the Generation
     */
    virtual void FinalizeMeshGeneration(ModelPart& rModelPart,
					MeshingParametersType& rMeshingVariables,
					ModelPart::IndexType MeshId=0);
  
    /**
     * Mesh Modeler :: Variables Transfer without remeshing
     */
    virtual void PerformTransferOnly(ModelPart& rModelPart,
				     MeshingParametersType& rMeshingVariables,
				     ModelPart::IndexType MeshId=0){};

    /**
     * Mesh Modeler :: Delaunay Tessellation
     */
    virtual void GenerateDT (ModelPart& rModelPart,
			     MeshingParametersType& rMeshingVariables,
			     ModelPart::IndexType MeshId=0){};
    
    /**
     * Mesh Modeler :: Constrained Delaunay Tessellation
     */
    virtual void GenerateCDT(ModelPart& rModelPart,
			     MeshingParametersType& rMeshingVariables,
			     ModelPart::IndexType MeshId=0){};
    
    /**
     * Mesh Modeler :: Refined Delaunay Tessellation
     */
    virtual void GenerateRDT(ModelPart& rModelPart,
			     MeshingParametersType& rMeshingVariables,
			     ModelPart::IndexType MeshId=0){};

    /**
     * Mesh Modeler :: Refined Constrained Delaunay Tessellation
     */
    virtual void GenerateRCDT(ModelPart& rModelPart,
			      MeshingParametersType& rMeshingVariables,
			      ModelPart::IndexType MeshId=0){};
        


    /**
     * Mesh Modeler :: Set Element Neighbours
     */
    virtual void SetElementNeighbours(ModelPart& rModelPart, 
				      MeshingParametersType& rMeshingVariables,
				      ModelPart::IndexType MeshId=0);


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

}; // Class MeshModeler

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MeshModeler& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MeshModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_MESH_MODELER_H_INCLUDED  defined 


