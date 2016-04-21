//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined (KRATOS_TETRAHEDRAL_MESH_3D_MODELER_H_INCLUDED)
#define  KRATOS_TETRAHEDRAL_MESH_3D_MODELER_H_INCLUDED

// External includes
#ifndef TETLIBRARY
#define TETLIBRARY
#endif

#if !defined(KRATOS_TETGEN_EXTERNAL_H_INCLUDED)
#define KRATOS_TETGEN_EXTERNAL_H_INCLUDED
#include "tetgen.h"
#endif

// System includes

// Project includes
#include "geometries/tetrahedra_3d_4.h"
#include "custom_modelers/mesh_modeler.hpp"

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
class TetrahedralMesh3DModeler
  : public MeshModeler
{
protected:

    enum TetgenErrors {INPUT_MEMORY_ERROR=1, INTERNAL_ERROR=2, INVALID_GEOMETRY_ERROR=3};

public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of TriGenCDT
    KRATOS_CLASS_POINTER_DEFINITION( TetrahedralMesh3DModeler );

    typedef ModelerUtilities::InfoParameters                     InfoParametersType;
    typedef ModelerUtilities::MeshingParameters               MeshingParametersType;
    typedef ModelerUtilities::RefiningParameters               RefineParametersType;

    typedef Node<3>                                                       PointType;
    typedef std::vector<PointType>                                      PointVector;
    typedef Node<3>::Pointer                                       PointPointerType;
    typedef std::vector<PointPointerType>                        PointPointerVector;
    typedef PointPointerVector::iterator                 PointPointerVectorIterator;
    typedef ModelPart::PropertiesType                                PropertiesType;
    typedef ModelPart::PropertiesContainerType              PropertiesContainerType;
    typedef ModelPart::MeshType                                            MeshType;
    typedef ModelPart::ElementsContainerType                  ElementsContainerType;
    typedef ModelPart::NodesContainerType                        NodesContainerType;
    typedef ModelPart::MeshType::GeometryType::PointsArrayType      PointsArrayType;

    //defintions for spatial search
    typedef std::vector<double>                                      DistanceVector;
    typedef std::vector<double>::iterator                          DistanceIterator;

    typedef Bucket<3, PointType, PointPointerVector, PointPointerType, PointPointerVectorIterator, DistanceIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> >                          KdtreeType; //Kdtree


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TetrahedralMesh3DModeler(): MeshModeler() {} //

    /// Copy constructor.
    TetrahedralMesh3DModeler(TetrahedralMesh3DModeler const& rOther): MeshModeler(rOther) {}

    /// Destructor.
    virtual ~TetrahedralMesh3DModeler() {}

    ///@}
    ///@name Operators
    ///@{


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
    virtual std::string Info() const
    {
	return "";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


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

  
    //*******************************************************************************************
    //*******************************************************************************************
    void PerformTransferOnly(ModelPart& rModelPart,
			     MeshingParametersType& rMeshingVariables,
			     ModelPart::IndexType MeshId=0);
    

  
    //*******************************************************************************************
    //*******************************************************************************************

    void GenerateDT (ModelPart& rModelPart,
		     MeshingParametersType& rMeshingVariables,
		     ModelPart::IndexType MeshId=0);
    



    //*******************************************************************************************
    //*******************************************************************************************
    void GenerateCDT(ModelPart& rModelPart,
		     MeshingParametersType& rMeshingVariables,
		     ModelPart::IndexType MeshId=0);
    


    //*******************************************************************************************
    //*******************************************************************************************

    void GenerateRDT(ModelPart& rModelPart,
		     MeshingParametersType& rMeshingVariables,
		     ModelPart::IndexType MeshId=0);
    


    //*******************************************************************************************
    //*******************************************************************************************
  
    void GenerateRCDT(ModelPart& rModelPart,
		      MeshingParametersType& rMeshingVariables,
		      ModelPart::IndexType MeshId=0);
        



    //*******************************************************************************************
    //*******************************************************************************************

    //Select elements after the Delaunay Tesselation
    void SelectMeshElements(ModelPart::NodesContainerType& rNodes,
			    MeshingParametersType& rMeshingVariables,
			    tetgenio& out);



    //Generate the Delaunay Tesselation
    int  GenerateTetrahedralization(Flags& MeshingOptions,
				    Flags& RefiningOptions,
				    tetgenio& in,
				    tetgenio& out);
  
    //Free memory of the mesher
    void ClearTetrahedraList  ( tetgenio& tr );
  
    void DeleteTetrahedraList ( tetgenio& tr );

    void DeletePointsList     ( tetgenio& tr );

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

    /// Assignment operator.
   TetrahedralMesh3DModeler& operator=(TetrahedralMesh3DModeler const& rOther);

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
    ///@name Unaccessible methods
    ///@{
  

    //METHODS CALLED BEFORE TESSELLATION


    //Set nodes in the triangulateio before the Delaunay Tesselation
    void SetTetrahedralizationNodes ( ModelPart &rModelPart,
				      MeshingParametersType & rMeshingVariables,
				      tetgenio &in,
				      tetgenio &out,
				      ModelPart::IndexType MeshId=0 );


    //recover the boundary position after an small offset when remeshing constrained
    void RecoverBoundaryPosition ( ModelPart &rModelPart,
				   MeshingParametersType & rMeshingVariables,
				   tetgenio &in,
				   tetgenio &out,
				   ModelPart::IndexType MeshId=0 );



    //METHODS CALLED AFTER TESSELLATION

    //Build elements in model_part after the Delaunay Tesselation
    void BuildMeshElements ( ModelPart& rModelPart,
			     MeshingParametersType& rMeshingVariables,
			     tetgenio& out,
			     ModelPart::IndexType MeshId=0 );
 

    //*******************************************************************************************
    //*******************************************************************************************


    //METHODS CALLED BEFORE REFINING THE TESSELLATION
    void SetDissipativeElements ( ModelPart& rModelPart, 
				  MeshingParametersType& rMeshingVariables,
				  ModelPart::IndexType MeshId=0 );
       

    //METHODS CALLED AFTER REFINING THE TESSELLATION   
    void RefineElements ( ModelPart &rModelPart,
			  MeshingParametersType& rMeshingVariables,
			  tetgenio &in,
			  tetgenio &out,
			  ModelPart::IndexType MeshId=0 );

    void GenerateNewParticles ( ModelPart& rModelPart,
				MeshingParametersType& rMeshingVariables,
				tetgenio &in, 
				tetgenio &out,
				ModelPart::IndexType MeshId=0 );

    //*******************************************************************************************
    //*******************************************************************************************

    void WriteTetrahedra     ( tetgenio& tr );
    void WritePoints         ( tetgenio& tr );
 
   ///@}

}; // Class TetrahedralMesh3DModeler

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
    inline std::istream& operator >> (std::istream& rIStream,
				      TetrahedralMesh3DModeler& rThis);

/// output stream function
    inline std::ostream& operator << (std::ostream& rOStream,
				      const TetrahedralMesh3DModeler& rThis)
    {
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
    }
///@}


}  // namespace Kratos.

#endif // KRATOS_TETRAHEDRAL_MESH_3D_MODELER_H_INCLUDED  defined 





//************************************************************************************
// Command line switches for TETGEN:


// -p Tetrahedralizes a piecewise linear complex (.poly or .smesh file).
// -q Quality mesh generation. A minimum radius-edge ratio may be specified (default 2.0).
// -a Applies a maximum tetrahedron volume constraint.
// -A Assigns attributes to identify tetrahedra in certain regions.
// -r Reconstructs and Refines a previously generated mesh.
// -Y Suppresses boundary facets/segments splitting.
// -i Inserts a list of additional points into mesh.
// -M Does not merge coplanar facets.
// -T Set a tolerance for coplanar test (default 1e-8).
// -d Detect intersections of PLC facets.
// -z Numbers all output items starting from zero.
// -o2 Generates second-order subparametric elements.
// -f Outputs faces (including non-boundary faces) to .face file.
// -e Outputs subsegments to .edge file.
// -n Outputs tetrahedra neighbors to .neigh file.
// -g Outputs mesh to .mesh file for viewing by Medit.
// -G Outputs mesh to .msh file for viewing by Gid.
// -O Outputs mesh to .off file for viewing by Geomview.
// -J No jettison of unused vertices from output .node file.
// -B Suppresses output of boundary information.
// -N Suppresses output of .node file.
// -E Suppresses output of .ele file.
// -F Suppresses output of .face file.
// -I Suppresses mesh iteration numbers.
// -C Checks the consistency of the final mesh.
// -Q Quiet: No terminal output except errors.
// -V Verbose: Detailed information, more terminal output.
// -v Prints the version information.
// -h Help: A brief instruction for using TetGen.

