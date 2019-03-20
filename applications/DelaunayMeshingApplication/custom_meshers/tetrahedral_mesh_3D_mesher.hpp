//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#if !defined (KRATOS_TETRAHEDRAL_MESH_3D_MESHER_H_INCLUDED)
#define  KRATOS_TETRAHEDRAL_MESH_3D_MESHER_H_INCLUDED

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
#include "custom_meshers/mesher.hpp"

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
class KRATOS_API(DELAUNAY_MESHING_APPLICATION) TetrahedralMesh3DMesher
  : public Mesher
{
protected:

    enum TetgenErrors {INPUT_MEMORY_ERROR=1, INTERNAL_ERROR=2, INVALID_GEOMETRY_ERROR=3};

public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of TriGenCDT
    KRATOS_CLASS_POINTER_DEFINITION( TetrahedralMesh3DMesher );

    typedef MesherUtilities::MeshingInfoParameters              InfoParametersType;
    typedef MesherUtilities::MeshingParameters               MeshingParametersType;
    typedef MesherUtilities::RefiningParameters               RefineParametersType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TetrahedralMesh3DMesher(): Mesher() {} //

    /// Copy constructor.
    TetrahedralMesh3DMesher(TetrahedralMesh3DMesher const& rOther): Mesher(rOther) {}

    /// Destructor.
    virtual ~TetrahedralMesh3DMesher() {}

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
    std::string Info() const override
    {
	return "";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override{}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override{}


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

    //generate :: meshing step call to delaunay tessellation
    void Generate (ModelPart& rModelPart, MeshingParametersType& rMeshingVariables) override;

    //generate the Delaunay Tesselation
    int  GenerateTessellation(MeshingParametersType& rMeshingVariables, tetgenio& in, tetgenio& out);


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
   TetrahedralMesh3DMesher& operator=(TetrahedralMesh3DMesher const& rOther);

    ///@}
    ///@name Private Operations
    ///@{

    //build the input for the mesher
    void BuildInput ( ModelPart &rModelPart,
		      MeshingParametersType & rMeshingVariables,
		      tetgenio& in );

    //set and get from mesh container in meshing variables
    void GetFromContainer(MesherUtilities::MeshContainer& rMesh, tetgenio& tr);

    void SetToContainer(MesherUtilities::MeshContainer& rMesh, tetgenio& tr);

    //delete in/out structures
    void DeleteInContainer(MesherUtilities::MeshContainer& rMesh, tetgenio& tr);

    void DeleteOutContainer(MesherUtilities::MeshContainer& rMesh, tetgenio& tr);

    //set faces in the tetgenio before the Delaunay Tesselation
    virtual void SetFaces( ModelPart &rModelPart, MeshingParametersType & rMeshingVariables, tetgenio &in );

    //check points
    void CheckInOutPoints    ( tetgenio& in, tetgenio& out );

    //print methods
    void WriteTetrahedra     ( tetgenio& tr );

    void WritePoints         ( tetgenio& tr );

    //free memory of the mesher
    void ClearTetgenIO ( tetgenio& tr );

    // void DeleteTetrahedraList ( tetgenio& tr );

    // void DeletePointsList ( tetgenio& tr );

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Unaccessible methods
    ///@{

    ///@}

}; // Class TetrahedralMesh3DMesher

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
    inline std::istream& operator >> (std::istream& rIStream,
				      TetrahedralMesh3DMesher& rThis);

/// output stream function
    inline std::ostream& operator << (std::ostream& rOStream,
				      const TetrahedralMesh3DMesher& rThis)
    {
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
    }
///@}


}  // namespace Kratos.

#endif // KRATOS_TETRAHEDRAL_MESH_3D_MESHER_H_INCLUDED  defined





//************************************************************************************
// Command line switches for TETGEN:

//1.5.0

// -p	Tetrahedralizes a piecewise linear complex (PLC).
// -Y	Preserves the input surface mesh (does not modify it).
// -r	Reconstructs a previously generated mesh.
// -q	Refines mesh (to improve mesh quality).
// -R	Mesh coarsening (to reduce the mesh elements).
// -A	Assigns attributes to tetrahedra in different regions.
// -a	Applies a maximum tetrahedron volume constraint.
// -m	Applies a mesh sizing function.
// -i	Inserts a list of additional points.
// -O	Specifies the level of mesh optimization.
// -S	Specifies maximum number of added points.
// -T	Sets a tolerance for coplanar test (default 10−8).
// -X	Suppresses use of exact arithmetic.
// -M	No merge of coplanar facets or very close vertices.
// -w	Generates weighted Delaunay (regular) triangulation.
// -c	Retains the convex hull of the PLC.
// -d	Detects self-intersections of facets of the PLC.
// -z	Numbers all output items starting from zero.
// -f	Outputs all faces to .face file.
// -e	Outputs all edges to .edge file.
// -n	Outputs tetrahedra neighbors to .neigh file.
// -v	Outputs Voronoi diagram to files.
// -g	Outputs mesh to .mesh file for viewing by Medit.
// -k	Outputs mesh to .vtk file for viewing by Paraview.
// -J	No jettison of unused vertices from output .node file.
// -B	Suppresses output of boundary information.
// -N	Suppresses output of .node file.
// -E	Suppresses output of .ele file.
// -F	Suppresses output of .face and .edge file.
// -I	Suppresses mesh iteration numbers.
// -C	Checks the consistency of the final mesh.
// -Q	Quiet: No terminal output except errors.
// -V	Verbose: Detailed information, more terminal output.
// -h	Help: A brief instruction for using TetGen.

//1.4.3

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

