//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_TRIANGULAR_MESH_2D_MODELER_H_INCLUDED )
#define  KRATOS_TRIANGULAR_MESH_2D_MODELER_H_INCLUDED

// If SINGLE is defined when triangle.o is compiled, it should also be defined here
// If SINGLE is NOT defined in compilation, it should not be defined here.         
// #define SINGLE

#ifdef   SINGLE
#define  REAL float
#else    // not SINGLE
#define  REAL double
#endif   // not SINGLE

// External includes
#ifndef TRILIBRARY
#define TRILIBRARY
#endif

#if !defined(KRATOS_TRIANGLE_EXTERNAL_H_INCLUDED) 
#define  KRATOS_TRIANGLE_EXTERNAL_H_INCLUDED 
#include "triangle.h"
#endif

// System includes

// Project includes
#include "geometries/triangle_2d_3.h"
#include "custom_modelers/mesh_modeler.hpp"

///VARIABLES used:
//Data:    
//StepData: 
//Flags:    (checked)  
//          (set)      
//          (modified)  
//          (reset)


namespace Kratos
{

    extern "C" {
	void triangulate(char *, struct triangulateio *, struct triangulateio *,struct triangulateio *);
	void trifree(void *);
    }


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
class KRATOS_API(PFEM_BASE_MECHANICS_APPLICATION) TriangularMesh2DModeler
  : public MeshModeler
{
protected:

    enum TriangleErrors {INPUT_MEMORY_ERROR=1, INTERNAL_ERROR=2, INVALID_GEOMETRY_ERROR=3};

public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of TriGenCDT
    KRATOS_CLASS_POINTER_DEFINITION( TriangularMesh2DModeler );
 
    typedef ModelerUtilities::InfoParameters                     InfoParametersType;
    typedef ModelerUtilities::MeshingParameters               MeshingParametersType;
    typedef ModelerUtilities::RefiningParameters               RefineParametersType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TriangularMesh2DModeler(): MeshModeler() {} //

    /// Copy constructor.
    TriangularMesh2DModeler(TriangularMesh2DModeler const& rOther): MeshModeler(rOther) {}

    /// Destructor.
    virtual ~TriangularMesh2DModeler() {}

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


    //generate :: meshing step call to delaunay tessellation
    void Generate(ModelPart& rModelPart, MeshingParametersType& rMeshingVariables);

    //generate the Delaunay Tesselation
    int  GenerateTessellation(MeshingParametersType& rMeshingVariables, struct triangulateio& in, struct triangulateio& out);
  

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
   TriangularMesh2DModeler& operator=(TriangularMesh2DModeler const& rOther);

    ///@}
    ///@name Private Operations
    ///@{

    //build the input for the mesher
    void BuildInput ( ModelPart &rModelPart, MeshingParametersType & rMeshingVariables, struct triangulateio &in);
 
    //set and get from mesh container in meshing variables
    void GetFromContainer(ModelerUtilities::MeshContainer& rMesh, struct triangulateio& tr);

    void SetToContainer(ModelerUtilities::MeshContainer& rMesh, struct triangulateio& tr);


    //delete in/out structures
    void DeleteInContainer(ModelerUtilities::MeshContainer& rMesh, struct triangulateio& tr);

    void DeleteOutContainer(ModelerUtilities::MeshContainer& rMesh, struct triangulateio& tr);


    //set faces in the triangulateio before the Delaunay Tesselation
    virtual void SetFaces ( ModelPart &rModelPart, MeshingParametersType & rMeshingVariables, struct triangulateio &in );
    

    //print methods
    void WriteTriangles      ( struct triangulateio& tr );

    void WritePoints         ( struct triangulateio& tr );


    //free memory of the mesher
    void ClearTrianglesList  ( struct triangulateio& tr );

    void DeleteTrianglesList ( struct triangulateio& tr );
    
    void DeletePointsList    ( struct triangulateio& tr );

    void FreeTrianglesList   ( struct triangulateio& tr );


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

}; // Class TriangularMesh2DModeler

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
    inline std::istream& operator >> (std::istream& rIStream,
				      TriangularMesh2DModeler& rThis);

/// output stream function
    inline std::ostream& operator << (std::ostream& rOStream,
				      const TriangularMesh2DModeler& rThis)
    {
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
    }
///@}


}  // namespace Kratos.

#endif // KRATOS_TRIANGULAR_MESH_2D_MODELER_H_INCLUDED  defined 




//************************************************************************************
// Command line switches for TRIANGLE:

// To run Triangle, the command line syntax is

// For some of the command line switches described below, you may click on the switch for detailed information.
// -p Triangulates a Planar Straight Line Graph (.poly file).
// -r Refines a previously generated mesh.
// -q Quality mesh generation with no angles smaller than 20 degrees. An alternate minimum angle may be specified after the `q'.
// -a Imposes a maximum triangle area constraint. A fixed area constraint (that applies to every triangle) may be specified after the `a', or varying area constraints may be read from a .poly file or .area file.
// -u Imposes a user-defined constraint on triangle size.
// -A Assigns a regional attribute to each triangle that identifies what segment-bounded region it belongs to.
// -c Encloses the convex hull with segments.
// -D Conforming Delaunay: use this switch if you want all triangles in the mesh to be Delaunay, and not just constrained Delaunay; or if you want to ensure that all Voronoi vertices lie within the triangulation.
// -j Jettisons vertices that are not part of the final triangulation from the output .node file (including duplicate input vertices and vertices ``eaten'' by holes).
// -e Outputs (to an .edge file) a list of edges of the triangulation.
// -v Outputs the Voronoi diagram associated with the triangulation. Does not attempt to detect degeneracies, so some Voronoi vertices may be duplicated.
// -n Outputs (to a .neigh file) a list of triangles neighboring each triangle.
// -g Outputs the mesh to an Object File Format (.off) file, suitable for viewing with the Geometry Center's Geomview package.
// -B Suppresses boundary markers in the output .node, .poly, and .edge output files.
// -P Suppresses the output .poly file. Saves disk space, but you lose the ability to maintain constraining segments on later refinements of the mesh.
// -N Suppresses the output .node file.
// -E Suppresses the output .ele file.
// -I Suppresses mesh iteration numbers.
// -O Suppresses holes: ignores the holes in the .poly file.
// -X Suppresses exact arithmetic.
// -z Numbers all items starting from zero (rather than one). Note that this switch is normally overrided by the value used to number the first vertex of the input .node or .poly file. However, this switch is useful when calling Triangle from another program.
// -o2 Generates second-order subparametric elements with six nodes each.
// -Y Prohibits the insertion of Steiner points on the mesh boundary. If specified twice (-YY), it prohibits the insertion of Steiner points on any segment, including internal segments.
// -S Specifies the maximum number of added Steiner points.
// -i Uses the incremental algorithm for Delaunay triangulation, rather than the divide-and-conquer algorithm.
// -F Uses Steven Fortune's sweepline algorithm for Delaunay triangulation, rather than the divide-and-conquer algorithm.
// -l Uses only vertical cuts in the divide-and-conquer algorithm. By default, Triangle uses alternating vertical and horizontal cuts, which usually improve the speed except with vertex sets that are small or short and wide. This switch is primarily of theoretical interest.
// -s Specifies that segments should be forced into the triangulation by recursively splitting them at their midpoints, rather than by generating a constrained Delaunay triangulation. Segment splitting is true to Ruppert's original algorithm, but can create needlessly small triangles. This switch is primarily of theoretical interest.
// -C Check the consistency of the final mesh. Uses exact arithmetic for checking, even if the -X switch is used. Useful if you suspect Triangle is buggy.
// -Q Quiet: Suppresses all explanation of what Triangle is doing, unless an error occurs.
// -V Verbose: Gives detailed information about what Triangle is doing. Add more `V's for increasing amount of detail. `-V' gives information on algorithmic progress and detailed statistics.
// -h Help: Displays complete instructions.
