//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_TRIANGLE_MESH_2D_MODELER_H_INCLUDED )
#define  KRATOS_TRIANGLE_MESH_2D_MODELER_H_INCLUDED


/* If SINGLE is defined when triangle.o is compiled, it should also be */
/*   defined here.  If not, it should not be defined here.             */

/* #define SINGLE */

#ifdef   SINGLE
#define  REAL float
#else    // not SINGLE
#define  REAL double
#endif   // not SINGLE

// External includes
#if !defined(KRATOS_TRIANGLE_EXTERNAL_H_INCLUDED)
#define  KRATOS_TRIANGLE_EXTERNAL_H_INCLUDED
#include "triangle.h"
#endif

// System includes

// Project includes
#include "geometries/triangle_2d_3.h"

#include "custom_processes/boundary_skin_build_process.hpp"

#include "custom_utilities/mesh_error_calculation_utilities.hpp"
#include "custom_utilities/change_tip_elements_utilities.hpp"

#include "custom_modelers/laplacian_smoothing.hpp"
#include "custom_modelers/modeler.hpp"

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
class TriangleMesh2DModeler
{
protected:

    enum TriangleErrors {INPUT_MEMORY_ERROR=1, INTERNAL_ERROR=2, INVALID_GEOMETRY_ERROR=3};

    struct RefineVariables
    {      
      int      NumberOfElements;     

      double   size_factor;       //nodal h factor

      double   critical_radius;   //critical area   size
      double   critical_side;     //critical length size

      double   critical_dissipation;
 
      double   reference_error;
    };


    struct BoundingBoxVariables
    {
      bool    is_set;
      double  Radius;
      Vector  Center;
      Vector  Velocity;
    };

    struct MeshingVariables
    {

	//Pointer variables
	const Element   *mpReferenceElement;
	const Condition *mpReferenceCondition;

	void SetReferenceElement   (const Element   & rElement)
	{
	    mpReferenceElement=&rElement;
	};
	void SetReferenceCondition (const Condition & rCondition)
	{
	    mpReferenceCondition=&rCondition;
	};

	Element const&    GetReferenceElement   ()
	{
	    return *mpReferenceElement;
	};
        Condition const&  GetReferenceCondition ()
	{
	    return *mpReferenceCondition;
	};

	//Other variables (general struct access)
	Flags   MeshingOptions;
	Flags   RefiningOptions;
      
        bool    remesh;
        bool    refine;
        bool    constrained;
        bool    mesh_smoothing;
        bool    jacobi_smoothing;
        bool    avoid_tip_elements;
      
        bool    idset;
	double  AlphaParameter;
        double  offset_factor;

	std::vector<int> PreIds;
	std::vector<int> NewIds;

        std::vector<int> PreservedElements;
        std::vector<std::vector<int> > NeighbourList; 

        RefineVariables  Refine;  
	// int     NumberElements;
        // double  h_factor;
	// double  critical_dissipation;
	// double  critical_radius;
	// double  critical_side;
	// double  reference_error;
      
        BoundingBoxVariables WallTip;
        BoundingBoxVariables BoundingBox;
 
    };


public:




  
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TriGenCDT
    KRATOS_CLASS_POINTER_DEFINITION(TriangleMesh2DModeler);


    typedef Node<3>                                                       PointType;
    typedef Node<3>::Pointer                                       PointPointerType;
    typedef std::vector<PointPointerType>                        PointPointerVector;
    typedef ModelPart::PropertiesType                                PropertiesType;
    typedef ModelPart::PropertiesContainerType              PropertiesContainerType;
    typedef ModelPart::MeshType                                            MeshType;
    typedef ModelPart::ElementsContainerType                  ElementsContainerType;
    typedef ModelPart::NodesContainerType                        NodesContainerType;
    typedef ModelPart::MeshType::GeometryType::PointsArrayType      PointsArrayType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TriangleMesh2DModeler() {} //

    /// Destructor.
    virtual ~TriangleMesh2DModeler() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    //*******************************************************************************************
    //*******************************************************************************************
    void SetInitialMeshData (int number_of_domains);
  
    void SetRemeshData (Element   const& rReferenceElement,
			Condition const& rReferenceCondition,
			bool remesh             = false,
			bool constrained        = false,
			bool mesh_smoothing     = false,
			bool jacobi_smoothing   = false,
			bool avoid_tip_elements = false,
			double alpha            = 1.4,
			double my_offset        = 1.0,
			int    MeshId           = 0);
 
    void SetRefineData (bool refine        = false,
			double h_factor    = 0.5,
			double dissipation = 40,
			double radius      = 0.00004,
			double error       = 2,
			int MeshId         = 0);

    void SetWallTip (double radius,
		     Vector center);

    void SetRefiningBox (double radius,
			 Vector center,
			 Vector velocity);


    void GenerateMesh(ModelPart& rModelPart);
    

  
    //*******************************************************************************************
    //*******************************************************************************************
    void PerformTransferOnly(ModelPart& rModelPart,
			     MeshingVariables & rMeshingVariables,
			     ModelPart::IndexType MeshId=0);
    

  
    //*******************************************************************************************
    //*******************************************************************************************

    void GenerateDT (ModelPart& rModelPart,
		     MeshingVariables &rMeshingVariables,
		     ModelPart::IndexType MeshId=0);
    



    //*******************************************************************************************
    //*******************************************************************************************
    void GenerateCDT(ModelPart& rModelPart,
		     MeshingVariables & rMeshingVariables,
		     ModelPart::IndexType MeshId=0);
    


    //*******************************************************************************************
    //*******************************************************************************************

    void GenerateWPDT(ModelPart& rModelPart,
		      MeshingVariables & rMeshingVariables,
		      ModelPart::IndexType MeshId=0);
    


    //*******************************************************************************************
    //*******************************************************************************************
  
    void GenerateWPCDT(ModelPart& rModelPart,
		       MeshingVariables & rMeshingVariables,
		       ModelPart::IndexType MeshId=0);
        

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
    std::vector<MeshingVariables> mVariables;

    ModelerUtilities mModelerUtilities;

    MeshDataTransferUtilities mMeshDataTransferUtilities;

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

protected:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
 

    ///@}
    ///@name Private Operators
    ///@{

    /// Assignment operator.
    TriangleMesh2DModeler& operator=(TriangleMesh2DModeler const& rOther);

    ///@}
    ///@name Private Operations
    ///@{

    void clean_triangulateio(struct triangulateio& tr );
 
    void clear_triangulateio(struct triangulateio& tr );

    void free_triangleio(struct triangulateio& tr );

    void free_pointio (struct triangulateio& tr );
 
    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Unaccessible methods
    ///@{

    void SetTriangulateNodes   (ModelPart &rModelPart,
				MeshingVariables & rMeshingVariables,
				struct triangulateio &in,
				struct triangulateio &out,
				ModelPart::IndexType MeshId=0);
    

    void RefineElements (ModelPart &rModelPart,MeshingVariables &rMeshingVariables,struct triangulateio &in,struct triangulateio &out,ModelPart::IndexType MeshId=0);
    
    void WriteTriangles  (struct triangulateio &in);
    
    void WritePoints  (struct triangulateio &in);
    

    int  GenerateTriangulation  (Flags &MeshingOptions,
				 Flags &RefiningOptions,
				 struct triangulateio &in,
				 struct triangulateio &out);
    
    void RecoverBoundaryPosition (ModelPart &rModelPart,
				  MeshingVariables & rMeshingVariables,
				  struct triangulateio &in,
				  struct triangulateio &out,
				  ModelPart::IndexType MeshId=0);
    
    //Select elements after the Delaunay Tesselation
    void SelectMeshElements (ModelPart::NodesContainerType& rNodes,
			     MeshingVariables & rMeshingVariables,
			     struct triangulateio & out);
    
  
   
    //Set elements in model_part after the Delaunay Tesselation
    void SetModelElements (ModelPart& rModelPart,
			   MeshingVariables& rMeshingVariables,
			   struct triangulateio &out,
			   ModelPart::IndexType MeshId=0);
    
    void RemoveCloseNodes (ModelPart& rModelPart, MeshingVariables& rMeshingVariables,ModelPart::IndexType MeshId=0);
    
    void RemoveNonConvexBoundary(ModelPart& rModelPart, MeshingVariables& rMeshingVariables,ModelPart::IndexType MeshId=0);
    
    void SetDissipativeElements (ModelPart& rModelPart, MeshingVariables& rMeshingVariables,ModelPart::IndexType MeshId=0);
    
    void RefineBoundary(ModelPart& rModelPart, MeshingVariables& rMeshingVariables,ModelPart::IndexType MeshId=0);
    

    void GenerateNewParticles (ModelPart& rModelPart,MeshingVariables &rMeshingVariables,struct triangulateio &in,struct triangulateio &out,ModelPart::IndexType MeshId=0);
    
  
						   
    
    void SetElementNeighbours(ModelPart& rModelPart, MeshingVariables & rMeshingVariables,ModelPart::IndexType MeshId=0);

    void SetConditionsBoundary(ModelPart& rModelPart, MeshingVariables& rMeshingVariables,ModelPart::IndexType MeshId=0);
    

   ///@}

}; // Class TriangleMesh2DModeler

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
    inline std::istream& operator >> (std::istream& rIStream,
				      TriangleMesh2DModeler& rThis);

/// output stream function
    inline std::ostream& operator << (std::ostream& rOStream,
				      const TriangleMesh2DModeler& rThis)
    {
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
    }
///@}


}  // namespace Kratos.

#endif // KRATOS_TRIANGLE_MESH_2D_MODELER_H_INCLUDED  defined 






//************************************************************************************
//************************************************************************************
//************************************************************************************
//************************************************************************************
//************************************************************************************





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
