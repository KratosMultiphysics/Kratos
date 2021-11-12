//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Angel Celigueta (maceli@cimne.upc.edu)
//


#if !defined(KRATOS_APPLY_LASER_PROCESS )
#define  KRATOS_APPLY_LASER_PROCESS

#include "includes/table.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "utilities/parallel_utilities.h"
#include "utilities/function_parser_utility.h"
#include "utilities/interval_utility.h"
#include "utilities/python_function_callback_utility.h"


#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
#include "geometries/tetrahedra_3d_4.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/timer.h"
#include "processes/node_erase_process.h"
#include "utilities/binbased_fast_point_locator.h"
#include "includes/deprecated_variables.h"

#include <boost/timer.hpp>
#include "utilities/timer.h"

namespace Kratos
{

class ApplyLaserProcess : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyLaserProcess);

    typedef Table<double, double> TableType;
    array_1d<int, 3> mPositionLaserTableId;
    double mradius;
    double mpower;

    /// Constructor
    ApplyLaserProcess(
        ModelPart& rModelPart,
        Parameters rParameters
        ) : Process() ,
            mrModelPart(rModelPart)
    {
        KRATOS_TRY

        Parameters default_parameters( R"(
            {
                "model_part_name":"MODEL_PART_NAME",
                "laser_profile" :{
            		"shape": "flat",
            		"radius": 0.0005,
            		"power": 5.0},
		"direction": [0.0, 0.0, -1.0],
		"path" : {
                	"table"    : [0.0, 0.0, 0.0,0.0]} 
		

            }  )" );


 //       Parameters default_parameters( R"(
 //           {
 //               //"model_part_name":"MODEL_PART_NAME"
 //                "x_coordinate": 0.0,
//		 "y_coordinate": 0.0,
//		 "z_coordinate": 0.0,
//		 "radius": 0.0,
//		 "face_heat_flux": 0.0
  //          }  )" );



        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        //mLaserProfile.clear();
        //mDirection.clear();
        mPositionLaserTable.clear();


	if(rParameters["laser_profile"]["radius"].IsNumber()) {
        	mradius=  rParameters["laser_profile"]["radius"].GetDouble();
	}

	if(rParameters["laser_profile"]["power"].IsNumber()) {
        	mpower=  rParameters["laser_profile"]["power"].GetDouble();
	}

	/*for(int i=0; i<3; i++) {	
	if(rParameters["path"]["table"][i].IsNull()) {
                mPositionLaserTableId[i] = 0;
            }
            else {
                mPositionLaserTableId[i] = rParameters["path"]["table"][i].GetInt();

		KRATOS_WATCH(mPositionLaserTableId[i])
            }
            mPositionLaserTable.push_back(rModelPart.pGetTable(mPositionLaserTableId[i])); 
        }*/

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ApplyLaserProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ApplyLaserProcess algorithms.
    void Execute() override
    {

    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY;
        KRATOS_CATCH("");
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override {
        KRATOS_TRY;

        double x= 0.57102;
	double y= 0.99997;
        double z= 0.5299;
        double radius= mradius;
        double face_heat_flux =mpower;
        const int TDim=3;

	//defintions for spatial search
	typedef Node < 3 > PointType;
        typedef Node < 3 > ::Pointer PointTypePointer;
        typedef std::vector<PointType::Pointer> PointVector;
        typedef std::vector<PointType::Pointer>::iterator PointIterator;
        typedef std::vector<double> DistanceVector;
        typedef std::vector<double>::iterator DistanceIterator;
	
        //creating an auxiliary list for the new nodes
        PointVector list_of_nodes;

        //  *************
        // Bucket types
        typedef Bucket< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
	
        typedef Tree< KDTreePartition<BucketType> > tree; //Kdtree;
	
	
        //starting calculating time of construction of the kdtree
        boost::timer kdtree_construction;
	
        for (ModelPart::NodesContainerType::iterator node_it = mrModelPart.NodesBegin();
	     node_it != mrModelPart.NodesEnd(); ++node_it)
	  {
            PointTypePointer pnode = *(node_it.base());
	    
            //putting the nodes of the destination_model part in an auxiliary list
            list_of_nodes.push_back(pnode);

	    (node_it)->FastGetSolutionStepValue(FACE_HEAT_FLUX)=0.0;

	  }
	
        std::cout << "kdt constructin time " << kdtree_construction.elapsed() << std::endl;
	
        //create a spatial database with the list of new nodes
        unsigned int bucket_size = 20;
        tree nodes_tree(list_of_nodes.begin(), list_of_nodes.end(), bucket_size);
	
        //work arrays
        Node < 3 > work_point(0, 0.0, 0.0, 0.0);
        unsigned int MaximumNumberOfResults = 10000;
        PointVector Results(MaximumNumberOfResults);
        DistanceVector SquaredResultsDistances(MaximumNumberOfResults);
	
        double sigma = 0.0;
        if (TDim == 2)
	  sigma = 10.0 / (7.0 * 3.1415926);
        else
	  sigma = 1.0 / 3.1415926;
	
        work_point.X() = x;
	work_point.Y() = y;
	work_point.Z() = z;

        //find all of the new nodes within the radius
        int number_of_points_in_radius;
		
        //look between the new nodes which of them is inside the radius of the circumscribed cyrcle
        number_of_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, Results.begin(), SquaredResultsDistances.begin(), MaximumNumberOfResults);
		
        if (number_of_points_in_radius > 0)
	{
		    
		double maximunweight = SPHCubicKernel(sigma, 0.0, radius);    
		for (int k = 0; k < number_of_points_in_radius; k++)
		{
			double distance = sqrt(*(SquaredResultsDistances.begin() + k));
	
			double weight = SPHCubicKernel(sigma, distance, radius);
			
			PointIterator it_found = Results.begin() + k;
			
			if((*it_found)->FastGetSolutionStepValue(IS_FREE_SURFACE)==1) //MATERIAL_VARIABLE
			  {
			
			    double& aux= (*it_found)->FastGetSolutionStepValue(FACE_HEAT_FLUX);
                            aux =  face_heat_flux * weight / maximunweight; 
			}
		}	     
	  }
        KRATOS_CATCH("")

    }


    /**
     * @brief This function will be executed at every time step AFTER performing the solve phase
     */
    void ExecuteFinalizeSolutionStep() override {

    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyLaserProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyLaserProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:



///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    ModelPart& mrModelPart;

    //mLaserProfile.clear();
    //mDirection.clear();
    std::vector<TableType::Pointer> mPositionLaserTable;


    /// Assignment operator.
    ApplyLaserProcess& operator=(ApplyLaserProcess const& rOther);

    /// Copy constructor.
    //ApplyLaserProcess(ApplyLaserProcess const& rOther);


	inline double SPHCubicKernel(const double sigma, const double r, const double hmax)
      {
        const int TDim=3;
        double h_half = 0.5 * hmax;
        const double s = r / h_half;
        const double coeff = sigma / pow(h_half, static_cast<int>(TDim));
	
        if (s <= 1.0)
	  return coeff * (1.0 - 1.5 * s * s + 0.75 * s * s * s);
        else if (s <= 2.0)
	  return 0.25 * coeff * pow(2.0 - s, 3);
        else
            return 0.0;
      }



}; // Class ApplyLaserProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyLaserProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyLaserProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_APPLY_LASER_PROCESS defined */
