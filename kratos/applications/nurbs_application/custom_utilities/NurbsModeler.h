//
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_NURBSMODELER_H_INCLUDED )
#define  KRATOS_NURBSMODELER_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <ctime>

// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "geometries/nurbs_3d.h"
#include "geometries/nurbs_2d.h"
#include "includes/variables.h"
#include "custom_elements/nurbs_poisson_2d.h"
#include "geometries/nurbs_base_geometry.h"
#include "spatial_containers/spatial_containers.h"

namespace Kratos
{

	///@name Kratos Globals
	///@{ 

	// Variables definition 



	///@} 
	///@name Type Definitions
	///@{ 
/**
 * Typedefs for search
 */




typedef PointerVectorSet< Node<3>, IndexedObject> NodesContainerType;
typedef Node<3> PointType;
typedef Node<3>::Pointer PointTypePointer;
typedef std::vector<PointType::Pointer >           PointVector;
typedef typename std::vector<PointType::Pointer >::iterator PointIterator;
typedef std::vector<double>               DistanceVector;
typedef typename std::vector<double>::iterator     DistanceIterator;
typedef Bins< 3, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > StaticBins;
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
    class NurbsModeler
	{
	public:
		///@name Type Definitions
		///@{
		

		/// Pointer definition of KratosNurbsTestcaseApplication
        KRATOS_CLASS_POINTER_DEFINITION(NurbsModeler);

		///@}
		///@name Life Cycle 
		///@{ 

        void ReadModelPart(ModelPart &NurbsModelPart,
                           int dimension,
                           int connectivities_start,
                           int connectivities_end,
                           int polynomial_degree_xi,
                           int polynomial_degree_eta,
                           Vector knots_xi,
                           Vector knots_eta,
                           int number_of_cp_xi,
                           int number_of_cp_eta,
                           Vector weights);

        void InterpolateDesignVariables(ModelPart &NurbsModelPart, ModelPart &TriangleModelPart);

        void ClosestPoint(ModelPart &FluidMesh, ModelPart &NurbsModelPart, double Max_search_radius);

		/// Default constructor.
        NurbsModeler();

		/// Destructor.
        virtual ~NurbsModeler();







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



		//       static const ApplicationCondition  msApplicationCondition; 

		///@} 
        ///@name Member Variables
		///@{ 


        std::vector<NurbsPatchGeometry< Node<3> >::Pointer> mNurbsPatchGeometry;




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
        NurbsModeler& operator=(NurbsModeler const& rOther);

		/// Copy constructor.
        NurbsModeler(NurbsModeler const& rOther);


		///@}    

	}; // Class KratosNurbsTestcaseApplication 

	///@} 


	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 

	///@} 


}  // namespace Kratos.

#endif // KRATOS_NURBSTESTCASE_APPLICATION_H_INCLUDED  defined 


