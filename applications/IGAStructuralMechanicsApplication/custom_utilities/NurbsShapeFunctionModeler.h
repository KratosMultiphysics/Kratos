//
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_NURBS_SHAPE_FUNCTION_MODELER_H_INCLUDED )
#define  KRATOS_NURBS_SHAPE_FUNCTION_MODELER_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <ctime>

// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "geometries/nurbs_3d.h"
#include "geometries/nurbs_2d.h"
//#include "custom_elements/meshless_base_element.h"
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

	//typedef PointerVectorSet< Node<3>, IndexedObject> NodesContainerType;
	//typedef Node<3> PointType;
	//typedef Node<3>::Pointer PointTypePointer;
	//typedef std::vector<PointType::Pointer >           PointVector;
	//typedef typename std::vector<PointType::Pointer >::iterator PointIterator;
	//typedef std::vector<double>               DistanceVector;
	//typedef typename std::vector<double>::iterator     DistanceIterator;
	//typedef Bins< 3, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > StaticBins;
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
    class NurbsShapeFunctionModeler
	{
	public:
		///@name Type Definitions
		///@{
		

		/// Pointer definition of KratosNurbsTestcaseApplication
        KRATOS_CLASS_POINTER_DEFINITION(NurbsShapeFunctionModeler);

		///@}
		///@name Life Cycle 
		///@{ 

		void SetUp(
			ModelPart &NurbsModelPart,
			int dimension,
			Vector ControlPoints, //integer IDs
			Vector Weights,
			int PolynomialDegreeXi,
			int PolynomialDegreeEta,
			Vector KnotsXi,
			Vector KnotsEta);
			
		void EvaluateShapeFunction(
			const array_1d<double, 3>& LocalCoordinatesOfQuadraturePoint,
			Vector &ShapeFunctionValues,
			Matrix &ShapeFunctionLocalDerivatives);

		void EvaluateShapeFunctionSecondOrder(
			const array_1d<double, 3>& LocalCoordinatesOfQuadraturePoint,
			Vector &ShapeFunctionValues,
			Matrix &ShapeFunctionLocalDerivatives);

		/// Default constructor.
		NurbsShapeFunctionModeler();

		/// Destructor.
        virtual ~NurbsShapeFunctionModeler();







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


		NurbsPatchGeometry2D< Node<3> >::Pointer mNurbsPatchGeometry;




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
		NurbsShapeFunctionModeler& operator=(NurbsShapeFunctionModeler const& rOther);

		/// Copy constructor.
		NurbsShapeFunctionModeler(NurbsShapeFunctionModeler const& rOther);


		///@}    

	}; // Class NurbsShapeFunctionModeler 

	///@} 


	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 

	///@} 


}  // namespace Kratos.

#endif // KRATOS_NURBS_SHAPE_FUNCTION_MODELER_APPLICATION_H_INCLUDED  defined 


