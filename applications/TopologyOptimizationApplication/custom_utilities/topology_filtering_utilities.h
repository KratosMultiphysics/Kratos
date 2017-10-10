// ==============================================================================
//  KratosTopologyOptimizationApplication
//
//  License:         BSD License
//                   license: TopologyOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//
// ==============================================================================

#if !defined(KRATOS_TOPOLOGY_FILTERING_UTILITIES_H_INCLUDED)
#define  KRATOS_TOPOLOGY_FILTERING_UTILITIES_H_INCLUDED

// System includes
#include <iostream>
#include <string>
#include <algorithm>
#include <iomanip>      // for std::setprecision

// External includes
#include <boost/python.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/process_info.h"
#include "containers/array_1d.h"
#include <vector>

// Application includes
#include "topology_optimization_application.h"
#include "spatial_containers/spatial_containers.h" // For kd-tree
#include "custom_utilities/filter_function.h"


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

/// Solution utility to filter results.
/** Detail class definition.

 */

class TopologyFilteringUtilities
{
public:

	// Define an auxiliary structure to hold a pointer to the element of interest and a position in space to be searched for.
	// This is needed to know the neighbours and where the element is for the filter.

	class ElementPositionItem: public Point
	{
	public:

		KRATOS_CLASS_POINTER_DEFINITION( ElementPositionItem );

		ElementPositionItem(): // This Constructor is used by the tree
			Point(),
			mpElement()
			{
			}

		ElementPositionItem(array_1d<double,3> Coords, Element::Pointer pElement):
			Point(Coords),
			mpElement(pElement)
			{}

		Element::Pointer GetOriginElement()
		{
			return mpElement;   // TEST FUNCTION
		}


	private:
		Element::Pointer mpElement;

	};

	///@name Type Definitions
	///@{

	// Element declarations to be used for the Filter (Determine elements location)
	typedef ModelPart::ElementsContainerType ElementsContainerType;
	typedef ModelPart::ElementsContainerType ElementsArrayType;
	typedef ModelPart::ElementIterator ElementIterator;

	typedef ModelPart::NodesContainerType NodesContainerType;
	typedef ModelPart::NodesContainerType NodesArrayType;
	typedef ModelPart::NodeIterator NodeIterator;
	typedef ModelPart::ConditionsContainerType ConditionsArrayType;

	// For ApplyFilter()
	typedef ElementPositionItem                         PointType;
	typedef ElementPositionItem::Pointer                PointTypePointer;
	typedef std::vector<PointType::Pointer>             ElementPositionVector;
	typedef std::vector<PointType::Pointer>::iterator   ElementPositionIterator;
	typedef std::vector<double>                         DistanceVector;
	typedef std::vector<double>::iterator               DistanceIterator;

	typedef Bucket< 3ul, PointType, ElementPositionVector, PointTypePointer, ElementPositionIterator, DistanceIterator > BucketType;
	typedef Tree< KDTreePartition<BucketType> > tree;

	/// Pointer definition of TopologyFilteringUtilities
	KRATOS_CLASS_POINTER_DEFINITION(TopologyFilteringUtilities);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	TopologyFilteringUtilities( ModelPart& model_part, const double SearchRadius, const int MaxElementsAffected )
	: mrModelPart(model_part),
	  mSearchRadius(SearchRadius),
	  mMaxElementsAffected(MaxElementsAffected)
	{
	}

	/// Destructor.
	virtual ~TopologyFilteringUtilities()
	{
	}


	///@}
	///@name Operators
	///@{


	///@}
	///@name Operations
	///@{

	// ---------------------------------------------------------------------------------------------------------------------------------------------
	// --------------------------------- FILTER SENSITIVITIES  -------------------------------------------------------------------------------------
	// ---------------------------------------------------------------------------------------------------------------------------------------------

	/// With the calculated sensitivities, an additional filter to avoid checkerboard effects is used here
	void ApplyFilter( char FilterType[], char FilterFunctionType[] )
	{

		KRATOS_TRY;

		// Create object of filter function
		FilterFunction FilterFunc(FilterFunctionType, mSearchRadius);

		// Function to Filter Sensitivities
		if ( strcmp( FilterType , "sensitivity" ) == 0 ){
			clock_t begin = clock();
			std::cout << "  Sensitivity filter chosen as filter for sensitivities" << std::endl;

			if ( strcmp( FilterFunctionType , "linear" ) == 0 )
				std::cout << "  Linear filter kernel selected" << std::endl;
			else
				KRATOS_ERROR << "No valid FilterFunction selected for the simulation. Selected one: " << FilterFunctionType << std::endl;

			// IMPORTANT: Tree data structure is re-created on each loop, although this has little impact in simulation time/computation,
			//            it is recommended to declare it in the constructor of the class (i.e. would be calculated just one time)

			// Create placeholder for any element in model. We can assign the center of the element as coordinates to the placeholder
			ElementPositionVector PositionList;
			for(ModelPart::ElementsContainerType::iterator elem_i = mrModelPart.ElementsBegin(); elem_i!=mrModelPart.ElementsEnd(); elem_i++)
			{
				// Find the center of the element
				array_1d<double,3> center_coord = ZeroVector(3);
				Geometry< Node<3> >& geom = elem_i->GetGeometry();
				for(unsigned int i=0; i<geom.size(); i++)
					noalias(center_coord) += geom[i].Coordinates();
				center_coord /= static_cast<double>(geom.size());

				// new "ElementPositionItem" for every element and assigns a pointer to the base element *it.base()
				PointTypePointer pGP = PointTypePointer(new ElementPositionItem( center_coord, *elem_i.base() ) );
				PositionList.push_back( pGP );
			}

			// Creates a tree space search structure - It will use a copy of mGaussPoinList (a std::vector which contains pointers)
			// Note that PositionList will be reordered by the tree for efficiency reasons
			const int BucketSize = 4;
			tree MyTree(PositionList.begin(),PositionList.end(),BucketSize);

			ElementPositionVector Results(mMaxElementsAffected);
			DistanceVector ResultingSquaredDistances(mMaxElementsAffected);

			clock_t tree_time = clock();
			std::cout << "  Filtered tree created                      [ spent time =  " << double(tree_time - begin) / CLOCKS_PER_SEC << " ] " << std::endl;

			// Compute filtered sensitivities
			Vector dcdx_filtered;
			dcdx_filtered.resize(mrModelPart.NumberOfElements());
			int i = 0;
			for(ModelPart::ElementsContainerType::iterator elem_i = mrModelPart.ElementsBegin(); elem_i!=mrModelPart.ElementsEnd(); elem_i++)
			{
				// Find the center of the element
				array_1d<double,3> center_coord = ZeroVector(3);
				Geometry< Node<3> >& geom = elem_i->GetGeometry();
				for(unsigned int i=0; i<geom.size(); i++)
					noalias(center_coord) += geom[i].Coordinates();
				center_coord /= static_cast<double>(geom.size());

				ElementPositionItem ElemPositionItem(center_coord,*elem_i.base());

				// This is broken. Bug found when using ResultingSquaredDistances, so we calculate our own distances
				int num_nodes_found;
				num_nodes_found = MyTree.SearchInRadius(ElemPositionItem,mSearchRadius,Results.begin(),ResultingSquaredDistances.begin(),mMaxElementsAffected);

				double Hxdc = 0.0;
				double Hxdc_sum = 0;
				double H = 0.0;
				double H_sum = 0;
				array_1d<double,3> elemental_distance;
				double distance = 0.0;

				for(int ElementPositionItem_j = 0; ElementPositionItem_j < num_nodes_found; ElementPositionItem_j++)
				{
					// Calculate distances
					elemental_distance = ZeroVector(3);
					elemental_distance = *Results[ElementPositionItem_j] - center_coord;
					distance = std::sqrt(inner_prod(elemental_distance,elemental_distance));

					// Creation of mesh independent convolution operator (weight factors)
					H  = FilterFunc.ComputeWeight(distance);
					Hxdc = H*(Results[ElementPositionItem_j]->GetOriginElement()->GetValue(DCDX))
						    *(Results[ElementPositionItem_j]->GetOriginElement()->GetValue(X_PHYS));
					H_sum    += H;
					Hxdc_sum += Hxdc;
				}

				// Calculate filtered sensitivities and assign to the elements
				dcdx_filtered[i++] = Hxdc_sum / (H_sum*std::max(0.001,elem_i->GetValue(X_PHYS)) );
			}

			// Overwrite sensitivities with filtered sensitivities
			i = 0;
			for(ModelPart::ElementsContainerType::iterator elem_i = mrModelPart.ElementsBegin();
					elem_i!=mrModelPart.ElementsEnd(); elem_i++)
				elem_i->SetValue(DCDX,dcdx_filtered[i++]);

			clock_t end = clock();
			std::cout << "  Filtered sensitivities calculated          [ spent time =  " << double(end - begin) / CLOCKS_PER_SEC << " ] " << std::endl;
		}
		else
			KRATOS_ERROR << "No valid FilterType selected for the simulation. Selected one: " << FilterType << std::endl;

		KRATOS_CATCH("");

	}


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
		return "TopologyFilteringUtilities";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream& rOStream) const
	{
		rOStream << "TopologyFilteringUtilities";
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


	///@}
	///@name Member Variables
	///@{

	ModelPart& mrModelPart;
	const double mSearchRadius;
	const int mMaxElementsAffected;

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
	//TopologyFilteringUtilities& operator=(TopologyFilteringUtilities const& rOther);

	/// Copy constructor.
	//TopologyFilteringUtilities(TopologyFilteringUtilities const& rOther);


	///@}

}; // Class TopologyFilteringUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif	/* KRATOS_TOPOLOGY_FILTERING_UTILITIES_H_INCLUDED */
