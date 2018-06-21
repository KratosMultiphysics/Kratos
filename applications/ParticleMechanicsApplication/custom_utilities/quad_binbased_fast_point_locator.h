//   Project Name:        KratosParticleMechanicsApplication $
//   Last modified by:    $Author:                 ilaria $
//   Date:                $Date:                April 2016 $
//   Revision:            $Revision:                  0.0 $
//
//


#if !defined(KRATOS_QUAD_BINBASED_FAST_POINT_LOCATOR_INCLUDED )
#define  KRATOS_QUAD_BINBASED_FAST_POINT_LOCATOR_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/node.h"

#include "spatial_containers/spatial_containers.h"
#include "spatial_containers/bounding_box.h"
#include "spatial_containers/cell.h"
#include "spatial_containers/bins_dynamic_objects.h"

//#include "utilities/spatial_containers_configure.h"

#include "custom_utilities/fast_quad_spatial_containers_configure.h"


namespace Kratos
{
///This class is designed to allow the fast location of MANY points on the top of a 3D mesh.
/** @author  Riccardo Rossi <rrossi@cimne.upc.edu>
*
*The utility relies on the creation of a Bin of objects that allows finding quikly a reduced number
*of element candidates for the location of a point.
*
*The basic idea is to allow finding the element in which a given spatial position sits
*
*The user should call the function "UpdateSearchDatabase" to mount the bin
*and subsequently locate the points as needed
*
*REMARK: the location function is threadsafe, and can be used in OpenMP loops
*/

template< unsigned int TDim, class ConfigureType = FastQuadSpatialContainersConfigure<TDim> >

class QuadBinBasedFastPointLocator
{
public:

    typedef ConfigureType Configure;

    typedef typename Configure::PointType PointType;
    //typedef PointType::CoordinatesArrayType           CoordinatesArrayType;
    typedef typename Configure::ContainerType ContainerType;
    //typedef Configure::PointerType                    PointerType;
    typedef typename Configure::IteratorType IteratorType;
    typedef typename Configure::ResultContainerType ResultContainerType;
    //typedef Configure::ResultPointerType              ResultPointerType;
    typedef typename Configure::ResultIteratorType ResultIteratorType;
    //typedef Configure::ContactPairType                ContactPairType;
    //typedef Configure::ContainerContactType           ContainerContactType;
    //typedef Configure::IteratorContactType            IteratorContactType;
    //typedef Configure::PointerContactType             PointerContactType;
    //typedef Configure::PointerTypeIterator            PointerTypeIterator;

    KRATOS_CLASS_POINTER_DEFINITION(QuadBinBasedFastPointLocator);

    QuadBinBasedFastPointLocator(ModelPart& model_part)
        : mr_model_part(model_part)
    {
    }

    ~QuadBinBasedFastPointLocator()
    {
    }

    ///function to construct or update the search database
    void UpdateSearchDatabase()
    {
        KRATOS_TRY
        //copy the elements to a new container, as the list will
        //be shuffled duringthe construction of the tree
        ContainerType& rElements = mr_model_part.ElementsArray();
        IteratorType it_begin = rElements.begin();
        IteratorType it_end = rElements.end();

        typename BinsObjectDynamic<Configure>::Pointer paux = typename BinsObjectDynamic<Configure>::Pointer(new BinsObjectDynamic<Configure > (it_begin, it_end));

        paux.swap(mpBinsObjectDynamic);

        KRATOS_CATCH("")
    }

    ///function to construct or update the search database
    void UpdateSearchDatabaseAssignedSize(double CellSize)
    {
        KRATOS_TRY

        //copy the elements to a new container, as the list will
        //be shuffled duringthe construction of the tree
        ContainerType& rElements = mr_model_part.ElementsArray();
        IteratorType it_begin = rElements.begin();
        IteratorType it_end = rElements.end();

        typename BinsObjectDynamic<Configure>::Pointer paux = typename BinsObjectDynamic<Configure>::Pointer(new BinsObjectDynamic<Configure > (it_begin, it_end,CellSize));
        paux.swap(mpBinsObjectDynamic);

        KRATOS_CATCH("")
    }
    ///this function should find the element into which a given node is located
    ///and return a pointer to the element
    ///if "false" is devolved the element is not found
    ///REMARK: this function is threadsafe and can be used within OpenMP loops
    bool FindPointOnMesh(const array_1d<double, 3 >& coords,
                        Element::Pointer& pelement,
                        ResultIteratorType result_begin,
                        const unsigned int MaxNumberOfResults = 1000)
    {
        typedef std::size_t SizeType;

        PointType  Result = coords;

        // Ask to the container for the list of candidate elements
        SizeType results_found = mpBinsObjectDynamic->SearchObjectsInCell(coords, result_begin, MaxNumberOfResults);

        if (results_found > 0)
        {
            // Loop over the candidate elements and check if the particle falls within
            for (SizeType i = 0; i < results_found; i++)
            {
                Geometry<Node < 3 > >& geom = (*(result_begin + i))->GetGeometry();
                
                bool is_found = geom.IsInside(coords, Result);
                
                if (is_found == true)
                {
                    pelement = (*(result_begin + i));
                    return true;
                }
            }
        }

        // Not found case
        return false;
    }


protected:


private:

    ModelPart& mr_model_part;

    typename BinsObjectDynamic<Configure>::Pointer mpBinsObjectDynamic;


};

} // namespace Kratos.

#endif // KRATOS_BINBASED_FAST_POINT_LOCATOR_INCLUDED  defined


