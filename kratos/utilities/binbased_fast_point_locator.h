//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                       license: license.txt
//
//  License:          BSD License
//  Main authors:     Riccardo Rossi
//                    Pablo Becker
//                    Carlos Roig
//                    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_BINBASED_FAST_POINT_LOCATOR_INCLUDED )
#define  KRATOS_BINBASED_FAST_POINT_LOCATOR_INCLUDED

// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "includes/node.h"

#include "spatial_containers/spatial_containers.h"
#include "spatial_containers/bounding_box.h"
#include "spatial_containers/cell.h"
#include "spatial_containers/bins_dynamic_objects.h"

#include "utilities/spatial_containers_configure.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// The size definition
    typedef std::size_t SizeType;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class BinBasedFastPointLocator
 * @ingroup KratosCore
 * @brief This class is designed to allow the fast location of MANY points on the top of a 3D mesh.
 * @details The utility relies on the creation of a Bin of objects that allows finding quikly a reduced number of element candidates for the location of a point.
 * The basic idea is to allow finding the element in which a given spatial position sits
 * The user should call the function "UpdateSearchDatabase" to mount the bin and subsequently locate the points as needed
 * @author  Riccardo Rossi <rrossi@cimne.upc.edu>
 * @note The location function is threadsafe, and can be used in OpenMP loops
 * @tparam TDim If we work in a 2D or 3D space
 * @tparam ConfigureType The spatial container
 */
template< SizeType TDim, class ConfigureType = SpatialContainersConfigure<TDim> >
class BinBasedFastPointLocator
{
public:
    ///@name Type Definitions
    ///@{

    /// The configure type
    typedef ConfigureType Configure;

    /// The definition of the different containers
    typedef typename Configure::PointType PointType;
    typedef typename Configure::ContainerType ContainerType;
    typedef typename Configure::IteratorType IteratorType;
    typedef typename Configure::ResultContainerType ResultContainerType;
    typedef typename Configure::ResultIteratorType ResultIteratorType;

    /// The definition of the node
    typedef Node<3> NodeType;

    /// The definition of the geometry
    typedef Geometry<NodeType> GeometryType;

    /// The index definition
    typedef std::size_t IndexType;

    /// Pointer definition of BinBasedFastPointLocator
    KRATOS_CLASS_POINTER_DEFINITION(BinBasedFastPointLocator);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief This is the default constructor
     * @param rModelPart The model part of the mesh used in the search
     */
    BinBasedFastPointLocator(ModelPart& rModelPart)
        : mrModelPart(rModelPart)
    {
    }

    /// Destructor.
    virtual ~BinBasedFastPointLocator() = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Function to construct or update the search database
     */
    void UpdateSearchDatabase()
    {
        KRATOS_TRY

        // Copy the elements to a new container, as the list will
        //be shuffled duringthe construction of the tree
        ContainerType& rElements = mrModelPart.ElementsArray();
        IteratorType it_begin = rElements.begin();
        IteratorType it_end = rElements.end();

        typename BinsObjectDynamic<Configure>::Pointer paux = typename BinsObjectDynamic<Configure>::Pointer(new BinsObjectDynamic<Configure > (it_begin, it_end));

        paux.swap(mpBinsObjectDynamic);

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to construct or update the search database
     * @param CellSize The current size of the cell used for search
     */
    void UpdateSearchDatabaseAssignedSize(double CellSize)
    {
        KRATOS_TRY

        //copy the elements to a new container, as the list will
        //be shuffled duringthe construction of the tree
        ContainerType& rElements = mrModelPart.ElementsArray();
        IteratorType it_begin = rElements.begin();
        IteratorType it_end = rElements.end();

        typename BinsObjectDynamic<Configure>::Pointer paux = typename BinsObjectDynamic<Configure>::Pointer(new BinsObjectDynamic<Configure > (it_begin, it_end,CellSize));
        paux.swap(mpBinsObjectDynamic);

        KRATOS_CATCH("")
    }

    /**
     * @brief This function should find the element into which a given node is located
     * and return a pointer to the element and the vector containing the
     * shape functions that define the postion within the element
     * @param rCoordinates The vector containign the coordinates of the point to be searched
     * @param rNShapeFunction The vector containing the shape function of the located point
     * @param pElement The pointer to the element containing the located point
     * @param ItResultBegin The iterator of the search
     * @param MaxNumberOfResults The max number of results to be considered
     * @param Tolerance The tolerance considered on the search
     * @return If "false" is devolved the element is not found
     * @note this function is threadsafe and can be used within OpenMP loops
     */
    bool FindPointOnMesh(
        const array_1d<double, 3 >& rCoordinates,
        array_1d<double, TDim + 1 > & rNShapeFunction,
        Element::Pointer& pElement,
        ResultIteratorType ItResultBegin,
        const SizeType MaxNumberOfResults = 1000,
        const double Tolerance = 1.0e-5
        )
    {
        // Ask to the container for the list of candidate elements
        SizeType results_found = mpBinsObjectDynamic->SearchObjectsInCell(rCoordinates, ItResultBegin, MaxNumberOfResults);

        if (results_found > 0) {
            // Loop over the candidate elements and check if the particle falls within
            for (IndexType i = 0; i < results_found; i++) {
                GeometryType& geom = (*(ItResultBegin + i))->GetGeometry();

                // Find local position
                bool is_found = CalculatePosition(geom, rCoordinates[0], rCoordinates[1], rCoordinates[2], rNShapeFunction, Tolerance);

                if (is_found) {
                    pElement = (*(ItResultBegin + i));
                    return true;
                }
            }
        }

        // Not found case
        return false;
    }

    /**
     * @brief Simplified (less efficient) function to find the element into which a given node is located and return a pointer to the element and the vector containing the shape functions that define the postion within the element
     * @param rCoordinates The vector containign the coordinates of the point to be searched
     * @param rNShapeFunction The vector containing the shape function of the located point
     * @param pElement The pointer to the element containing the located point
     * @param MaxNumberOfResults The max number of results to be considered
     * @param Tolerance The tolerance considered on the search
     * @return If "false" is devolved the element is not found
     * @note this function is threadsafe and can be used within OpenMP loops
     */
    bool FindPointOnMeshSimplified(
        const array_1d<double, 3 >& rCoordinates,
        Vector& rNShapeFunction,
        Element::Pointer& pElement,
        const SizeType MaxNumberOfResults = 1000,
        const double Tolerance = 1.0e-5
        )
    {
        ResultContainerType results(MaxNumberOfResults);

        if(rNShapeFunction.size() != TDim+1) {
            rNShapeFunction.resize(TDim+1,false);
        }

        array_1d<double,TDim+1> aux;

        const bool is_found = FindPointOnMesh(rCoordinates, aux, pElement, results.begin(), MaxNumberOfResults, Tolerance);

        if(is_found) {
            noalias(rNShapeFunction) = aux;
        }

        return is_found;
    }
    
    //***************************************
    //***************************************

    inline bool CalculatePosition(
        GeometryType& geom,
        const double xc, 
        const double yc, 
        const double zc,
        array_1d<double, 3 > & N,
        const double Tolerance = 1.0e-5
        )
    {
        const double x0 = geom[0].X();
        const double y0 = geom[0].Y();
        const double x1 = geom[1].X();
        const double y1 = geom[1].Y();
        const double x2 = geom[2].X();
        const double y2 = geom[2].Y();

        const double area = CalculateVol(x0, y0, x1, y1, x2, y2);

        double inv_area = 0.0;
        if (area == 0.0) {
            KRATOS_ERROR << "Element with zero area found with the current geometry " << geom << std::endl;
        } else {
            inv_area = 1.0 / area;
        }

        N[0] = CalculateVol(x1, y1, x2, y2, xc, yc) * inv_area;
        N[1] = CalculateVol(x2, y2, x0, y0, xc, yc) * inv_area;
        N[2] = CalculateVol(x0, y0, x1, y1, xc, yc) * inv_area;

        if ((N[0] >= -Tolerance) && (N[1] >= -Tolerance) && (N[2] >= -Tolerance) &&
            (N[0] <= (1.0 + Tolerance)) && (N[1] <= (1.0 + Tolerance)) && (N[2] <= (1.0 + Tolerance))) //if the xc yc is inside the triangle return true
        {
            return true;
        }

        return false;
    }

    //***************************************
    //***************************************

    inline bool CalculatePosition(
        GeometryType& geom,
        const double xc, 
        const double yc, 
        const double zc,
        array_1d<double, 4 > & N,
        const double Tolerance = 1.0e-5
        )
    {
        const double x0 = geom[0].X();
        const double y0 = geom[0].Y();
        const double z0 = geom[0].Z();
        const double x1 = geom[1].X();
        const double y1 = geom[1].Y();
        const double z1 = geom[1].Z();
        const double x2 = geom[2].X();
        const double y2 = geom[2].Y();
        const double z2 = geom[2].Z();
        const double x3 = geom[3].X();
        const double y3 = geom[3].Y();
        const double z3 = geom[3].Z();

        const double vol = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);

        double inv_vol = 0.0;
        if (vol == 0.0) {
            KRATOS_ERROR << "Element with zero area found with the current geometry " << geom << std::endl;
        } else {
            inv_vol = 1.0 / vol;
        }

        N[0] = CalculateVol(x1,y1,z1,x3,y3,z3,x2,y2,z2,xc,yc,zc) * inv_vol;
        N[1] = CalculateVol(x3,y3,z3,x0,y0,z0,x2,y2,z2,xc,yc,zc) * inv_vol;
        N[2] = CalculateVol(x3,y3,z3,x1,y1,z1,x0,y0,z0,xc,yc,zc) * inv_vol;
        N[3] = CalculateVol(x0,y0,z0,x1,y1,z1,x2,y2,z2,xc,yc,zc) * inv_vol;

        if ((N[0] >= -Tolerance) && (N[1] >= -Tolerance) && (N[2] >= -Tolerance) && (N[3] >= -Tolerance) &&
            (N[0] <= (1.0 + Tolerance)) && (N[1] <= (1.0 + Tolerance)) && (N[2] <= (1.0 + Tolerance)) && (N[3] <= (1.0 + Tolerance)))
            //if the xc yc zc is inside the tetrahedron return true
        {
            return true;
        }

        return false;
    }

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

    ModelPart& mrModelPart; /// The model part containing the mesh for the search

    typename BinsObjectDynamic<Configure>::Pointer mpBinsObjectDynamic; /// The pointer of the bins used for the search

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    inline double CalculateVol(const double x0, const double y0,
                               const double x1, const double y1,
                               const double x2, const double y2
                              )
    {
        return 0.5 * ((x1 - x0)*(y2 - y0)- (y1 - y0)*(x2 - x0));
    }

    //***************************************
    //***************************************

    inline double CalculateVol(const double x0, const double y0, const double z0,
                               const double x1, const double y1, const double z1,
                               const double x2, const double y2, const double z2,
                               const double x3, const double y3, const double z3
                              )
    {
        const double x10 = x1 - x0;
        const double y10 = y1 - y0;
        const double z10 = z1 - z0;

        const double x20 = x2 - x0;
        const double y20 = y2 - y0;
        const double z20 = z2 - z0;

        const double x30 = x3 - x0;
        const double y30 = y3 - y0;
        const double z30 = z3 - z0;

        double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;
        return detJ * 0.1666666666666666666667;
    }

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
};

} // namespace Kratos.

#endif // KRATOS_BINBASED_FAST_POINT_LOCATOR_INCLUDED  defined


