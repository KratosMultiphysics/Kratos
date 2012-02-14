/*
==============================================================================
KratosIncompressibleFluidApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2011
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pbecker $
//   Date:                $Date: 2011-09-21 12:30:32 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_BINBASED_FAST_POINT_LOCATOR_INCLUDED )
#define  KRATOS_BINBASED_FAST_POINT_LOCATOR_INCLUDED



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

#include "utilities/spatial_containers_configure.h"



namespace Kratos
{
    ///this class is designed to allow the fast location of MANY points on the top of a 3D mesh.
    ///the utility relies on the creation of a Bin of objects that allows finding quikly a reduced number
    ///of element candidates for the location of a point.
    ///The user should call the function "UpdateSearchDatabase" to mount the bin
    ///and subsequently locate the points as needed
    ///REMARK: the location function is threadsafe, and can be used in OpenMP loops
    template< unsigned int TDim>
    class BinBasedFastPointLocator
    {
    public:

        typedef SpatialContainersConfigure<TDim> Configure;
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

        KRATOS_CLASS_POINTER_DEFINITION(BinBasedFastPointLocator);

        BinBasedFastPointLocator(ModelPart& model_part)
        : mr_model_part(model_part)
        {
        }

        ~BinBasedFastPointLocator()
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
        ///and return a pointer to the element and the vector containing the
        ///shape functions that define the postion within the element
        ///if "false" is devolved the element is not found
        ///REMARK: this function is threadsafe and can be used within OpenMP loops
        bool FindPointOnMesh(const array_1d<double, 3 >& coords,
                array_1d<double, TDim + 1 > & N,
                Element::Pointer& pelement,
                ResultIteratorType result_begin,
                const unsigned int MaxNumberOfResults)
        {
            typedef std::size_t SizeType;

            //ask to the container for the list of candidate elements
            SizeType results_found = mpBinsObjectDynamic->SearchObjectsInCell(coords, result_begin, MaxNumberOfResults);

            if (results_found > 0)
            {
                //loop over the candidate elements and check if the particle falls within
                for (SizeType i = 0; i < results_found; i++)
                {
                    Geometry<Node < 3 > >& geom = (*(result_begin + i))->GetGeometry();

                    //find local position
                    bool is_found = CalculatePosition(geom, coords[0], coords[1], coords[2], N);

                    if (is_found == true)
                    {
                        pelement = (*(result_begin + i));
                        return true;
                    }
                }
            }

            //not found case
            return false;
        }


        ///simplified (less efficient) function to find the element into which a given node is located
        ///and return a pointer to the element and the vector containing the
        ///shape functions that define the postion within the element
        ///if "false" is devolved the element is not found
        ///REMARK: this function is threadsafe and can be used within OpenMP loops
        bool FindPointOnMeshSimplified(const array_1d<double, 3 >& coords,
                Vector& N,
                Element::Pointer& pelement)
        {
            const int max_results = 1000;
	    ResultContainerType results(max_results);
            if(N.size() != TDim+1)
                N.resize(TDim+1,false);

            array_1d<double,TDim+1> aux;

            bool is_found = FindPointOnMesh(coords, aux, pelement, results.begin(), max_results);

            if(is_found == true)
                noalias(N) = aux;

            return is_found;
        }

    protected:


    private:

        //***************************************
        //***************************************
        inline bool CalculatePosition(Geometry<Node < 3 > >&geom,
                const double xc, const double yc, const double zc,
                array_1d<double, 3 > & N
                )
        {
            double x0 = geom[0].X();
            double y0 = geom[0].Y();
            double x1 = geom[1].X();
            double y1 = geom[1].Y();
            double x2 = geom[2].X();
            double y2 = geom[2].Y();

            double area = CalculateVol(x0, y0, x1, y1, x2, y2);
            double inv_area = 0.0;
            if (area == 0.0)
            {
                KRATOS_ERROR(std::logic_error, "element with zero area found", "");
            } else
            {
                inv_area = 1.0 / area;
            }


            N[0] = CalculateVol(x1, y1, x2, y2, xc, yc) * inv_area;
            N[1] = CalculateVol(x2, y2, x0, y0, xc, yc) * inv_area;
            N[2] = CalculateVol(x0, y0, x1, y1, xc, yc) * inv_area;


            if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0) //if the xc yc is inside the triangle return true
                return true;

            return false;
        }

        //***************************************
        //***************************************

        inline bool CalculatePosition(Geometry<Node < 3 > >&geom,
                const double xc, const double yc, const double zc,
                array_1d<double, 4 > & N
                )
        {

            double x0 = geom[0].X();
            double y0 = geom[0].Y();
            double z0 = geom[0].Z();
            double x1 = geom[1].X();
            double y1 = geom[1].Y();
            double z1 = geom[1].Z();
            double x2 = geom[2].X();
            double y2 = geom[2].Y();
            double z2 = geom[2].Z();
            double x3 = geom[3].X();
            double y3 = geom[3].Y();
            double z3 = geom[3].Z();

            double vol = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);

            double inv_vol = 0.0;
            if (vol < 0.0000000000001)
            {
                KRATOS_ERROR(std::logic_error, "element with zero vol found", "");
            } else
            {
                inv_vol = 1.0 / vol;
            }

            N[0] = CalculateVol(x1, y1, z1, x3, y3, z3, x2, y2, z2, xc, yc, zc) * inv_vol;
            N[1] = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, xc, yc, zc) * inv_vol;
            N[2] = CalculateVol(x3, y3, z3, x1, y1, z1, x0, y0, z0, xc, yc, zc) * inv_vol;
            N[3] = CalculateVol(x3, y3, z3, x0, y0, z0, x2, y2, z2, xc, yc, zc) * inv_vol;


            if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[3] >= 0.0 &&
                    N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0 && N[3] <= 1.0)
                //if the xc yc zc is inside the tetrahedron return true
                return true;

            return false;
        }

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
            double x10 = x1 - x0;
            double y10 = y1 - y0;
            double z10 = z1 - z0;

            double x20 = x2 - x0;
            double y20 = y2 - y0;
            double z20 = z2 - z0;

            double x30 = x3 - x0;
            double y30 = y3 - y0;
            double z30 = z3 - z0;

            double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;
            return detJ * 0.1666666666666666666667;
        }

        ModelPart& mr_model_part;

        typename BinsObjectDynamic<Configure>::Pointer mpBinsObjectDynamic;


    };

} // namespace Kratos.

#endif // KRATOS_BINBASED_FAST_POINT_LOCATOR_INCLUDED  defined


