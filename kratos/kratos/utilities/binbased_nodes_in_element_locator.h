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


#if !defined(KRATOS_BINBASED_NODES_IN_ELEMENT_LOCATOR_INCLUDED )
#define  KRATOS_BINBASED_NODES_IN_ELEMENT_LOCATOR_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/node.h"

#include "spatial_containers/spatial_containers.h"



namespace Kratos
{
///This class is designed to allow the fast location of the nodes of a fixed mesh with respect to a moving mesh.
///The utility relies on the creation of a static Bin that allows finding quikly the nodes of a fixed mesh that
///are inside each of the elements of the moving one.
///After the creation of the "BinBasedNodesInElementLocator",
///the user should call the function "UpdateSearchDatabase" or"UpdateSearchDatabaseAssignedSize(hmin)"  to mount the bin
///and subsequently locate the points as needed
///An application of this utility can be found in 

///REMARK: the location function is threadsafe, and can be used in OpenMP loops
template< unsigned int TDim>
class BinBasedNodesInElementLocator
{
public:

    typedef Node<3> PointType;
    typedef Node<3>::Pointer PointTypePointer;
    typedef std::vector<PointType::Pointer >           PointVector;
    typedef typename std::vector<PointType::Pointer >::iterator PointIterator;
    typedef std::vector<double>               DistanceVector;
    typedef typename std::vector<double>::iterator     DistanceIterator;

    // Bucket types
    typedef Bins< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > StaticBins;
//     typedef Tree< StaticBins > tree; 		     		//Binstree;

    KRATOS_CLASS_POINTER_DEFINITION(BinBasedNodesInElementLocator);

    BinBasedNodesInElementLocator(ModelPart& model_part)
            : mr_model_part(model_part)
    {
    }

    ~BinBasedNodesInElementLocator()
    {
    }

    ///Function to construct or update the search database
    void UpdateSearchDatabase()
    {
        KRATOS_TRY

        mlist_of_new_nodes.clear();
        for (ModelPart::NodesContainerType::iterator node_it = mr_model_part.NodesBegin();
                node_it != mr_model_part.NodesEnd(); ++node_it)
        {
            //PointType::Pointer pnode(new PointType(*node_it));
            Node<3>::Pointer pnode = *(node_it.base());

            //putting the nodes of the destination_model part in an auxiliary list
            mlist_of_new_nodes.push_back( pnode );
        }

        //create a spatial database with the list of new nodes
        typename StaticBins::SizeType bucket_size = 20;

        typename StaticBins::Pointer paux = typename StaticBins::Pointer(new StaticBins(mlist_of_new_nodes.begin(),mlist_of_new_nodes.end(),bucket_size) );
        paux.swap(mp_search_structure); //the class mp_search_structure remains in the memory, while the paux just constucted is destryed.

        KRATOS_CATCH("")
    }

    /// Function to construct or update the search database
    /// The cell size is requested as input parameter. This action reduces consideribly the searching procedure.
    /// One possibility is to give the dimension of the mesh.
    /**
    * @param CellSize: the bin cell size
    */
    void UpdateSearchDatabaseAssignedSize(double CellSize)
    {
        KRATOS_TRY

        mlist_of_new_nodes.clear();
        for (ModelPart::NodesContainerType::iterator node_it = mr_model_part.NodesBegin();
                node_it != mr_model_part.NodesEnd(); ++node_it)
        {
            //PointType::Pointer pnode(new PointType(*node_it));
            Node<3>::Pointer pnode = *(node_it.base());

            //putting the nodes of the destination_model part in an auxiliary list
            mlist_of_new_nodes.push_back( pnode );
        }

        typename StaticBins::Pointer paux = typename StaticBins::Pointer(new StaticBins(mlist_of_new_nodes.begin(),mlist_of_new_nodes.end(),CellSize) );
        paux.swap(mp_search_structure); //the class mp_search_structure remains in the memory, while the paux just constucted is destryed.

        KRATOS_CATCH("")
    }

    ///function to find all teh nodes of a fixed mesh contained in the elements of a moving mesh
    /**
    *	It is called each time a projection from a moving to a fixed mesh is performed using the function of @see binbased_projection.h
    */
    unsigned int FindNodesInElement(Element::Pointer& pelement,
                                    boost::numeric::ublas::vector<int>& positions,
                                    Matrix& Nmat,
                                    const unsigned int& max_results,
				    PointIterator work_results,
				    DistanceIterator work_distances,
				    Node<3>& work_point  // is it needed from outside???????????????
                                   )
    {
	Geometry<Node<3> >&geom = pelement->GetGeometry();
	array_1d<double,TDim+1> N;
				
	double xc, yc, zc,  radius;	
	CalculateCenterAndSearchRadius( geom, xc,yc,zc, radius,N);  
	work_point.X() = xc; work_point.Y() = yc; work_point.Z() = zc;
	
	//find all of the new nodes within the radius
	int number_of_points_in_radius;

	//look between the new nodes which of them is inside the radius of the circumscribed cyrcle
	number_of_points_in_radius = mp_search_structure->SearchInRadius(work_point, radius, work_results, work_distances, max_results);

	//check if inside 
	unsigned int counter = 0;
	
	for( PointIterator it_found = work_results; it_found != work_results + number_of_points_in_radius; it_found++)			
	{	
 		bool is_inside = false;
		is_inside = CalculatePosition(geom,	(*it_found)->X(),(*it_found)->Y(),(*it_found)->Z(),N);

		//if the node falls inside the element interpolate
		if(is_inside == true)
		{
		    positions[counter] = it_found - work_results;
		    for(unsigned int k=0; k<TDim+1;k++)
		      Nmat(counter,k) = N[k];
		    counter++;
		}
	}
	return counter;
    }

protected:


private:

	inline void CalculateCenterAndSearchRadius(Geometry<Node<3> >&geom,
					double& xc, double& yc, double& zc, double& R, array_1d<double,3>& N		
					)
	{
		  double x0 = geom[0].X();double  y0 = geom[0].Y(); 
		  double x1 = geom[1].X();double  y1 = geom[1].Y(); 
		  double x2 = geom[2].X();double  y2 = geom[2].Y(); 
		  

		xc = 0.3333333333333333333*(x0+x1+x2);
		yc = 0.3333333333333333333*(y0+y1+y2);
		zc = 0.0;

		double R1 = (xc-x0)*(xc-x0) + (yc-y0)*(yc-y0);
		double R2 = (xc-x1)*(xc-x1) + (yc-y1)*(yc-y1);
		double R3 = (xc-x2)*(xc-x2) + (yc-y2)*(yc-y2);
		
		R = R1;
		if(R2 > R) R = R2;
		if(R3 > R) R = R3;
		
		R = 1.01 * sqrt(R);
	}
	//***************************************
	//***************************************
	inline void CalculateCenterAndSearchRadius(Geometry<Node<3> >&geom,
					double& xc, double& yc, double& zc, double& R, array_1d<double,4>& N		

					)
	{
		  double x0 = geom[0].X();double  y0 = geom[0].Y();double  z0 = geom[0].Z();
		  double x1 = geom[1].X();double  y1 = geom[1].Y();double  z1 = geom[1].Z();
		  double x2 = geom[2].X();double  y2 = geom[2].Y();double  z2 = geom[2].Z();
		  double x3 = geom[3].X();double  y3 = geom[3].Y();double  z3 = geom[3].Z();	


		xc = 0.25*(x0+x1+x2+x3);
		yc = 0.25*(y0+y1+y2+y3);
		zc = 0.25*(z0+z1+z2+z3);			 

		double R1 = (xc-x0)*(xc-x0) + (yc-y0)*(yc-y0) + (zc-z0)*(zc-z0);
		double R2 = (xc-x1)*(xc-x1) + (yc-y1)*(yc-y1) + (zc-z1)*(zc-z1);
		double R3 = (xc-x2)*(xc-x2) + (yc-y2)*(yc-y2) + (zc-z2)*(zc-z2);
		double R4 = (xc-x3)*(xc-x3) + (yc-y3)*(yc-y3) + (zc-z3)*(zc-z3);
		
		R = R1;
		if(R2 > R) R = R2;
		if(R3 > R) R = R3;
		if(R4 > R) R = R4;
		  
		R = sqrt(R);
	}  
  
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
        }
        else
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
        }
        else
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


    PointVector mlist_of_new_nodes;

    ModelPart& mr_model_part;

    typename StaticBins::Pointer mp_search_structure;


};

} // namespace Kratos.

#endif // KRATOS_BINBASED_NODES_IN_ELEMENT_LOCATOR_INCLUDED  defined


