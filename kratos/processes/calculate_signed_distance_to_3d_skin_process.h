/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2007-10-31 17:51:34 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_CALCULATE_DISTANCE_PROCESS_H_INCLUDED )
#define  KRATOS_CALCULATE_DISTANCE_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "spatial_containers/octree_binary.h"
#include "utilities/spatial_containers_configure.h"
#include "utilities/timer.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "geometries/triangle_3d_3.h"
#include "utilities/body_normal_calculation_utils.h"
#include "utilities/math_utils.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#ifdef _OPENMP
#include "omp.h"
#endif

using namespace boost::numeric::ublas;


namespace Kratos
{

class DistanceSpatialContainersConfigure
{

    class CellNodeData
    {
        double mDistance;
        double mCoordinates[3];
        std::size_t mId;
    public:
        double& Distance(){return mDistance;}
        double& X() {return mCoordinates[0];}
        double& Y() {return mCoordinates[1];}
        double& Z() {return mCoordinates[2];}
        double& Coordinate(int i) {return mCoordinates[i-1];}
        std::size_t& Id(){return mId;}
    };

public:


    ///@name Type Definitions
    ///@{

    enum { Dimension = 3,
           DIMENSION = 3,
           MAX_LEVEL = 12,
           MIN_LEVEL = 2    // this cannot be less than 2!!!
         };

    typedef Point<3, double>                                PointType;  /// always the point 3D
    typedef std::vector<double>::iterator                   DistanceIteratorType;
    typedef ModelPart::ElementsContainerType::ContainerType ContainerType;
    typedef ContainerType::value_type                       PointerType;
    typedef ContainerType::iterator                         IteratorType;
    typedef ModelPart::ElementsContainerType::ContainerType ResultContainerType;
    typedef ResultContainerType::value_type                 ResultPointerType;
    typedef ResultContainerType::iterator                   ResultIteratorType;

    typedef Element::Pointer                                        pointer_type;
    typedef CellNodeData                cell_node_data_type;
    typedef std::vector<CellNodeData*> data_type;

    typedef std::vector<PointerType>::iterator             PointerTypeIterator;




    /// Pointer definition of DistanceSpatialContainersConfigure
    KRATOS_CLASS_POINTER_DEFINITION(DistanceSpatialContainersConfigure);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DistanceSpatialContainersConfigure() {}

    /// Destructor.
    virtual ~DistanceSpatialContainersConfigure() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    static data_type* AllocateData() {
        return new data_type(27, NULL);
    }

    static void CopyData(data_type* source, data_type* destination) {
        *destination = *source;
    }

    static void DeleteData(data_type* data) {
        delete data;
    }

    static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint)
    {
        rHighPoint = rObject->GetGeometry().GetPoint(0);
        rLowPoint  = rObject->GetGeometry().GetPoint(0);
        for (unsigned int point = 0; point<rObject->GetGeometry().PointsNumber(); point++)
        {
            for(std::size_t i = 0; i<3; i++)
            {
                rLowPoint[i]  =  (rLowPoint[i]  >  rObject->GetGeometry().GetPoint(point)[i] ) ?  rObject->GetGeometry().GetPoint(point)[i] : rLowPoint[i];
                rHighPoint[i] =  (rHighPoint[i] <  rObject->GetGeometry().GetPoint(point)[i] ) ?  rObject->GetGeometry().GetPoint(point)[i] : rHighPoint[i];
            }
        }
    }

    static inline void GetBoundingBox(const PointerType rObject, double* rLowPoint, double* rHighPoint)
    {

        for(std::size_t i = 0; i<3; i++)
        {
            rLowPoint[i]  =  rObject->GetGeometry().GetPoint(0)[i];
            rHighPoint[i] =  rObject->GetGeometry().GetPoint(0)[i];
        }

        for (unsigned int point = 0; point<rObject->GetGeometry().PointsNumber(); point++)
        {
            for(std::size_t i = 0; i<3; i++)
            {
                rLowPoint[i]  =  (rLowPoint[i]  >  rObject->GetGeometry().GetPoint(point)[i] ) ?  rObject->GetGeometry().GetPoint(point)[i] : rLowPoint[i];
                rHighPoint[i] =  (rHighPoint[i] <  rObject->GetGeometry().GetPoint(point)[i] ) ?  rObject->GetGeometry().GetPoint(point)[i] : rHighPoint[i];
            }
        }
    }

    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2)
    {
        Element::GeometryType& geom_1 = rObj_1->GetGeometry();
        Element::GeometryType& geom_2 = rObj_2->GetGeometry();
        return  geom_1.HasIntersection(geom_2);

    }


    static inline bool  IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint)
    {
        return rObject->GetGeometry().HasIntersection(rLowPoint, rHighPoint);
    }


    static  inline bool  IsIntersected(const Element::Pointer rObject, double Tolerance, const double* rLowPoint, const double* rHighPoint)
    {
        Point<3,double> low_point(rLowPoint[0] - Tolerance, rLowPoint[1] - Tolerance, rLowPoint[2] - Tolerance);
        Point<3,double> high_point(rHighPoint[0] + Tolerance, rHighPoint[1] + Tolerance, rHighPoint[2] + Tolerance);

        KRATOS_ERROR(std::logic_error, "Not Implemented method", "")
                //return HasIntersection(rObject->GetGeometry(), low_point, high_point);
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
        return " Spatial Containers Configure";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


    ///@}

protected:

private:

    /// Assignment operator.
    DistanceSpatialContainersConfigure& operator=(DistanceSpatialContainersConfigure const& rOther);

    /// Copy constructor.
    DistanceSpatialContainersConfigure(DistanceSpatialContainersConfigure const& rOther);


}; // Class DistanceSpatialContainersConfigure

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
class CalculateSignedDistanceTo3DSkinProcess
        : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CalculateSignedDistanceTo3DSkinProcess
    KRATOS_CLASS_POINTER_DEFINITION(CalculateSignedDistanceTo3DSkinProcess);

    typedef DistanceSpatialContainersConfigure ConfigurationType;
    typedef OctreeBinaryCell<ConfigurationType> CellType;
    typedef OctreeBinary<CellType> OctreeType;
    typedef DistanceSpatialContainersConfigure::cell_node_data_type CellNodeDataType;
    typedef Point<3, double> PointType;  /// always the point 3D
    typedef OctreeType::cell_type::object_container_type object_container_type;
    typedef struct{
        array_1d<double,3>  Coordinates;
        array_1d<double,3>  StructElemNormal;
        unsigned int EdgeNode1;
        unsigned int EdgeNode2;
    }IntersectionNodeStruct;
    typedef struct{
        std::vector<IntersectionNodeStruct> IntNodes;
    }TetEdgeStruct;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    CalculateSignedDistanceTo3DSkinProcess(ModelPart& rThisModelPartStruc, ModelPart& rThisModelPartFluid)
        : mrSkinModelPart(rThisModelPartStruc), mrBodyModelPart(rThisModelPartStruc), mrFluidModelPart(rThisModelPartFluid)
    {
    }

    /// Destructor.
    virtual ~CalculateSignedDistanceTo3DSkinProcess()
    {
    }


    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    virtual void Execute()
    {
        KRATOS_TRY

                GenerateOctree();

        DistanceFluidStructure();

        //        for( ModelPart::ElementIterator i_fluidElement = mrFluidModelPart.ElementsBegin();
        //                                        i_fluidElement != mrFluidModelPart.ElementsEnd();
        //                                        i_fluidElement++)
        //        {
        //            Geometry< Node<3> >& fluidGeom = i_fluidElement->GetGeometry();
        //            fluidGeom[0].GetSolutionStepValue(VELOCITY,0) = ZeroVector(3);
        //            fluidGeom[1].GetSolutionStepValue(VELOCITY,0) = ZeroVector(3);
        //            fluidGeom[2].GetSolutionStepValue(VELOCITY,0) = ZeroVector(3);
        //            fluidGeom[3].GetSolutionStepValue(VELOCITY,0) = ZeroVector(3);
        //        }

        //        for( ModelPart::ElementIterator i_fluidElement = mrFluidModelPart.ElementsBegin();
        //                                        i_fluidElement != mrFluidModelPart.ElementsEnd();
        //                                        i_fluidElement++)
        //        {
        //            array_1d<double,3> vel_elem = i_fluidElement->GetValue(EMBEDDED_VELOCITY);

        //            if(i_fluidElement->GetValue(SPLIT_ELEMENT) == true)
        //            {
        //                KRATOS_WATCH(vel_elem);

        //                Geometry< Node<3> >& fluidGeom = i_fluidElement->GetGeometry();
        //                fluidGeom[0].GetSolutionStepValue(VELOCITY,0) = vel_elem;
        //                fluidGeom[1].GetSolutionStepValue(VELOCITY,0) = vel_elem;
        //                fluidGeom[2].GetSolutionStepValue(VELOCITY,0) = vel_elem;
        //                fluidGeom[3].GetSolutionStepValue(VELOCITY,0) = vel_elem;
        //            }
        //        }

        //          ------------------------------------------------------------------
        //          GenerateNodes();
        //          CalculateDistance2(); // I have to change this. Pooyan.
        //          mrSkinModelPart.GetCommunicator().AssembleCurrentData(DISTANCE);
        //          std::ofstream mesh_file1("octree1.post.msh");
        //          std::ofstream res_file("octree1.post.res");
        //          Timer::Start("Writing Gid conform Mesh");
        //          PrintGiDMesh(mesh_file1);
        //          PrintGiDResults(res_file);
        //          octree.PrintGiDMeshNew(mesh_file2);
        //          Timer::Stop("Writing Gid conform Mesh");
        //          delete octree. TODO: Carlos
        //          ------------------------------------------------------------------

        KRATOS_CATCH("");
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************
    /**
       * This function maps the nodal pressure values computed in the CFD analysis to the respective
       * structural nodes, i.e. for each structural node inside a fluid tetrahedra positive and negative
       * face pressure is computed by mapping between the nodal values of the tetrahedra. Afterwards
       * the resulting delta is applied as new nodal pressure.
       */
    void MappingPressureToStructure(const double conditionVariable)
    {
        unsigned int aux_var_pos = 0; // to be deleted, just to test with some outputs
        unsigned int aux_var_neg = 0;

        // loop over all fluid elements
        // this loop is parallelized using openmp
#ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
#else
        int number_of_threads = 1;
#endif

        ModelPart::ElementsContainerType& pElements = mrFluidModelPart.Elements();

        vector<unsigned int> element_partition;
        CreatePartition(number_of_threads, pElements.size(), element_partition);
        KRATOS_WATCH(number_of_threads);
        KRATOS_WATCH(element_partition);

#pragma omp parallel for
        for (int k = 0; k < number_of_threads; k++)
        {

            ModelPart::ElementsContainerType::iterator it_begin = pElements.ptr_begin() + element_partition[k];
            ModelPart::ElementsContainerType::iterator it_end = pElements.ptr_begin() + element_partition[k+1];

            for(ModelPart::ElementIterator i_fluidElement = it_begin; i_fluidElement != it_end; ++i_fluidElement)
            {
                const array_1d<double,4> elementalDistances = i_fluidElement->GetValue(ELEMENTAL_DISTANCES);
                array_1d<double,4> Npos;
                array_1d<double,4> Nneg;
                Geometry< Node<3> >& fluidGeom = i_fluidElement->GetGeometry();
                array_1d<double,4> nodalPressures;
                std::vector<OctreeType::cell_type*> leaves;

                // Get leaves of octree intersecting with fluid element
                mOctree.GetIntersectedLeaves(*(i_fluidElement).base(),leaves);

                // loop over all octree cells which are intersected by the fluid element
                for(unsigned int i_cell = 0 ; i_cell < leaves.size() ; i_cell++)
                {
                    // Structural element contained in one cell of the octree
                    object_container_type* struct_elem = (leaves[i_cell]->pGetObjects());

                    // loop over all structural elements within each octree cell
                    for(object_container_type::iterator i_StructElement = struct_elem->begin(); i_StructElement != struct_elem->end(); i_StructElement++)
                    {
                        // check if a node of the structure lies within the tetrahedra
                        // --> if yes, then do the pressure mapping
                        Geometry< Node<3> >& structGeom = (*i_StructElement)->GetGeometry();

                        // loop over all nodes of the structural element
                        for(unsigned int i_structElementNode = 0; i_structElementNode < 3 ; i_structElementNode++)
                        {
                            // Comparison of volumes to identify whether a structure points lies inside a tetrahedra or not
                            double Vol = fluidGeom.Volume();
                            double Vol1;
                            double Vol2;
                            double Vol3;
                            double Vol4;

                            Node<3>& S = structGeom[i_structElementNode];
                            Node<3>& F1 = fluidGeom.GetPoint(0);
                            Node<3>& F2 = fluidGeom.GetPoint(1);
                            Node<3>& F3 = fluidGeom.GetPoint(2);
                            Node<3>& F4 = fluidGeom.GetPoint(3);

                            CalculateJacobian(F1,F2,F3,S,Vol1);
                            CalculateJacobian(F1,F2,S,F4,Vol2);
                            CalculateJacobian(F1,S,F3,F4,Vol3);
                            CalculateJacobian(S,F2,F3,F4,Vol4);

                            double SumSingleVolumes = fabs(Vol1) + fabs(Vol2) + fabs(Vol3) + fabs(Vol4);

                            // When computed volumes are equal structure point lies inside
                            if(fabs(SumSingleVolumes-Vol) < epsilon)
                            {
                                // Extract nodal pressures
                                nodalPressures[0] = F1.FastGetSolutionStepValue(PRESSURE);
                                nodalPressures[1] = F2.FastGetSolutionStepValue(PRESSURE);
                                nodalPressures[2] = F3.FastGetSolutionStepValue(PRESSURE);
                                nodalPressures[3] = F4.FastGetSolutionStepValue(PRESSURE);

                                // Do mapping
                                ComputeDiscontinuousInterpolation(S,fluidGeom,elementalDistances,Npos,Nneg,conditionVariable);

                                // Compute face pressure
                                double p_positive_structure = inner_prod(nodalPressures,Npos);
                                double p_negative_structure = inner_prod(nodalPressures,Nneg);

                                double sum_pos = Npos[0]+Npos[1]+Npos[2]+Npos[3];
                                double sum_neg = Nneg[0]+Nneg[1]+Nneg[2]+Nneg[3];

                                // Can be removed after implementing a robust radial shape function ////////////////
                                if((fabs(sum_pos-1) > 0.0001) & (aux_var_pos == 0))
                                {
                                    std::cout << "Radial shape functions in the pressure mapping step are violating the partition of unity condition!" << std::endl;
                                    KRATOS_WATCH(sum_pos);
                                    aux_var_pos = 1;
                                }

                                if((fabs(sum_neg-1) > 0.0001) & (aux_var_neg == 0))
                                {
                                    std::cout << "Radial shape functions in the pressure mapping step are violating the partition of unity condition!" << std::endl;
                                    KRATOS_WATCH(sum_neg);
                                    aux_var_neg = 1;
                                }
                                ///////////////////////////////////////////////////////////////////////////////////

                                // Assign ModelPart::ElementIteratorface pressure to structure node
                                S.FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE,0) = p_positive_structure;
                                S.FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE,0) = p_negative_structure;
                            }
                        }
                    }
                }
            }
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************
    /**
       *
       */
    void AveragePressureToStructure()
    {
        // loop over all fluid elements
        // this loop is parallelized using openmp
#ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
#else
        int number_of_threads = 1;
#endif

        ModelPart::ElementsContainerType& pElements = mrFluidModelPart.Elements();

        vector<unsigned int> element_partition;
        CreatePartition(number_of_threads, pElements.size(), element_partition);
        KRATOS_WATCH(number_of_threads);
        KRATOS_WATCH(element_partition);

#pragma omp parallel for
        for (int k = 0; k < number_of_threads; k++)
        {
            ModelPart::ElementsContainerType::iterator it_begin = pElements.ptr_begin() + element_partition[k];
            ModelPart::ElementsContainerType::iterator it_end = pElements.ptr_begin() + element_partition[k+1];

            for(ModelPart::ElementIterator i_fluidElement = it_begin; i_fluidElement != it_end; ++i_fluidElement)
            {
                const array_1d<double,4> elementalDistances = i_fluidElement->GetValue(ELEMENTAL_DISTANCES);

                Geometry< Node<3> >& fluidGeom = i_fluidElement->GetGeometry();
                array_1d<double,4> nodalPressures;
                std::vector<OctreeType::cell_type*> leaves;

                // Get leaves of octree intersecting with fluid element
                mOctree.GetIntersectedLeaves(*(i_fluidElement).base(),leaves);

                // loop over all octree cells which are intersected by the fluid element
                for(unsigned int i_cell = 0 ; i_cell < leaves.size() ; i_cell++)
                {
                    // Structural element contained in one cell of the octree
                    object_container_type* struct_elem = (leaves[i_cell]->pGetObjects());

                    // loop over all structural elements within each octree cell
                    for(object_container_type::iterator i_StructElement = struct_elem->begin(); i_StructElement != struct_elem->end(); i_StructElement++)
                    {
                        // check if a node of the structure lies within the tetrahedra
                        // --> if yes, then do the pressure mapping
                        Geometry< Node<3> >& structGeom = (*i_StructElement)->GetGeometry();

                        // loop over all nodes of the structural element
                        for(unsigned int i_structElementNode = 0; i_structElementNode < 3 ; i_structElementNode++)
                        {
                            // Comparison of volumes to identify whether a structure points lies inside a tetrahedra or not
                            double Vol = fluidGeom.Volume();
                            double Vol1;
                            double Vol2;
                            double Vol3;
                            double Vol4;

                            Node<3>& S = structGeom[i_structElementNode];
                            Node<3>& F1 = fluidGeom.GetPoint(0);
                            Node<3>& F2 = fluidGeom.GetPoint(1);
                            Node<3>& F3 = fluidGeom.GetPoint(2);
                            Node<3>& F4 = fluidGeom.GetPoint(3);

                            CalculateJacobian(F1,F2,F3,S,Vol1);
                            CalculateJacobian(F1,F2,S,F4,Vol2);
                            CalculateJacobian(F1,S,F3,F4,Vol3);
                            CalculateJacobian(S,F2,F3,F4,Vol4);

                            double SumSingleVolumes = fabs(Vol1) + fabs(Vol2) + fabs(Vol3) + fabs(Vol4);

                            // When computed volumes are equal structure point lies inside
                            if(fabs(SumSingleVolumes-Vol) < epsilon)
                            {
                                // Extract nodal pressures
                                nodalPressures[0] = F1.GetSolutionStepValue(PRESSURE);
                                nodalPressures[1] = F2.GetSolutionStepValue(PRESSURE);
                                nodalPressures[2] = F3.GetSolutionStepValue(PRESSURE);
                                nodalPressures[3] = F4.GetSolutionStepValue(PRESSURE);

                                // Compute average of all positive and all negative values
                                double positiveAverage = 0;
                                int nPos = 0;
                                double negativeAverage = 0;
                                int nNeg = 0;

                                // for output of
                                for(unsigned int i = 0 ; i<4 ; i++)
                                {
                                    if(elementalDistances[i]>=0)
                                    {
                                        positiveAverage += nodalPressures[i];
                                        nPos++;
                                    }
                                    else
                                    {
                                        negativeAverage += nodalPressures[i];
                                        nNeg++;
                                    }
                                }

                                positiveAverage /=nPos;
                                negativeAverage /=nNeg;

                                // Assign Pressures
                                S.SetLock();
                                S.GetSolutionStepValue(POSITIVE_FACE_PRESSURE,0) = positiveAverage;
                                S.GetSolutionStepValue(NEGATIVE_FACE_PRESSURE,0) = negativeAverage;
                                S.UnSetLock();
                            }
                        }
                    }
                }
            }
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void CalculateJacobian( Node<3> node0,
                            Node<3> node1,
                            Node<3> node2,
                            Node<3> node3,
                            double& Volume)
    {
        double x10 = node1.X() - node0.X();
        double y10 = node1.Y() - node0.Y();
        double z10 = node1.Z() - node0.Z();

        double x20 = node2.X() - node0.X();
        double y20 = node2.Y() - node0.Y();
        double z20 = node2.Z() - node0.Z();

        double x30 = node3.X() - node0.X();
        double y30 = node3.Y() - node0.Y();
        double z30 = node3.Z() - node0.Z();

        double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;
        Volume = detJ*0.1666666666666666666667;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void ComputeDiscontinuousInterpolation( const Node<3> pNode,
                                            const Geometry< Node<3> >& geom,
                                            const array_1d<double,4>& distances,
                                            array_1d<double,4>& Npos,
                                            array_1d<double,4>& Nneg,
                                            const double conditionVariable)
    {
        //count positives
        int n_positives = 0;
        for(unsigned int i=0; i<distances.size(); i++)
            if(distances[i]>0) n_positives++;

        //generate the points on the edges at the zero of the distance function
        //generate "father nodes", defined as the end nodes of the edge on which the local point is located
        std::vector< Point<3> > edge_points;
        edge_points.reserve(4);
        array_1d<unsigned int, 4> positive_fathers,  negative_fathers;	//there are at most 4 cut edges
        unsigned int k=0;
        unsigned int l=0;

        for(unsigned int i=0; i<3; i++)
        {
            for(unsigned int j=i+1; j<4; j++) // go through the edges 01, 02, 03, 12, 13, 23
            {
                double di = distances[i];
                double dj = distances[j];

                if(di*dj < 0) //edge is cut
                {
                    //generate point on edge by linear interpolation
                    double Ni = fabs(dj) / ( fabs(di) + fabs(dj) );
                    double Nj = 1.0 - Ni;
                    Point<3> edge_point(Ni * geom[i] + Nj * geom[j]);
                    edge_points.push_back(edge_point);

                    //store the id of the positive and negative fathers
                    if(di > 0.0)
                    {
                        positive_fathers[k++] = i;
                        negative_fathers[l++] = j;
                    }
                    else
                    {
                        positive_fathers[k++] = j;
                        negative_fathers[l++] = i;
                    }
                }
            }
        }

        //detect a h radius as the maximum distance to another node
        array_1d<double,3> xmin = pNode.Coordinates();
        array_1d<double,3> xmax = pNode.Coordinates();
        for(unsigned int i=0; i<edge_points.size(); i++)
        {
            for(unsigned int k=0; k<3; k++)
            {
                if(xmin[k] < edge_points[i][k]) xmin[k] = edge_points[i][k];
                if(xmax[k] > edge_points[i][k]) xmax[k] = edge_points[i][k];
            }
        }
        noalias(xmax) -= xmin; //here i overwrite xmax with the difference
        double h = norm_2(xmax); // h corresponds to parameter value in the gaussian like shape function for the radial interpolation

        //compute local shape functions by sph interpolation
        unsigned int ncut_nodes = edge_points.size();
        Matrix interpolants(ncut_nodes,ncut_nodes,0.0);
        for(unsigned int i=0; i<ncut_nodes; i++)
            for(unsigned int j=0; j<ncut_nodes; j++)
            {
                interpolants(i,j) = RadialShapeFunctionEvaluation(edge_points[j], edge_points[i], h);

                if(i == j)
                {
                    interpolants(i,j) += conditionVariable;
                }
            }

        // pressure for the computation of the weights is assumed to be one at the intersection nodes
        Vector f(ncut_nodes);
        for(unsigned int i = 0; i<ncut_nodes; i++)
            f[i] = 1;

        // Solve equation system
        Vector weights(ncut_nodes);

        // create a permutation matrix for the LU-factorization
        typedef permutation_matrix<std::size_t> pmatrix;
        pmatrix pm(ncut_nodes);
        weights = f;

        // perform LU-factorization to solve LSE
        lu_factorize(interpolants, pm);
        lu_substitute(interpolants, pm, weights);

        //                                KRATOS_WATCH(interpolants)
        //                                KRATOS_WATCH(weights)
        //                                Vector test(ncut_nodes);
        //                                test = prod(interpolants, weights);
        //                                KRATOS_WATCH(test);

        //compute local shape functions (tell how to interpolate from the edge nodes)
        Vector Nlocal(ncut_nodes);

        // Projection of structure point to surface
        for(unsigned int i=0; i<ncut_nodes; i++)
        {
            Nlocal[i] = weights[i]*RadialShapeFunctionEvaluation(pNode, edge_points[i], h);
        }

        noalias(Npos) = ZeroVector(4);
        noalias(Nneg) = ZeroVector(4);
        for(unsigned int i=0; i<ncut_nodes; i++)
        {
            Npos[ positive_fathers[i] ] += Nlocal[i];
            Nneg[ negative_fathers[i] ] += Nlocal[i];
        }

    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    double RadialShapeFunctionEvaluation(const Point<3> edge_point_j,
                                         const Point<3> edge_point_i,
                                         const double h)
    {
        // choose Gaussian radial base function
        array_1d<double,3> radial_vector;

        radial_vector[0] = edge_point_i[0] - edge_point_j[0];
        radial_vector[1] = edge_point_i[1] - edge_point_j[1];
        radial_vector[2] = edge_point_i[2] - edge_point_j[2];

        const double r = norm_2(radial_vector);

        double phi = exp(-(r*r)/(2*h*h));
        return phi;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void DistanceFluidStructure()
    {
        // Initialize elemental distances in the domain
        Initialize();

        // Initialize index table that defines line Edges of fluid element
        bounded_matrix<unsigned int,6,2> TetEdgeIndexTable;
        SetIndexTable(TetEdgeIndexTable);

        // loop over all fluid elements
        // this loop is parallelized using openmp
#ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
#else
        int number_of_threads = 1;
#endif

        ModelPart::ElementsContainerType& pElements = mrFluidModelPart.Elements();

        vector<unsigned int> element_partition;
        CreatePartition(number_of_threads, pElements.size(), element_partition);
        KRATOS_WATCH(number_of_threads);
        KRATOS_WATCH(element_partition);

#pragma omp parallel for
        for (int k = 0; k < number_of_threads; k++)
        {
            ModelPart::ElementsContainerType::iterator it_begin = pElements.ptr_begin() + element_partition[k];
            ModelPart::ElementsContainerType::iterator it_end = pElements.ptr_begin() + element_partition[k+1];

            // assemble all elements
            for (ModelPart::ElementIterator it = it_begin; it != it_end; ++it)
            {
                CalcElementDistances( it , TetEdgeIndexTable );
            }
        }

        // Finally, each tetrahedral element has 4 distance values. But each node belongs to
        // several elements, such that it is assigned several distance values
        // --> now synchronize these values by finding the minimal distance and assign to each node a minimal nodal distance
        AssignMinimalNodalDistance();
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void GenerateSkinModelPart( ModelPart& mrNewSkinModelPart )
    {

        std::vector< array_1d<double,3> >       nodesMatrix;
        std::vector< array_1d<unsigned int,3> > elementsMatrix;

        ReproduceStructuralSkinMesh(nodesMatrix,elementsMatrix);

        ModelPart& new_model_part = mrNewSkinModelPart;

        //new_model_part.GetNodalSolutionStepVariablesList() = mrSkinModelPart.GetNodalSolutionStepVariablesList();

        // ######## ADDING NEW NODES #########

        const unsigned int size_nodes = nodesMatrix.size();
        new_model_part.Nodes().reserve(size_nodes);

        for (unsigned int nodesIndex=0; nodesIndex!=size_nodes; ++nodesIndex)
        {
            int id = nodesIndex + 1;
            double node_X = nodesMatrix[nodesIndex][0];
            double node_Y = nodesMatrix[nodesIndex][1];
            double node_Z = nodesMatrix[nodesIndex][2];
            Node < 3 >::Pointer pnode = Node < 3 > ::Pointer (new Node < 3 >(id,
                                                                             node_X,
                                                                             node_Y,
                                                                             node_Z));
            //pnode->SetBufferSize(mrSkinModelPart.NodesBegin()->GetBufferSize());
            //pnode->X0() = node_X;
            //pnode->Y0() = node_Y;
            //pnode->Z0() = node_Z;
            //pnode->FastGetSolutionStepValue(PARTITION_INDEX)=local_Partition_ov[index];
            new_model_part.Nodes().push_back(pnode);
        }

        // ######## ADDING NEW CONDITIONS #########
        // required information for a new condition: Id, Geometry, Property

        const unsigned int size_conditions = elementsMatrix.size();

        Condition const& rReferenceCondition = KratosComponents<Condition>::Get("Condition3D");
        Properties::Pointer properties = mrSkinModelPart.GetMesh().pGetProperties(1);
        new_model_part.Conditions().reserve(size_conditions);

        for (unsigned int condIndex=0; condIndex!=size_conditions; ++condIndex)
        {
            unsigned int id = condIndex + 1;
            unsigned int node1_ID = elementsMatrix[condIndex][0];
            unsigned int node2_ID = elementsMatrix[condIndex][1];
            unsigned int node3_ID = elementsMatrix[condIndex][2];

            Triangle3D3<Node < 3 > > triangle(new_model_part.Nodes()(node1_ID),
                                              new_model_part.Nodes()(node2_ID),
                                              new_model_part.Nodes()(node3_ID));
            // Index must take care of the IDs of the other processors --> triangle_id_int + 1 + elements_before + total_existing_elements
            Condition::Pointer p_condition = rReferenceCondition.Create(id,
                                                                        triangle,
                                                                        properties);
            new_model_part.Conditions().push_back(p_condition);
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void Initialize()
    {
        const double initial_distance = 1.0;

        ModelPart::NodesContainerType::ContainerType& nodes = mrFluidModelPart.NodesArray();

        // reset the node distance to 1.0 which is the maximum distance in our normalized space.
        unsigned int nodesSize = nodes.size();
        for(unsigned int i = 0 ; i < nodesSize ; i++)
            nodes[i]->GetSolutionStepValue(DISTANCE) = initial_distance;

        ModelPart::ElementsContainerType::ContainerType& fluid_elements = mrFluidModelPart.ElementsArray();

        array_1d<double,4> ElementalDistances;
        ElementalDistances[0] = initial_distance;
        ElementalDistances[1] = initial_distance;
        ElementalDistances[2] = initial_distance;
        ElementalDistances[3] = initial_distance;

        // reset the elemental distance to 1.0 which is the maximum distance in our normalized space.
        // also initialize the embedded velocity of the fluid element
        for(unsigned int i = 0 ; i < fluid_elements.size() ; i++)
        {
            fluid_elements[i]->GetValue(ELEMENTAL_DISTANCES) = ElementalDistances;
            fluid_elements[i]->GetValue(SPLIT_ELEMENT) = false;
            fluid_elements[i]->GetValue(EMBEDDED_VELOCITY)=ZeroVector(3);
        }

    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void SetIndexTable( bounded_matrix<unsigned int,6,2>& TetEdgeIndexTable )
    {
        // Initialize index table to define line Edges of fluid element
        TetEdgeIndexTable(0,0) = 0;
        TetEdgeIndexTable(0,1) = 1;
        TetEdgeIndexTable(1,0) = 0;
        TetEdgeIndexTable(1,1) = 2;
        TetEdgeIndexTable(2,0) = 0;
        TetEdgeIndexTable(2,1) = 3;
        TetEdgeIndexTable(3,0) = 1;
        TetEdgeIndexTable(3,1) = 2;
        TetEdgeIndexTable(4,0) = 1;
        TetEdgeIndexTable(4,1) = 3;
        TetEdgeIndexTable(5,0) = 2;
        TetEdgeIndexTable(5,1) = 3;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void CalcElementDistances( ModelPart::ElementsContainerType::iterator& i_fluidElement,
                               bounded_matrix<unsigned int,6,2>            TetEdgeIndexTable )
    {
        std::vector<OctreeType::cell_type*> leaves;
        std::vector<TetEdgeStruct>          IntersectedTetEdges;
        unsigned int NumberIntersectionsOnTetCorner = 0;

        // Get leaves of octree intersecting with fluid element
        mOctree.GetIntersectedLeaves(*(i_fluidElement).base(),leaves);

        int intersection_counter = 0;

        // Loop over all 6 line Edges of the tetrahedra
        for(unsigned int i_tetEdge = 0;
            i_tetEdge < 6;
            i_tetEdge++)
        {
            IdentifyIntersectionNodes( i_fluidElement, i_tetEdge, leaves, IntersectedTetEdges, NumberIntersectionsOnTetCorner, TetEdgeIndexTable, intersection_counter );
        }

        if (intersection_counter!=0)
        {
            i_fluidElement->GetValue(EMBEDDED_VELOCITY) /= intersection_counter;
        }

        if(IntersectedTetEdges.size() > 0)
            CalcDistanceTo3DSkin( IntersectedTetEdges , i_fluidElement , NumberIntersectionsOnTetCorner );
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void IdentifyIntersectionNodes( ModelPart::ElementsContainerType::iterator&   i_fluidElement,
                                    unsigned int                                  i_tetEdge,
                                    std::vector<OctreeType::cell_type*>&          leaves,
                                    std::vector<TetEdgeStruct>&                   IntersectedTetEdges,
                                    unsigned int&                                 NumberIntersectionsOnTetCorner,
                                    bounded_matrix<unsigned int,6,2>              TetEdgeIndexTable,
                                    int&                                          intersection_counter )
    {
        std::vector<unsigned int> IntersectingStructElemID;
        TetEdgeStruct             NewTetEdge;

        // Get nodes of line Edge
        unsigned int EdgeStartIndex = TetEdgeIndexTable(i_tetEdge,0);
        unsigned int EdgeEndIndex   = TetEdgeIndexTable(i_tetEdge,1);

        PointType& P1 = i_fluidElement->GetGeometry()[EdgeStartIndex];
        PointType& P2 = i_fluidElement->GetGeometry()[EdgeEndIndex];

        double EdgeNode1[3] = {P1.X() , P1.Y() , P1.Z()};
        double EdgeNode2[3] = {P2.X() , P2.Y() , P2.Z()};

        // loop over all octree cells which are intersected by the fluid element
        for(unsigned int i_cell = 0 ; i_cell < leaves.size() ; i_cell++)
        {
            // Structural element contained in one cell of the octree
            object_container_type* struct_elem = (leaves[i_cell]->pGetObjects());

            // loop over all structural elements within each octree cell
            for(object_container_type::iterator i_StructElement = struct_elem->begin(); i_StructElement != struct_elem->end(); i_StructElement++)
            {

                if( StructuralElementNotYetConsidered( (*i_StructElement)->Id() , IntersectingStructElemID ) )
                {

                    // Calculate and associate intersection point to the current fluid element
                    double IntersectionPoint[3] = {0.0 , 0.0 , 0.0};
                    int TetEdgeHasIntersections = IntersectionTriangleSegment( (*i_StructElement)->GetGeometry() , EdgeNode1 , EdgeNode2 , IntersectionPoint );

                    if( TetEdgeHasIntersections == 1 )
                    {
                        IntersectionNodeStruct NewIntersectionNode;

                        // Assign information to the intersection node
                        NewIntersectionNode.Coordinates[0] = IntersectionPoint[0];
                        NewIntersectionNode.Coordinates[1] = IntersectionPoint[1];
                        NewIntersectionNode.Coordinates[2] = IntersectionPoint[2];

                        if ( IsNewIntersectionNode( NewIntersectionNode , IntersectedTetEdges ) )
                        {

                            if( IsIntersectionNodeOnTetEdge( IntersectionPoint , EdgeNode1 , EdgeNode2 ) )
                            {

                                // Calculate normal of the structural element at the position of the intersection point
                                CalculateNormal3D((*i_StructElement)->GetGeometry()[0],
                                                  (*i_StructElement)->GetGeometry()[1],
                                                  (*i_StructElement)->GetGeometry()[2],
                                                  NewIntersectionNode.StructElemNormal);

                                // check, how many intersection nodes are located on corner points of the tetrahedra
                                if ( IsIntersectionOnCorner( NewIntersectionNode , EdgeNode1 , EdgeNode2) )
                                {
                                    NumberIntersectionsOnTetCorner++;
                                }
                                else // CONSIDER ONLY THE EDGES THAT ARE CUT "NOT AT THE VERTEX" -->
                                {
                                    array_1d<double,3> emb_vel = (*i_StructElement)->GetGeometry()[0].GetSolutionStepValue(VELOCITY);
                                    emb_vel += (*i_StructElement)->GetGeometry()[1].GetSolutionStepValue(VELOCITY);
                                    emb_vel += (*i_StructElement)->GetGeometry()[2].GetSolutionStepValue(VELOCITY);

                                    i_fluidElement->GetValue(EMBEDDED_VELOCITY) += emb_vel/3;
                                    intersection_counter++;
                                }

                                if(NumberIntersectionsOnTetCorner < 2)
                                {
                                    // add the new intersection point to the list of intersection points of the fluid element
                                    NewIntersectionNode.EdgeNode1 = EdgeStartIndex;
                                    NewIntersectionNode.EdgeNode2 = EdgeEndIndex;
                                    NewTetEdge.IntNodes.push_back(NewIntersectionNode);

                                    // if tet edge belonging to this intersection point is not already marked as "IntersectedTetEdge" --> put it into the respective container
                                    // when a second intersection node is found, then it is not necessary to push_back again
                                    if( NewTetEdge.IntNodes.size() == 1 )
                                        IntersectedTetEdges.push_back(NewTetEdge);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    bool StructuralElementNotYetConsidered( unsigned int                IDCurrentStructElem,
                                            std::vector<unsigned int>&  IntersectingStructElemID )
    {
        // check if the structural element was already considered as intersecting element
        for(unsigned int k = 0 ; k < IntersectingStructElemID.size() ; k++)
        {
            if( IDCurrentStructElem == IntersectingStructElemID[k] )
                return false;
        }

        // if structural element has not been considered in another octree, which also intersects the fluid element
        // add the new object ID to the vector
        IntersectingStructElemID.push_back( IDCurrentStructElem );
        return true;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    bool IsIntersectionNodeOnTetEdge( double* IntersectionPoint,
                                      double* EdgeNode1,
                                      double* EdgeNode2 )
    {
        // check, if intersection point is located on any edge of the fluid element
        array_1d<double,3> ConnectVectTetNodeIntNode1;
        array_1d<double,3> ConnectVectTetNodeIntNode2;
        array_1d<double,3> EdgeVector;

        ConnectVectTetNodeIntNode1[0] = IntersectionPoint[0] - EdgeNode1[0];
        ConnectVectTetNodeIntNode1[1] = IntersectionPoint[1] - EdgeNode1[1];
        ConnectVectTetNodeIntNode1[2] = IntersectionPoint[2] - EdgeNode1[2];

        ConnectVectTetNodeIntNode2[0] = IntersectionPoint[0] - EdgeNode2[0];
        ConnectVectTetNodeIntNode2[1] = IntersectionPoint[1] - EdgeNode2[1];
        ConnectVectTetNodeIntNode2[2] = IntersectionPoint[2] - EdgeNode2[2];

        double LengthConnectVect1 = norm_2( ConnectVectTetNodeIntNode1 );
        double LengthConnectVect2 = norm_2( ConnectVectTetNodeIntNode2 );

        EdgeVector[0] = EdgeNode2[0] - EdgeNode1[0];
        EdgeVector[1] = EdgeNode2[1] - EdgeNode1[1];
        EdgeVector[2] = EdgeNode2[2] - EdgeNode1[2];

        double MaxEdgeLength = norm_2( EdgeVector );

        // if both connection vectors (corner point --> intersection point)
        // are smaller or equal to the edge length of tetrahedra,
        // then intersection point is located on the edge
        if( (LengthConnectVect1 <= (MaxEdgeLength)) && (LengthConnectVect2 <= (MaxEdgeLength)) )
            return true;
        else
            return false;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    bool IsNewIntersectionNode( IntersectionNodeStruct&    NewIntersectionNode,
                                std::vector<TetEdgeStruct> IntersectedTetEdges )
    {
        array_1d<double,3> DiffVector;
        double NormDiffVector = 0;
        unsigned int NumberIntNodes = 0;

        for( unsigned int i_TetEdge = 0 ; i_TetEdge < IntersectedTetEdges.size() ; i_TetEdge++ )
        {
            NumberIntNodes = IntersectedTetEdges[i_TetEdge].IntNodes.size();
            for( unsigned int i_IntNode = 0 ; i_IntNode < NumberIntNodes ; i_IntNode++ )
            {
                DiffVector[0] = NewIntersectionNode.Coordinates[0] - IntersectedTetEdges[i_TetEdge].IntNodes[i_IntNode].Coordinates[0];
                DiffVector[1] = NewIntersectionNode.Coordinates[1] - IntersectedTetEdges[i_TetEdge].IntNodes[i_IntNode].Coordinates[1];
                DiffVector[2] = NewIntersectionNode.Coordinates[2] - IntersectedTetEdges[i_TetEdge].IntNodes[i_IntNode].Coordinates[2];

                NormDiffVector = norm_2(DiffVector);

                if( NormDiffVector < epsilon )
                    return false;
            }
        }

        // if the new intersection node is not existing (as intersection with a corner point), then return false
        return true;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    bool IsIntersectionOnCorner( IntersectionNodeStruct& NewIntersectionNode,
                                 double*                 EdgeNode1,
                                 double*                 EdgeNode2 )
    {
        array_1d<double,3> DiffVector;
        double NormDiffVector;

        DiffVector[0] = EdgeNode1[0] - NewIntersectionNode.Coordinates[0];
        DiffVector[1] = EdgeNode1[1] - NewIntersectionNode.Coordinates[1];
        DiffVector[2] = EdgeNode1[2] - NewIntersectionNode.Coordinates[2];
        NormDiffVector = norm_2(DiffVector);

        if( NormDiffVector < epsilon )
            return true;

        DiffVector[0] = EdgeNode2[0] - NewIntersectionNode.Coordinates[0];
        DiffVector[1] = EdgeNode2[1] - NewIntersectionNode.Coordinates[1];
        DiffVector[2] = EdgeNode2[2] - NewIntersectionNode.Coordinates[2];
        NormDiffVector = norm_2(DiffVector);

        if( NormDiffVector < epsilon )
            return true;
        else
            return false;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void CalculateNormal3D( Point<3>&       Point1,
                            Point<3>&       Point2,
                            Point<3>&       Point3,
                            array_1d<double,3>&   rResultNormal )
    {
        array_1d<double,3> v1 = Point2 - Point1;
        array_1d<double,3> v2 = Point3 - Point1;

        MathUtils<double>::CrossProduct(rResultNormal,v1,v2);
        rResultNormal *= 0.5;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void CalcDistanceTo3DSkin( std::vector<TetEdgeStruct>&                 IntersectedTetEdges,
                               ModelPart::ElementsContainerType::iterator& i_fluid_element,
                               unsigned int                                NumberIntersectionsOnTetCorner )
    {
        std::vector<IntersectionNodeStruct> NodesOfApproximatedStructure;
        array_1d<double,4> ElementalDistances;

        // Reduce all found intersection nodes located on each tetdrahedra edge to just one intersection node by averaging
        ComputeApproximationNodes(IntersectedTetEdges,NodesOfApproximatedStructure);

        // Intersection with one corner point
        if( NodesOfApproximatedStructure.size() == 1 && NumberIntersectionsOnTetCorner == 1 )
        {
            CalcSignedDistancesToOneIntNode(i_fluid_element,NodesOfApproximatedStructure,ElementalDistances);
            i_fluid_element->GetValue(SPLIT_ELEMENT) = true;
        }

        // Intersection with two corner points / one tetrahedra edge
        if( NodesOfApproximatedStructure.size() == 2 && NumberIntersectionsOnTetCorner == 2 )
        {
            CalcSignedDistancesToTwoIntNodes(i_fluid_element,NodesOfApproximatedStructure,ElementalDistances);
            i_fluid_element->GetValue(SPLIT_ELEMENT) = true;
        }

        // Intersection with three tetrahedra edges
        if( NodesOfApproximatedStructure.size() == 3 )
        {
            CalcSignedDistancesToThreeIntNodes(i_fluid_element,NodesOfApproximatedStructure,IntersectedTetEdges,ElementalDistances);
            i_fluid_element->GetValue(SPLIT_ELEMENT) = true;
        }

        // Intersection with more than three tetrahedra edges
        if( NodesOfApproximatedStructure.size() > 3 )
        {
            CalcSignedDistancesToMoreThanThreeIntNodes(i_fluid_element,NodesOfApproximatedStructure,IntersectedTetEdges,ElementalDistances);
            i_fluid_element->GetValue(SPLIT_ELEMENT) = true;
        }

        // Postprocessing treatment of elemental distances
        if( i_fluid_element->GetValue(SPLIT_ELEMENT) == true )
            AvoidZeroDistances(i_fluid_element, ElementalDistances);

        // In case there is intersection with fluid element: assign distances to the element
        if( i_fluid_element->GetValue(SPLIT_ELEMENT) == true )
            i_fluid_element->GetValue(ELEMENTAL_DISTANCES) = ElementalDistances;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void ComputeApproximationNodes( std::vector<TetEdgeStruct>           IntersectedTetEdges,
                                    std::vector<IntersectionNodeStruct>& NodesOfApproximatedStructure )
    {
        unsigned int NumberIntNodes = 0;
        double sum_X;
        double sum_Y;
        double sum_Z;

        unsigned int NumberCutEdges = IntersectedTetEdges.size();
        std::vector<unsigned int> IndicesOfDoubleCutEdges;
        IndicesOfDoubleCutEdges.reserve(NumberCutEdges);

        // calculate average of all intersection nodes of each tetrahedra edge
        for(unsigned int i_TetEdge = 0 ; i_TetEdge < NumberCutEdges ; i_TetEdge++)
        {
            NumberIntNodes = IntersectedTetEdges[i_TetEdge].IntNodes.size();
            if(NumberIntNodes == 2)
            {
                IndicesOfDoubleCutEdges.push_back(i_TetEdge);
            }
        }

        const unsigned int NumberEdgesDoubleCut = IndicesOfDoubleCutEdges.size();
        if(NumberEdgesDoubleCut == 1)
        {
            // kick these intersection nodes on this edge out of the container
            unsigned int index = IndicesOfDoubleCutEdges[0];
            IntersectedTetEdges.erase(IntersectedTetEdges.begin() + index);
            NumberCutEdges -= 1;
        }

        for(unsigned int i_TetEdge = 0 ; i_TetEdge < NumberCutEdges ; i_TetEdge++)
        {
            sum_X = 0;
            sum_Y = 0;
            sum_Z = 0;

            for( unsigned int i_IntNode = 0 ; i_IntNode < NumberIntNodes ; i_IntNode++ )
            {
                sum_X += IntersectedTetEdges[i_TetEdge].IntNodes[i_IntNode].Coordinates[0];
                sum_Y += IntersectedTetEdges[i_TetEdge].IntNodes[i_IntNode].Coordinates[1];
                sum_Z += IntersectedTetEdges[i_TetEdge].IntNodes[i_IntNode].Coordinates[2];
            }

            IntersectionNodeStruct NewApproximationNode;
            NewApproximationNode.Coordinates[0] = sum_X / NumberIntNodes;
            NewApproximationNode.Coordinates[1] = sum_Y / NumberIntNodes;
            NewApproximationNode.Coordinates[2] = sum_Z / NumberIntNodes;

            if(IntersectedTetEdges.size() <= 2)
                NewApproximationNode.StructElemNormal = IntersectedTetEdges[i_TetEdge].IntNodes[0].StructElemNormal;

            NewApproximationNode.EdgeNode1 = IntersectedTetEdges[i_TetEdge].IntNodes[0].EdgeNode1;
            NewApproximationNode.EdgeNode2 = IntersectedTetEdges[i_TetEdge].IntNodes[0].EdgeNode2;

            NodesOfApproximatedStructure.push_back(NewApproximationNode);
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void CalcSignedDistancesToOneIntNode( ModelPart::ElementsContainerType::iterator& i_fluid_element,
                                          std::vector<IntersectionNodeStruct>         NodesOfApproximatedStructure,
                                          array_1d<double,4>&                         ElementalDistances )
    {
        const array_1d<double,3>& IntersectionNodeCoord = NodesOfApproximatedStructure[0].Coordinates;
        array_1d<double,3> DistVecTetNode;
        array_1d<double,3> TetNode;
        array_1d<double,3> NormalAtIntersectionNode;
        double             NormDistTetNode;
        double             InnerProduct;

        Geometry< Node<3> >& rFluidGeom = i_fluid_element->GetGeometry();

        for(unsigned int i_TetNode = 0 ; i_TetNode < 4 ; i_TetNode++)
        {
            // Get coordinates of the fluid elmenent nodes
            TetNode  = rFluidGeom[i_TetNode].Coordinates();

            // Compute unsigned distance
            DistVecTetNode[0] = TetNode[0] - IntersectionNodeCoord[0];
            DistVecTetNode[1] = TetNode[1] - IntersectionNodeCoord[1];
            DistVecTetNode[2] = TetNode[2] - IntersectionNodeCoord[2];
            NormDistTetNode = norm_2( DistVecTetNode );

            // Get normal at intersection
            NormalAtIntersectionNode = NodesOfApproximatedStructure[0].StructElemNormal;
            InnerProduct = inner_prod(DistVecTetNode,NormalAtIntersectionNode);

            // Assign distances as nodal solution values
            if(InnerProduct>epsilon)
                ElementalDistances[i_TetNode] = NormDistTetNode;
            else if(InnerProduct>-epsilon)
                ElementalDistances[i_TetNode] = 0;
            else
                ElementalDistances[i_TetNode] = -NormDistTetNode;
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void CalcSignedDistancesToTwoIntNodes( ModelPart::ElementsContainerType::iterator& i_fluid_element,
                                           std::vector<IntersectionNodeStruct>         NodesOfApproximatedStructure,
                                           array_1d<double,4>&                         ElementalDistances )
    {
        const array_1d<double,3>& IntersectionNode1Coord = NodesOfApproximatedStructure[0].Coordinates;
        const array_1d<double,3>& IntersectionNode2Coord = NodesOfApproximatedStructure[1].Coordinates;
        array_1d<double,3> TetNode;
        array_1d<double,3> DistVecTetNode;
        array_1d<double,3> NormalAtIntersectionNode1;
        array_1d<double,3> NormalAtIntersectionNode2;
        array_1d<double,3> ResNormal;
        double             InnerProduct;
        double             NormDistTetNode;

        const Point<3> LinePoint1 = Point<3>(IntersectionNode1Coord[0] , IntersectionNode1Coord[1] , IntersectionNode1Coord[2]);
        const Point<3> LinePoint2 = Point<3>(IntersectionNode2Coord[0] , IntersectionNode2Coord[1] , IntersectionNode2Coord[2]);

        Geometry< Node<3> >& rFluidGeom = i_fluid_element->GetGeometry();

        for(unsigned int i_TetNode = 0 ; i_TetNode < 4 ; i_TetNode++)
        {
            // Get coordinates of the fluid element nodes
            TetNode  = rFluidGeom(i_TetNode)->Coordinates();

            // Compute distance to point
            NormDistTetNode = GeometryUtils::PointDistanceToLineSegment3D(LinePoint1, LinePoint2 , Point<3>(TetNode[0],TetNode[1],TetNode[2]));

            // Compute unsigned distance vector by assuming the mean position vector of the two intersection points
            DistVecTetNode[0] = TetNode[0] - IntersectionNode1Coord[0];
            DistVecTetNode[1] = TetNode[1] - IntersectionNode1Coord[1];
            DistVecTetNode[2] = TetNode[2] - IntersectionNode1Coord[2];

            // Get normal at intersections, average them and check direction of distances
            NormalAtIntersectionNode1 = NodesOfApproximatedStructure[0].StructElemNormal;
            NormalAtIntersectionNode2 = NodesOfApproximatedStructure[1].StructElemNormal;

            // Compute unsigned distance
            ResNormal[0] = 0.5*(NormalAtIntersectionNode1[0] + NormalAtIntersectionNode2[0]);
            ResNormal[1] = 0.5*(NormalAtIntersectionNode1[1] + NormalAtIntersectionNode2[1]);
            ResNormal[2] = 0.5*(NormalAtIntersectionNode1[2] + NormalAtIntersectionNode2[2]);
            InnerProduct = inner_prod(DistVecTetNode,ResNormal);

            // Assign distances as nodal solution values
            if(InnerProduct>epsilon)
                ElementalDistances[i_TetNode] = NormDistTetNode;
            else if(InnerProduct>-epsilon)
                ElementalDistances[i_TetNode] = 0;
            else
                ElementalDistances[i_TetNode] = -NormDistTetNode;
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void CalcSignedDistancesToThreeIntNodes( ModelPart::ElementsContainerType::iterator& i_fluid_element,
                                             std::vector<IntersectionNodeStruct>         NodesOfApproximatedStructure,
                                             std::vector<TetEdgeStruct>                  IntersectedTetEdges,
                                             array_1d<double,4>&                         ElementalDistances )
    {
        Geometry< Node<3> >& rFluidGeom = i_fluid_element->GetGeometry();

        Point<3> P1;
        Point<3> P2;
        Point<3> P3;

        P1.Coordinates() = NodesOfApproximatedStructure[0].Coordinates;
        P2.Coordinates() = NodesOfApproximatedStructure[1].Coordinates;
        P3.Coordinates() = NodesOfApproximatedStructure[2].Coordinates;
        array_1d<double,3> Normal;
        CalculateNormal3D(P1,P2,P3,Normal);

        for(unsigned int i_TetNode = 0 ; i_TetNode < 4 ; i_TetNode++)
        {
            array_1d<double,3> TetNode;
            array_1d<double,3> IntersectionNode1Coord;
            array_1d<double,3> IntersectionNode2Coord;
            array_1d<double,3> IntersectionNode3Coord;
            Point<3>           ApproxTrianglePoint1;
            Point<3>           ApproxTrianglePoint2;
            Point<3>           ApproxTrianglePoint3;
            double             UnsignedDistance;
            double             InnerProduct;

            // Get coordinates of the fluid element nodes
            TetNode = rFluidGeom(i_TetNode)->Coordinates();

            IntersectionNode1Coord = NodesOfApproximatedStructure[0].Coordinates;
            IntersectionNode2Coord = NodesOfApproximatedStructure[1].Coordinates;
            IntersectionNode3Coord = NodesOfApproximatedStructure[2].Coordinates;

            ApproxTrianglePoint1 = Point<3>(IntersectionNode1Coord[0] , IntersectionNode1Coord[1] , IntersectionNode1Coord[2]);
            ApproxTrianglePoint2 = Point<3>(IntersectionNode2Coord[0] , IntersectionNode2Coord[1] , IntersectionNode2Coord[2]);
            ApproxTrianglePoint3 = Point<3>(IntersectionNode3Coord[0] , IntersectionNode3Coord[1] , IntersectionNode3Coord[2]);

            // Compute distance from tet node to current triangle
            UnsignedDistance = PointDistanceToTriangle3DPlane(ApproxTrianglePoint1, ApproxTrianglePoint2 , ApproxTrianglePoint3 , rFluidGeom[i_TetNode]);

            bool TetNodeIsInsideStructure = true;
            bool TetNodeIsOnStructure = true;
            array_1d <double,3> DistVec;

            DistVec[0] = TetNode[0] - IntersectionNode1Coord[0];
            DistVec[1] = TetNode[1] - IntersectionNode1Coord[1];
            DistVec[2] = TetNode[2] - IntersectionNode1Coord[2];

            InnerProduct = inner_prod(DistVec,Normal);

            if(InnerProduct > epsilon)
            {
                TetNodeIsInsideStructure = false;
                TetNodeIsOnStructure = false;
            }
            else if (InnerProduct < -epsilon)
                TetNodeIsOnStructure = false;

            // Assign distances as nodal solution values ( + = outside of structure, - = inside structure)
            if( TetNodeIsInsideStructure == true )
                ElementalDistances[i_TetNode] = -UnsignedDistance;
            else if( TetNodeIsOnStructure == true )
                ElementalDistances[i_TetNode] = 0;
            else
                ElementalDistances[i_TetNode] = +UnsignedDistance;
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void CalcSignedDistancesToMoreThanThreeIntNodes(  ModelPart::ElementsContainerType::iterator& i_fluid_element,
                                                      std::vector<IntersectionNodeStruct>         NodesOfApproximatedStructure,
                                                      std::vector<TetEdgeStruct>                  IntersectedTetEdges,
                                                      array_1d<double,4>&                         ElementalDistances )
    {
        unsigned int numberCutEdges = NodesOfApproximatedStructure.size();

        // Compute average of the intersection nodes which is a node on the plane we look for
        Point<3> r_mean;
        for(unsigned int k=0; k<numberCutEdges; k++)
            for(unsigned int i=0; i<3; i++)
                r_mean.Coordinates()[i] += NodesOfApproximatedStructure[k].Coordinates[i];

        for(unsigned int i=0; i<3; i++)
            r_mean.Coordinates()[i] /= numberCutEdges;

        // Compute average of normals and areas
        array_1d<double,3> N_mean;

        Matrix coordinates(numberCutEdges,3);
        for(unsigned int i=0; i<numberCutEdges; i++)
            for(unsigned int j=0; j<3; j++)
                coordinates(i,j) = NodesOfApproximatedStructure[i].Coordinates[j] - r_mean[j];

        Matrix A = prod(trans(coordinates),coordinates);
        Matrix V(3,3);
        Vector lambda(3);

        // Calculate the eigenvectors V and the corresponding eigenvalues lambda
        EigenVectors(A, V, lambda);

        // Look for the minimal eigenvalue all lambdas
        unsigned int min_pos = 0;
        double min_lambda = lambda[min_pos];
        for(unsigned int i=1;i<3; i++)
            if(min_lambda > lambda[i])
            {
                min_lambda = lambda[i];
                min_pos = i;
            }

        // the normal equals to the eigenvector which corresponds to the minimal eigenvalue
        for(unsigned int i=0;i<3; i++) N_mean[i] = V(min_pos,i);
        N_mean /= norm_2(N_mean);

        // Determine about the minimal distance by considering the distances to both triangles
        double UnsignedDistance;
        array_1d<double,3> TetNode;

        for(unsigned int i_TetNode = 0 ; i_TetNode < 4 ; i_TetNode++)
        {
            TetNode = i_fluid_element->GetGeometry()[i_TetNode].Coordinates();

            // Calculate Distance
            UnsignedDistance = PointDistanceTo3DPlane(r_mean, N_mean, i_fluid_element->GetGeometry()[i_TetNode] );

            // Assign sings according to closest normals
            bool TetNodeIsInsideStructure = true;
            bool TetNodeIsOnStructure = true;
            array_1d <double,3> DistVec;

            DistVec[0] = TetNode[0] - r_mean[0];
            DistVec[1] = TetNode[1] - r_mean[1];
            DistVec[2] = TetNode[2] - r_mean[2];

            double InnerProduct = inner_prod(DistVec,N_mean);
            if(InnerProduct > epsilon)
            {
                TetNodeIsInsideStructure = false;
                TetNodeIsOnStructure = false;
            }
            else if (InnerProduct < -epsilon)
                TetNodeIsOnStructure = false;

            // Assign distances as nodal solution values ( + = outside of structure, - = inside structure)
            if( TetNodeIsInsideStructure == true )
                ElementalDistances[i_TetNode] = -UnsignedDistance;
            else if( TetNodeIsOnStructure == true )
                ElementalDistances[i_TetNode] = 0;
            else
                ElementalDistances[i_TetNode] = +UnsignedDistance;
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    /**
       * This function calculates the distance of a 3D point to a plane spanned by a 3D triangle
       * Copyright 2001 softSurfer, 2012 Dan Sunday; Source: http://geomalgorithms.com/a04-_planes.html
       * @param TrianglePoint1 First point of triangle
       * @param TrianglePoint2 Second point of triangle
       * @param TrianglePoint3 Third point of triangle
       * @param ToPoint The point which distance is required
       * @return The distance between the point and the plane spanned by the 3D triangle
       */
    double PointDistanceToTriangle3DPlane( Point<3>& TrianglePoint1,
                                           Point<3>& TrianglePoint2,
                                           Point<3>& TrianglePoint3,
                                           Point<3>& ToPoint)
    {
        // Compute normal of triangle
        array_1d<double,3> normalTriangle;
        CalculateNormal3D(TrianglePoint1,
                          TrianglePoint2,
                          TrianglePoint3,
                          normalTriangle);

        // calculate vector pointing from a node in the plane (e.g. triangle point 1) to the considered node ToPoint
        array_1d<double,3> planeToPointVec = ToPoint - TrianglePoint1;

        // projection of node on the plane
        const double sn = -inner_prod(normalTriangle,planeToPointVec);
        const double sd = inner_prod(normalTriangle,normalTriangle);
        const double sb = sn / sd;

        // project node onto the plane
        //const array_1d<double,3> PointProjectedOnPlane = ToPoint + sb * normalTriangle;
        const array_1d<double,3> ProjectedPointToPoint = sb * normalTriangle;
        const double DistanceToPlane = norm_2(ProjectedPointToPoint);

        return DistanceToPlane;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    /**
       * This function calculates the distance of a 3D point to a plane spanned by a 3D triangle
       * @param Plane base point
       * @param planeNormal
       * @param ToPoint The point which distance is required
       * @return The distance between the point and the plane spanned by the 3D triangle
       */
    double PointDistanceTo3DPlane( Point<3>&            planeBasePoint,
                                   array_1d<double, 3>& planeNormal,
                                   Point<3>&            ToPoint)
    {
        // calculate vector pointing from a node in the plane (e.g. triangle point 1) to the considered node ToPoint
        array_1d<double,3> planeToPointVec = ToPoint - planeBasePoint;

        // projection of node on the plane
        const double sn = -inner_prod(planeNormal,planeToPointVec);
        const double sd = inner_prod(planeNormal,planeNormal);
        const double sb = sn / sd;

        // project node onto the plane
        const array_1d<double,3> ProjectedPointToPoint = sb * planeNormal;
        const double DistanceToPlane = norm_2(ProjectedPointToPoint);

        return DistanceToPlane;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    /**
       * This function calculates the minimal distances of a node considering all neighboring elements
       */
    void AssignMinimalNodalDistance()
    {
        // loop over all fluid elements
        for( ModelPart::ElementIterator i_fluid_element = mrFluidModelPart.ElementsBegin();
             i_fluid_element != mrFluidModelPart.ElementsEnd();
             i_fluid_element++)
        {
            Geometry< Node<3> >& geom = i_fluid_element->GetGeometry();
            array_1d<double,4> ElementalDistances = i_fluid_element->GetValue(ELEMENTAL_DISTANCES);

            // Assign distances to the single nodes, if a smaller value is found
            for( unsigned int i_TetNode = 0; i_TetNode < 4; i_TetNode++ )
            {
                double currentNodeDist = geom[i_TetNode].GetSolutionStepValue(DISTANCE);
                double nodeInElemDist  = ElementalDistances[i_TetNode];

                if( fabs( nodeInElemDist ) < fabs( currentNodeDist ) )
                    geom[i_TetNode].GetSolutionStepValue(DISTANCE) = nodeInElemDist; // overwrite nodal distance (which is global)
            } // loop i_TetNode
        } // loop i_fluidElement
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    /**
       * If structure directly passes through the corner point of a tetrahedra (leading to zero distances
       * in the respective node), then a small distance value (different from zero) will be stored for
       * that point. This is necessary since the embedded solver cannot handle zero distances.
       * @param element            current element which was cut by the structure (flag SPLIT_ELEMENT is set to one)
       * @param ElementalDistances elemental distances calculated by the intersection pattern
       */
    void AvoidZeroDistances( ModelPart::ElementsContainerType::iterator& element,
                             array_1d<double,4>&                         ElementalDistances)
    {
        // Assign a distance limit
        double dist_limit = 1e-6;
        bool distChangedToLimit = false; //variable to indicate that a distance value < tolerance is set to a limit distance = tolerance

        for(unsigned int i_node = 0; i_node < 4; i_node++)
        {
            if(fabs(ElementalDistances[i_node]) < dist_limit)
            {
                ElementalDistances[i_node] = dist_limit;
                distChangedToLimit = true;
            }
        }

        // Check, if this approach changes the split-flag (might be, that element is not cut anymore if node with zero distance gets a positive limit distance value
        unsigned int numberNodesPositiveDistance = 0;
        for(unsigned int i_node = 0; i_node < 4; i_node++)
        {
            if((ElementalDistances[i_node]) > 0)
                numberNodesPositiveDistance++;
        }

        // element is not set
        if(numberNodesPositiveDistance == 4 && distChangedToLimit == true)
            element->GetValue(SPLIT_ELEMENT) = false;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void ReproduceStructuralSkinMesh( std::vector< array_1d<double,3> >&       nodesMatrix,
                                      std::vector< array_1d<unsigned int,3> >& elementsMatrix )
    {
        // Extract distance information of intersected elements and generate triangular geometry
        // in order to visualize in GiD as a mdpa-file

        // loop over all fluid elements
        for(ModelPart::ElementIterator i_fluid_element = mrFluidModelPart.ElementsBegin();
            i_fluid_element != mrFluidModelPart.ElementsEnd();
            i_fluid_element++)
        {
            bool split_element = i_fluid_element->GetValue(SPLIT_ELEMENT);
            if (split_element == true)
            {
                array_1d<double,4> distances = i_fluid_element->GetValue(ELEMENTAL_DISTANCES);

                // Initialize index table to define line Edges of fluid element
                bounded_matrix<unsigned int,6,2> TetEdgeIndexTable;
                SetIndexTable(TetEdgeIndexTable);

                unsigned int numberCutEdges = 0;

                // First, compute all intersection nodes along the tetrahedra edges
                for(unsigned int i_tetEdge = 0; i_tetEdge < 6; i_tetEdge++)
                {
                    // Get nodes of tet edge
                    unsigned int EdgeStartIndex = TetEdgeIndexTable(i_tetEdge,0);
                    unsigned int EdgeEndIndex   = TetEdgeIndexTable(i_tetEdge,1);

                    // Edge is cut, if distance values on both nodes have opposite sign
                    if(distances[EdgeStartIndex]*distances[EdgeEndIndex] < 0)
                        numberCutEdges++;
                }

                /*
                  if numberCutEdges == 0: tet is intersected, BUT all nodes are nevertheless located outside of structure
                                          (think of very thin structures)
                  if numberCutEdges == 1: 1 node  = 0 (intersected tet corner node) --> NEGLECTED!!
                  if numberCutEdges == 2: 2 nodes = 0 (intersected tet edge)        --> NEGLECTED!!
                  if numberCutEdges == 3: 3 nodes with same sign, 1 node with opposite sign  --> 3 intersection nodes
                  if numberCutEdges == 4: 2 nodes with same sign, 2 nodes with opposite sign --> 4 intersection nodes
                  */

                if(numberCutEdges == 3)
                {
                    ReproduceThreeIntNodes(i_fluid_element,distances,nodesMatrix,elementsMatrix);
                }
                else if(numberCutEdges == 4)
                {
                    ReproduceFourIntNodes(i_fluid_element,distances,nodesMatrix,elementsMatrix);
                }
            }
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void ReproduceThreeIntNodes( ModelPart::ElementsContainerType::iterator& i_fluid_element,
                                 array_1d<double,4>                          distances,
                                 std::vector< array_1d<double,3> >&          nodesMatrix,
                                 std::vector< array_1d<unsigned int,3> >&    elementsMatrix)
    {
        // case: 3 intersection nodes --> write 1 triangle directly
        unsigned int nextNodeID = nodesMatrix.size()+1;
        unsigned int ID_Node_A = nextNodeID;
        unsigned int ID_Node_B = nextNodeID + 1;
        unsigned int ID_Node_C = nextNodeID + 2;

        array_1d<double,3> intNode_A;
        array_1d<double,3> intNode_B;
        array_1d<double,3> intNode_C;

        array_1d<unsigned int,3> triangle;

        // ********* ONE NODE NEGATIVE; REMAINING NODES POSITIVE **********
        if( distances[0] < 0 && distances[1] > 0 && distances[2] > 0 && distances[3] > 0 )
        { // node 0 negative, remaining nodes positive
            CalcIntersectionNode(i_fluid_element,distances,0,1,intNode_A); // Intersection Node A
            CalcIntersectionNode(i_fluid_element,distances,0,2,intNode_B); // Intersection Node B
            CalcIntersectionNode(i_fluid_element,distances,0,3,intNode_C); // Intersection Node C
        }
        else if( distances[0] > 0 && distances[1] < 0 && distances[2] > 0 && distances[3] > 0 )
        { // node 1 negative, remaining nodes positive
            CalcIntersectionNode(i_fluid_element,distances,1,0,intNode_A); // Intersection Node A
            CalcIntersectionNode(i_fluid_element,distances,1,2,intNode_B); // Intersection Node B
            CalcIntersectionNode(i_fluid_element,distances,1,3,intNode_C); // Intersection Node C
        }
        else if( distances[0] > 0 && distances[1] > 0 && distances[2] < 0 && distances[3] > 0 )
        { // node 2 negative, remaining nodes positive
            CalcIntersectionNode(i_fluid_element,distances,2,0,intNode_A); // Intersection Node A
            CalcIntersectionNode(i_fluid_element,distances,2,1,intNode_B); // Intersection Node B
            CalcIntersectionNode(i_fluid_element,distances,2,3,intNode_C); // Intersection Node C
        }
        else if( distances[0] > 0 && distances[1] > 0 && distances[2] > 0 && distances[3] < 0 )
        { // node 3 negative, remaining nodes positive
            CalcIntersectionNode(i_fluid_element,distances,3,0,intNode_A); // Intersection Node A
            CalcIntersectionNode(i_fluid_element,distances,3,1,intNode_B); // Intersection Node B
            CalcIntersectionNode(i_fluid_element,distances,3,2,intNode_C); // Intersection Node C
        }
        // ********* ONE NODE POSITIVE; REMAINING NODES NEGATIVE **********
        else if( distances[0] > 0 && distances[1] < 0 && distances[2] < 0 && distances[3] < 0 )
        { // node 0 positive, remaining nodes negative
            CalcIntersectionNode(i_fluid_element,distances,1,0,intNode_A); // Intersection Node A
            CalcIntersectionNode(i_fluid_element,distances,2,0,intNode_B); // Intersection Node B
            CalcIntersectionNode(i_fluid_element,distances,3,0,intNode_C); // Intersection Node C
        }
        else if( distances[0] < 0 && distances[1] > 0 && distances[2] < 0 && distances[3] < 0 )
        { // node 1 positive, remaining nodes negative
            CalcIntersectionNode(i_fluid_element,distances,0,1,intNode_A); // Intersection Node A
            CalcIntersectionNode(i_fluid_element,distances,2,1,intNode_B); // Intersection Node B
            CalcIntersectionNode(i_fluid_element,distances,3,1,intNode_C); // Intersection Node C
        }
        else if( distances[0] < 0 && distances[1] < 0 && distances[2] > 0 && distances[3] < 0 )
        { // node 2 positive, remaining nodes negative
            CalcIntersectionNode(i_fluid_element,distances,0,2,intNode_A); // Intersection Node A
            CalcIntersectionNode(i_fluid_element,distances,1,2,intNode_B); // Intersection Node B
            CalcIntersectionNode(i_fluid_element,distances,3,2,intNode_C); // Intersection Node C
        }
        else if( distances[0] < 0 && distances[1] < 0 && distances[2] < 0 && distances[3] > 0 )
        { // node 3 positive, remaining nodes negative
            CalcIntersectionNode(i_fluid_element,distances,0,3,intNode_A); // Intersection Node A
            CalcIntersectionNode(i_fluid_element,distances,1,3,intNode_B); // Intersection Node B
            CalcIntersectionNode(i_fluid_element,distances,2,3,intNode_C); // Intersection Node C
        }

        nodesMatrix.push_back(intNode_A);
        nodesMatrix.push_back(intNode_B);
        nodesMatrix.push_back(intNode_C);

        // triangle: ABC
        triangle[0] = ID_Node_A;
        triangle[1] = ID_Node_B;
        triangle[2] = ID_Node_C;

        AdaptIndexingByNormal(i_fluid_element,distances,triangle,nodesMatrix);

        elementsMatrix.push_back(triangle);
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void ReproduceFourIntNodes(  ModelPart::ElementsContainerType::iterator& i_fluid_element,
                                 array_1d<double,4>                          distances,
                                 std::vector< array_1d<double,3> >&          nodesMatrix,
                                 std::vector< array_1d<unsigned int,3> >&    elementsMatrix)
    {
        // case: 4 intersection nodes --> write 2 triangles
        // or more abstractly spoken: search the nodes which form the 2 triangles

        unsigned int nextNodeID = nodesMatrix.size()+1;
        unsigned int ID_Node_A = nextNodeID;
        unsigned int ID_Node_B = nextNodeID + 1;
        unsigned int ID_Node_C = nextNodeID + 2;
        unsigned int ID_Node_D = nextNodeID + 3;

        array_1d<double,3> intNode_A;
        array_1d<double,3> intNode_B;
        array_1d<double,3> intNode_C;
        array_1d<double,3> intNode_D;

        array_1d<unsigned int,3> triangle1;
        array_1d<unsigned int,3> triangle2;

        if( distances[0] < 0 )
        {
            if( distances[1] < 0 ) // nodes [0,1] negative
            {
                CalcIntersectionNode(i_fluid_element,distances,0,2,intNode_A); // Intersection Node A
                CalcIntersectionNode(i_fluid_element,distances,0,3,intNode_B); // Intersection Node B
                CalcIntersectionNode(i_fluid_element,distances,1,2,intNode_C); // Intersection Node C
                CalcIntersectionNode(i_fluid_element,distances,1,3,intNode_D); // Intersection Node D
            }
            else if( distances[2] < 0 ) // nodes [0,2] negative
            {
                CalcIntersectionNode(i_fluid_element,distances,0,1,intNode_A); // Intersection Node A
                CalcIntersectionNode(i_fluid_element,distances,0,3,intNode_B); // Intersection Node B
                CalcIntersectionNode(i_fluid_element,distances,2,1,intNode_C); // Intersection Node C
                CalcIntersectionNode(i_fluid_element,distances,2,3,intNode_D); // Intersection Node D
            }
            else if( distances[3] < 0 ) // nodes [0,3] negative
            {
                CalcIntersectionNode(i_fluid_element,distances,0,1,intNode_A); // Intersection Node A
                CalcIntersectionNode(i_fluid_element,distances,0,2,intNode_B); // Intersection Node B
                CalcIntersectionNode(i_fluid_element,distances,3,1,intNode_C); // Intersection Node C
                CalcIntersectionNode(i_fluid_element,distances,3,2,intNode_D); // Intersection Node D
            }
        }
        else if( distances[0] > 0 )
        {
            if( distances[1] > 0 ) // nodes [2,3] negative
            {
                CalcIntersectionNode(i_fluid_element,distances,2,0,intNode_A); // Intersection Node A
                CalcIntersectionNode(i_fluid_element,distances,2,1,intNode_B); // Intersection Node B
                CalcIntersectionNode(i_fluid_element,distances,3,0,intNode_C); // Intersection Node C
                CalcIntersectionNode(i_fluid_element,distances,3,1,intNode_D); // Intersection Node D
            }
            else if( distances[2] > 0 ) // nodes [1,3] negative
            {
                CalcIntersectionNode(i_fluid_element,distances,1,0,intNode_A); // Intersection Node A
                CalcIntersectionNode(i_fluid_element,distances,1,2,intNode_B); // Intersection Node B
                CalcIntersectionNode(i_fluid_element,distances,3,0,intNode_C); // Intersection Node C
                CalcIntersectionNode(i_fluid_element,distances,3,2,intNode_D); // Intersection Node D
            }
            else if( distances[3] > 0 ) // nodes [1,2] negative
            {
                CalcIntersectionNode(i_fluid_element,distances,1,0,intNode_A); // Intersection Node A
                CalcIntersectionNode(i_fluid_element,distances,1,3,intNode_B); // Intersection Node B
                CalcIntersectionNode(i_fluid_element,distances,2,0,intNode_C); // Intersection Node C
                CalcIntersectionNode(i_fluid_element,distances,2,3,intNode_D); // Intersection Node D
            }
        }

        // now combine these intersection nodes to 2 triangles
        // Applying the ordered numbering scheme, the 2 triangles are always
        // formed out of nodes ABC and BCD

        // triangle 1: ABC
        triangle1[0] = ID_Node_A;
        triangle1[1] = ID_Node_B;
        triangle1[2] = ID_Node_C;

        // triangle 2: BCD
        triangle2[0] = ID_Node_B;
        triangle2[1] = ID_Node_C;
        triangle2[2] = ID_Node_D;

        nodesMatrix.push_back(intNode_A);
        nodesMatrix.push_back(intNode_B);
        nodesMatrix.push_back(intNode_C);
        nodesMatrix.push_back(intNode_D);

        AdaptIndexingByNormal(i_fluid_element,distances,triangle1,nodesMatrix);
        AdaptIndexingByNormal(i_fluid_element,distances,triangle2,nodesMatrix);

        elementsMatrix.push_back(triangle1);
        elementsMatrix.push_back(triangle2);
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void CalcIntersectionNode( ModelPart::ElementsContainerType::iterator&  i_fluid_element,
                               array_1d<double,4>                           distances,
                               unsigned int                                 Index1,
                               unsigned int                                 Index2,
                               array_1d<double,3> &                         intNodeCoord )
    {
        // Interpolate the intersection node given two tet nodes and the distance values

        // tet_edge is intersected --> calculate one triangle node by interpolation
        Element::GeometryType& geom = i_fluid_element->GetGeometry(); // Nodos del elemento
        array_1d<double,3> node1 = geom[Index1].Coordinates();
        array_1d<double,3> node2 = geom[Index2].Coordinates();

        // Compute intersection node by means of linear interpolation
        double scale = (-distances[Index1] / (distances[Index2] - distances[Index1]));
        double x_Int = node1[0] + (node2[0] - node1[0]) * scale;
        double y_Int = node1[1] + (node2[1] - node1[1]) * scale;
        double z_Int = node1[2] + (node2[2] - node1[2]) * scale;

        intNodeCoord[0] = x_Int;
        intNodeCoord[1] = y_Int;
        intNodeCoord[2] = z_Int;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void AdaptIndexingByNormal( ModelPart::ElementsContainerType::iterator& i_fluid_element,
                                array_1d<double,4>                          distances,
                                array_1d<unsigned int,3>&                   IDNodesInTri,
                                std::vector< array_1d<double,3> >           nodesMatrix )
    {
        // Check, if normal of reproduced triangular element points into the direction of a tet node with positive dist,
        // else adapt nodal indexes of the element

        // First, calculate normal of triangle
        array_1d<double,3> triEdge12; // vector connecting node 1 and 2 of the reproduced triangle
        array_1d<double,3> triEdge13; // vector connecting node 1 and 3 of the reproduced triangle

        triEdge12[0] = nodesMatrix[IDNodesInTri[1]-1][0] - nodesMatrix[IDNodesInTri[0]-1][0];
        triEdge12[1] = nodesMatrix[IDNodesInTri[1]-1][1] - nodesMatrix[IDNodesInTri[0]-1][1];
        triEdge12[2] = nodesMatrix[IDNodesInTri[1]-1][2] - nodesMatrix[IDNodesInTri[0]-1][2];

        triEdge13[0] = nodesMatrix[IDNodesInTri[2]-1][0] - nodesMatrix[IDNodesInTri[0]-1][0];
        triEdge13[1] = nodesMatrix[IDNodesInTri[2]-1][1] - nodesMatrix[IDNodesInTri[0]-1][1];
        triEdge13[2] = nodesMatrix[IDNodesInTri[2]-1][2] - nodesMatrix[IDNodesInTri[0]-1][2];

        // Compute normal by cross product of the triangle edges
        array_1d<double,3> normalVec;
        MathUtils<double>::CrossProduct(normalVec,triEdge12,triEdge13);
        normalVec *= 0.5;

        // Check, if normal vector is pointing into a tetnode with a positive distance value
        Element::GeometryType& geom = i_fluid_element->GetGeometry();
        array_1d<double,3> tetNode = geom[0].Coordinates();
        array_1d<double,3> intNodeTetNodeVec; // Vector connecting the tet node and the first intersection node
        intNodeTetNodeVec[0] = tetNode[0] - nodesMatrix[IDNodesInTri[0]-1][0];
        intNodeTetNodeVec[1] = tetNode[1] - nodesMatrix[IDNodesInTri[0]-1][1];
        intNodeTetNodeVec[2] = tetNode[2] - nodesMatrix[IDNodesInTri[0]-1][2];

        // Determine inner product as measure for the orientation of normal vector relative to tet node 0
        double InnerProduct = inner_prod(intNodeTetNodeVec,normalVec);

        if((InnerProduct < 0.0 and distances[0] > 0) or (InnerProduct > 0.0 and distances[0] < 0))
        {
            // in this case, the normal is NOT pointing towards a tet node with a positive value
            // --> switch two arbitrary nodes of the triangle
            // IDNodesInTri[0] = IDNodesInTri[0];
            unsigned int getIndex1 = IDNodesInTri[1];
            IDNodesInTri[1] = IDNodesInTri[2];
            IDNodesInTri[2] = getIndex1;
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void GenerateOctree()
    {
        Timer::Start("Generating Octree");
        std::cout << "Generating the Octree..." << std::endl;

        //mOctree.RefineWithUniformSize(0.0625);

        // loop over all fluid elements
        // this loop is parallelized using openmp

#ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
#else
        int number_of_threads = 1;
#endif

        ModelPart::NodesContainerType& pNodes = mrSkinModelPart.Nodes();

        vector<unsigned int> node_partition;
        CreatePartition(number_of_threads, pNodes.size(), node_partition);
        KRATOS_WATCH(number_of_threads);
        KRATOS_WATCH(node_partition);

#pragma omp parallel for
        for (int k = 0; k < number_of_threads; k++)
        {
            ModelPart::NodesContainerType::iterator it_begin = pNodes.ptr_begin() + node_partition[k];
            ModelPart::NodesContainerType::iterator it_end = pNodes.ptr_begin() + node_partition[k+1];

            // assemble all nodes
            for (ModelPart::NodeIterator it = it_begin; it != it_end; ++it)
            {
                double temp_point[3];
                temp_point[0] = it->Coordinate(1);
                temp_point[1] = it->Coordinate(2);
                temp_point[2] = it->Coordinate(3);
                mOctree.Insert(temp_point);
            }
        }

        //mOctree.Constrain2To1(); // To be removed. Pooyan.

        ModelPart::ElementsContainerType& pElements = mrSkinModelPart.Elements();

        vector<unsigned int> element_partition;
        CreatePartition(number_of_threads, pElements.size(), element_partition);
        KRATOS_WATCH(number_of_threads);
        KRATOS_WATCH(element_partition);

#pragma omp parallel for
        for (int k = 0; k < number_of_threads; k++)
        {
            ModelPart::ElementsContainerType::iterator it_begin = pElements.ptr_begin() + element_partition[k];
            ModelPart::ElementsContainerType::iterator it_end = pElements.ptr_begin() + element_partition[k+1];

            for (ModelPart::ElementIterator it = it_begin; it != it_end; ++it)
            {
                mOctree.Insert(*(it).base());
            }
        }

        Timer::Stop("Generating Octree");

        //KRATOS_WATCH(mOctree);

        //          std::cout << "######## WRITING OCTREE MESH #########" << std::endl;
        //          std::ofstream myfile;
        //          myfile.open ("octree.post.msh");
        //          mOctree.PrintGiDMesh(myfile);
        //          myfile.close();

        std::cout << "Generating the Octree finished" << std::endl;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void GenerateNodes()
    {
        Timer::Start("Generating Nodes");
        std::vector<OctreeType::cell_type*> all_leaves;
        mOctree.GetAllLeavesVector(all_leaves);

#pragma omp parallel for
        for (unsigned int i = 0; i < all_leaves.size(); i++)
        {
            *(all_leaves[i]->pGetDataPointer()) = ConfigurationType::AllocateData();
        }


        std::size_t last_id = mrBodyModelPart.NumberOfNodes() + 1;
        KRATOS_WATCH(all_leaves.size());
        for (std::size_t i = 0; i < all_leaves.size(); i++)
        {
            KRATOS_WATCH(i)
                    CellType* cell = all_leaves[i];
            GenerateCellNode(cell, last_id);
        }

        Timer::Stop("Generating Nodes");

    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void GenerateCellNode(CellType* pCell, std::size_t& LastId)
    {
        for (int i_pos=0; i_pos < 8; i_pos++) // position 8 is for center
        {
            DistanceSpatialContainersConfigure::cell_node_data_type* p_node = (*(pCell->pGetData()))[i_pos];
            if(p_node == 0)
            {
                (*(pCell->pGetData()))[i_pos] = new DistanceSpatialContainersConfigure::cell_node_data_type;

                (*(pCell->pGetData()))[i_pos]->Id() = LastId++;
                KRATOS_WATCH(LastId)

                        mOctreeNodes.push_back((*(pCell->pGetData()))[i_pos]);

                SetNodeInNeighbours(pCell,i_pos,(*(pCell->pGetData()))[i_pos]);
            }

        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void SetNodeInNeighbours(CellType* pCell, int Position, CellNodeDataType* pNode)
    {
        CellType::key_type point_key[3];
        pCell->GetKey(Position, point_key);

        for (std::size_t i_direction = 0; i_direction < 8; i_direction++) {
            CellType::key_type neighbour_key[3];
            if (pCell->GetNeighbourKey(Position, i_direction, neighbour_key)) {
                CellType* neighbour_cell = mOctree.pGetCell(neighbour_key);
                if (!neighbour_cell || (neighbour_cell == pCell))
                    continue;

                std::size_t position = neighbour_cell->GetLocalPosition(point_key);
                if((*neighbour_cell->pGetData())[position])
                {
                    std::cout << "ERROR!! Bad Position calculated!!!!!!!!!!! position :" << position << std::endl;
                    continue;
                }

                (*neighbour_cell->pGetData())[position] = pNode;
            }
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void CalculateDistance2()
    {
        Timer::Start("Calculate Distances2");
        ModelPart::NodesContainerType::ContainerType& nodes = mrFluidModelPart.NodesArray();
        int nodes_size = nodes.size();
        //         // first of all we reset the node distance to 1.00 which is the maximum distnace in our normalized space.
        //#pragma omp parallel for firstprivate(nodes_size)
        //         for(int i = 0 ; i < nodes_size ; i++)
        //             nodes[i]->GetSolutionStepValue(DISTANCE) = 1.00;

        std::vector<CellType*> leaves;

        mOctree.GetAllLeavesVector(leaves);
        //int leaves_size = leaves.size();

        //         for(int i = 0 ; i < leaves_size ; i++)
        //             CalculateNotEmptyLeavesDistance(leaves[i]);

#pragma omp parallel for firstprivate(nodes_size)
        for(int i = 0 ; i < nodes_size ; i++)
        {
            CalculateNodeDistance(*(nodes[i]));
        }
        Timer::Stop("Calculate Distances2");

    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    //     void CalculateDistance3()
    //     {
    //         Timer::Start("Calculate Distances2");
    //         ModelPart::NodesContainerType::ContainerType& nodes = mrFluidModelPart.NodesArray();
    //         int nodes_size = nodes.size();
    ////         // first of all we reset the node distance to 1.00 which is the maximum distnace in our normalized space.
    //#pragma omp parallel for firstprivate(nodes_size)
    //         for(int i = 0 ; i < nodes_size ; i++)
    //             nodes[i]->GetSolutionStepValue(DISTANCE) = 1.00;

    //           std::vector<CellType*> leaves;

    //         mOctree.GetAllLeavesVector(leaves);
    //         int leaves_size = leaves.size();

    //         for(int i = 0 ; i < leaves_size ; i++)
    //             CalculateNotEmptyLeavesDistance(leaves[i]);

    //#pragma omp parallel for firstprivate(nodes_size)
    //         for(int i = 0 ; i < nodes_size ; i++)
    //         {
    //             CalculateNodeDistance(*(nodes[i]));
    //         }
    //         Timer::Stop("Calculate Distances2");

    //     }
    //     void CalculateDistance4()
    //     {
    //         Timer::Start("Calculate Distances3");
    //         ModelPart::NodesContainerType::ContainerType& nodes = mrFluidModelPart.NodesArray();
    //         int nodes_size = nodes.size();
    //           std::vector<CellType*> leaves;

    //         mOctree.GetAllLeavesVector(leaves);
    //         int leaves_size = leaves.size();

    //#pragma omp parallel for firstprivate(nodes_size)
    //         for(int i = 0 ; i < nodes_size ; i++)
    //         {
    //             CalculateNodeDistanceFromCell(*(nodes[i]));
    //         }
    //         Timer::Stop("Calculate Distances3");

    //     }

    void CalculateDistance()
    {
        Timer::Start("Calculate Distances");
        DistanceSpatialContainersConfigure::data_type& nodes = mOctreeNodes;
        int nodes_size = nodes.size();
        // first of all we reste the node distance to 1.00 which is the maximum distnace in our normalized space.
#pragma omp parallel for firstprivate(nodes_size)
        for(int i = 0 ; i < nodes_size ; i++)
            nodes[i]->Distance() = 1.00;


        std::vector<CellType*> leaves;

        mOctree.GetAllLeavesVector(leaves);
        int leaves_size = leaves.size();

        for(int i = 0 ; i < leaves_size ; i++)
            CalculateNotEmptyLeavesDistance(leaves[i]);

        for(int i_direction = 0 ; i_direction < 1 ; i_direction++)
        {

            //#pragma omp parallel for firstprivate(nodes_size)
            for(int i = 0 ; i < nodes_size ; i++)
            {
                if(nodes[i]->X() < 1.00 && nodes[i]->Y() < 1.00 && nodes[i]->Z() < 1.00)
                    // if((*nodes[i])[i_direction] == 0.00)
                    CalculateDistance(*(nodes[i]), i_direction);
            }
        }
        Timer::Stop("Calculate Distances");

    }

    void CalculateDistance(CellNodeDataType& rNode, int i_direction)
    {
        double coords[3] = {rNode.X(), rNode.Y(), rNode.Z()};
        // KRATOS_WATCH_3(coords);

        //This function must color the positions in space defined by 'coords'.
        //coords is of dimension (3) normalized in (0,1)^3 space

        typedef Element::GeometryType triangle_type;
        typedef std::vector<std::pair<double, triangle_type*> > intersections_container_type;

        intersections_container_type intersections;
        DistanceSpatialContainersConfigure::data_type nodes_array;


        const double epsilon = 1e-12;

        double distance = 1.0;

        // Creating the ray
        double ray[3] = {coords[0], coords[1], coords[2]};
        ray[i_direction] = 0; // starting from the lower extreme

        //            KRATOS_WATCH_3(ray)
        GetIntersectionsAndNodes(ray, i_direction, intersections, nodes_array);
        //            KRATOS_WATCH(nodes_array.size())
        for (std::size_t i_node = 0; i_node < nodes_array.size() ; i_node++)
        {
            double coord = nodes_array[i_node]->Coordinate(i_direction+1);
            //             KRATOS_WATCH(intersections.size());

            int ray_color= 1;
            std::vector<std::pair<double, Element::GeometryType*> >::iterator i_intersection = intersections.begin();
            while (i_intersection != intersections.end()) {
                double d = coord - i_intersection->first;
                if (d > epsilon) {

                    ray_color = -ray_color;
                    distance = d;
                } else if (d > -epsilon) {//interface
                    distance = 0.00;
                    break;
                } else {
                    if(distance > -d)
                        distance = -d;
                    break;
                }

                i_intersection++;
            }

            distance *= ray_color;

            double& node_distance = nodes_array[i_node]->Distance();
            if(fabs(distance) < fabs(node_distance))
                node_distance = distance;
            else if (distance*node_distance < 0.00) // assigning the correct sign
                node_distance = -node_distance;


        }
    }

    void CalculateNotEmptyLeavesDistance(CellType* pCell)
    {
        typedef Element::GeometryType triangle_type;
        typedef OctreeType::cell_type::object_container_type object_container_type;

        object_container_type* objects = (pCell->pGetObjects());

        // There are no intersection in empty cells
        if (objects->empty())
            return;


        for (int i_pos=0; i_pos < 8; i_pos++) // position 8 is for center
        {
            double distance = 1.00; // maximum distance is 1.00

            for(object_container_type::iterator i_object = objects->begin(); i_object != objects->end(); i_object++)
            {
                CellType::key_type keys[3];
                pCell->GetKey(i_pos,keys);

                double cell_point[3];
                mOctree.CalculateCoordinates(keys,cell_point);

                //                cell_point[0] = pCell->GetCoordinate(keys[0]);
                //                cell_point[1] = pCell->GetCoordinate(keys[1]);
                //                cell_point[2] = pCell->GetCoordinate(keys[2]);

                double d = GeometryUtils::PointDistanceToTriangle3D((*i_object)->GetGeometry()[0], (*i_object)->GetGeometry()[1], (*i_object)->GetGeometry()[2], Point<3>(cell_point[0], cell_point[1], cell_point[2]));

                if(d < distance)
                    distance = d;
            }

            double& node_distance = (*(pCell->pGetData()))[i_pos]->Distance();
            if(distance < node_distance)
                node_distance = distance;

        }

    }

    void CalculateNodeDistance(Node<3>& rNode)
    {
        double coord[3] = {rNode.X(), rNode.Y(), rNode.Z()};
        double distance = DistancePositionInSpace(coord);
        double& node_distance =  rNode.GetSolutionStepValue(DISTANCE);

        //const double epsilon = 1.00e-12;
        if(fabs(node_distance) > fabs(distance))
            node_distance = distance;
        else if (distance*node_distance < 0.00) // assigning the correct sign
            node_distance = -node_distance;
    }

    //      void CalculateNodeDistanceFromCell(Node<3>& rNode)
    //      {
    //          OctreeType::key_type node_key[3] = {octree->CalcKeyNormalized(rNode.X()), octree->CalcKeyNormalized(rNode.Y()), octree->CalcKeyNormalized(rNode.Z())};
    //          OctreeType::cell_type* pcell = octree->pGetCell(node_key);

    //          object_container_type* objects = (pCell->pGetObjects());

    //          // We interpolate the cell distances for the node in empty cells
    //          if (objects->empty())
    //          {

    //          }

    //          double distance = DistancePositionInSpace(coord);
    //          double& node_distance =  rNode.GetSolutionStepValue(DISTANCE);

    //          //const double epsilon = 1.00e-12;
    //          if(fabs(node_distance) > fabs(distance))
    //            node_distance = distance;
    //          else if (distance*node_distance < 0.00) // assigning the correct sign
    //              node_distance = -node_distance;
    //      }

    double DistancePositionInSpace(double* coords)
    {
        //This function must color the positions in space defined by 'coords'.
        //coords is of dimension (3) normalized in (0,1)^3 space

        typedef Element::GeometryType triangle_type;
        typedef std::vector<std::pair<double, triangle_type*> > intersections_container_type;

        intersections_container_type intersections;

        const int dimension = 3;
        const double epsilon = 1e-12;

        double distances[3] = {1.0, 1.0, 1.0};

        for (int i_direction = 0; i_direction < dimension; i_direction++)
        {
            // Creating the ray
            double ray[3] = {coords[0], coords[1], coords[2]};
            ray[i_direction] = 0; // starting from the lower extreme

            GetIntersections(ray, i_direction, intersections);

            //                if(intersections.size() == 1)
            //                    KRATOS_WATCH_3(ray)

            //             KRATOS_WATCH(intersections.size());

            int ray_color= 1;
            std::vector<std::pair<double, Element::GeometryType*> >::iterator i_intersection = intersections.begin();
            while (i_intersection != intersections.end()) {
                double d = coords[i_direction] - i_intersection->first;
                if (d > epsilon) {

                    ray_color = -ray_color;
                    distances[i_direction] = d;
                    //                        if(distances[i_direction] > d) // I think this is redundunt. Pooyan.
                    //                        {
                    //                            if(ray_color > 0.00)
                    //                                    distances[i_direction] = d;
                    //                            else
                    //                                distances[i_direction] = -d;
                    //                        }
                } else if (d > -epsilon) {//interface
                    distances[i_direction] = 0.00;
                    break;
                } else {
                    if(distances[i_direction] > -d)
                        distances[i_direction] = -d;
                    break;
                }

                i_intersection++;
            }

            distances[i_direction] *= ray_color;
        }

        //            if(distances[0]*distances[1] < 0.00 || distances[2]*distances[1] < 0.00)
        //                KRATOS_WATCH_3(distances);

#ifdef _DEBUG
        std::cout << "colors : " << colors[0] << ", " << colors[1] << ", " << colors[2] << std::endl;
#endif
        double distance = (fabs(distances[0]) > fabs(distances[1])) ? distances[1] : distances[0];
        distance = (fabs(distance) > fabs(distances[2])) ? distances[2] : distance;

        return distance;

    }

    void GetIntersectionsAndNodes(double* ray, int direction, std::vector<std::pair<double,Element::GeometryType*> >& intersections, DistanceSpatialContainersConfigure::data_type& rNodesArray)
    {
        //This function passes the ray through the model and gives the hit point to all objects in its way
        //ray is of dimension (3) normalized in (0,1)^3 space
        // direction can be 0,1,2 which are x,y and z respectively

        const double epsilon = 1.00e-12;

        // first clearing the intersections points vector
        intersections.clear();

        OctreeType* octree = &mOctree;

        OctreeType::key_type ray_key[3] = {octree->CalcKeyNormalized(ray[0]), octree->CalcKeyNormalized(ray[1]), octree->CalcKeyNormalized(ray[2])};
        OctreeType::key_type cell_key[3];

        // getting the entrance cell from lower extreme
        ray_key[direction] = 0;
        OctreeType::cell_type* cell = octree->pGetCell(ray_key);

        while (cell) {
            std::size_t position = cell->GetLocalPosition(ray_key); // Is this the local position!?!?!?!
            OctreeType::key_type node_key[3];
            cell->GetKey(position, node_key);
            if((node_key[0] == ray_key[0]) && (node_key[1] == ray_key[1]) && (node_key[2] == ray_key[2]))
            {
                if(cell->pGetData())
                {
                    if(cell->pGetData()->size() > position)
                    {
                        CellNodeDataType* p_node = (*cell->pGetData())[position];
                        if(p_node)
                        {
                            //KRATOS_WATCH(p_node->Id())
                            rNodesArray.push_back(p_node);
                        }
                    }
                    else
                        KRATOS_WATCH(cell->pGetData()->size())
                }
            }


            //        std::cout << ".";
            GetCellIntersections(cell, ray, ray_key, direction, intersections);

            // Add the cell's middle node if existed
            //      cell->GetKey(8, cell_key); // 8 is the central position
            //      ray_key[direction]=cell_key[direction]; // positioning the ray in the middle of cell in its direction

            //      position = cell->GetLocalPosition(ray_key);
            //      if(position < 27) // principal nodes
            //      {
            //          if(cell->pGetData())
            //          {
            //              if(cell->pGetData()->size() > position)
            //              {
            //                  Node<3>* p_node = (*cell->pGetData())[position];
            //                  if(p_node)
            //                  {
            //                      //KRATOS_WATCH(p_node->Id())
            //                      rNodesArray.push_back(p_node);
            //                  }
            //              }
            //              else
            //                  KRATOS_WATCH(cell->pGetData()->size())
            //          }
            //      }
            //      else
            //      {
            //          KRATOS_WATCH(position);
            //          KRATOS_WATCH(*cell);
            //      }


            // go to the next cell
            if (cell->GetNeighbourKey(1 + direction * 2, cell_key)) {
                ray_key[direction] = cell_key[direction];
                cell = octree->pGetCell(ray_key);
                ray_key[direction] -= 1 ;//the key returned by GetNeighbourKey is inside the cell (minkey +1), to ensure that the corresponding
                //cell get in pGetCell is the right one.
#ifdef _DEBUG
                Octree_Pooyan::key_type min_key[3];
                cell->GetMinKey(min_key[0],min_key[1],min_key[2]);
                Octree_Pooyan::key_type tmp;
                tmp= min_key[direction];
                assert(ray_key[direction]==tmp);
#endif
            } else
                cell = NULL;
        }



        //   KRATOS_WATCH(rNodesArray.size());
        // now eliminating the repeated objects
        if (!intersections.empty()) {
            //sort
            std::sort(intersections.begin(), intersections.end());
            // unique
            std::vector<std::pair<double, Element::GeometryType*> >::iterator i_begin = intersections.begin();
            std::vector<std::pair<double, Element::GeometryType*> >::iterator i_intersection = intersections.begin();
            while (++i_begin != intersections.end()) {
                // considering the very near points as the same points
                if (fabs(i_begin->first - i_intersection->first) > epsilon) // if the hit points are far enough they are not the same
                    *(++i_intersection) = *i_begin;
            }
            intersections.resize((++i_intersection) - intersections.begin());

        }
    }

    void GetIntersections(double* ray, int direction, std::vector<std::pair<double,Element::GeometryType*> >& intersections)
    {
        //This function passes the ray through the model and gives the hit point to all objects in its way
        //ray is of dimension (3) normalized in (0,1)^3 space
        // direction can be 0,1,2 which are x,y and z respectively

        const double epsilon = 1.00e-12;

        // first clearing the intersections points vector
        intersections.clear();

        OctreeType* octree = &mOctree;

        OctreeType::key_type ray_key[3] = {octree->CalcKeyNormalized(ray[0]), octree->CalcKeyNormalized(ray[1]), octree->CalcKeyNormalized(ray[2])};
        OctreeType::key_type cell_key[3];

        // getting the entrance cell from lower extreme
        OctreeType::cell_type* cell = octree->pGetCell(ray_key);

        while (cell) {
            //        std::cout << ".";
            GetCellIntersections(cell, ray, ray_key, direction, intersections);
            // go to the next cell
            if (cell->GetNeighbourKey(1 + direction * 2, cell_key)) {
                ray_key[direction] = cell_key[direction];
                cell = octree->pGetCell(ray_key);
                ray_key[direction] -= 1 ;//the key returned by GetNeighbourKey is inside the cell (minkey +1), to ensure that the corresponding
                //cell get in pGetCell is the right one.
#ifdef _DEBUG
                Octree_Pooyan::key_type min_key[3];
                cell->GetMinKey(min_key[0],min_key[1],min_key[2]);
                Octree_Pooyan::key_type tmp;
                tmp= min_key[direction];
                assert(ray_key[direction]==tmp);
#endif
            } else
                cell = NULL;
        }


        // now eliminating the repeated objects
        if (!intersections.empty()) {
            //sort
            std::sort(intersections.begin(), intersections.end());
            // unique
            std::vector<std::pair<double, Element::GeometryType*> >::iterator i_begin = intersections.begin();
            std::vector<std::pair<double, Element::GeometryType*> >::iterator i_intersection = intersections.begin();
            while (++i_begin != intersections.end()) {
                // considering the very near points as the same points
                if (fabs(i_begin->first - i_intersection->first) > epsilon) // if the hit points are far enough they are not the same
                    *(++i_intersection) = *i_begin;
            }
            intersections.resize((++i_intersection) - intersections.begin());

        }
    }

    int GetCellIntersections(OctreeType::cell_type* cell, double* ray,
                             OctreeType::key_type* ray_key, int direction,
                             std::vector<std::pair<double, Element::GeometryType*> >& intersections)  {
        //This function passes the ray through the cell and gives the hit point to all objects in its way
        //ray is of dimension (3) normalized in (0,1)^3 space
        // direction can be 0,1,2 which are x,y and z respectively

        typedef Element::GeometryType triangle_type;
        typedef OctreeType::cell_type::object_container_type object_container_type;

        object_container_type* objects = (cell->pGetObjects());

        // There are no intersection in empty cells
        if (objects->empty())
            return 0;

        //      std::cout << "X";
        // calculating the two extreme of the ray segment inside the cell
        double ray_point1[3] = {ray[0], ray[1], ray[2]};
        double ray_point2[3] = {ray[0], ray[1], ray[2]};
        double normalized_coordinate;
        mOctree.CalculateCoordinateNormalized(ray_key[direction], normalized_coordinate);
        ray_point1[direction] = normalized_coordinate;
        ray_point2[direction] = ray_point1[direction] + mOctree.CalcSizeNormalized(cell);

        mOctree.ScaleBackToOriginalCoordinate(ray_point1);
        mOctree.ScaleBackToOriginalCoordinate(ray_point2);

        for (object_container_type::iterator i_object = objects->begin(); i_object != objects->end(); i_object++) {
            double intersection[3]={0.00,0.00,0.00};

            int is_intersected = IntersectionTriangleSegment((*i_object)->GetGeometry(), ray_point1, ray_point2, intersection); // This intersection has to be optimized for axis aligned rays

            if (is_intersected == 1) // There is an intersection but not coplanar
                intersections.push_back(std::pair<double, Element::GeometryType*>(intersection[direction], &((*i_object)->GetGeometry())));
            //else if(is_intersected == 2) // coplanar case
        }

        return 0;
    }

    int IntersectionTriangleSegment(Element::GeometryType& rGeometry, double* RayPoint1, double* RayPoint2, double* IntersectionPoint)
    {
        // This is the adaption of the implemnetation provided in:
        // http://www.softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm#intersect_RayTriangle()

        const double epsilon = 1.00e-12;

        array_1d<double,3>    u, v, n;             // triangle vectors
        array_1d<double,3>    dir, w0, w;          // ray vectors
        double     r, a, b;             // params to calc ray-plane intersect


        // get triangle edge vectors and plane normal
        u = rGeometry[1] - rGeometry[0];
        v = rGeometry[2] - rGeometry[0];

        MathUtils<double>::CrossProduct(n, u, v);             // cross product

        if (norm_2(n) == 0)            // triangle is degenerate
            return -1;                 // do not deal with this case

        for(int i = 0 ; i < 3 ; i++)
        {
            dir[i] = RayPoint2[i] - RayPoint1[i];             // ray direction vector
            w0[i] = RayPoint1[i] - rGeometry[0][i];
        }

        a = -inner_prod(n,w0);
        b = inner_prod(n,dir);

        if (fabs(b) < epsilon) {     // ray is parallel to triangle plane
            if (a == 0)                // ray lies in triangle plane
                return 2;
            else return 0;             // ray disjoint from plane
        }

        // get intersect point of ray with triangle plane
        r = a / b;
        if (r < 0.0)                   // ray goes away from triangle
            return 0;                  // => no intersect
        // for a segment, also test if (r > 1.0) => no intersect

        for(int i = 0 ; i < 3 ; i++)
            IntersectionPoint[i]  = RayPoint1[i] + r * dir[i];           // intersect point of ray and plane

        // is I inside T?
        double    uu, uv, vv, wu, wv, D;
        uu = inner_prod(u,u);
        uv = inner_prod(u,v);
        vv = inner_prod(v,v);


        for(int i = 0 ; i < 3 ; i++)
            w[i] = IntersectionPoint[i] - rGeometry[0][i];


        wu = inner_prod(w,u);
        wv = inner_prod(w,v);
        D = uv * uv - uu * vv;

        // get and test parametric coords
        double s, t;
        s = (uv * wv - vv * wu) / D;
        if (s < 0.0 - epsilon || s > 1.0 + epsilon)        // I is outside T
            return 0;
        t = (uv * wu - uu * wv) / D;
        if (t < 0.0 - epsilon || (s + t) > 1.0 + epsilon)  // I is outside T
            return 0;

        return 1;                      // I is in T

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
        return "CalculateSignedDistanceTo3DSkinProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "CalculateSignedDistanceTo3DSkinProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    void PrintGiDMesh(std::ostream & rOStream) const {
        std::vector<CellType*> leaves;

        mOctree.GetAllLeavesVector(leaves);

        std::cout << "writing " << leaves.size() << " leaves" << std::endl;
        rOStream << "MESH \"leaves\" dimension 3 ElemType Hexahedra Nnode 8" << std::endl;
        rOStream << "# color 96 96 96" << std::endl;
        rOStream << "Coordinates" << std::endl;
        rOStream << "# node number coordinate_x coordinate_y coordinate_z  " << std::endl;

        for(DistanceSpatialContainersConfigure::data_type::const_iterator i_node = mOctreeNodes.begin() ; i_node != mOctreeNodes.end() ; i_node++)
        {
            rOStream << (*i_node)->Id() << "  " << (*i_node)->Coordinate(1) << "  " << (*i_node)->Coordinate(2) << "  " << (*i_node)->Coordinate(3) << std::endl;
            //mOctree.Insert(temp_point);
        }
        std::cout << "Nodes written..." << std::endl;
        rOStream << "end coordinates" << std::endl;
        rOStream << "Elements" << std::endl;
        rOStream << "# element node_1 node_2 node_3 material_number" << std::endl;

        for (std::size_t i = 0; i < leaves.size(); i++) {
            if ((leaves[i]->pGetData()))
            {
                DistanceSpatialContainersConfigure::data_type& nodes = (*(leaves[i]->pGetData()));

                rOStream << i + 1;
                for(int j = 0 ; j < 8 ; j++)
                    rOStream << "  " << nodes[j]->Id();
                rOStream << std::endl;
            }
        }
        rOStream << "end elements" << std::endl;

    }

    void PrintGiDResults(std::ostream & rOStream) const {
        std::vector<CellType*> leaves;

        mOctree.GetAllLeavesVector(leaves);

        rOStream << "GiD Post Results File 1.0" << std::endl << std::endl;

        rOStream << "Result \"Distance\" \"Kratos\" 1 Scalar OnNodes" << std::endl;

        rOStream << "Values" << std::endl;

        for(DistanceSpatialContainersConfigure::data_type::const_iterator i_node = mOctreeNodes.begin() ; i_node != mOctreeNodes.end() ; i_node++)
        {
            rOStream << (*i_node)->Id() << "  " << (*i_node)->Distance() << std::endl;
        }
        rOStream << "End Values" << std::endl;

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
    ModelPart& mrSkinModelPart;
    ModelPart& mrBodyModelPart;
    ModelPart& mrFluidModelPart;

    DistanceSpatialContainersConfigure::data_type mOctreeNodes;

    OctreeType mOctree;

    static const double epsilon;

    /**
         * @}
         */
    /**
         * calculates the eigenvectors and eigenvalues of given symmetric matrix A.
         * The eigenvectors and eigenvalues are calculated using the iterative
         * Gauss-Seidel-method
         * @param A the given symmetric matrix the eigenvectors are to be calculated.
         * :WARNING: Matrix A will be overwritten and has to be symmetric
         * @param V the result matrix (will be overwritten with the eigenvectors)
         * @param zero_tolerance the largest value considered to be zero
         */

    static inline void EigenVectors(const Matrix& A, Matrix& vectors, Vector& lambda, double zero_tolerance =1e-9, int max_iterations = 10)
    {
        Matrix Help= A;

        for(int i=0; i<3; i++)
            for(int j=0; j<3; j++)
                Help(i,j)= Help(i,j);


        vectors.resize(Help.size1(),Help.size2(),false);

        lambda.resize(Help.size1(),false);

        Matrix HelpDummy(Help.size1(),Help.size2());

        bool is_converged = false;

        Matrix unity=ZeroMatrix(Help.size1(),Help.size2());

        for(unsigned int i=0; i< Help.size1(); i++)
            unity(i,i)= 1.0;

        Matrix V= unity;

        Matrix VDummy(Help.size1(),Help.size2());

        Matrix Rotation(Help.size1(),Help.size2());


        for(int iterations=0; iterations<max_iterations; iterations++)
        {

            is_converged= true;

            double a= 0.0;

            unsigned int index1= 0;

            unsigned int index2= 1;

            for(unsigned int i=0; i< Help.size1(); i++)
            {
                for(unsigned int j=(i+1); j< Help.size2(); j++)
                {
                    if((fabs(Help(i,j)) > a ) && (fabs(Help(i,j)) > zero_tolerance))
                    {
                        a= fabs(Help(i,j));

                        index1= i;
                        index2= j;

                        is_converged= false;
                    }
                }
            }

            //                 KRATOS_WATCH(Help);

            if(is_converged)
                break;

            //Calculation of Rotationangle

            double gamma= (Help(index2,index2)-Help(index1,index1))/(2*Help(index1,index2));

            double u=1.0;

            if(fabs(gamma) > zero_tolerance && fabs(gamma)< (1/zero_tolerance))
            {
                u= gamma/fabs(gamma)*1.0/(fabs(gamma)+sqrt(1.0+gamma*gamma));
            }
            else
            {
                if  (fabs(gamma)>= (1.0/zero_tolerance))
                    u= 0.5/gamma;
            }

            double c= 1.0/(sqrt(1.0+u*u));

            double s= c*u;

            double teta= s/(1.0+c);

            //Ratotion of the Matrix
            HelpDummy= Help;

            HelpDummy(index2,index2)= Help(index2,index2)+u*Help(index1,index2);
            HelpDummy(index1,index1)= Help(index1,index1)-u*Help(index1,index2);
            HelpDummy(index1,index2)= 0.0;
            HelpDummy(index2,index1)= 0.0;

            for(unsigned int i=0; i<Help.size1(); i++)
            {
                if((i!= index1) && (i!= index2))
                {
                    HelpDummy(index2,i)=Help(index2,i)+s*(Help(index1,i)- teta*Help(index2,i));
                    HelpDummy(i,index2)=Help(index2,i)+s*(Help(index1,i)- teta*Help(index2,i));

                    HelpDummy(index1,i)=Help(index1,i)-s*(Help(index2,i)+ teta*Help(index1,i));
                    HelpDummy(i,index1)=Help(index1,i)-s*(Help(index2,i)+ teta*Help(index1,i));
                }
            }


            Help= HelpDummy;

            //Calculation of the eigenvectors V
            Rotation =unity;
            Rotation(index2,index1)=-s;
            Rotation(index1,index2)=s;
            Rotation(index1,index1)=c;
            Rotation(index2,index2)=c;

            //                 Help=ZeroMatrix(A.size1(),A.size1());

            VDummy = ZeroMatrix(Help.size1(), Help.size2());

            for(unsigned int i=0; i< Help.size1(); i++)
            {
                for(unsigned int j=0; j< Help.size1(); j++)
                {
                    for(unsigned int k=0; k< Help.size1(); k++)
                    {
                        VDummy(i,j) += V(i,k)*Rotation(k,j);
                    }
                }
            }
            V= VDummy;
        }

        if(!(is_converged))
        {
            std::cout<<"########################################################"<<std::endl;
            std::cout<<"Max_Iterations exceed in Jacobi-Seidel-Iteration (eigenvectors)"<<std::endl;
            std::cout<<"########################################################"<<std::endl;
        }

        for(unsigned int i=0; i< Help.size1(); i++)
        {
            for(unsigned int j=0; j< Help.size1(); j++)
            {
                vectors(i,j)= V(j,i);
            }
        }

        for(unsigned int i=0; i<Help.size1(); i++)
            lambda(i)= Help(i,i);

        return;
    }


    inline void CreatePartition(unsigned int number_of_threads, const int number_of_rows, vector<unsigned int>& partitions)
    {
        partitions.resize(number_of_threads + 1);
        int partition_size = number_of_rows / number_of_threads;
        partitions[0] = 0;
        partitions[number_of_threads] = number_of_rows;
        for (unsigned int i = 1; i < number_of_threads; i++)
            partitions[i] = partitions[i - 1] + partition_size;
    }

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
    CalculateSignedDistanceTo3DSkinProcess& operator=(CalculateSignedDistanceTo3DSkinProcess const& rOther);

    /// Copy constructor.
    //CalculateSignedDistanceTo3DSkinProcess(CalculateSignedDistanceTo3DSkinProcess const& rOther);


    ///@}

}; // Class CalculateSignedDistanceTo3DSkinProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  CalculateSignedDistanceTo3DSkinProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const CalculateSignedDistanceTo3DSkinProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

const double CalculateSignedDistanceTo3DSkinProcess::epsilon = 1e-18;


}  // namespace Kratos.

#endif // KRATOS_CALCULATE_DISTANCE_PROCESS_H_INCLUDED  defined 


