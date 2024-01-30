//
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    hbui
//

// System includes

// External includes

// Project includes
#include "spatial_containers/bounding_volume_tree.h"

namespace Kratos
{

const double kDOP::msDirection[][3] = {{1, 1, 1}};
const double _6DOP::msDirection[][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
const double _8DOP::msDirection[][3] = {{1, 1, 1}, {-1, 1, 1}, {1, -1, 1}, {1, 1, -1}};
const double _12DOP::msDirection[][3] = {{1, 1, 0}, {1, -1, 0}, {1, 0, 1}, {1, 0, -1}, {0, 1, 1}, {0, 1, -1}};
const double _14DOP::msDirection[][3] = {  {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
                                           {1, 1, 1}, {-1, 1, 1}, {1, -1, 1}, {1, 1, -1} };
const double _18DOP::msDirection[][3] = {  {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
                                           {1, 1, 0}, {1, -1, 0}, {1, 0, 1}, {1, 0, -1}, {0, 1, 1}, {0, 1, -1} };
const double _20DOP::msDirection[][3] = {  {1, 1, 1}, {-1, 1, 1}, {1, -1, 1}, {1, 1, -1},
                                           {1, 1, 0}, {1, -1, 0}, {1, 0, 1}, {1, 0, -1}, {0, 1, 1}, {0, 1, -1} };
const double _26DOP::msDirection[][3] = {  {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
                                           {1, 1, 1}, {-1, 1, 1}, {1, -1, 1}, {1, 1, -1},
                                           {1, 1, 0}, {1, -1, 0}, {1, 0, 1}, {1, 0, -1}, {0, 1, 1}, {0, 1, -1} };

kDOP::Array2DType kDOP::Direction() const {return msDirection;}
kDOP::Array2DType _6DOP::Direction() const {return msDirection;}
kDOP::Array2DType _8DOP::Direction() const {return msDirection;}
kDOP::Array2DType _12DOP::Direction() const {return msDirection;}
kDOP::Array2DType _14DOP::Direction() const {return msDirection;}
kDOP::Array2DType _18DOP::Direction() const {return msDirection;}
kDOP::Array2DType _20DOP::Direction() const {return msDirection;}
kDOP::Array2DType _26DOP::Direction() const {return msDirection;}

void LineRegressionVolumePartitioner::Partition(ConditionsContainerType& rAllConditions,
                                                const kDOP& rBoundingVolume,
                                                ConditionsContainerType& rOutputSet1,
                                                ConditionsContainerType& rOutputSet2)
{
    // firstly extract all the nodes of the providing geometry
    std::set<std::size_t> NodeIds;
    for(ConditionsContainerType::ptr_iterator it = rAllConditions.ptr_begin(); it != rAllConditions.ptr_end(); ++it)
        for(std::size_t j = 0; j < (*it)->GetGeometry().size(); ++j)
            NodeIds.insert((*it)->GetGeometry()[j].Id());

    // build the 3D regression line of the point set
    // Here I used a simple approach:
    //  +   I perform linear fit for (x,z) data, then I got the surface perpendicular to Oxz. Similarly, I fit (y,z) data to get the surface perpendicular to Oyz. Intersection of 2 surfaces gives me an approximation of best line fit of my 3D data
    //  +   Theoretically, it's best to use PCA (Principal Component Analysis) to reduce the 3D cloud data to a line. But it would be an overkill due to SVD.
    //  REF: http://stackoverflow.com/questions/24747643/3d-linear-regression
    //       http://stackoverflow.com/questions/2298390/fitting-a-line-in-3d
    // TODO
}

/// REF: https://github.com/brandonpelfrey/Fast-BVH
void SimpleBoundingVolumePartitioner::Partition(ConditionsContainerType& rAllConditions,
                                                const kDOP& rBoundingVolume,
                                                ConditionsContainerType& rOutputSet1,
                                                ConditionsContainerType& rOutputSet2)
{
    // for simplest case, if there is only one segment, then output 1 is set straightforward
    if(rAllConditions.size() == 1)
    {
        for(ConditionsContainerType::ptr_iterator it = rAllConditions.ptr_begin(); it != rAllConditions.ptr_end(); ++it)
            rOutputSet1.push_back(*it);
    }

    // for another simple case, if there are only 2 segments then divide by half and quit
    if(rAllConditions.size() == 2)
    {
        int cnt = 0;
        for(ConditionsContainerType::ptr_iterator it = rAllConditions.ptr_begin(); it != rAllConditions.ptr_end(); ++it)
        {
            if(cnt == 0)
                rOutputSet1.push_back(*it);
            else
                rOutputSet2.push_back(*it);
            ++cnt;
        }
        return;
    }

    // extract the longest axis of the bounding volume
    std::size_t longest_axis = rBoundingVolume.GetLongestAxis();
    const double* longest_dir = rBoundingVolume.Direction(longest_axis);

    // compute the centroid of the geometries underlying bounding volume
    double C[3];
    double Median[3];
    Median[0] = 0.0;
    Median[1] = 0.0;
    Median[2] = 0.0;
    for(ConditionsContainerType::ptr_iterator it = rAllConditions.ptr_begin(); it != rAllConditions.ptr_end(); ++it)
    {
        this->ComputeCentroid((*it)->GetGeometry(), C);
        Median[0] += C[0];
        Median[1] += C[1];
        Median[2] += C[2];
    }
    Median[0] /= rAllConditions.size();
    Median[1] /= rAllConditions.size();
    Median[2] /= rAllConditions.size();

    // divide each geometry in the InputSet to each of OutputSet
    for(ConditionsContainerType::ptr_iterator it = rAllConditions.ptr_begin(); it != rAllConditions.ptr_end(); ++it)
    {
        this->ComputeCentroid((*it)->GetGeometry(), C);
        double v = 0.0;
        for(std::size_t j = 0; j < 3; ++j)
            v += (C[j] - Median[j]) * longest_dir[j];
        if(v < 0)
            rOutputSet1.push_back(*it);
        else
            rOutputSet2.push_back(*it);
    }

    // make a size check
    if(rOutputSet1.size() == 0 || rOutputSet2.size() == 0)
        KRATOS_THROW_ERROR(std::logic_error, "The partitioning is wrong. One set of output cannot be empty", "")
}

}  // namespace Kratos.

