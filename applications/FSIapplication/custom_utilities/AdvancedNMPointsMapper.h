/*
 * File:   AdvancedNMPointsMapper.h
 * Author: jcotela
 *
 * Created on 19 January 2010, 10:20
 */

#if !defined(KRATOS_ADVANCED_NM_POINTS_MAPPER_H_INCLUDED )
#define  KRATOS_ADVANCED_NM_POINTS_MAPPER_H_INCLUDED

//#include "includes/define.h"
#include <iostream>
#include <vector>
#include "fsi_application.h"
#include "includes/model_part.h"
#include "containers/array_1d.h"
#include "spatial_containers/spatial_containers.h" //kd-tree
#include "utilities/math_utils.h" // Cross Product

namespace Kratos
{

    // Custom Gauss Point container to be used by the mapper
    class GaussPointItem: public Point<3>
    {
    public:

        typedef boost::numeric::ublas::bounded_matrix<double,3,3> MatrixVar;
        typedef boost::shared_ptr< GaussPointItem > Pointer;

        GaussPointItem(): // This Constructor is used by the tree
                Point<3>(),
                mArea(0),
                mProjStatus(0)
        {
            mNormal[0] = mNormal[1] = mNormal[2] = 0.0;
        }

        GaussPointItem(array_1d<double,3> Coords, double Area, array_1d<double,3> Normal):
                Point<3>(Coords),
                mArea(Area),
                mNormal(Normal),
                mProjStatus(0)
                {}

        GaussPointItem(const GaussPointItem& rhs): // Copy Constructor (not really requred)
                Point<3>(rhs),
                mArea(rhs.mArea),
                mNormal(rhs.mNormal),
                mProjStatus(rhs.mProjStatus),
                mpOriginCond(rhs.mpOriginCond),
                mpOriginNode(rhs.mpOriginNode)
        {
            mOriginCoords[0] = rhs.mOriginCoords[0];
            mOriginCoords[1] = rhs.mOriginCoords[1];
        }

//        void GetNormal(array_1d<double,3>& Normal) {Normal = mNormal;} // Unused accessors
//        void GetOriginLocalCoords(array_1d<double,2>& Coords ) {Coords = mOriginCoords;}
//        void GetOriginCond(Condition::WeakPointer& pCond) {pCond = mpOriginCond;}
//        void GetOriginNode(Node<3>::WeakPointer& pNode) {pNode = mpOriginNode;}
        void GetArea(double& Area) {Area = mArea;}
        void GetDist(double& Dist) {Dist = mDist;}
        void GetProjStatus(int& Proj) {Proj = mProjStatus;}

        boost::weak_ptr<Condition> GetOriginCond() {return mpOriginCond;} // TEST FUNCTION

        void SetProjection(Condition::WeakPointer Cond,
                array_1d<double,2> Coords, double Dist)
        {
            mpOriginCond = Cond;
            mOriginCoords = Coords;
            mDist = Dist;
            mProjStatus = 1;
        }

        void SetProjection(Node<3>::WeakPointer pNode,const double SqDist)
        {
            mpOriginNode = pNode;
            mDist = SqDist;
            mProjStatus = 2;
            mOriginCoords[0] = 0;
            mOriginCoords[1] = 0;
        }

        void Project(Condition::Pointer pOriginCond,
            array_1d<double,2> & Coords, double & Dist);

        void GetProjectedValue(const Variable<double> & rOriginVar, double& Value); //Scalar variables
        void GetProjectedValue(const Variable< array_1d<double,3> > & rOriginVar, array_1d<double,3>& Value); //Vector variables

    private:
        double mArea;
        array_1d<double,3> mNormal; // Destinationn condition's normal
        double mDist; // For GP projected to a Condition, Distance along Normal from Gauss Point to Condition
                      // For GP projected to a Node, SQUARED Distance to Node
        int mProjStatus;    // 0: Not Projected
                            // 1: Projected to a condition
                            // 2: Couldn't be projected, but a value can be obtained from a nearby node
        Condition::WeakPointer mpOriginCond;
        array_1d<double,2> mOriginCoords; // For GP projected to a condition
        Node<3>::WeakPointer mpOriginNode; // For GP projected to a node
    };


    // Mapper Class
    class AdvancedNMPointsMapper
    {
        // Type definitions for the tree
        typedef GaussPointItem                              PointType;
        typedef GaussPointItem::Pointer                     PointTypePointer;
        typedef std::vector<PointType::Pointer>             GaussPointVector;
        typedef std::vector<PointType::Pointer>::iterator   GaussPointIterator;
        typedef std::vector<double>                         DistanceVector;
        typedef std::vector<double>::iterator               DistanceIterator;

        // 3 if for dimension
        typedef Bucket< 3ul, PointType, GaussPointVector, PointTypePointer, GaussPointIterator, DistanceIterator > BucketType;
        typedef Tree< KDTreePartition<BucketType> > tree;

        typedef boost::numeric::ublas::bounded_matrix<double,3,3> MatrixVar;

    public:
        // Class Constructor
        // WARNING: Input ModelParts are expected to contain interface nodes and conditions ONLY
        // Use an InterfacePreprocess object to create such a model part from a regular one:
        // InterfaceMapper = InterfacePreprocess()
        // InterfacePart = InterfaceMapper.GenerateInterfacePart(Complete_Model_Part)
        AdvancedNMPointsMapper(const ModelPart & rOriginModelPart,ModelPart & rDestinationModelPart);

        // Class Destructor
        //~AdvancedNMPointsMapper();

        void ScalarMap(const Variable<double> & rOriginVar, Variable<double> & rDestVar,
                const int MaxIter, const double TolIter); // Scalar Version
        void VectorMap(const Variable< array_1d<double,3> > & rOriginVar,Variable< array_1d<double,3> > & rDestVar,
                const int MaxIter,const double TolIter); // Vector Version

        void FindNeighbours(double SearchRadiusFactor);

    private:
        void CalcNormalAndArea(Condition::Pointer pCond, array_1d<double,3>& Normal,
                double& Area);
        void TriangleCenterAndRadius(const Condition::Pointer pCond,
                Point<3>& Center, double& Radius);
        void SetProjectionToCond(GaussPointItem& GaussPoint, Condition::Pointer pCandidateCond); // Desired outcome

        void SetProjectionToNode(GaussPointItem& GaussPoint,
                Node<3>::Pointer pCandidateNode, double& Dist); // Alternative when no condition is available

        const ModelPart& mrOriginModelPart;
        ModelPart& mrDestinationModelPart;
        unsigned int mBucketSize; //Bucket size for kd-tree
        GaussPointVector mGaussPointList;

        void DistanceCheck(); // TEST FUNCTION
    };
}

#endif