/*
 * File:   AdvancedNMPointsMapper.cpp
 * Author: jcotela
 * Co-author: vmataix, rzorrilla
 *
 * Created on 19 January 2010, 10:20
 * Last update on 31 August 2016, 10:28
 */

// System includes

// External includes

// Project includes
#include "AdvancedNMPointsMapper.hpp"

// TODO: Clean code, use internal functions from Kratos, as Geometry.Area()....

namespace Kratos
{
// GaussPointItem Methods
void GaussPointItem::Project(
        Condition::Pointer pOriginCond,
        array_1d<double,2>& Coords,
        double& Dist,
        const int dimension
        )
{
    
    if (dimension == 2)
    {
        Point<3> point_to_project;
        point_to_project.Coordinate(1) = this->Coordinate(1);
        point_to_project.Coordinate(2) = this->Coordinate(2);
        point_to_project.Coordinate(3) = this->Coordinate(3);
        
        Point<3> point_projected;
        ProjectPointToPlane(pOriginCond->GetGeometry()[0],point_to_project,point_projected,Dist,mNormal);
                
        array_1d<double, 3> point_projected_local_coor;
        point_projected_local_coor = pOriginCond->GetGeometry().PointLocalCoordinates(point_projected_local_coor,point_projected);
        
        Coords[0] = point_projected_local_coor[0];
        Coords[1] = 0.0;
    }
    else
    {
        // xi,yi,zi are Nodal Coordinates, n is the destination condition's unit normal
        // and d is the distance along n from the point to its projection in the condition
        // | DestX-x3 |   | x1-x3 x2-x3 nx |   | Chi |
        // | DestY-y3 | = | y1-y3 y2-y3 ny | . | Eta |
        // | DestZ-z3 |   | z1-z3 z2-z3 nz |   |  d  |

        Matrix ChangeMatrix(3, 3, false);
        Matrix InvChange(3, 3, false);
        double det;

        array_1d<double, 3> RHS, Res;

        RHS[0] = this->Coordinate(1) - pOriginCond->GetGeometry()[2].X();
        RHS[1] = this->Coordinate(2) - pOriginCond->GetGeometry()[2].Y();
        RHS[2] = this->Coordinate(3) - pOriginCond->GetGeometry()[2].Z();

        ChangeMatrix(0, 0) = pOriginCond->GetGeometry()[0].X() - pOriginCond->GetGeometry()[2].X();
        ChangeMatrix(1, 0) = pOriginCond->GetGeometry()[0].Y() - pOriginCond->GetGeometry()[2].Y();
        ChangeMatrix(2, 0) = pOriginCond->GetGeometry()[0].Z() - pOriginCond->GetGeometry()[2].Z();

        ChangeMatrix(0, 1) = pOriginCond->GetGeometry()[1].X() - pOriginCond->GetGeometry()[2].X();
        ChangeMatrix(1, 1) = pOriginCond->GetGeometry()[1].Y() - pOriginCond->GetGeometry()[2].Y();
        ChangeMatrix(2, 1) = pOriginCond->GetGeometry()[1].Z() - pOriginCond->GetGeometry()[2].Z();

        ChangeMatrix(0, 2) = mNormal[0];
        ChangeMatrix(1, 2) = mNormal[1];
        ChangeMatrix(2, 2) = mNormal[2];

        MathUtils<double>::InvertMatrix3(ChangeMatrix,InvChange,det);
        noalias(Res) = prod(InvChange, RHS);

        Coords[0] = Res[0];
        Coords[1] = Res[1];
        // Keep distance positive, regardless of normal orientation
        Dist = (Res[2] < 0)? -Res[2] : Res[2] ;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void GaussPointItem::ProjectPointToPlane(
        const Point<3>& PointInPlane,
        const Point<3>& PointToBeProjected,
        Point<3>& PointProjected,
        double& dist,
        const array_1d<double,3>& Normal
        )
    {
         array_1d<double,3> vector_points;
         noalias(vector_points) = PointToBeProjected.Coordinates() - PointInPlane.Coordinates();

         dist = inner_prod(vector_points, Normal);

         PointProjected.Coordinates() = PointToBeProjected.Coordinates() - Normal * dist;
    }

/***********************************************************************************/
/***********************************************************************************/

void GaussPointItem::GetProjectedValue(
        const Variable<double> & rOriginVar,
        double& Value,
        const int dimension
        )
{
    if (mProjStatus == 1) // Get Interpolated value from origin condition
    {
        if (dimension == 2)
        {
            Point<3> GPloccoords;
            GPloccoords.Coordinate(1) = mOriginCoords[0];
            GPloccoords.Coordinate(2) = 0.0;
            GPloccoords.Coordinate(3) = 0.0;
            
            Vector shfunc_values;
            mpOriginCond.lock()->GetGeometry().ShapeFunctionsValues(shfunc_values,GPloccoords);

            Value =  (mpOriginCond.lock())->GetGeometry()[0].FastGetSolutionStepValue(rOriginVar) * shfunc_values[0]
                   + (mpOriginCond.lock())->GetGeometry()[1].FastGetSolutionStepValue(rOriginVar) * shfunc_values[1];
        }
        else
        {
            Value = ( (mpOriginCond.lock())->GetGeometry()[0].FastGetSolutionStepValue(rOriginVar) * mOriginCoords[0]
                    + (mpOriginCond.lock())->GetGeometry()[1].FastGetSolutionStepValue(rOriginVar) * mOriginCoords[1]
                    + (mpOriginCond.lock())->GetGeometry()[2].FastGetSolutionStepValue(rOriginVar) * (1.0 - mOriginCoords[0] - mOriginCoords[1]) );
        }
    }
    else if (mProjStatus == 2)   // Get Value from origin node
    {
        Value = (mpOriginNode.lock())->FastGetSolutionStepValue(rOriginVar);
    }
    else   // mProjStatus == 0: Return 0
    {
        Value = 0.0;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void GaussPointItem::GetProjectedValue(
        const Variable<array_1d<double,3> >& rOriginVar,
        array_1d<double,3>& Value,
        const int dimension
        )
{

    if (mProjStatus == 1) // Get Interpolated value from origin condition
    {
        if (dimension == 2)
        {
            for (unsigned int i = 0; i < 2; i++)
            {
                Point<3> GPloccoords;
                GPloccoords.Coordinate(1) = mOriginCoords[0];
                GPloccoords.Coordinate(2) = 0.0;
                GPloccoords.Coordinate(3) = 0.0;

                Vector shfunc_values;
                mpOriginCond.lock()->GetGeometry().ShapeFunctionsValues(shfunc_values,GPloccoords);

                Value[i] = (mpOriginCond.lock())->GetGeometry()[0].FastGetSolutionStepValue(rOriginVar)[i] * shfunc_values[0]
                         + (mpOriginCond.lock())->GetGeometry()[1].FastGetSolutionStepValue(rOriginVar)[i] * shfunc_values[1];
            }
            
            Value[2] = 0.0;

        }
        else
        {
            for (unsigned int i = 0; i < 3; i++)
            {
                Value[i] = ((mpOriginCond.lock())->GetGeometry()[0].FastGetSolutionStepValue(rOriginVar)[i] * mOriginCoords[0]
                          + (mpOriginCond.lock())->GetGeometry()[1].FastGetSolutionStepValue(rOriginVar)[i] * mOriginCoords[1]
                          + (mpOriginCond.lock())->GetGeometry()[2].FastGetSolutionStepValue(rOriginVar)[i] * (1.0 - mOriginCoords[0] - mOriginCoords[1]) );
            }
        }
    }
    else if (mProjStatus == 2)   // Get Value from origin node
    {
        Value = (mpOriginNode.lock())->FastGetSolutionStepValue(rOriginVar);
    }
    else   // mProjStatus == 0: Return 0
    {
        Value = ZeroVector(3);
    }
}

/***********************************************************************************/
/***********************************************************************************/

// Mapper Methods
AdvancedNMPointsMapper::AdvancedNMPointsMapper(
        const ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart
        ):

    mrOriginModelPart(rOriginModelPart),
    mrDestinationModelPart(rDestinationModelPart),
    mBucketSize(4)
{
    const unsigned int dimension = mrDestinationModelPart.ConditionsBegin()->GetGeometry().WorkingSpaceDimension();

    array_1d<double, 3> GPCoord; // Will store the coordinates of a condition's Gauss Points
    GPCoord = ZeroVector(3);
    array_1d<double, 3> Normal;

    if (dimension == 2) // 2D case
    {
        boost::numeric::ublas::bounded_matrix<double,2,3> Nodes;
        boost::numeric::ublas::bounded_matrix<double,2,3> GaussPoints;
        boost::numeric::ublas::bounded_matrix<double,2,2> GPPos;
        double Length;
        
        // 2 Gauss-Legendre point quadrature
        // eps = +-(1/3)^0.5
        // N1 = (1-eps)/2
        // N2 = (1+eps)/2
        
        GPPos(0, 0) = (1.0+(-std::sqrt(1.0/3.0)))/2.0;
        GPPos(0, 1) = (1.0-(-std::sqrt(1.0/3.0)))/2.0;
        GPPos(1, 0) = (1.0+(+std::sqrt(1.0/3.0)))/2.0;
        GPPos(1, 1) = (1.0-(+std::sqrt(1.0/3.0)))/2.0;

        for (ModelPart::ConditionsContainerType::iterator cond_it = rDestinationModelPart.ConditionsBegin();
            cond_it != rDestinationModelPart.ConditionsEnd();
            cond_it++)
        {
            CalcNormalAndArea(*cond_it.base(), Normal, Length,  dimension);

            for(unsigned int i = 0; i < 2; i++)
            {
                for(unsigned int j = 0; j < 3; j++)
                {
                    Nodes(i,j) = cond_it->GetGeometry()[i].Coordinate(j + 1);
                }
            }

            noalias(GaussPoints) = prod(GPPos, Nodes);
            
            for(unsigned int k = 0; k < 2; k++)
            {
                for(unsigned int l = 0; l < 3; l++)
                {
                    GPCoord[l] = GaussPoints(k,l);
                }
                GaussPointItem::Pointer pGP = GaussPointItem::Pointer(new GaussPointItem(GPCoord, Length, Normal));
                mGaussPointList.push_back( pGP );
            }
        }
    }
    else // 3D case
    {
        // The constructor defines the Gauss Points in the destination interface, stores them as
        // GaussPointItem instances and builds a std::vector of pointers to their position in memory
        // The Gauss Point Coordinates (Gix,Giy,Gyz) are obtained from the nodes' (A,B,C) coordinates
        // as follows:
        // | G1x G1y G1z |   | 0.6 0.2 0.2 |   | Ax Ay Az |
        // | G2x G2y G2z | = | 0.2 0.6 0.2 | . | Bx By Bz |
        // | G3x G3y G3z |   | 0.2 0.2 0.6 |   | Cx Cy Cz |
        
        double Area;
        MatrixVar Nodes, GaussPoints, GPPos;
        GPPos(0,0) = 0.6;
        GPPos(0,1) = 0.2;
        GPPos(0,2) = 0.2;
        GPPos(1,0) = 0.2;
        GPPos(1,1) = 0.6;
        GPPos(1,2) = 0.2;
        GPPos(2,0) = 0.2;
        GPPos(2,1) = 0.2;
        GPPos(2,2) = 0.6;

        for (
            ModelPart::ConditionsContainerType::iterator cond_it = rDestinationModelPart.ConditionsBegin();
            cond_it != rDestinationModelPart.ConditionsEnd();
            cond_it++)
        {
            CalcNormalAndArea(*cond_it.base(), Normal, Area, dimension);

            for(unsigned int i = 0; i < 3; i++)
            {
                for(unsigned int j = 0; j < 3; j++)
                {
                    Nodes(i,j) = cond_it->GetGeometry()[i].Coordinate(j + 1);
                }
            }

            noalias(GaussPoints) = prod(GPPos,Nodes);

            for(unsigned int k = 0; k < 3; k++)
            {
                for(unsigned int l = 0; l < 3; l++)
                {
                    GPCoord[l] = GaussPoints(k,l);
                }
                GaussPointItem::Pointer pGP = GaussPointItem::Pointer(new GaussPointItem(GPCoord,Area,Normal) );
                mGaussPointList.push_back( pGP );
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AdvancedNMPointsMapper::FindNeighbours(double SearchRadiusFactor)
{
    // Initialize some values
    const unsigned int MaxResults = 5000; // Maximum number of points to find in a single tree search
    GaussPointVector Results(MaxResults);
    std::vector<double> ResultDistances(MaxResults);
    array_1d<double,3> ZeroVect = ZeroVector(3);

    const unsigned int dimension = mrDestinationModelPart.ConditionsBegin()->GetGeometry().WorkingSpaceDimension();

    // Create a tree
    // It will use a copy of mGaussPoinList (a std::vector which contains pointers)
    // Copying the list is required because the tree will reorder it for efficiency
    GaussPointVector pGaussPoints = mGaussPointList;
    tree Tree_conds(pGaussPoints.begin(), pGaussPoints.end(), mBucketSize);

    double Radius, SearchRadius;
    double MaxRadius = 0.0;

    for(ModelPart::ConditionsContainerType::const_iterator cond_it = mrOriginModelPart.ConditionsBegin();
            cond_it != mrOriginModelPart.ConditionsEnd();
            cond_it++)
    {
        Point<3> Center;
        if (dimension == 2) // 2D case
        {
            LineCenterAndRadius(*cond_it.base(), Center, Radius);
        }
        else // 3D case
        {
            TriangleCenterAndRadius(*cond_it.base(), Center, Radius);
        }

        // Create a fake GaussPoint with the same coordinates
        // (to ensure that the object used in the kd-tree search is of the same type as the tree contents)
        GaussPointItem CenterGP(Center.Coordinates(), 0, ZeroVect);

        if (Radius > MaxRadius)
        {
            MaxRadius = Radius;
        }

        // Count Gauss Points within SearchRadius
        SearchRadius = SearchRadiusFactor * Radius;
        unsigned int Found = Tree_conds.SearchInRadius(CenterGP, SearchRadius, Results.begin(), ResultDistances.begin(), MaxResults);
    
        for(unsigned int i = 0; i < Found; i++)
        {
            SetProjectionToCond(*Results[i], *cond_it.base(), dimension);
        }
    }

    // Try to find reasonable origin values in nearby nodes for Gauss Points that couldn't be projected to a condition
    GaussPointVector ProjectionlessGP;

    for (GaussPointIterator gauss_it = mGaussPointList.begin(); gauss_it != mGaussPointList.end(); gauss_it++)
    {
        int Status = 0;
        (*gauss_it)->GetProjStatus(Status);

        if( Status != 1)
        {
            ProjectionlessGP.push_back( *gauss_it );
        }
    }
    
    // Try to use a nearby node as reference for points that couldn't be projected to a condition
    if (ProjectionlessGP.size() != 0)
    {
        std::cout << "AdvancedNMPointsMapper: " << ProjectionlessGP.size()
                  << " Gauss points could not be projected to a condition.\n"
                  << "    Attempring to get a reasonable value for them from a node..." << std::endl;

        GaussPointVector NodeResults(MaxResults);
        std::vector<double> NodeResultDistances(MaxResults);

        tree kdtree_nodes(ProjectionlessGP.begin(), ProjectionlessGP.end(), mBucketSize);
        SearchRadius = SearchRadiusFactor * MaxRadius;

        for (ModelPart::NodesContainerType::const_iterator node_it = mrOriginModelPart.NodesBegin();
                node_it != mrOriginModelPart.NodesEnd();
                node_it++)
        {
            // Our tree uses GaussPointItem objects as input to sort
            GaussPointItem NodePos(node_it->Coordinates(), 0, ZeroVect);
            unsigned int Found = kdtree_nodes.SearchInRadius(NodePos, SearchRadius, NodeResults.begin(), NodeResultDistances.begin(), MaxResults);

            for (unsigned int i = 0; i < Found; i++)
            {
                SetProjectionToNode(*NodeResults[i],*node_it.base(), NodeResultDistances[i]);
            }
        }
        //DistanceCheck(); // Test function

        // Count how many points without a projection remain
        unsigned int counter = 0;
        for (GaussPointIterator gauss_it = mGaussPointList.begin(); gauss_it != mGaussPointList.end(); gauss_it++)
        {
            int Status = 0;
            (*gauss_it)->GetProjStatus(Status);

            if(Status == 0)
            {
                counter++;
            }
        }

        std::cout << "   ... " << counter << " Gauss Points without a reference value remain." << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AdvancedNMPointsMapper::CalcNormalAndArea(
        const Condition::Pointer Cond,
        array_1d<double,3> & Normal,
        double & Area,
        const int dimension
        )
{
    array_1d<double,3> v1,v2;

    v1[0] = Cond->GetGeometry()[1].X() - Cond->GetGeometry()[0].X();
    v1[1] = Cond->GetGeometry()[1].Y() - Cond->GetGeometry()[0].Y();
    v1[2] = Cond->GetGeometry()[1].Z() - Cond->GetGeometry()[0].Z();

    if (dimension == 3)
    {
        v2[0] = Cond->GetGeometry()[2].X() - Cond->GetGeometry()[0].X();
        v2[1] = Cond->GetGeometry()[2].Y() - Cond->GetGeometry()[0].Y();
        v2[2] = Cond->GetGeometry()[2].Z() - Cond->GetGeometry()[0].Z();
    }
    else // Assuming plane X-Y
    {
        v2[0] = 0.0;
        v2[1] = 0.0;
        v2[2] = 1.0;
    }

    // Compute the condition normal
    MathUtils<double>::CrossProduct(Normal,v1,v2);

    double NNorm = std::sqrt(Normal[0]*Normal[0] + Normal[1]*Normal[1] + Normal[2]*Normal[2]);

    Normal /= NNorm;
    
    if (dimension ==2)
    {
        Area = std::sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]); // Note that in 2D, the Area variable represents the condition length
    }
    else
    {
        Area = 0.5 * NNorm;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AdvancedNMPointsMapper::LineCenterAndRadius(
        const Condition::Pointer Cond,
        Point<3>& Center,
        double& Radius
        )
{
    Radius = 0.0;
    Center.Coordinate(1) = 0.5 * (Cond->GetGeometry()[0].X() + Cond->GetGeometry()[1].X());
    Center.Coordinate(2) = 0.5 * (Cond->GetGeometry()[0].Y() + Cond->GetGeometry()[1].Y());
    Center.Coordinate(3) = 0.5 * (Cond->GetGeometry()[0].Z() + Cond->GetGeometry()[1].Z());

    for(int i = 0; i < 2; i++)
    {
        double dx = Center.Coordinate(1) - Cond->GetGeometry()[i].X();
        double dy = Center.Coordinate(2) - Cond->GetGeometry()[i].Y();
        double dz = Center.Coordinate(3) - Cond->GetGeometry()[i].Z();

        double tmp = dx*dx + dy*dy + dz*dz;

        if(tmp > Radius)
        {
            Radius = tmp;
        }
    }

    Radius = std::sqrt(Radius);
}

/***********************************************************************************/
/***********************************************************************************/

void AdvancedNMPointsMapper::TriangleCenterAndRadius(
        const Condition::Pointer Cond,
        Point<3>& Center,
        double& Radius
        )
{
    Radius = 0.0;
    Center.Coordinate(1) = 0.33333333333333333333*(Cond->GetGeometry()[0].X() + Cond->GetGeometry()[1].X() + Cond->GetGeometry()[2].X());
    Center.Coordinate(2) = 0.33333333333333333333*(Cond->GetGeometry()[0].Y() + Cond->GetGeometry()[1].Y() + Cond->GetGeometry()[2].Y());
    Center.Coordinate(3) = 0.33333333333333333333*(Cond->GetGeometry()[0].Z() + Cond->GetGeometry()[1].Z() + Cond->GetGeometry()[2].Z());

    for(unsigned int i = 0; i < 3; i++)
    {
        double dx = Center.Coordinate(1) - Cond->GetGeometry()[i].X();
        double dy = Center.Coordinate(2) - Cond->GetGeometry()[i].Y();
        double dz = Center.Coordinate(3) - Cond->GetGeometry()[i].Z();

        double tmp = dx*dx + dy*dy + dz*dz;

        if(tmp > Radius)
        {
            Radius = tmp;
        }
    }
    Radius = std::sqrt(Radius);
}

/***********************************************************************************/
/***********************************************************************************/

void AdvancedNMPointsMapper::SetProjectionToCond(
        GaussPointItem& GaussPoint,
        Condition::Pointer pCandidateCond,
        const int dimension
        )
{    
    array_1d<double,2> LocalCoords;
    double Dist;
    GaussPoint.Project(pCandidateCond, LocalCoords, Dist, dimension);

    if (dimension == 2)
    {
        // Point belongs to Condition (-1 <= chi <= 1)
        if (LocalCoords[0] >= -1.0 && LocalCoords[0] <= 1.0)
        {
            int Status = 0;
            GaussPoint.GetProjStatus(Status);

            if ( Status != 1) // No good projection found (yet)
            {
                GaussPoint.SetProjection(pCandidateCond, LocalCoords, Dist);
            }
            else // A good projection was found in a previous iteration
            {
                double CurrentDist = 0.0;
                GaussPoint.GetDist(CurrentDist);

                if (Dist < CurrentDist)
                {
                    GaussPoint.SetProjection(pCandidateCond, LocalCoords, Dist);
                }
            }
        }
    }
    else
    {
        // Point belongs to Condition
        if (LocalCoords[0] >= 0.0 && LocalCoords[1] >= 0.0 && (1.0 - LocalCoords[0] - LocalCoords[1]) >= 0.0)
        {
            int Status = 0;
            GaussPoint.GetProjStatus(Status);

            if ( Status != 1) // No good projection found (yet)
            {
                GaussPoint.SetProjection(pCandidateCond ,LocalCoords, Dist);
            }
            else // A good projection was found in a previous iteration
            {
                double CurrentDist = 0.0;
                GaussPoint.GetDist(CurrentDist);

                if (Dist < CurrentDist)
                {
                    GaussPoint.SetProjection(pCandidateCond, LocalCoords, Dist);
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AdvancedNMPointsMapper::SetProjectionToNode(
        GaussPointItem& GaussPoint,
        Node<3>::Pointer pCandidateNode,
        const double& Dist
        )
{
    int Status = 0;
    GaussPoint.GetProjStatus(Status);

    if ( Status == 0)
    {
        GaussPoint.SetProjection(pCandidateNode, Dist);
    }
    else
    {
        double CurrentDist = 0.0;
        GaussPoint.GetDist(CurrentDist);

        if (Dist < CurrentDist)
        {
            GaussPoint.SetProjection(pCandidateNode, Dist);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AdvancedNMPointsMapper::ScalarToNormalVectorMap( 
        const Variable<double> & rOriginVar,
        Variable<array_1d<double,3> >& rDestVar,
        const int MaxIter,
        const double TolIter,
        const bool sign_pos
        )
{
    array_1d<double,3> ZeroVect = ZeroVector(3);
    
    // Define if the mapping swaps the sign of the variable values
    double sign = 1.0;
    if (sign_pos == false)
    {
        sign = -1.0;
    }

    // Initialize results and NODAL_MAUX
    for ( ModelPart::NodesContainerType::iterator node_it = mrDestinationModelPart.NodesBegin();
            node_it != mrDestinationModelPart.NodesEnd();
            node_it++)
    {
        node_it->GetValue(NODAL_MAUX) = 0.0;
        node_it->FastGetSolutionStepValue(rDestVar) = ZeroVect;
    }
    
    // Compute nodal lengths/areas in both origin and destination modelparts
    ComputeNodalLengthArea();
    
    const unsigned int dimension = mrDestinationModelPart.ConditionsBegin()->GetGeometry().WorkingSpaceDimension();
    
    if (dimension == 2) // 2D case
    {
        // Interpolation matrix obtention
        boost::numeric::ublas::bounded_matrix<double,2,2> MCons; // Elemental Consistent Mass Matrix = L/2 * MCons
        boost::numeric::ublas::bounded_matrix<double,2,2> MInterp;
        boost::numeric::ublas::bounded_matrix<double,2,2> aux_GPPos;
        boost::numeric::ublas::bounded_matrix<double,2,2> inv_aux_GPPos;

        MCons(0, 0) = 2.0/3.0; 
        MCons(0, 1) = 1.0/3.0;
        MCons(1, 0) = 1.0/3.0;
        MCons(1, 1) = 2.0/3.0;

        aux_GPPos(0, 0) = (1.0+(-std::sqrt(1.0/3.0)))/2.0;
        aux_GPPos(0, 1) = (1.0-(-std::sqrt(1.0/3.0)))/2.0;
        aux_GPPos(1, 0) = (1.0+(+std::sqrt(1.0/3.0)))/2.0;
        aux_GPPos(1, 1) = (1.0-(+std::sqrt(1.0/3.0)))/2.0;

        double det_aux_GPPos;
        det_aux_GPPos = (aux_GPPos(0,0)*aux_GPPos(1,1))-(aux_GPPos(0,1)*aux_GPPos(1,0));

        inv_aux_GPPos(0, 0) = aux_GPPos(1, 1)/det_aux_GPPos;
        inv_aux_GPPos(0, 1) = -aux_GPPos(0, 1)/det_aux_GPPos;
        inv_aux_GPPos(1, 0) = -aux_GPPos(1, 0)/det_aux_GPPos;
        inv_aux_GPPos(1, 1) = aux_GPPos(0, 0)/det_aux_GPPos;

        MInterp = prod(MCons, inv_aux_GPPos); // Interpolation matrix (NodalValues = (L/6)*MInterp*GaussValues)

        std::vector< boost::shared_ptr<array_1d<double,2> > > pInterpValues;

        // Here we loop both the Destination Model Part and the Gauss Point List, using the
        // fact that the GP list was created by a loop over the Dest. Model Part's conditions
        int GPi = 0;
        for (ModelPart::ConditionIterator cond_it = mrDestinationModelPart.ConditionsBegin();
                cond_it != mrDestinationModelPart.ConditionsEnd();
                cond_it++)
        {
            array_1d<double, 2> GPValues;
            array_1d<double, 2> NodalValues;

            // Get the length of the parent condition of the two GP considered
            double CondLength = 0.0;
            mGaussPointList[GPi]->GetArea(CondLength); 
            
            // Get the condition normal
            array_1d<double,3> NormalVector = ZeroVect;
            mGaussPointList[GPi]->GetNormal(NormalVector);

            // Store the condition normal in each node to compute the nodal normal
            cond_it->GetGeometry()[0].GetValue(NORMAL) += NormalVector; 
            cond_it->GetGeometry()[1].GetValue(NORMAL) += NormalVector; 
            
            // Store in GPValues the projected value in the destiny condition Gauss Points
            // Note that currently the implementation is valid for only 2 GP
            double TempValue = 0.0;
            
            mGaussPointList[GPi]->GetProjectedValue(rOriginVar, TempValue, dimension);
            GPValues[0] = TempValue;

            mGaussPointList[GPi + 1]->GetProjectedValue(rOriginVar, TempValue, dimension);
            GPValues[1] = TempValue;

            const double K = CondLength/2.0;
            
            NodalValues[0] = K*(MInterp(0,0)*GPValues[0] + MInterp(0,1)*GPValues[1]);
            NodalValues[1] = K*(MInterp(1,0)*GPValues[0] + MInterp(1,1)*GPValues[1]);

            boost::shared_ptr< array_1d<double,2> > pNodalValues(new array_1d<double,2>(NodalValues) );
            pInterpValues.push_back(pNodalValues); // This is computed here because it is the constant part of RHS

            GPi += 2; // 1 Condition = 2 Gauss Points
        }

        // Compute the nodal normals (normalised sum of the conditions normals) in the destination modelpart
        for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
                node_it != mrDestinationModelPart.NodesEnd();
                node_it++)
        {
            const array_1d<double,3> NormalVector = node_it->GetValue(NORMAL);
            node_it->GetValue(NORMAL) = NormalVector/norm_2(NormalVector);
        }

        // Iteration
        for (int k = 0; k < MaxIter; k++)
        {
            // At the begining of each iteration initialize the variable containing the assembled RHS as 0
            for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
                    node_it != mrDestinationModelPart.NodesEnd();
                    node_it++)
            {
                node_it->GetValue(MAPPER_SCALAR_PROJECTION_RHS) = 0.0;
            }

            double LocalRHS0, LocalRHS1; // Local RHS for each node
            array_1d<double, 2> LastSolution;
            int IV_iter = 0;

            for (ModelPart::ConditionIterator cond_it = mrDestinationModelPart.ConditionsBegin();
                    cond_it != mrDestinationModelPart.ConditionsEnd();
                    cond_it++)
            {
                double CondLength = 0.0;
                mGaussPointList[2 * IV_iter]->GetArea(CondLength);
                
                LocalRHS0 = 0.0;
                LocalRHS1 = 0.0;
                
                const array_1d<double,3> NormalVector0 = cond_it->GetGeometry()[0].GetValue(NORMAL);
                const array_1d<double,3> NormalVector1 = cond_it->GetGeometry()[1].GetValue(NORMAL);

                LastSolution[0] = sign * inner_prod(cond_it->GetGeometry()[0].FastGetSolutionStepValue(rDestVar), NormalVector0);
                LastSolution[1] = sign * inner_prod(cond_it->GetGeometry()[1].FastGetSolutionStepValue(rDestVar), NormalVector1);             

                const double K = CondLength/6.0;
                array_1d<double,2> CondValues = *pInterpValues[IV_iter];
                
                LocalRHS0 = CondValues[0] - K * (2.0*LastSolution[0] + 1.0*LastSolution[1]);
                LocalRHS1 = CondValues[1] - K * (1.0*LastSolution[0] + 2.0*LastSolution[1]);
                // We are taking advantage of 1 Condition = 2 Gauss Points to iterate the interpolation results and the Gauss Point Vector Again

                cond_it->GetGeometry()[0].GetValue(MAPPER_SCALAR_PROJECTION_RHS) += LocalRHS0;
                cond_it->GetGeometry()[1].GetValue(MAPPER_SCALAR_PROJECTION_RHS) += LocalRHS1;

                IV_iter++;
            }

            // Solve
            double dValNorm      = 0.0;
            double ValNorm       = 0.0;
            double RelativeError = 0.0;
            unsigned int NodeNum = 0;

            for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
                    node_it != mrDestinationModelPart.NodesEnd();
                    node_it++)
            {
                const array_1d<double,3> NormalVector = node_it->GetValue(NORMAL);
                const double & NodeLength = node_it->GetValue(NODAL_MAUX);
                //node_it->FastGetSolutionStepValue(NODAL_MAUX) = NodeArea; // TEST: store nodal area so GiD can paint it later
                double dVal = node_it->GetValue(MAPPER_SCALAR_PROJECTION_RHS)/NodeLength;
                node_it->FastGetSolutionStepValue(rDestVar) += sign * dVal * NormalVector;

                // Variables for convergence check
                dValNorm += dVal * dVal;
                ValNorm  += std::pow(inner_prod(node_it->FastGetSolutionStepValue(rDestVar), NormalVector), 2);

                NodeNum++;
            }
            //std::cout << "ScalarToNormalVectorMap iteration: " << k+1 << std::endl;

            // Check Convergence
            if(ValNorm > 10e-15)
            {
                RelativeError = dValNorm / ValNorm;
            }
            if( (ValNorm/NodeNum < 0.00001 * std::pow(TolIter, 2.00)) || (RelativeError < std::pow(TolIter, 2.00)) )
            {
                std::cout << "ScalarToNormalVectorMap converged in " << k + 1 << " iterations." << std::endl;
                break;
            }
            else if ((k + 1) == MaxIter)
            {
                std::cout << "WARNING: ScalarToNormalVectorMap did not converge in " << k + 1 << " iterations." << std::endl;
            }
        }
    }
    else // 3D case
    {
        // Define some variables that will be used in the iteration
        MatrixVar MCons; // Elemental Consistent Mass Matrix = Aelem/12 * MCons
        MCons(0, 0) = 2.0;
        MCons(0, 1) = 1.0;
        MCons(0, 2) = 1.0;
        MCons(1, 0) = 1.0;
        MCons(1, 1) = 2.0;
        MCons(1, 2) = 1.0;
        MCons(2, 0) = 1.0;
        MCons(2, 1) = 1.0;
        MCons(2, 2) = 2.0;

        MatrixVar MInterp; // Interpolation Matrix (NodalValues = (A/24)*MInterp*GaussValues)
        MInterp(0, 0) = 6.0;
        MInterp(0, 1) = 1.0;
        MInterp(0, 2) = 1.0;
        MInterp(1, 0) = 1.0;
        MInterp(1, 1) = 6.0;
        MInterp(1, 2) = 1.0;
        MInterp(2, 0) = 1.0;
        MInterp(2, 1) = 1.0;
        MInterp(2, 2) = 6.0;

        std::vector< boost::shared_ptr< array_1d<double, 3> > > pInterpValues;

        // Here we loop both the Destination Model Part and the Gauss Point List, using the
        // fact that the GP list was created by a loop over the Dest. Model Part's conditions
        int GPi = 0;
        for (ModelPart::ConditionsContainerType::iterator cond_it = mrDestinationModelPart.ConditionsBegin();
                cond_it != mrDestinationModelPart.ConditionsEnd();
                cond_it++)
        {
            array_1d<double,3> GPValues;
            array_1d<double,3> NodalValues;

            double CondArea = 0.0;
            mGaussPointList[GPi]->GetArea(CondArea);

            array_1d<double,3> NormalVector = ZeroVect;
            mGaussPointList[GPi]->GetNormal(NormalVector);

            for (unsigned int i = 0; i < 3; i++)
            {
                //~ cond_it->GetGeometry()[i].GetValue(NODAL_MAUX) += 0.333333333333333333 * CondArea;
                mGaussPointList[GPi + i]->GetProjectedValue(rOriginVar, GPValues[i], 3);
                
                cond_it->GetGeometry()[i].GetValue(NORMAL) += NormalVector; // Store the condition normal in each node to compute the nodal normal
            }

            noalias(NodalValues) = (CondArea/24.0) * prod(MInterp, GPValues);

            boost::shared_ptr< array_1d<double, 3> > pNodalValues(new array_1d<double, 3>(NodalValues) );
            pInterpValues.push_back(pNodalValues); // This is computed here because it is the constant part of RHS

            GPi += 3; // 1 Condition = 3 Gauss Points
        }

        // Compute the nodal normal (normalised sum of the conditions normals)
        for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
                node_it != mrDestinationModelPart.NodesEnd();
                node_it++)
        {
            const array_1d<double,3> NormalVector = node_it->GetValue(NORMAL);
            node_it->GetValue(NORMAL) = NormalVector/norm_2(NormalVector);
        }

        // Iteration
        for (int k = 0; k < MaxIter; k++)
        {
            // At the begining of each iteration initialize the variable containing the assembled RHS as 0
            for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
                    node_it != mrDestinationModelPart.NodesEnd();
                    node_it++)
            {
                node_it->GetValue(MAPPER_SCALAR_PROJECTION_RHS) = 0.0;
            }

            array_1d<double, 3> LocalRHS;
            array_1d<double, 3> LastSolution;
            int IV_iter = 0;

            for (ModelPart::ConditionsContainerType::iterator cond_it = mrDestinationModelPart.ConditionsBegin();
                    cond_it != mrDestinationModelPart.ConditionsEnd();
                    cond_it++)
            {
                for(unsigned int j = 0; j < 3; j++)
                {
                    const array_1d<double,3> NormalVector = cond_it->GetGeometry()[j].GetValue(NORMAL);
                    LocalRHS[j] = 0.0;
                    LastSolution[j] = sign * inner_prod(cond_it->GetGeometry()[j].FastGetSolutionStepValue(rDestVar), NormalVector);
                }

                double CondArea = 0.0;
                mGaussPointList[3 * IV_iter]->GetArea(CondArea);

                noalias(LocalRHS) = *pInterpValues[IV_iter] - (CondArea/12.0) * prod(MCons, LastSolution);
                // We are taking advantage of 1 Condition = 3 Gauss Points to iterate the interpolation results and the Gauss Point Vector Again

                for (unsigned int j = 0; j < 3 ; j++)
                {
                    cond_it->GetGeometry()[j].GetValue(MAPPER_SCALAR_PROJECTION_RHS) += LocalRHS[j];
                }

                IV_iter++;
            }

            // Solve
            double dValNorm      = 0.0;
            double ValNorm       = 0.0;
            double RelativeError = 0.0;
            unsigned int NodeNum = 0;

            for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
                    node_it != mrDestinationModelPart.NodesEnd();
                    node_it++)
            {
                const array_1d<double,3> NormalVector = node_it->GetValue(NORMAL);
                const double & NodeArea = node_it->GetValue(NODAL_MAUX);
                //node_it->FastGetSolutionStepValue(NODAL_MAUX) = NodeArea; // TEST: store nodal area so GiD can paint it later
                double dVal = node_it->GetValue(MAPPER_SCALAR_PROJECTION_RHS)/NodeArea;
                node_it->FastGetSolutionStepValue(rDestVar) += sign * dVal * NormalVector;

                // Variables for convergence check
                dValNorm += dVal * dVal;
                ValNorm  += std::pow(inner_prod(node_it->FastGetSolutionStepValue(rDestVar), NormalVector), 2);

                NodeNum++;
            }
            //std::cout << "ScalarToNormalVectorMap iteration: " << k+1 << std::endl;

            // Check Convergence
            if(ValNorm > 10e-15)
            {
                RelativeError = dValNorm / ValNorm;
            }
            if( (ValNorm/NodeNum < 0.00001 * std::pow(TolIter, 2.00)) || (RelativeError < std::pow(TolIter, 2.00)) )
            {
                std::cout << "ScalarToNormalVectorMap converged in " << k + 1 << " iterations." << std::endl;
                break;
            }
            else if ((k + 1) == MaxIter)
            {
                std::cout << "WARNING: ScalarToNormalVectorMap did not converge in " << k + 1 << " iterations." << std::endl;
            }
        } // End of Iteration
    }
} // End of Map (scalar to normal vector version)

/***********************************************************************************/
/***********************************************************************************/

void AdvancedNMPointsMapper::NormalVectorToScalarMap( // Note: JUST 3D!!!!
        const Variable<array_1d<double,3> >& rOriginVar,
        Variable<double> & rDestVar,
        const int MaxIter,
        const double TolIter,
        const bool sign_pos
        )
{
    double sign = 1.0;
    if (sign_pos == false)
    {
        sign = -1.0;
    }

    array_1d<double,3> ZeroVect = ZeroVector(3);

    // Initialize results and NODAL_MAUX
    for ( ModelPart::NodesContainerType::iterator node_it = mrDestinationModelPart.NodesBegin();
            node_it != mrDestinationModelPart.NodesEnd();
            node_it++)
    {
        node_it->GetValue(NODAL_MAUX) = 0.0;
        node_it->FastGetSolutionStepValue(rDestVar) = 0.0;
    }
    
    // Compute nodal lengths/areas in both origin and destination modelparts
    ComputeNodalLengthArea();
    
    const unsigned int dimension = mrDestinationModelPart.ConditionsBegin()->GetGeometry().WorkingSpaceDimension();
    
    if (dimension == 2)
    {
        // Interpolation matrix obtention
        boost::numeric::ublas::bounded_matrix<double,2,2> MCons; // Elemental Consistent Mass Matrix = L/2 * MCons
        boost::numeric::ublas::bounded_matrix<double,2,2> MInterp;
        boost::numeric::ublas::bounded_matrix<double,2,2> aux_GPPos;
        boost::numeric::ublas::bounded_matrix<double,2,2> inv_aux_GPPos;

        MCons(0, 0) = 2.0/3.0; 
        MCons(0, 1) = 1.0/3.0;
        MCons(1, 0) = 1.0/3.0;
        MCons(1, 1) = 2.0/3.0;

        aux_GPPos(0, 0) = (1.0+(-std::sqrt(1.0/3.0)))/2.0;
        aux_GPPos(0, 1) = (1.0-(-std::sqrt(1.0/3.0)))/2.0;
        aux_GPPos(1, 0) = (1.0+(+std::sqrt(1.0/3.0)))/2.0;
        aux_GPPos(1, 1) = (1.0-(+std::sqrt(1.0/3.0)))/2.0;

        double det_aux_GPPos;
        det_aux_GPPos = (aux_GPPos(0,0)*aux_GPPos(1,1))-(aux_GPPos(0,1)*aux_GPPos(1,0));

        inv_aux_GPPos(0, 0) = aux_GPPos(1, 1)/det_aux_GPPos;
        inv_aux_GPPos(0, 1) = -aux_GPPos(0, 1)/det_aux_GPPos;
        inv_aux_GPPos(1, 0) = -aux_GPPos(1, 0)/det_aux_GPPos;
        inv_aux_GPPos(1, 1) = aux_GPPos(0, 0)/det_aux_GPPos;

        MInterp = prod(MCons, inv_aux_GPPos); // Interpolation matrix (NodalValues = (L/6)*MInterp*GaussValues)

        std::vector< boost::shared_ptr<array_1d<double,6> > > pInterpValues;

        // Here we loop both the Destination Model Part and the Gauss Point List, using the
        // fact that the GP list was created by a loop over the Dest. Model Part's conditions
        int GPi = 0;
        for (ModelPart::ConditionIterator cond_it = mrDestinationModelPart.ConditionsBegin();
                cond_it != mrDestinationModelPart.ConditionsEnd();
                cond_it++)
        {
            array_1d<double, 6> GPValues;
            array_1d<double, 6> NodalValues;

            // Get the length of the parent condition of the two GP considered
            double CondLength = 0.0;
            mGaussPointList[GPi]->GetArea(CondLength); 
            
            // Get the condition normal
            array_1d<double,3> NormalVector = ZeroVect;
            mGaussPointList[GPi]->GetNormal(NormalVector);

            // Store the condition normal in each node to compute the nodal normal
            cond_it->GetGeometry()[0].GetValue(NORMAL) += NormalVector; 
            cond_it->GetGeometry()[1].GetValue(NORMAL) += NormalVector; 
            
            // Store in GPValues the projected value in the destiny condition Gauss Points
            // Note that currently the implementation is valid for only 2 GP
            array_1d<double, 3> TempValue = ZeroVect;
            
            for (unsigned int j = 0; j < 3; j++)
            {
                mGaussPointList[GPi]->GetProjectedValue(rOriginVar, TempValue, dimension);
                GPValues[j] = TempValue[j];

                mGaussPointList[GPi + 1]->GetProjectedValue(rOriginVar, TempValue, dimension);
                GPValues[j + 3] = TempValue[j];
            }

            const double K = CondLength/2.0;
            
            for (unsigned int j = 0; j < 3; j++)
            {
                NodalValues[j] = K*(MInterp(0,0)*GPValues[j] + MInterp(0,1)*GPValues[3 + j]);
                NodalValues[3 + j] = K*(MInterp(1,0)*GPValues[j] + MInterp(1,1)*GPValues[3 + j]);
            }

            boost::shared_ptr< array_1d<double,6> > pNodalValues(new array_1d<double,6>(NodalValues) );
            pInterpValues.push_back(pNodalValues); // This is computed here because it is the constant part of RHS

            GPi += 2; // 1 Condition = 2 Gauss Points
        }

        // Compute the nodal normals (normalised sum of the conditions normals) in the destination modelpart
        for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
                node_it != mrDestinationModelPart.NodesEnd();
                node_it++)
        {
            const array_1d<double,3> NormalVector = node_it->GetValue(NORMAL);
            node_it->GetValue(NORMAL) = NormalVector/norm_2(NormalVector);
        }

        // Iteration
        for (int k = 0; k < MaxIter; k++)
        {
            // At the begining of each iteration initialize the variable containing the assembled RHS as 0
            for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
                    node_it != mrDestinationModelPart.NodesEnd();
                    node_it++)
            {
                node_it->GetValue(MAPPER_VECTOR_PROJECTION_RHS) = ZeroVect;
            }

            array_1d<double, 3> LocalRHS0, LocalRHS1; // Local RHS for each node
            array_1d<double, 6> LastSolution;
            int IV_iter = 0;

            for (ModelPart::ConditionIterator cond_it = mrDestinationModelPart.ConditionsBegin();
                    cond_it != mrDestinationModelPart.ConditionsEnd();
                    cond_it++)
            {
                double CondLength = 0.0;
                mGaussPointList[2 * IV_iter]->GetArea(CondLength);
                
                LocalRHS0 = ZeroVect;
                LocalRHS1 = ZeroVect;
                
                array_1d<double, 3> NormalVector0 = cond_it->GetGeometry()[0].GetValue(NORMAL);
                array_1d<double, 3> NormalVector1 = cond_it->GetGeometry()[1].GetValue(NORMAL);

                for(unsigned int j = 0; j < 3; j++)
                {
                    LastSolution[j]     = sign * cond_it->GetGeometry()[0].FastGetSolutionStepValue(rDestVar) * NormalVector0[j];
                    LastSolution[3 + j] = sign * cond_it->GetGeometry()[1].FastGetSolutionStepValue(rDestVar) * NormalVector1[j];
                }            

                const double K = CondLength/6.0;
                array_1d<double,6> CondValues = *pInterpValues[IV_iter];
                
                for(unsigned int j = 0; j < 3; j++)
                {
                    LocalRHS0[j] = CondValues[j] - K * (2.0*LastSolution[j] + 1.0*LastSolution[3 + j]);
                    LocalRHS1[j] = CondValues[3 + j] - K * (1.0*LastSolution[j] + 2.0*LastSolution[3 + j]);
                }
                // We are taking advantage of 1 Condition = 2 Gauss Points to iterate the interpolation results and the Gauss Point Vector Again

                cond_it->GetGeometry()[0].GetValue(MAPPER_VECTOR_PROJECTION_RHS) += LocalRHS0;
                cond_it->GetGeometry()[1].GetValue(MAPPER_VECTOR_PROJECTION_RHS) += LocalRHS1;

                IV_iter++;
            }

            // Solve
            array_1d<double,3> dVal = ZeroVect;
            double dValNorm      = 0.0;
            double ValNorm       = 0.0;
            double RelativeError = 0.0;
            unsigned int NodeNum = 0;

            for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
                    node_it != mrDestinationModelPart.NodesEnd();
                    node_it++)
            {
                const array_1d<double,3> NormalVector = node_it->GetValue(NORMAL);
                const double & NodeLength = node_it->GetValue(NODAL_MAUX);
                dVal = node_it->GetValue(MAPPER_VECTOR_PROJECTION_RHS)/NodeLength;
                node_it->FastGetSolutionStepValue(rDestVar) += sign * inner_prod(dVal, NormalVector);

                // Variables for convergence check
                for (unsigned int j = 0; j < 3; j++)
                {
                    dValNorm += dVal[j] * dVal[j];
                    ValNorm += std::pow(node_it->FastGetSolutionStepValue(rDestVar) * NormalVector[j], 2);
                }

                NodeNum++;
            }
            //std::cout << "NormalVectorToScalarMap iteration: " << k+1 << std::endl;

            // Check Convergence
            if(ValNorm > 10e-15)
            {
                RelativeError = dValNorm / ValNorm;
            }
            if( (ValNorm/NodeNum < 0.00001 * std::pow(TolIter, 2.00)) || (RelativeError < std::pow(TolIter, 2.00)) )
            {
                std::cout << "NormalVectorToScalarMap converged in " << k + 1 << " iterations." << std::endl;
                break;
            }
            else if ((k + 1) == MaxIter)
            {
                std::cout << "WARNING: NormalVectorToScalarMap did not converge in " << k + 1 << " iterations." << std::endl;
            }
        }
    }
    else
    {

        std::vector< boost::shared_ptr<array_1d<double,9> > > pInterpValues;

        // Here we loop both the Destination Model Part and the Gauss Point List, using the
        // fact that the GP list was created by a loop over the Dest. Model Part's conditions
        int GPi = 0;
        for (ModelPart::ConditionIterator cond_it = mrDestinationModelPart.ConditionsBegin();
                cond_it != mrDestinationModelPart.ConditionsEnd();
                cond_it++)
        {
            array_1d<double, 9> GPValues;
            array_1d<double, 9> NodalValues;

            double CondArea = 0.0;
            mGaussPointList[GPi]->GetArea(CondArea);
            const double K = CondArea/24.0;

            array_1d<double,3> NormalVector = ZeroVect;
            mGaussPointList[GPi]->GetNormal(NormalVector);

            for (unsigned int i = 0; i < 3; i++)
            {
                //~ cond_it->GetGeometry()[i].GetValue(NODAL_MAUX) += 0.333333333333333333 * CondArea;
                cond_it->GetGeometry()[i].GetValue(NORMAL) += NormalVector;
            }

            array_1d<double,3> TempValues = ZeroVect;

            // Gauss point 1
            mGaussPointList[GPi]->GetProjectedValue(rOriginVar, TempValues, 3);

            for (unsigned int i = 0; i < 3; i++)
            {
                GPValues[i] = TempValues[i];
            }

            // Gauss point 2
            mGaussPointList[GPi + 1]->GetProjectedValue(rOriginVar, TempValues, 3);

            for (unsigned int i = 0; i < 3; i++)
            {
                GPValues[3 + i] = TempValues[i];
            }

            // Gauss point 3
            mGaussPointList[GPi + 2]->GetProjectedValue(rOriginVar, TempValues, 3);

            for (unsigned int i = 0; i < 3; i++)
            {
                GPValues[6 + i] = TempValues[i];
            }

            // Nodal values from GP projections
            for (unsigned int i = 0; i < 3; i++)
            {
                NodalValues[i]     = K*(6.0 * GPValues[i] + 1.0 * GPValues[3 + i] + 1.0 * GPValues[6 + i]);
                NodalValues[3 + i] = K*(1.0 * GPValues[i] + 6.0 * GPValues[3 + i] + 1.0 * GPValues[6 + i]);
                NodalValues[6 + i] = K*(1.0 * GPValues[i] + 1.0 * GPValues[3 + i] + 6.0 * GPValues[6 + i]);
            }

            boost::shared_ptr< array_1d<double,9> > pNodalValues(new array_1d<double,9>(NodalValues) );
            pInterpValues.push_back(pNodalValues); // This is computed here because it is the constant part of RHS

            GPi += 3; // 1 Condition = 3 Gauss Points
        }

        // Compute the destination interface nodal normals
        for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
                node_it != mrDestinationModelPart.NodesEnd();
                node_it++)
        {
            const array_1d<double,3> NormalVector = node_it->GetValue(NORMAL);
            node_it->GetValue(NORMAL) = NormalVector/norm_2(NormalVector);
        }

        // Iteration
        for (int k = 0; k < MaxIter; k++)
        {
            // At the begining of each iteration initialize the variable containing the assembled RHS as 0
            for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
                    node_it != mrDestinationModelPart.NodesEnd();
                    node_it++)
            {
                node_it->GetValue(MAPPER_VECTOR_PROJECTION_RHS) = ZeroVect;
            }

            array_1d<double, 3> LocalRHS0, LocalRHS1, LocalRHS2; // Local RHS for each node
            array_1d<double, 9> LastSolution;
            int IV_iter = 0;

            for (ModelPart::ConditionIterator cond_it = mrDestinationModelPart.ConditionsBegin();
                    cond_it != mrDestinationModelPart.ConditionsEnd();
                    cond_it++)
            {
                double CondArea = 0.0;
                mGaussPointList[3 * IV_iter]->GetArea(CondArea);
                const double K = CondArea/12.0;

                LocalRHS0 = ZeroVect;
                LocalRHS1 = ZeroVect;
                LocalRHS2 = ZeroVect;

                for(unsigned int j = 0; j < 3; j++)
                {
                    array_1d<double,3> NormalVector = cond_it->GetGeometry()[0].GetValue(NORMAL);
                    LastSolution[j    ] = sign * cond_it->GetGeometry()[0].FastGetSolutionStepValue(rDestVar) * NormalVector[j];
                    NormalVector = cond_it->GetGeometry()[1].GetValue(NORMAL);
                    LastSolution[j + 3] = sign * cond_it->GetGeometry()[1].FastGetSolutionStepValue(rDestVar) * NormalVector[j];
                    NormalVector = cond_it->GetGeometry()[2].GetValue(NORMAL);
                    LastSolution[j + 6] = sign * cond_it->GetGeometry()[2].FastGetSolutionStepValue(rDestVar) * NormalVector[j];
                }

                array_1d<double,9> CondValues = *pInterpValues[IV_iter];

                for(unsigned int j = 0; j < 3; j++)
                {
                    LocalRHS0[j] = CondValues[j]     - K * (2.0 * LastSolution[j] + 1.0 * LastSolution[j + 3] + 1.0 * LastSolution[j + 6]);
                    LocalRHS1[j] = CondValues[j + 3] - K * (1.0 * LastSolution[j] + 2.0 * LastSolution[j + 3] + 1.0 * LastSolution[j + 6]);
                    LocalRHS2[j] = CondValues[j + 6] - K * (1.0 * LastSolution[j] + 1.0 * LastSolution[j + 3] + 2.0 * LastSolution[j + 6]);
                }
                // We are taking advantage of 1 Condition = 3 Gauss Points to iterate the interpolation results and the Gauss Point Vector Again

                cond_it->GetGeometry()[0].GetValue(MAPPER_VECTOR_PROJECTION_RHS) += LocalRHS0;
                cond_it->GetGeometry()[1].GetValue(MAPPER_VECTOR_PROJECTION_RHS) += LocalRHS1;
                cond_it->GetGeometry()[2].GetValue(MAPPER_VECTOR_PROJECTION_RHS) += LocalRHS2;

                IV_iter++;
            }

            // Solve
            array_1d<double,3> dVal = ZeroVect;
            double dValNorm      = 0.0;
            double ValNorm       = 0.0;
            double RelativeError = 0.0;
            unsigned int NodeNum = 0;

            for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
                    node_it != mrDestinationModelPart.NodesEnd();
                    node_it++)
            {
                const array_1d<double,3> NormalVector = node_it->GetValue(NORMAL);
                const double & NodeArea = node_it->GetValue(NODAL_MAUX);
                noalias(dVal) = node_it->GetValue(MAPPER_VECTOR_PROJECTION_RHS)/NodeArea;
                node_it->FastGetSolutionStepValue(rDestVar) += sign * inner_prod(dVal, NormalVector);

                // Variables for convergence check
                for (unsigned int j = 0; j < 3; j++)
                {
                    dValNorm += dVal[j] * dVal[j];
                    ValNorm += std::pow(node_it->FastGetSolutionStepValue(rDestVar) * NormalVector[j], 2);
                }

                NodeNum++;
            }
            //std::cout << "NormalVectorToScalarMap iteration: " << k+1 << std::endl;

            // Check Convergence
            if(ValNorm > 10e-15)
            {
                RelativeError = dValNorm / ValNorm;
            }
            if( (ValNorm/NodeNum < 0.00001 * TolIter*TolIter) || RelativeError < TolIter * TolIter)
            {
                std::cout << "NormalVectorToScalarMap converged in " << k + 1 << " iterations." << std::endl;
                break;
            }
            else if ( (k + 1) == MaxIter)
            {
                std::cout << "WARNING: NormalVectorToScalarMap did not converge in " << k + 1 << " iterations." << std::endl;
            }
        } // End of Iteration
    }
} // End of Map (normal vector to scalar version)

/***********************************************************************************/
/***********************************************************************************/

void AdvancedNMPointsMapper::ScalarMap(
        const Variable<double> & rOriginVar,
        Variable<double> & rDestVar,
        const int MaxIter,
        const double TolIter,
        const bool sign_pos
        )
{

    // Define if the mapping swaps the sign of the variable values
    double sign = 1.0;
    if (sign_pos == false)
    {
        sign = -1.0;
    }

    // Initialize results and NODAL_MAUX
    for ( ModelPart::NodesContainerType::iterator node_it = mrDestinationModelPart.NodesBegin();
            node_it != mrDestinationModelPart.NodesEnd();
            node_it++)
    {
        node_it->GetValue(NODAL_MAUX) = 0.0;
        node_it->FastGetSolutionStepValue(rDestVar) = 0.0;
    }

    const unsigned int dimension = mrDestinationModelPart.ConditionsBegin()->GetGeometry().WorkingSpaceDimension();
    
    // Compute nodal lengths/areas in both origin and destination modelparts
    ComputeNodalLengthArea();

    if (dimension == 2) // 2D case
    {

        // Interpolation matrix obtention
        boost::numeric::ublas::bounded_matrix<double,2,2> MCons; // Elemental Consistent Mass Matrix = L/6 * MCons
        boost::numeric::ublas::bounded_matrix<double,2,2> MInterp;
        boost::numeric::ublas::bounded_matrix<double,2,2> aux_GPPos;
        boost::numeric::ublas::bounded_matrix<double,2,2> inv_aux_GPPos;
        
        MCons(0, 0) = 2.0/3.0; 
        MCons(0, 1) = 1.0/3.0;
        MCons(1, 0) = 1.0/3.0;
        MCons(1, 1) = 2.0/3.0;
        
        aux_GPPos(0, 0) = (1.0+(-std::sqrt(1.0/3.0)))/2.0;
        aux_GPPos(0, 1) = (1.0-(-std::sqrt(1.0/3.0)))/2.0;
        aux_GPPos(1, 0) = (1.0+(+std::sqrt(1.0/3.0)))/2.0;
        aux_GPPos(1, 1) = (1.0-(+std::sqrt(1.0/3.0)))/2.0;
        
        double det_aux_GPPos;
        det_aux_GPPos = (aux_GPPos(0,0)*aux_GPPos(1,1))-(aux_GPPos(0,1)*aux_GPPos(1,0));
        
        inv_aux_GPPos(0, 0) = aux_GPPos(1, 1)/det_aux_GPPos;
        inv_aux_GPPos(0, 1) = -aux_GPPos(0, 1)/det_aux_GPPos;
        inv_aux_GPPos(1, 0) = -aux_GPPos(1, 0)/det_aux_GPPos;
        inv_aux_GPPos(1, 1) = aux_GPPos(0, 0)/det_aux_GPPos;

        MInterp = prod(MCons, inv_aux_GPPos); // Interpolation Matrix (NodalValues = (L/6)*MInterp*GaussValues)
        
        std::vector< boost::shared_ptr<array_1d<double,2> > > pInterpValues;

        // Here we loop both the Destination Model Part and the Gauss Point List, using the
        // fact that the GP list was created by a loop over the Dest. Model Part's conditions
        int GPi = 0;
        for (ModelPart::ConditionIterator cond_it = mrDestinationModelPart.ConditionsBegin();
                cond_it != mrDestinationModelPart.ConditionsEnd();
                cond_it++)
        {
            array_1d<double, 2> GPValues;
            array_1d<double, 2> NodalValues;

            double CondLength = 0.0;
            mGaussPointList[GPi]->GetArea(CondLength); // Gets the length of the parent condition of the two points considered
            const double K = CondLength/2.0;
    
            // Store in GPValues the projected value in the destiny condition Gauss Points
            // Note that currently the implementation is valid for only 2 GP
            double TempValue = 0.0;

            mGaussPointList[GPi]->GetProjectedValue(rOriginVar, TempValue, dimension);
            GPValues[0] = TempValue;

            mGaussPointList[GPi + 1]->GetProjectedValue(rOriginVar, TempValue, dimension);
            GPValues[1] = TempValue;

            NodalValues[0] = K*(MInterp(0,0)*GPValues[0] + MInterp(0,1)*GPValues[1]);
            NodalValues[1] = K*(MInterp(1,0)*GPValues[0] + MInterp(1,1)*GPValues[1]);

            boost::shared_ptr< array_1d<double,2> > pNodalValues(new array_1d<double,2>(NodalValues) );
            pInterpValues.push_back(pNodalValues); // This is computed here because it is the constant part of RHS

            GPi += 2; // 1 Condition = 2 Gauss Points
        }
        
        // Iteration
        for (int k = 0; k < MaxIter; k++)
        {
            // At the begining of each iteration initialize the variable containing the assembled RHS as 0
            for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
                    node_it != mrDestinationModelPart.NodesEnd();
                    node_it++)
            {
                node_it->GetValue(MAPPER_SCALAR_PROJECTION_RHS) = 0.0;
            }

            double LocalRHS0, LocalRHS1; // Local RHS for each node
            array_1d<double, 2> LastSolution;
            int IV_iter = 0;

            for (ModelPart::ConditionIterator cond_it = mrDestinationModelPart.ConditionsBegin();
                    cond_it != mrDestinationModelPart.ConditionsEnd();
                    cond_it++)
            {
                double CondLength = 0.0;
                mGaussPointList[2 * IV_iter]->GetArea(CondLength);
                const double K = CondLength/6.0;
                
                LocalRHS0 = 0.0;
                LocalRHS1 = 0.0;
                
                LastSolution[0] = sign * cond_it->GetGeometry()[0].FastGetSolutionStepValue(rDestVar);
                LastSolution[1] = sign * cond_it->GetGeometry()[1].FastGetSolutionStepValue(rDestVar);

                array_1d<double,2> CondValues = *pInterpValues[IV_iter];

                LocalRHS0 = CondValues[0] - K * (2.0*LastSolution[0] + 1.0*LastSolution[1]);
                LocalRHS1 = CondValues[1] - K * (1.0*LastSolution[0] + 2.0*LastSolution[1]);

                // We are taking advantage of 1 Condition = 2 Gauss Points to iterate the interpolation results and the Gauss Point Vector Again
                cond_it->GetGeometry()[0].GetValue(MAPPER_SCALAR_PROJECTION_RHS) += LocalRHS0;
                cond_it->GetGeometry()[1].GetValue(MAPPER_SCALAR_PROJECTION_RHS) += LocalRHS1;

                IV_iter++;
            }

            // Solve
            double dVal          = 0.0;
            double dValNorm      = 0.0;
            double ValNorm       = 0.0;
            double RelativeError = 0.0;
            unsigned int NodeNum = 0;

            for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
                    node_it != mrDestinationModelPart.NodesEnd();
                    node_it++)
            {
                const double & NodeLength = node_it->GetValue(NODAL_MAUX);
                dVal = node_it->GetValue(MAPPER_SCALAR_PROJECTION_RHS)/NodeLength;
                node_it->FastGetSolutionStepValue(rDestVar) += sign * dVal;

                // Variables for convergence check
                dValNorm += dVal;
                ValNorm += node_it->FastGetSolutionStepValue(rDestVar);

                NodeNum++;
            }
            //std::cout << "ScalarMap iteration: " << k+1 << std::endl;

            // Check Convergence
            if(ValNorm > 10e-15)
            {
                RelativeError = dValNorm / ValNorm;
            }
            if( (ValNorm/NodeNum < 0.00001 * TolIter * TolIter) || RelativeError < TolIter * TolIter)
            {
                std::cout << "ScalarMap converged in " << k + 1 << " iterations." << std::endl;
                break;
            }
            else if ( (k + 1) == MaxIter)
            {
                std::cout << "WARNING: VectorMap did not converge in " << k + 1 << " iterations." << std::endl;
            }
        } // End of Iteration
    }
    
    else // 3D case
    {
        // Define some variables that will be used in the iteration
        MatrixVar MCons; // Elemental Consistent Mass Matrix = Aelem/12 * MCons
        MCons(0, 0) = 2.0;
        MCons(0, 1) = 1.0;
        MCons(0, 2) = 1.0;
        MCons(1, 0) = 1.0;
        MCons(1, 1) = 2.0;
        MCons(1, 2) = 1.0;
        MCons(2, 0) = 1.0;
        MCons(2, 1) = 1.0;
        MCons(2, 2) = 2.0;

        MatrixVar MInterp; // Interpolation Matrix (NodalValues = (A/24)*MInterp*GaussValues)
        MInterp(0, 0) = 6.0;
        MInterp(0, 1) = 1.0;
        MInterp(0, 2) = 1.0;
        MInterp(1, 0) = 1.0;
        MInterp(1, 1) = 6.0;
        MInterp(1, 2) = 1.0;
        MInterp(2, 0) = 1.0;
        MInterp(2, 1) = 1.0;
        MInterp(2, 2) = 6.0;

        std::vector< boost::shared_ptr< array_1d<double, 3> > > pInterpValues;

        // Here we loop both the Destination Model Part and the Gauss Point List, using the
        // fact that the GP list was created by a loop over the Dest. Model Part's conditions
        int GPi = 0;
        for (ModelPart::ConditionsContainerType::iterator cond_it = mrDestinationModelPart.ConditionsBegin();
                cond_it != mrDestinationModelPart.ConditionsEnd();
                cond_it++)
        {
            array_1d<double,3> GPValues;
            array_1d<double,3> NodalValues;

            double CondArea = 0.0;
            mGaussPointList[GPi]->GetArea(CondArea);

            for (unsigned int i = 0; i < 3; i++)
            {
                mGaussPointList[GPi + i]->GetProjectedValue(rOriginVar, GPValues[i], dimension);
            }

            noalias(NodalValues) = (CondArea/24.0) * prod(MInterp, GPValues);

            boost::shared_ptr< array_1d<double, 3> > pNodalValues(new array_1d<double, 3>(NodalValues) );
            pInterpValues.push_back(pNodalValues); // This is computed here because it is the constant part of RHS

            GPi += 3; // 1 Condition = 3 Gauss Points
        }

        // Iteration
        for (int k = 0; k < MaxIter; k++)
        {
            // At the begining of each iteration initialize the variable containing the assembled RHS as 0
            for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
                    node_it != mrDestinationModelPart.NodesEnd();
                    node_it++)
            {
                node_it->GetValue(MAPPER_SCALAR_PROJECTION_RHS) = 0.0;
            }

            array_1d<double, 3> LocalRHS;
            array_1d<double, 3> LastSolution;
            int IV_iter = 0;

            for (ModelPart::ConditionsContainerType::iterator cond_it = mrDestinationModelPart.ConditionsBegin();
                    cond_it != mrDestinationModelPart.ConditionsEnd();
                    cond_it++)
            {
                for(unsigned int j = 0; j < 3; j++)
                {
                    LocalRHS[j] = 0.0;
                    LastSolution[j] = sign * cond_it->GetGeometry()[j].FastGetSolutionStepValue(rDestVar);
                }

                double CondArea = 0.0;
                mGaussPointList[3 * IV_iter]->GetArea(CondArea);

                noalias(LocalRHS) = *pInterpValues[IV_iter] - (CondArea/12.0) * prod(MCons, LastSolution);
                // We are taking advantage of 1 Condition = 3 Gauss Points to iterate the interpolation results and the Gauss Point Vector Again

                for (unsigned int j = 0; j < 3 ; j++)
                {
                    cond_it->GetGeometry()[j].GetValue(MAPPER_SCALAR_PROJECTION_RHS) += LocalRHS[j];
                }

                IV_iter++;
            }

            // Solve
            double dValNorm      = 0.0;
            double ValNorm       = 0.0;
            double RelativeError = 0.0;
            unsigned int NodeNum = 0;

            for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
                    node_it != mrDestinationModelPart.NodesEnd();
                    node_it++)
            {
                const double & NodeArea = node_it->GetValue(NODAL_MAUX);
                //node_it->FastGetSolutionStepValue(NODAL_MAUX) = NodeArea; // TEST: store nodal area so GiD can paint it later
                double dVal = node_it->GetValue(MAPPER_SCALAR_PROJECTION_RHS)/NodeArea;
                node_it->FastGetSolutionStepValue(rDestVar) += sign * dVal;

                // Variables for convergence check
                dValNorm += dVal * dVal;
                ValNorm  += std::pow(node_it->FastGetSolutionStepValue(rDestVar), 2);

                NodeNum++;
            }
            //std::cout << "ScalarMap iteration: " << k+1 << std::endl;

            // Check Convergence
            if(ValNorm > 10e-15)
            {
                RelativeError = dValNorm / ValNorm;
            }
            if( (ValNorm/NodeNum < 0.00001 * std::pow(TolIter, 2.00)) || (RelativeError < std::pow(TolIter, 2.00)) )
            {
                std::cout << "ScalarMap converged in " << k + 1 << " iterations." << std::endl;
                break;
            }
            else if ((k + 1) == MaxIter)
            {
                std::cout << "WARNING: ScalarMap did not converge in " << k + 1 << " iterations." << std::endl;
            }
        } // End of Iteration
    }

} // End of Map (scalar version)

/***********************************************************************************/
/***********************************************************************************/

void AdvancedNMPointsMapper::VectorMap(
        const Variable<array_1d<double,3> >& rOriginVar,
        Variable<array_1d<double,3> >& rDestVar,
        const int MaxIter,
        const double TolIter,
        const bool sign_pos,
        const bool distributed
        )
{
    array_1d<double,3> ZeroVect = ZeroVector(3);
    
    // Define if the mapping swaps the sign of the variable values
    double sign = 1.0;
    if (sign_pos == false)
    {
        sign = -1.0;
    }

    // Initialize results and NODAL_MAUX
    for ( ModelPart::NodesContainerType::iterator node_it = mrDestinationModelPart.NodesBegin();
            node_it != mrDestinationModelPart.NodesEnd();
            node_it++)
    {
        node_it->GetValue(NODAL_MAUX) = 0.0;
        node_it->FastGetSolutionStepValue(rDestVar) = ZeroVect;
    }
    
    // Compute nodal lengths/areas in both origin and destination modelparts
    ComputeNodalLengthArea();
    
    // If dealing with punctual loads, obtain their equivalent tractions
    if (distributed == true)
    {
        int MaxIterTractions = MaxIter;
        double TolIterTractions = TolIter;
        ComputeEquivalentTractions(rOriginVar, MaxIterTractions, TolIterTractions);
    }

    const unsigned int dimension = mrDestinationModelPart.ConditionsBegin()->GetGeometry().WorkingSpaceDimension();
    
    if (dimension == 2) // 2D case
    {
        // Interpolation matrix obtention
        boost::numeric::ublas::bounded_matrix<double,2,2> MCons; // Elemental Consistent Mass Matrix = L/2 * MCons
        boost::numeric::ublas::bounded_matrix<double,2,2> MInterp;
        boost::numeric::ublas::bounded_matrix<double,2,2> aux_GPPos;
        boost::numeric::ublas::bounded_matrix<double,2,2> inv_aux_GPPos;
        
        MCons(0, 0) = 2.0/3.0; 
        MCons(0, 1) = 1.0/3.0;
        MCons(1, 0) = 1.0/3.0;
        MCons(1, 1) = 2.0/3.0;
        
        aux_GPPos(0, 0) = (1.0+(-std::sqrt(1.0/3.0)))/2.0;
        aux_GPPos(0, 1) = (1.0-(-std::sqrt(1.0/3.0)))/2.0;
        aux_GPPos(1, 0) = (1.0+(+std::sqrt(1.0/3.0)))/2.0;
        aux_GPPos(1, 1) = (1.0-(+std::sqrt(1.0/3.0)))/2.0;
        
        double det_aux_GPPos;
        det_aux_GPPos = (aux_GPPos(0,0)*aux_GPPos(1,1))-(aux_GPPos(0,1)*aux_GPPos(1,0));
        
        inv_aux_GPPos(0, 0) = aux_GPPos(1, 1)/det_aux_GPPos;
        inv_aux_GPPos(0, 1) = -aux_GPPos(0, 1)/det_aux_GPPos;
        inv_aux_GPPos(1, 0) = -aux_GPPos(1, 0)/det_aux_GPPos;
        inv_aux_GPPos(1, 1) = aux_GPPos(0, 0)/det_aux_GPPos;

        MInterp = prod(MCons, inv_aux_GPPos); // Interpolation matrix (NodalValues = (L/6)*MInterp*GaussValues)
        
        std::vector< boost::shared_ptr<array_1d<double,6> > > pInterpValues;

        // Here we loop both the Destination Model Part and the Gauss Point List, using the
        // fact that the GP list was created by a loop over the Dest. Model Part's conditions
        int GPi = 0;
        for (ModelPart::ConditionIterator cond_it = mrDestinationModelPart.ConditionsBegin();
                cond_it != mrDestinationModelPart.ConditionsEnd();
                cond_it++)
        {
            
            double CondLength = 0.0;
            mGaussPointList[GPi]->GetArea(CondLength); // Gets the length of the parent condition of the two points considered
            const double K = CondLength/2.0;
    
            // Store in GPValues the projected value in the destiny condition Gauss Points
            // Note that currently the implementation is valid for only 2 GP           
            array_1d<double, 6> GPValues = ZeroVector(6);
            array_1d<double, 6> NodalValues = ZeroVector(6);
            array_1d<double, 3> TempValues = ZeroVect;
            
            // Gauss point 1
            if (distributed == false)
            {
                mGaussPointList[GPi]->GetProjectedValue(rOriginVar, TempValues, dimension);
            }
            else
            {   
                mGaussPointList[GPi]->GetProjectedValue(VAUX_EQ_TRACTION, TempValues, dimension);
            }

            for (unsigned int i = 0; i < 3; i++)
            {
                GPValues[i] = TempValues[i];
            }
            
            // Gauss point 2
            if (distributed == false)
            {
                mGaussPointList[GPi + 1]->GetProjectedValue(rOriginVar, TempValues, dimension);
            }
            else
            {
                mGaussPointList[GPi + 1]->GetProjectedValue(VAUX_EQ_TRACTION, TempValues, dimension);
            }

            for (unsigned int i = 0; i < 3; i++)
            {
                GPValues[3 + i] = TempValues[i];
            }

            // Compute the nodal values from the projected ones using the interpolation matrix
            NodalValues[0] = K*(MInterp(0,0)*GPValues[0] + MInterp(0,1)*GPValues[3]);
            NodalValues[1] = K*(MInterp(0,0)*GPValues[1] + MInterp(0,1)*GPValues[4]);
            NodalValues[2] = 0.0;
            NodalValues[3] = K*(MInterp(1,0)*GPValues[0] + MInterp(1,1)*GPValues[3]);
            NodalValues[4] = K*(MInterp(1,0)*GPValues[1] + MInterp(1,1)*GPValues[4]);
            NodalValues[5] = 0.0;

            boost::shared_ptr< array_1d<double,6> > pNodalValues(new array_1d<double,6>(NodalValues) );
            pInterpValues.push_back(pNodalValues); // This is computed here because it is the constant part of RHS

            GPi += 2; // 1 Condition = 2 Gauss Points
        }
        
        // Iteration
        for (int k = 0; k < MaxIter; k++)
        {
            // At the begining of each iteration initialize the variable containing the assembled RHS as 0
            for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
                    node_it != mrDestinationModelPart.NodesEnd();
                    node_it++)
            {
                node_it->GetValue(MAPPER_VECTOR_PROJECTION_RHS) = ZeroVect;
            }

            array_1d<double, 3> LocalRHS0, LocalRHS1; // Local RHS for each node
            array_1d<double, 6> LastSolution;
            int IV_iter = 0;

            for (ModelPart::ConditionIterator cond_it = mrDestinationModelPart.ConditionsBegin();
                    cond_it != mrDestinationModelPart.ConditionsEnd();
                    cond_it++)
            {
                double CondLength = 0.0;
                mGaussPointList[2 * IV_iter]->GetArea(CondLength);
                const double K = CondLength/6.0;

                LocalRHS0 = ZeroVector(3);
                LocalRHS1 = ZeroVector(3);
    
                if (distributed == false)
                {
                    for(unsigned int j = 0; j < 3; j++)
                    {
                        LastSolution[j]     = sign * cond_it->GetGeometry()[0].FastGetSolutionStepValue(rDestVar)[j];
                        LastSolution[3 + j] = sign * cond_it->GetGeometry()[1].FastGetSolutionStepValue(rDestVar)[j];
                    }
                }
                else
                {
                    for(unsigned int j = 0; j < 3; j++)
                    {
                        LastSolution[j]     = sign * cond_it->GetGeometry()[0].FastGetSolutionStepValue(VAUX_EQ_TRACTION)[j];
                        LastSolution[3 + j] = sign * cond_it->GetGeometry()[1].FastGetSolutionStepValue(VAUX_EQ_TRACTION)[j];
                    }
                }
    
                array_1d<double,6> CondValues = *pInterpValues[IV_iter];

                for(unsigned int j = 0; j < 3; j++)
                {
                    LocalRHS0[j] = CondValues[j]     - K * (2.0*LastSolution[j] + 1.0*LastSolution[j + 3]);
                    LocalRHS1[j] = CondValues[j + 3] - K * (1.0*LastSolution[j] + 2.0*LastSolution[j + 3]);
                }
                // We are taking advantage of 1 Condition = 2 Gauss Points to iterate the interpolation results and the Gauss Point Vector Again

                cond_it->GetGeometry()[0].GetValue(MAPPER_VECTOR_PROJECTION_RHS) += LocalRHS0;
                cond_it->GetGeometry()[1].GetValue(MAPPER_VECTOR_PROJECTION_RHS) += LocalRHS1;

                IV_iter++;
            }

            // Solve
            array_1d<double,3> dVal = ZeroVect;
            double dValNorm      = 0.0;
            double ValNorm       = 0.0;
            double RelativeError = 0.0;
            unsigned int NodeNum = 0;

            if (distributed == false)
            {
                for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
                        node_it != mrDestinationModelPart.NodesEnd();
                        node_it++)
                {
                    const double & NodeLength = node_it->GetValue(NODAL_MAUX);
                    noalias(dVal) = node_it->GetValue(MAPPER_VECTOR_PROJECTION_RHS)/NodeLength;
                    node_it->FastGetSolutionStepValue(rDestVar) += sign * dVal;

                    // Variables for convergence check
                    for (unsigned int j = 0; j < 3; j++)
                    {
                        dValNorm += dVal[j] * dVal[j];
                        ValNorm += std::pow(node_it->FastGetSolutionStepValue(rDestVar)[j], 2);
                    }

                    NodeNum++;
                }
            }
            else
            {
                for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
                        node_it != mrDestinationModelPart.NodesEnd();
                        node_it++)
                {
                    const double & NodeLength = node_it->GetValue(NODAL_MAUX);
                    noalias(dVal) = node_it->GetValue(MAPPER_VECTOR_PROJECTION_RHS)/NodeLength;
                    node_it->FastGetSolutionStepValue(VAUX_EQ_TRACTION) += sign * dVal;

                    // Variables for convergence check
                    for (unsigned int j = 0; j < 3; j++)
                    {
                        dValNorm += dVal[j] * dVal[j];
                        ValNorm += std::pow(node_it->FastGetSolutionStepValue(VAUX_EQ_TRACTION)[j], 2);
                    }

                    NodeNum++;
                }
            }
            //~ std::cout << "VectorMap iteration: " << k+1 << "dValNorm: " << dValNorm << std::endl;

            // Check Convergence
            if(ValNorm > 10e-15)
            {
                RelativeError = dValNorm / ValNorm;
            }
            if( (ValNorm/NodeNum < 0.00001 * TolIter * TolIter) || RelativeError < TolIter * TolIter)
            {
                std::cout << "VectorMap converged in " << k + 1 << " iterations." << std::endl;
                break;
            }
            else if ( (k + 1) == MaxIter)
            {
                std::cout << "WARNING: VectorMap did not converge in " << k + 1 << " iterations." << std::endl;
            }
        } 
        // End of Iteration
        
    }
    else // 3D case
    {
        std::vector< boost::shared_ptr<array_1d<double,9> > > pInterpValues;

        // Here we loop both the Destination Model Part and the Gauss Point List, using the
        // fact that the GP list was created by a loop over the Dest. Model Part's conditions
        int GPi = 0;
        for (ModelPart::ConditionIterator cond_it = mrDestinationModelPart.ConditionsBegin();
                cond_it != mrDestinationModelPart.ConditionsEnd();
                cond_it++)
        {
            array_1d<double, 9> GPValues;
            array_1d<double, 9> NodalValues;

            double CondArea = 0.0;
            mGaussPointList[GPi]->GetArea(CondArea);
            const double K = CondArea/24.0;

            array_1d<double,3> TempValues = ZeroVect;
            
            if (distributed == false)
            {
                // Gauss point 1
                mGaussPointList[GPi]->GetProjectedValue(rOriginVar, TempValues, dimension);

                for (unsigned int i = 0; i < 3; i++)
                {
                    GPValues[i] = TempValues[i];
                }

                // Gauss point 2
                mGaussPointList[GPi + 1]->GetProjectedValue(rOriginVar, TempValues, dimension);

                for (unsigned int i = 0; i < 3; i++)
                {
                    GPValues[3 + i] = TempValues[i];
                }

                // Gauss point 3
                mGaussPointList[GPi + 2]->GetProjectedValue(rOriginVar, TempValues, dimension);

                for (unsigned int i = 0; i < 3; i++)
                {
                    GPValues[6 + i] = TempValues[i];
                }
            }
            else
            {
                // Gauss point 1
                mGaussPointList[GPi]->GetProjectedValue(VAUX_EQ_TRACTION, TempValues, dimension);

                for (unsigned int i = 0; i < 3; i++)
                {
                    GPValues[i] = TempValues[i];
                }

                // Gauss point 2
                mGaussPointList[GPi + 1]->GetProjectedValue(VAUX_EQ_TRACTION, TempValues, dimension);

                for (unsigned int i = 0; i < 3; i++)
                {
                    GPValues[3 + i] = TempValues[i];
                }

                // Gauss point 3
                mGaussPointList[GPi + 2]->GetProjectedValue(VAUX_EQ_TRACTION, TempValues, dimension);

                for (unsigned int i = 0; i < 3; i++)
                {
                    GPValues[6 + i] = TempValues[i];
                }
            }

            // Compute the nodal values from the projected ones using the interpolation matrix
            for (unsigned int i = 0; i < 3; i++)
            {
                NodalValues[i]     = K*(6.0 * GPValues[i] + 1.0 * GPValues[3 + i] + 1.0 * GPValues[6 + i]);
                NodalValues[3 + i] = K*(1.0 * GPValues[i] + 6.0 * GPValues[3 + i] + 1.0 * GPValues[6 + i]);
                NodalValues[6 + i] = K*(1.0 * GPValues[i] + 1.0 * GPValues[3 + i] + 6.0 * GPValues[6 + i]);
            }

            boost::shared_ptr< array_1d<double,9> > pNodalValues(new array_1d<double,9>(NodalValues) );
            pInterpValues.push_back(pNodalValues); // This is computed here because it is the constant part of RHS

            GPi += 3; // 1 Condition = 3 Gauss Points
        }

        // Iteration
        for (int k = 0; k < MaxIter; k++)
        {
            // At the begining of each iteration initialize the variable containing the assembled RHS as 0
            for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
                    node_it != mrDestinationModelPart.NodesEnd();
                    node_it++)
            {
                node_it->GetValue(MAPPER_VECTOR_PROJECTION_RHS) = ZeroVect;
            }

            array_1d<double, 3> LocalRHS0, LocalRHS1, LocalRHS2; // Local RHS for each node
            array_1d<double, 9> LastSolution;
            int IV_iter = 0;

            for (ModelPart::ConditionIterator cond_it = mrDestinationModelPart.ConditionsBegin();
                    cond_it != mrDestinationModelPart.ConditionsEnd();
                    cond_it++)
            {
                double CondArea = 0.0;
                mGaussPointList[3 * IV_iter]->GetArea(CondArea);
                const double K = CondArea/12.0;

                LocalRHS0 = ZeroVect;
                LocalRHS1 = ZeroVect;
                LocalRHS2 = ZeroVect;

                if (distributed == false)
                {
                    for(unsigned int j = 0; j < 3; j++)
                    {
                        LastSolution[j]     = sign * cond_it->GetGeometry()[0].FastGetSolutionStepValue(rDestVar)[j];
                        LastSolution[3 + j] = sign * cond_it->GetGeometry()[1].FastGetSolutionStepValue(rDestVar)[j];
                        LastSolution[6 + j] = sign * cond_it->GetGeometry()[2].FastGetSolutionStepValue(rDestVar)[j];
                    }
                }
                else
                {
                    for(unsigned int j = 0; j < 3; j++)
                    {
                        LastSolution[j]     = sign * cond_it->GetGeometry()[0].FastGetSolutionStepValue(VAUX_EQ_TRACTION)[j];
                        LastSolution[3 + j] = sign * cond_it->GetGeometry()[1].FastGetSolutionStepValue(VAUX_EQ_TRACTION)[j];
                        LastSolution[6 + j] = sign * cond_it->GetGeometry()[2].FastGetSolutionStepValue(VAUX_EQ_TRACTION)[j];
                    }
                }

                array_1d<double,9> CondValues = *pInterpValues[IV_iter];

                for(unsigned int j = 0; j < 3; j++)
                {
                    LocalRHS0[j] = CondValues[j]     - K * (2.0 * LastSolution[j] + 1.0 * LastSolution[j + 3] + 1.0 * LastSolution[j + 6]);
                    LocalRHS1[j] = CondValues[j + 3] - K * (1.0 * LastSolution[j] + 2.0 * LastSolution[j + 3] + 1.0 * LastSolution[j + 6]);
                    LocalRHS2[j] = CondValues[j + 6] - K * (1.0 * LastSolution[j] + 1.0 * LastSolution[j + 3] + 2.0 * LastSolution[j + 6]);
                }
                // We are taking advantage of 1 Condition = 3 Gauss Points to iterate the interpolation results and the Gauss Point Vector Again

                cond_it->GetGeometry()[0].GetValue(MAPPER_VECTOR_PROJECTION_RHS) += LocalRHS0;
                cond_it->GetGeometry()[1].GetValue(MAPPER_VECTOR_PROJECTION_RHS) += LocalRHS1;
                cond_it->GetGeometry()[2].GetValue(MAPPER_VECTOR_PROJECTION_RHS) += LocalRHS2;

                IV_iter++;
            }

            // Solve
            array_1d<double,3> dVal = ZeroVect;
            double dValNorm      = 0.0;
            double ValNorm       = 0.0;
            double RelativeError = 0.0;
            unsigned int NodeNum = 0;

            if (distributed == false)
            {
                for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
                        node_it != mrDestinationModelPart.NodesEnd();
                        node_it++)
                {
                    const double & NodeArea = node_it->GetValue(NODAL_MAUX);
                    noalias(dVal) = node_it->GetValue(MAPPER_VECTOR_PROJECTION_RHS)/NodeArea;
                    node_it->FastGetSolutionStepValue(rDestVar) += sign * dVal;

                    // Variables for convergence check
                    for (unsigned int j = 0; j < 3; j++)
                    {
                        dValNorm += dVal[j] * dVal[j];
                        ValNorm += std::pow(node_it->FastGetSolutionStepValue(rDestVar)[j], 2);
                    }

                    NodeNum++;
                }
            }
            else
            {
                for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
                        node_it != mrDestinationModelPart.NodesEnd();
                        node_it++)
                {
                    const double & NodeArea = node_it->GetValue(NODAL_MAUX);
                    noalias(dVal) = node_it->GetValue(MAPPER_VECTOR_PROJECTION_RHS)/NodeArea;
                    node_it->FastGetSolutionStepValue(VAUX_EQ_TRACTION) += sign * dVal;

                    // Variables for convergence check
                    for (unsigned int j = 0; j < 3; j++)
                    {
                        dValNorm += dVal[j] * dVal[j];
                        ValNorm += std::pow(node_it->FastGetSolutionStepValue(VAUX_EQ_TRACTION)[j], 2);
                    }

                    NodeNum++;
                }
            }
            //std::cout << "VectorMap iteration: " << k+1 << std::endl;

            // Check Convergence
            if(ValNorm > 10e-15)
            {
                RelativeError = dValNorm / ValNorm;
            }
            if( (ValNorm/NodeNum < 0.00001 * TolIter * TolIter) || RelativeError < TolIter * TolIter)
            {
                std::cout << "VectorMap converged in " << k + 1 << " iterations." << std::endl;
                break;
            }
            else if ( (k + 1) == MaxIter)
            {
                std::cout << "WARNING: VectorMap did not converge in " << k + 1 << " iterations." << std::endl;
            }
        } 
        // End of Iteration       
    }

    // Convert the computed equivalent tractions to punctual nodal loads
    if (distributed == true)
    {
        ComputeNodalLoadsFromTractions(rDestVar);
    }
    
} // End of Map (vector version)

/***********************************************************************************/
/***********************************************************************************/

void AdvancedNMPointsMapper::DistanceCheck()
{
    unsigned int GPiter = 0;
    for (ModelPart::ConditionsContainerType::iterator cond_it = mrDestinationModelPart.ConditionsBegin();
            cond_it != mrDestinationModelPart.ConditionsEnd();
            cond_it++)
    {
        for (unsigned int i = 0; i < 3; i++)
        {
            double dist = 0.0;
            int Status  = 0;
            mGaussPointList[GPiter + i]->GetProjStatus(Status);
            std::cout << GPiter + i << " " << Status << std::endl;

            if (Status != 0)
            {
                mGaussPointList[GPiter+i]->GetDist(dist);

                if (Status == 2)
                {
                    dist = -std::sqrt(dist);
                }
            }
            else
            {
                dist = -10000;
            }

            cond_it->GetGeometry()[i].FastGetSolutionStepValue(FICTITIOUS_FLUID_DENSITY) = dist;
        }

        GPiter += 3;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AdvancedNMPointsMapper::ComputeNodalLengthArea()
{
    const unsigned int dimension = mrDestinationModelPart.ConditionsBegin()->GetGeometry().WorkingSpaceDimension();
    
    // NODAL_MAUX initialization
    for (ModelPart::NodesContainerType::const_iterator node_it = mrOriginModelPart.NodesBegin();
            node_it != mrOriginModelPart.NodesEnd();
            node_it++)
        {
            node_it->GetValue(NODAL_MAUX) = 0.0;
        }
        
    for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
            node_it != mrDestinationModelPart.NodesEnd();
            node_it++)
        {
            node_it->GetValue(NODAL_MAUX) = 0.0;
        }
    
    // DestinationModelPart
    if (dimension == 2)
    {
        double CondLength = 0.0;
        
        for (ModelPart::ConditionIterator cond_it = mrDestinationModelPart.ConditionsBegin();
                cond_it != mrDestinationModelPart.ConditionsEnd();
                cond_it++)
            {
                for (unsigned int i = 0; i < 2; i++)
                {
                    CondLength = cond_it->GetGeometry().Length();
                    cond_it->GetGeometry()[i].GetValue(NODAL_MAUX) += 0.5 * CondLength;
                }
            }
    }
    else
    {
        double CondArea = 0.0;
        
        for (ModelPart::ConditionIterator cond_it = mrDestinationModelPart.ConditionsBegin();
                cond_it != mrDestinationModelPart.ConditionsEnd();
                cond_it++)
            {                
                for (unsigned int i = 0; i < 3; i++)
                {
                    CondArea = cond_it->GetGeometry().Area();
                    cond_it->GetGeometry()[i].GetValue(NODAL_MAUX) += (1.0/3.0) * CondArea;
                }
            }
    }
    
    // TEST: store nodal area so GiD can paint it later
    for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
        node_it != mrDestinationModelPart.NodesEnd();
        node_it++)
        {   
            node_it->FastGetSolutionStepValue(NODAL_MAUX) = node_it->GetValue(NODAL_MAUX); 
        }
    
    // OriginModelPart
    if (dimension == 2)
    {
        double CondLength = 0.0;
            
        for(ModelPart::ConditionsContainerType::const_iterator cond_it = mrOriginModelPart.ConditionsBegin();
            cond_it != mrOriginModelPart.ConditionsEnd();
            cond_it++)
            {
                for (unsigned int i = 0; i < 2; i++)
                {
                    CondLength = cond_it->GetGeometry().Length();
                    cond_it->GetGeometry()[i].GetValue(NODAL_MAUX) += 0.5 * CondLength;
                }
            }
    }
    else
    {
        double CondArea = 0.0;
        
        for(ModelPart::ConditionsContainerType::const_iterator cond_it = mrOriginModelPart.ConditionsBegin();
            cond_it != mrOriginModelPart.ConditionsEnd();
            cond_it++)
            {
                for (unsigned int i = 0; i < 3; i++)
                {
                    CondArea = cond_it->GetGeometry().Area();
                    cond_it->GetGeometry()[i].GetValue(NODAL_MAUX) += (1.0/3.0) * CondArea;
                }
            }
    }
    
    // TEST: store nodal area so GiD can paint it later
    for (ModelPart::NodesContainerType::const_iterator node_it = mrOriginModelPart.NodesBegin();
        node_it != mrOriginModelPart.NodesEnd();
        node_it++)
        {   
            node_it->FastGetSolutionStepValue(NODAL_MAUX) = node_it->GetValue(NODAL_MAUX); 
        }

}

/***********************************************************************************/
/***********************************************************************************/

void AdvancedNMPointsMapper::ComputeEquivalentTractions(
     const Variable<array_1d<double,3> >& rOriginVar,
     const int MaxIter,
     const double TolIter)
{
    const unsigned int dimension = mrOriginModelPart.ConditionsBegin()->GetGeometry().WorkingSpaceDimension();

    if (dimension == 2)
    {

        boost::numeric::ublas::bounded_matrix<double,2,2> MassMat; // Elemental consistent mass matrix 2 Gauss points in a 1D element
            
        MassMat(0, 0) = 2.0/3.0; 
        MassMat(0, 1) = 1.0/3.0;
        MassMat(1, 0) = 1.0/3.0;
        MassMat(1, 1) = 2.0/3.0;
        
        // Initialization of equivalent tractions
        for (ModelPart::NodesContainerType::const_iterator node_it = mrOriginModelPart.NodesBegin();
             node_it != mrOriginModelPart.NodesEnd();
             node_it++)
        {
            node_it->FastGetSolutionStepValue(VAUX_EQ_TRACTION) = ZeroVector(3);
        }
        
        // Store the initial guess        
        for (ModelPart::NodesContainerType::const_iterator node_it = mrOriginModelPart.NodesBegin();
             node_it != mrOriginModelPart.NodesEnd();
             node_it++)
        {
            for(unsigned int j = 0; j < 2; j++)
            {
                const double & NodeLength = node_it->GetValue(NODAL_MAUX);
                node_it->FastGetSolutionStepValue(VAUX_EQ_TRACTION)[j] = (node_it->FastGetSolutionStepValue(rOriginVar)[j])/NodeLength;
            }
        }                   
        
        // Iteration
        for (int k = 0; k < MaxIter; k++)
        {
            // At the begining of each iteration initialize the variable containing the assembled RHS as 0
            for (ModelPart::NodesContainerType::const_iterator node_it = mrOriginModelPart.NodesBegin();
                 node_it != mrOriginModelPart.NodesEnd();
                 node_it++)
            {
                node_it->GetValue(MAPPER_VECTOR_PROJECTION_RHS) = ZeroVector(3);
            }
            
            array_1d<double, 3> LocalRHS0, LocalRHS1; // Local RHS for each node
            array_1d<double, 6> LastSolution;
            array_1d<double, 6> OriginNodalValues;
        
            for(ModelPart::ConditionsContainerType::const_iterator cond_it = mrOriginModelPart.ConditionsBegin();
                cond_it != mrOriginModelPart.ConditionsEnd();
                cond_it++)
            {
                
                double CondLength = 0.0;
                CondLength = cond_it->GetGeometry().Length();
                const double Jac = 0.5*CondLength;

                LocalRHS0 = ZeroVector(3);
                LocalRHS1 = ZeroVector(3);
                OriginNodalValues = ZeroVector(6);
                
                // Original nodal values (assembled)
                for(unsigned int j = 0; j < 3; j++)
                {
                    OriginNodalValues[j]     = cond_it->GetGeometry()[0].FastGetSolutionStepValue(rOriginVar)[j];
                    OriginNodalValues[3 + j] = cond_it->GetGeometry()[1].FastGetSolutionStepValue(rOriginVar)[j];
                }
                
                // Unassemble the nodal values to obtain their elemental contributions
                double aux0 = Jac/(cond_it->GetGeometry()[0].GetValue(NODAL_MAUX));
                double aux1 = Jac/(cond_it->GetGeometry()[1].GetValue(NODAL_MAUX));
                                
                for(unsigned int j = 0; j < 3; j++)
                {
                    OriginNodalValues[j]     *= aux0;
                    OriginNodalValues[3 + j] *= aux1;
                }
                
                // Previous iteration solution
                for(unsigned int j = 0; j < 3; j++)
                {
                    LastSolution[j]     = cond_it->GetGeometry()[0].FastGetSolutionStepValue(VAUX_EQ_TRACTION)[j];
                    LastSolution[3 + j] = cond_it->GetGeometry()[1].FastGetSolutionStepValue(VAUX_EQ_TRACTION)[j];
                }

                // Compute the nodal RHS
                for(unsigned int j = 0; j < 3; j++)
                {
                    LocalRHS0[j] = OriginNodalValues[j]     - Jac * (MassMat(0,0)*LastSolution[j] + MassMat(0,1)*LastSolution[j + 3]);
                    LocalRHS1[j] = OriginNodalValues[j + 3] - Jac * (MassMat(1,0)*LastSolution[j] + MassMat(1,1)*LastSolution[j + 3]);
                }

                // Accumulate the nodal RHS (localRHS)
                cond_it->GetGeometry()[0].GetValue(MAPPER_VECTOR_PROJECTION_RHS) += LocalRHS0;
                cond_it->GetGeometry()[1].GetValue(MAPPER_VECTOR_PROJECTION_RHS) += LocalRHS1;
            }

            // Solve
            array_1d<double,3> dVal = ZeroVector(3);
            double dValNorm      = 0.0;
            double ValNorm       = 0.0;
            double RelativeError = 0.0;
            unsigned int NodeNum = 0;

            for (ModelPart::NodesContainerType::const_iterator node_it = mrOriginModelPart.NodesBegin();
                 node_it != mrOriginModelPart.NodesEnd();
                 node_it++)
            {
                const double & NodeLength = node_it->GetValue(NODAL_MAUX);
                noalias(dVal) = node_it->GetValue(MAPPER_VECTOR_PROJECTION_RHS)/NodeLength;              
                node_it->FastGetSolutionStepValue(VAUX_EQ_TRACTION) += dVal;

                // Variables for convergence check
                for (unsigned int j = 0; j < 3; j++)
                {
                    dValNorm += dVal[j] * dVal[j];
                    ValNorm += std::pow(node_it->FastGetSolutionStepValue(VAUX_EQ_TRACTION)[j], 2);
                }

                NodeNum++;
            }
            
            //~ std::cout << "Compute equivalent tractions iteration: " << k+1 << " dValNorm " << dValNorm << std::endl;

            // Check Convergence
            if(ValNorm > 10e-15)
            {
                RelativeError = dValNorm / ValNorm;
            }
            if( (ValNorm/NodeNum < 0.00001 * TolIter * TolIter) || RelativeError < TolIter * TolIter)
            {
                std::cout << "ComputeEquivalentTractions converged in " << k + 1 << " iterations." << std::endl;
                break;
            }
            else if ( (k + 1) == MaxIter)
            {
                std::cout << "WARNING: ComputeEquivalentTractions did not converge in " << k + 1 << " iterations." << std::endl;
            }
        } // End of Iteration        
    }
    else
    {
        boost::numeric::ublas::bounded_matrix<double,3,3> MassMat; // Elemental consistent mass matrix 3 Gauss points in a 2D triangular element
            
        MassMat(0, 0) = 1.0/12.0; 
        MassMat(0, 1) = 1.0/24.0;
        MassMat(0, 2) = 1.0/24.0;
        MassMat(1, 0) = 1.0/24.0;
        MassMat(1, 1) = 1.0/12.0;
        MassMat(1, 2) = 1.0/24.0;
        MassMat(2, 0) = 1.0/24.0;
        MassMat(2, 1) = 1.0/24.0;
        MassMat(2, 2) = 1.0/12.0;
        
        // Initialization of equivalent tractions
        for (ModelPart::NodesContainerType::const_iterator node_it = mrOriginModelPart.NodesBegin();
             node_it != mrOriginModelPart.NodesEnd();
             node_it++)
        {
            node_it->FastGetSolutionStepValue(VAUX_EQ_TRACTION) = ZeroVector(3);
        }
        
        // Store the initial guess        
        for (ModelPart::NodesContainerType::const_iterator node_it = mrOriginModelPart.NodesBegin();
             node_it != mrOriginModelPart.NodesEnd();
             node_it++)
        {
            for(unsigned int j = 0; j < 3; j++)
            {
                const double & NodeArea = node_it->GetValue(NODAL_MAUX);
                node_it->FastGetSolutionStepValue(VAUX_EQ_TRACTION)[j] = (node_it->FastGetSolutionStepValue(rOriginVar)[j])/NodeArea;
            }
        }
        
        // Iteration
        for (int k = 0; k < MaxIter; k++)
        {
            // At the begining of each iteration initialize the variable containing the assembled RHS as 0
            for (ModelPart::NodesContainerType::const_iterator node_it = mrOriginModelPart.NodesBegin();
                 node_it != mrOriginModelPart.NodesEnd();
                 node_it++)
            {
                node_it->GetValue(MAPPER_VECTOR_PROJECTION_RHS) = ZeroVector(3);
            }
            
            array_1d<double, 3> LocalRHS0, LocalRHS1, LocalRHS2; // Local RHS for each node
            array_1d<double, 9> LastSolution;
            array_1d<double, 9> OriginNodalValues;
        
            for(ModelPart::ConditionsContainerType::const_iterator cond_it = mrOriginModelPart.ConditionsBegin();
                cond_it != mrOriginModelPart.ConditionsEnd();
                cond_it++)
            {
                double CondArea = 0.0;
                CondArea = cond_it->GetGeometry().Area();
                const double Jac = 2.0*CondArea;

                LocalRHS0 = ZeroVector(3);
                LocalRHS1 = ZeroVector(3);
                LocalRHS2 = ZeroVector(3);
                OriginNodalValues = ZeroVector(9);
                
                // Original nodal values (assembled)
                for(unsigned int j = 0; j < 3; j++)
                {
                    OriginNodalValues[j]     = cond_it->GetGeometry()[0].FastGetSolutionStepValue(rOriginVar)[j];
                    OriginNodalValues[3 + j] = cond_it->GetGeometry()[1].FastGetSolutionStepValue(rOriginVar)[j];
                    OriginNodalValues[6 + j] = cond_it->GetGeometry()[2].FastGetSolutionStepValue(rOriginVar)[j];
                }
                
                // Unassemble the nodal values to obtain their elemental contributions
                double aux0 = (CondArea/3.0)/(cond_it->GetGeometry()[0].GetValue(NODAL_MAUX));
                double aux1 = (CondArea/3.0)/(cond_it->GetGeometry()[1].GetValue(NODAL_MAUX));
                double aux2 = (CondArea/3.0)/(cond_it->GetGeometry()[2].GetValue(NODAL_MAUX));
                                
                for(unsigned int j = 0; j < 3; j++)
                {
                    OriginNodalValues[j]     *= aux0;
                    OriginNodalValues[3 + j] *= aux1;
                    OriginNodalValues[6 + j] *= aux2;
                }

                // Previous iteration solution
                for(unsigned int j = 0; j < 3; j++)
                {
                    LastSolution[j]     = cond_it->GetGeometry()[0].FastGetSolutionStepValue(VAUX_EQ_TRACTION)[j];
                    LastSolution[3 + j] = cond_it->GetGeometry()[1].FastGetSolutionStepValue(VAUX_EQ_TRACTION)[j];
                    LastSolution[6 + j] = cond_it->GetGeometry()[2].FastGetSolutionStepValue(VAUX_EQ_TRACTION)[j];
                }
                
                // Compute the nodal RHS
                for(unsigned int j = 0; j < 3; j++)
                {
                    LocalRHS0[j] = OriginNodalValues[j]     - Jac * (MassMat(0,0)*LastSolution[j] + MassMat(0,1)*LastSolution[j + 3] + MassMat(0,2)*LastSolution[j + 6]);
                    LocalRHS1[j] = OriginNodalValues[j + 3] - Jac * (MassMat(1,0)*LastSolution[j] + MassMat(1,1)*LastSolution[j + 3] + MassMat(1,2)*LastSolution[j + 6]);
                    LocalRHS2[j] = OriginNodalValues[j + 6] - Jac * (MassMat(2,0)*LastSolution[j] + MassMat(2,1)*LastSolution[j + 3] + MassMat(2,2)*LastSolution[j + 6]);
                }

                // Accumulate the nodal RHS (localRHS)
                cond_it->GetGeometry()[0].GetValue(MAPPER_VECTOR_PROJECTION_RHS) += LocalRHS0;
                cond_it->GetGeometry()[1].GetValue(MAPPER_VECTOR_PROJECTION_RHS) += LocalRHS1;
                cond_it->GetGeometry()[2].GetValue(MAPPER_VECTOR_PROJECTION_RHS) += LocalRHS2;
            }

            // Solve
            array_1d<double,3> dVal = ZeroVector(3);
            double dValNorm      = 0.0;
            double ValNorm       = 0.0;
            double RelativeError = 0.0;
            unsigned int NodeNum = 0;

            for (ModelPart::NodesContainerType::const_iterator node_it = mrOriginModelPart.NodesBegin();
                 node_it != mrOriginModelPart.NodesEnd();
                 node_it++)
            {
                const double & NodeArea = node_it->GetValue(NODAL_MAUX);
                noalias(dVal) = node_it->GetValue(MAPPER_VECTOR_PROJECTION_RHS)/NodeArea;
                node_it->FastGetSolutionStepValue(VAUX_EQ_TRACTION) += dVal;
                
                // Variables for convergence check
                for (unsigned int j = 0; j < 3; j++)
                {
                    dValNorm += dVal[j] * dVal[j];
                    ValNorm += std::pow(node_it->FastGetSolutionStepValue(VAUX_EQ_TRACTION)[j], 2);
                }

                NodeNum++;
            }
            
            //~ std::cout << "Compute equivalent tractions iteration: " << k+1 << " dValNorm " << dValNorm << std::endl;

            // Check Convergence
            if(ValNorm > 10e-15)
            {
                RelativeError = dValNorm / ValNorm;
            }
            if( (ValNorm/NodeNum < 0.00001 * TolIter * TolIter) || RelativeError < TolIter * TolIter)
            {
                std::cout << "Compute equivalent tractions converged in " << k + 1 << " iterations." << std::endl;
                break;
            }
            else if ( (k + 1) == MaxIter)
            {
                std::cout << "WARNING: Compute equivalent tractions did not converge in " << k + 1 << " iterations." << std::endl;
            }
        } // End of Iteration
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AdvancedNMPointsMapper::ComputeNodalLoadsFromTractions(
     const Variable<array_1d<double,3> >& rDestVar)
{
    const unsigned int dimension = mrDestinationModelPart.ConditionsBegin()->GetGeometry().WorkingSpaceDimension();
    array_1d<double,3> ZeroVect = ZeroVector(3);

    if (dimension == 2)
    {

        boost::numeric::ublas::bounded_matrix<double,2,2> MassMat; // Elemental consistent mass matrix 2 Gauss points in a 1D element
            
        MassMat(0, 0) = 2.0/3.0; 
        MassMat(0, 1) = 1.0/3.0;
        MassMat(1, 0) = 1.0/3.0;
        MassMat(1, 1) = 2.0/3.0;
        
        array_1d<double, 3> Node0Values;
        array_1d<double, 3> Node1Values;
        array_1d<double, 6> NodalTractions;
    
        for(ModelPart::ConditionsContainerType::const_iterator cond_it = mrDestinationModelPart.ConditionsBegin();
            cond_it != mrDestinationModelPart.ConditionsEnd();
            cond_it++)
        {
            double CondLength = 0.0;
            CondLength = cond_it->GetGeometry().Length();
            const double Jac = 0.5*CondLength;
            
            Node0Values = ZeroVector(3);
            Node1Values = ZeroVector(3);
            NodalTractions = ZeroVector(6);
            
            for(unsigned int j = 0; j < 3; j++)
            {
                NodalTractions[j]     = cond_it->GetGeometry()[0].FastGetSolutionStepValue(VAUX_EQ_TRACTION)[j];
                NodalTractions[3 + j] = cond_it->GetGeometry()[1].FastGetSolutionStepValue(VAUX_EQ_TRACTION)[j];
            }
                        
            for(unsigned int j = 0; j < 3; j++)
            {
                Node0Values[j] = Jac*(MassMat(0,0)*NodalTractions[j]+MassMat(0,1)*NodalTractions[j+3]);
                Node1Values[j] = Jac*(MassMat(1,0)*NodalTractions[j]+MassMat(1,1)*NodalTractions[j+3]);
            }

            // We are taking advantage of 1 Condition = 2 Gauss Points to iterate the interpolation results and the Gauss Point Vector Again
            cond_it->GetGeometry()[0].FastGetSolutionStepValue(rDestVar) += Node0Values;
            cond_it->GetGeometry()[1].FastGetSolutionStepValue(rDestVar) += Node1Values;
        }
    }
    else
    {
        boost::numeric::ublas::bounded_matrix<double,3,3> MassMat; // Elemental consistent mass matrix 3 Gauss points in a 2D triangular element
            
        MassMat(0, 0) = 1.0/12.0; 
        MassMat(0, 1) = 1.0/24.0;
        MassMat(0, 2) = 1.0/24.0;
        MassMat(1, 0) = 1.0/24.0;
        MassMat(1, 1) = 1.0/12.0;
        MassMat(1, 2) = 1.0/24.0;
        MassMat(2, 0) = 1.0/24.0;
        MassMat(2, 1) = 1.0/24.0;
        MassMat(2, 2) = 1.0/12.0;
        
        array_1d<double, 3> Node0Values;
        array_1d<double, 3> Node1Values;
        array_1d<double, 3> Node2Values;
        array_1d<double, 9> NodalTractions;
        int IV_iter = 0;
    
        for(ModelPart::ConditionsContainerType::const_iterator cond_it = mrDestinationModelPart.ConditionsBegin();
            cond_it != mrDestinationModelPart.ConditionsEnd();
            cond_it++)
        {
            double CondArea = 0.0;
            CondArea = cond_it->GetGeometry().Area();
            const double Jac = 2.0*CondArea;
            
            Node0Values = ZeroVector(3);
            Node1Values = ZeroVector(3);
            Node2Values = ZeroVector(3);
            NodalTractions = ZeroVector(9);
            
            for(unsigned int j = 0; j < 3; j++)
            {
                NodalTractions[j]     = cond_it->GetGeometry()[0].FastGetSolutionStepValue(VAUX_EQ_TRACTION)[j];
                NodalTractions[3 + j] = cond_it->GetGeometry()[1].FastGetSolutionStepValue(VAUX_EQ_TRACTION)[j];
                NodalTractions[6 + j] = cond_it->GetGeometry()[2].FastGetSolutionStepValue(VAUX_EQ_TRACTION)[j];
            }
                        
            for(unsigned int j = 0; j < 3; j++)
            {
                Node0Values[j] = Jac*(MassMat(0,0)*NodalTractions[j] + MassMat(0,1)*NodalTractions[j+3] + MassMat(0,2)*NodalTractions[j+6]);
                Node1Values[j] = Jac*(MassMat(1,0)*NodalTractions[j] + MassMat(1,1)*NodalTractions[j+3] + MassMat(1,2)*NodalTractions[j+6]);
                Node2Values[j] = Jac*(MassMat(2,0)*NodalTractions[j] + MassMat(2,1)*NodalTractions[j+3] + MassMat(2,2)*NodalTractions[j+6]);
            }

            // We are taking advantage of 1 Condition = 2 Gauss Points to iterate the interpolation results and the Gauss Point Vector Again
            cond_it->GetGeometry()[0].FastGetSolutionStepValue(rDestVar) += Node0Values;
            cond_it->GetGeometry()[1].FastGetSolutionStepValue(rDestVar) += Node1Values;
            cond_it->GetGeometry()[2].FastGetSolutionStepValue(rDestVar) += Node2Values;
            
            IV_iter++;
        }
    }
}

} // Namespace Kratos.
