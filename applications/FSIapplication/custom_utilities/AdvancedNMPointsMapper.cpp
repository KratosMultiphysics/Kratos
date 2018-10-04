//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela, Vicente Mataix and Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "AdvancedNMPointsMapper.hpp"

namespace Kratos
{
// GaussPointItem Methods

/**
 * It projects in 2D/3D for a line/triangle a returns the local coordinates and distance
 * @param pOriginCond: Pointer to the Gauss point origin condition
 * @return ProjectedLocalCoords: Projection local coordinates
 * @return dist: The distance between the point and the plane
 */
void GaussPointItem::Project(Condition::Pointer pOriginCond,
                             array_1d<double,2>& ProjectedLocalCoords,
                             double& Dist)
{
    GeometryType& rOriginGeom = pOriginCond->GetGeometry();
    const unsigned int dimension = rOriginGeom.WorkingSpaceDimension();

    if (dimension == 2)
    {
        Point point_projected;
        Point point_to_project = Point(this->X(), this->Y(), this->Z());
        ProjectPointToLine(rOriginGeom[0], point_to_project, point_projected, Dist);

        array_1d<double, 3> point_projected_local_coor;
        point_projected_local_coor = rOriginGeom.PointLocalCoordinates(point_projected_local_coor, point_projected);

        ProjectedLocalCoords[0] = point_projected_local_coor[0];
        ProjectedLocalCoords[1] = 0.0;
    }
    else
    {
        // xi,yi,zi are Nodal Coordinates, n is the destination condition's unit normal
        // and d is the distance along n from the point to its projection in the condition
        // | DestX-x0 |   | x1-x0 x2-x0 nx |   | Chi |
        // | DestY-y0 | = | y1-y0 y2-y0 ny | . | Eta |
        // | DestZ-z0 |   | z1-z0 z2-z0 nz |   |  d  |

        Matrix ChangeMatrix(3, 3, false);
        Matrix InvChange(3, 3, false);
        double det;

        array_1d<double, 3> RHS, Res;

        noalias(RHS) = this->Coordinates() - rOriginGeom[0].Coordinates();

        ChangeMatrix(0, 0) = rOriginGeom[1].X() - rOriginGeom[0].X();
        ChangeMatrix(1, 0) = rOriginGeom[1].Y() - rOriginGeom[0].Y();
        ChangeMatrix(2, 0) = rOriginGeom[1].Z() - rOriginGeom[0].Z();

        ChangeMatrix(0, 1) = rOriginGeom[2].X() - rOriginGeom[0].X();
        ChangeMatrix(1, 1) = rOriginGeom[2].Y() - rOriginGeom[0].Y();
        ChangeMatrix(2, 1) = rOriginGeom[2].Z() - rOriginGeom[0].Z();

        ChangeMatrix(0, 2) = mNormal[0];
        ChangeMatrix(1, 2) = mNormal[1];
        ChangeMatrix(2, 2) = mNormal[2];

        MathUtils<double>::InvertMatrix3(ChangeMatrix,InvChange,det);
        noalias(Res) = prod(InvChange, RHS);

        ProjectedLocalCoords[0] = Res[0];
        ProjectedLocalCoords[1] = Res[1];
        // Keep distance positive, regardless of normal orientation
        Dist = (Res[2] < 0)? -Res[2] : Res[2] ;
    }
}

/**
 * Projects a point over a line
 * @param PointInPlane: A point in the plane
 * @param PointToBeProjected: The point to be projected
 * @return PointProjected: The point pojected over the plane
 * @return dist: The distance between the point and the plane
 */
void GaussPointItem::ProjectPointToLine(const Point& PointInPlane,
                                        const Point& PointToBeProjected,
                                        Point& PointProjected,
                                        double& dist)
{
     array_1d<double,3> vector_points;
     noalias(vector_points) = PointToBeProjected.Coordinates() - PointInPlane.Coordinates();

     dist = inner_prod(vector_points, mNormal);
     PointProjected.Coordinates() = PointToBeProjected.Coordinates() - mNormal * dist;
}

/**
 * It gets the projected value for scalar variables
 * @param rOriginVar: The variable (scalar) in the original condition
 * @return Value: The projected value (scalar)
 */
void GaussPointItem::GetProjectedValue(const Variable<double> & rOriginVar,
                                       double& Value)
{
    Value = 0.0;    // Value initialization (if mProjStatus == 2 it will remain as 0.0)

    if (mProjStatus == 1) // Get Interpolated value from origin condition
    {
        GeometryType& rOriginGeom = (mpOriginCond.lock())->GetGeometry();
        const unsigned int dimension = rOriginGeom.WorkingSpaceDimension();

        // Shape functions values in the projected Gauss pt.
        Point GPloccoords = (dimension == 2) ? Point(mOriginCoords[0], 0.0, 0.0) : Point(mOriginCoords[0], mOriginCoords[1], 0.0);
        Vector shfunc_values;
        rOriginGeom.ShapeFunctionsValues(shfunc_values, GPloccoords);

        // Interpolate the nodal values in the projected Gauss pt.
        for (unsigned int i=0; i<rOriginGeom.size(); ++i)
        {
            Value += rOriginGeom[i].FastGetSolutionStepValue(rOriginVar)*shfunc_values[i];
        }
    }
    else if (mProjStatus == 2)   // Get Value from origin node
    {
        Value = (mpOriginNode.lock())->FastGetSolutionStepValue(rOriginVar);
    }
}

/**
 * It gets the projected value for vector variables
 * @param rOriginVar: The variable (vector) in the original condition
 * @return Value: The projected value (vector)
 */
void GaussPointItem::GetProjectedValue(const Variable<array_1d<double,3> >& rOriginVar,
                                       array_1d<double,3>& Value)
{
    Value = ZeroVector(3);    // Value initialization (if mProjStatus == 2 it will remain as 0.0)
    GeometryType& rOriginGeom = (mpOriginCond.lock())->GetGeometry();
    const unsigned int dimension = rOriginGeom.WorkingSpaceDimension();

    if (mProjStatus == 1) // Get Interpolated value from origin condition
    {
        // Shape functions values in the projected Gauss pt.
        Point GPloccoords = (dimension == 2) ? Point(mOriginCoords[0], 0.0, 0.0) : Point(mOriginCoords[0], mOriginCoords[1], 0.0);
        Vector shfunc_values;
        rOriginGeom.ShapeFunctionsValues(shfunc_values, GPloccoords);

        // Interpolate the nodal values in the projected Gauss pt.
        for (unsigned int comp=0; comp<dimension; ++comp)
        {
            for (unsigned int i=0; i<rOriginGeom.size(); ++i)
            {
                Value[comp] += rOriginGeom[i].FastGetSolutionStepValue(rOriginVar)[comp]*shfunc_values[i];
            }
        }
    }
    else if (mProjStatus == 2)   // Get Value from origin node
    {
        Value = (mpOriginNode.lock())->FastGetSolutionStepValue(rOriginVar);
    }
}

/***********************************************************************************/
/***********************************************************************************/

// Mapper Methods
AdvancedNMPointsMapper::AdvancedNMPointsMapper(ModelPart& rOriginModelPart,
                                               ModelPart& rDestinationModelPart):
    mrOriginModelPart(rOriginModelPart),
    mrDestinationModelPart(rDestinationModelPart),
    mBucketSize(4)
{
    const unsigned int dimension = mrDestinationModelPart.ConditionsBegin()->GetGeometry().WorkingSpaceDimension();
    const unsigned int ngauss = (dimension == 2) ? 2 : 3 ; //TODO: Get them from the quadratures
    const unsigned int nnodes = (rDestinationModelPart.ConditionsBegin()->GetGeometry()).size();

    // Get the OpenMP partitioned vectors
    const int number_of_threads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::PartitionVector conditions_partition;
    OpenMPUtils::DivideInPartitions(static_cast<int>(rDestinationModelPart.NumberOfConditions()), number_of_threads, conditions_partition);

    // Get the shape functions values at the Gauss pts.
    boost::numeric::ublas::matrix<double>         GPPos(ngauss, nnodes);
    if (dimension == 2) // 2D case
    {
        // 2 Gauss-Legendre point quadrature
        // eps = +-(1/3)^0.5
        // N1 = (1-eps)/2
        // N2 = (1+eps)/2

        GPPos(0, 0) = (1.0+(-std::sqrt(1.0/3.0)))/2.0; GPPos(0, 1) = (1.0-(-std::sqrt(1.0/3.0)))/2.0;
        GPPos(1, 0) = (1.0+(+std::sqrt(1.0/3.0)))/2.0; GPPos(1, 1) = (1.0-(+std::sqrt(1.0/3.0)))/2.0;
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

        GPPos(0,0) = 0.6; GPPos(0,1) = 0.2; GPPos(0,2) = 0.2;
        GPPos(1,0) = 0.2; GPPos(1,1) = 0.6; GPPos(1,2) = 0.2;
        GPPos(2,0) = 0.2; GPPos(2,1) = 0.2; GPPos(2,2) = 0.6;
    }

    std::vector<GaussPointVector> AuxGaussPointVector(number_of_threads);

    #pragma omp parallel for firstprivate(GPPos) shared(AuxGaussPointVector)
    for (int l=0; l<number_of_threads; ++l)
    {
        GaussPointVector PrivateGaussPointList;
        array_1d<double, 3> GPCoord, Normal;
        boost::numeric::ublas::matrix<double>              Nodes(nnodes, 3); // Matrix that stores the nodal coordinates of the condition
        boost::numeric::ublas::matrix<double>        GaussPoints(ngauss, 3); // Matrix to store the obtained coordinates of the Gauss pts.

        for (int k=conditions_partition[l]; k<conditions_partition[l+1]; ++k)
        {
            const ModelPart::ConditionsContainerType::iterator cond_it = rDestinationModelPart.ConditionsBegin() + k;
            const GeometryType& rGeom = cond_it->GetGeometry();
            const double DomainSize = rGeom.DomainSize(); // Length in 2D and area in 3D
            ComputeConditionNormal(*(cond_it.base()), Normal);

            for(unsigned int i = 0; i < nnodes; i++)
            {
                const array_1d<double,3>& r_coordinates = rGeom[i].Coordinates();
                for(unsigned int j = 0; j < 3; j++)
                {
                    Nodes(i,j) = r_coordinates[j];
                }
            }

            noalias(GaussPoints) = prod(GPPos, Nodes);

            for(unsigned int i = 0; i < ngauss; i++)
            {
                for(unsigned int j = 0; j < 3; j++)
                {
                    GPCoord[j] = GaussPoints(i,j);
                }
                GaussPointItem::Pointer pGP = GaussPointItem::Pointer(new GaussPointItem(GPCoord, DomainSize, Normal));
                PrivateGaussPointList.push_back( pGP );
            }
        }

        const unsigned int thread_id = OpenMPUtils::ThisThread();
        AuxGaussPointVector[thread_id] = PrivateGaussPointList;
    }

    // This is deliverately done outside the previous loop to insert the Gauss pts. preserving the conditions order
    for (int i=0; i<number_of_threads; ++i)
    {
        mGaussPointList.insert(mGaussPointList.end(), AuxGaussPointVector[i].begin(), AuxGaussPointVector[i].end());
    }

}

/**
 * It searches neighbours nodes in a specific radius
 * @param SearchRadiusFactor: The radius of search
 */
void AdvancedNMPointsMapper::FindNeighbours(double SearchRadiusFactor)
{
    // Initialize some values
    const unsigned int MaxResults = 5000; // Maximum number of points to find in a single tree search
    GaussPointVector Results(MaxResults);
    std::vector<double> ResultDistances(MaxResults);
    const array_1d<double,3> ZeroVect = ZeroVector(3);

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
        Point Center;
        ComputeGeometryCenterAndRadius(*cond_it.base(), Center, Radius);
        MaxRadius = std::max(MaxRadius, Radius);

        // Create a fake GaussPoint with the condition center coordinates
        // (to ensure that the object used in the kd-tree search is of the same type as the tree contents)
        GaussPointItem CenterGP(Center.Coordinates(), 0, ZeroVect);

        // Count Gauss Points within SearchRadius
        SearchRadius = SearchRadiusFactor * Radius;
        unsigned int Found = Tree_conds.SearchInRadius(CenterGP, SearchRadius, Results.begin(), ResultDistances.begin(), MaxResults);

        for(unsigned int i = 0; i < Found; i++)
        {
            SetProjectionToCond(*Results[i], *cond_it.base());
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
        KRATOS_WARNING("AdvancedNMPointsMapper") << ProjectionlessGP.size()
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

        KRATOS_WARNING("AdvancedNMPointsMapper") << counter << " Gauss Points without a reference value remain." << std::endl;
    }
}

/**
 * It calculates the normal and area of a condition
 * @param pCond: The pointer to the condition
 * @return Normal: The normal of the condition
 */
void AdvancedNMPointsMapper::ComputeConditionNormal(const Condition::Pointer Cond,
                                                    array_1d<double,3> & Normal)
{
    const unsigned int dimension = Cond->WorkingSpaceDimension();

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

    // Compute the condition unit normal vector
    MathUtils<double>::CrossProduct(Normal,v1,v2);
    double NNorm = std::sqrt(Normal[0]*Normal[0] + Normal[1]*Normal[1] + Normal[2]*Normal[2]);
    Normal /= NNorm;
}

/**
 * It calculates the the center and radius of the condition
 * @param pCond: The pointer to the condition
 * @return Center: Center point (3D)
 * @return Radius: Geometry radius
 */
void AdvancedNMPointsMapper::ComputeGeometryCenterAndRadius(const Condition::Pointer Cond,
                                                            Point& Center,
                                                            double& Radius)
{
    Radius = 0.0;
    const GeometryType& rGeom = Cond->GetGeometry();
    const unsigned int nnodes = rGeom.size();

    Center = rGeom.Center();

    for(unsigned int i = 0; i < nnodes; i++)
    {
        const Node<3>& r_node = Cond->GetGeometry()[i];
        double dx = Center.X() - r_node.X();
        double dy = Center.Y() - r_node.Y();
        double dz = Center.Z() - r_node.Z();

        double tmp = dx*dx + dy*dy + dz*dz;
        Radius = std::max(Radius,tmp);
    }

    Radius = std::sqrt(Radius);
}

/**
 * Desired outcome: It sets the projection of a Gauss node to a condition
 * @param GaussPoint: The origin Gauss Point
 * @param pCandidateCond: The candidate condition
 */
void AdvancedNMPointsMapper::SetProjectionToCond(GaussPointItem& GaussPoint,
                                                 Condition::Pointer pCandidateCond)
{
    array_1d<double,2> LocalCoords;
    double Dist;
    GaussPoint.Project(pCandidateCond, LocalCoords, Dist);
    const unsigned int dimension = (pCandidateCond->GetGeometry()).WorkingSpaceDimension();

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

/**
 * Alternative when no condition is available: It sets the projection of a Gauss point to a node
 * @param GaussPoint: The origin Gauss Point
 * @param pCandidateNode: The candidate node
 * @param Dist: The distance between the node and the Gauss Point
 */
void AdvancedNMPointsMapper::SetProjectionToNode(GaussPointItem& GaussPoint,
                                                 Node<3>::Pointer pCandidateNode,
                                                 const double& Dist)
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

/**
 * It maps a scalar variable to a normal vector from a model part to other
 * @param rOriginVar: The original value (scalar) of the variable
 * @param rDestVar: The variable (normal vector) in the destiny modelpart
 * @param MaxIter: Maximum number of iteration allowed
 * @param TolIter: Tolerance accepted in the iteration
 * @param sign_pos: Positive or negative projection
 */
void AdvancedNMPointsMapper::ScalarToNormalVectorMap(const Variable<double> & rOriginVar,
                                                     const Variable<array_1d<double,3> >& rDestVar,
                                                     const int MaxIter,
                                                     const double TolIter,
                                                     const bool sign_pos)
{
    const unsigned int dimension = mrDestinationModelPart.ConditionsBegin()->GetGeometry().WorkingSpaceDimension();

    // Define if the mapping swaps the sign of the variable values
    const double sign = (sign_pos == false) ? -1.0 : 1.0;

    // Initialize destination variable
    VariableUtils().SetToZero_VectorVar(rDestVar, mrDestinationModelPart.Nodes());

    // Compute nodal lengths/areas (NODAL_MAUX) in both origin and destination modelparts
    ComputeNodalLengthArea();

    // Compute the nodal unit normal vector in the destination modelpart
    #pragma omp parallel for
    for (int k=0; k<static_cast<int>(mrDestinationModelPart.NumberOfNodes()); ++k)
    {
        ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin() + k;
        array_1d<double, 3>& node_normal = node_it->GetValue(NORMAL);
        node_normal /= norm_2(node_normal);
    }

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

        std::vector< Kratos::shared_ptr<array_1d<double,2> > > pInterpValues;

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
            array_1d<double,3> NormalVector = ZeroVector(3);
            mGaussPointList[GPi]->GetNormal(NormalVector);

            // Store the condition normal in each node to compute the nodal normal
            cond_it->GetGeometry()[0].GetValue(NORMAL) += NormalVector;
            cond_it->GetGeometry()[1].GetValue(NORMAL) += NormalVector;

            // Store in GPValues the projected value in the destiny condition Gauss Points
            // Note that currently the implementation is valid for only 2 GP
            double TempValue = 0.0;

            mGaussPointList[GPi]->GetProjectedValue(rOriginVar, TempValue);
            GPValues[0] = TempValue;

            mGaussPointList[GPi + 1]->GetProjectedValue(rOriginVar, TempValue);
            GPValues[1] = TempValue;

            const double K = CondLength/2.0;

            NodalValues[0] = K*(MInterp(0,0)*GPValues[0] + MInterp(0,1)*GPValues[1]);
            NodalValues[1] = K*(MInterp(1,0)*GPValues[0] + MInterp(1,1)*GPValues[1]);

            Kratos::shared_ptr< array_1d<double,2> > pNodalValues(new array_1d<double,2>(NodalValues) );
            pInterpValues.push_back(pNodalValues); // This is computed here because it is the constant part of RHS

            GPi += 2; // 1 Condition = 2 Gauss Points
        }

        // Iteration
        for (int k = 0; k < MaxIter; k++)
        {
            // At the begining of each iteration initialize the variable containing the assembled RHS as 0
            #pragma omp parallel for
            for (int k=0; k < static_cast<int>(mrDestinationModelPart.NumberOfNodes()); ++k)
            {
                ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin() + k;
                node_it->SetValue(MAPPER_SCALAR_PROJECTION_RHS, 0.0);
            }

            double LocalRHS0, LocalRHS1; // Local RHS for each node
            array_1d<double, 2> LastSolution;
            int IV_iter = 0;

            for (ModelPart::ConditionIterator cond_it = mrDestinationModelPart.ConditionsBegin();
                    cond_it != mrDestinationModelPart.ConditionsEnd();
                    cond_it++)
            {
                double CondLength = 0.0;
                GeometryType& rGeom = cond_it->GetGeometry();
                mGaussPointList[2 * IV_iter]->GetArea(CondLength);

                LocalRHS0 = 0.0;
                LocalRHS1 = 0.0;

                const array_1d<double,3> NormalVector0 = rGeom[0].GetValue(NORMAL);
                const array_1d<double,3> NormalVector1 = rGeom[1].GetValue(NORMAL);

                LastSolution[0] = sign * inner_prod(rGeom[0].FastGetSolutionStepValue(rDestVar), NormalVector0);
                LastSolution[1] = sign * inner_prod(rGeom[1].FastGetSolutionStepValue(rDestVar), NormalVector1);

                const double K = CondLength/6.0;
                array_1d<double,2> CondValues = *pInterpValues[IV_iter];

                LocalRHS0 = CondValues[0] - K * (2.0*LastSolution[0] + 1.0*LastSolution[1]);
                LocalRHS1 = CondValues[1] - K * (1.0*LastSolution[0] + 2.0*LastSolution[1]);
                // We are taking advantage of 1 Condition = 2 Gauss Points to iterate the interpolation results and the Gauss Point Vector Again

                rGeom[0].GetValue(MAPPER_SCALAR_PROJECTION_RHS) += LocalRHS0;
                rGeom[1].GetValue(MAPPER_SCALAR_PROJECTION_RHS) += LocalRHS1;

                IV_iter++;
            }

            // Solve
            double dValNorm      = 0.0;
            double ValNorm       = 0.0;
            const unsigned int NodeNum = mrDestinationModelPart.NumberOfNodes();

            #pragma omp parallel for reduction(+ : dValNorm, ValNorm)
            for (int i=0; i<static_cast<int>(mrDestinationModelPart.NumberOfNodes()); ++i)
            {
                ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin() + i;
                const array_1d<double,3> NormalVector = node_it->GetValue(NORMAL);
                const double & NodeLength = node_it->GetValue(NODAL_MAUX);
                double dVal = node_it->GetValue(MAPPER_SCALAR_PROJECTION_RHS)/NodeLength;
                node_it->FastGetSolutionStepValue(rDestVar) += sign * dVal * NormalVector;

                // Variables for convergence check
                dValNorm += dVal * dVal;
                ValNorm  += std::pow(inner_prod(node_it->FastGetSolutionStepValue(rDestVar), NormalVector), 2);
            }
            //std::cout << "ScalarToNormalVectorMap iteration: " << k+1 << std::endl;

            // Check Convergence
            const double RelativeError = (ValNorm > 10e-15) ? dValNorm / ValNorm : 0.0;
            if( (ValNorm/NodeNum < 0.00001 * std::pow(TolIter, 2.00)) || (RelativeError < std::pow(TolIter, 2.00)) )
            {
                // std::cout << "ScalarToNormalVectorMap converged in " << k + 1 << " iterations." << std::endl;
                break;
            }
            else if ((k + 1) == MaxIter)
            {
                KRATOS_WARNING("AdvancedNMPointsMapper") << "ScalarToNormalVectorMap did not converge in " << k + 1 << " iterations." << std::endl;
            }
        }
    }
    else // 3D case
    {
        // Define some variables that will be used in the iteration
        MatrixVar MCons; // Elemental Consistent Mass Matrix = Aelem/12 * MCons
        MCons(0, 0) = 2.0; MCons(0, 1) = 1.0; MCons(0, 2) = 1.0;
        MCons(1, 0) = 1.0; MCons(1, 1) = 2.0; MCons(1, 2) = 1.0;
        MCons(2, 0) = 1.0; MCons(2, 1) = 1.0; MCons(2, 2) = 2.0;

        MatrixVar MInterp; // Interpolation Matrix (NodalValues = (A/24)*MInterp*GaussValues)
        MInterp(0, 0) = 6.0; MInterp(0, 1) = 1.0; MInterp(0, 2) = 1.0;
        MInterp(1, 0) = 1.0; MInterp(1, 1) = 6.0; MInterp(1, 2) = 1.0;
        MInterp(2, 0) = 1.0; MInterp(2, 1) = 1.0; MInterp(2, 2) = 6.0;

        std::vector< Kratos::shared_ptr< array_1d<double, 3> > > pInterpValues;

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

            array_1d<double,3> NormalVector = ZeroVector(3);
            mGaussPointList[GPi]->GetNormal(NormalVector);

            for (unsigned int i = 0; i < 3; i++)
            {
                //~ cond_it->GetGeometry()[i].GetValue(NODAL_MAUX) += 0.333333333333333333 * CondArea;
                mGaussPointList[GPi + i]->GetProjectedValue(rOriginVar, GPValues[i]);

                cond_it->GetGeometry()[i].GetValue(NORMAL) += NormalVector; // Store the condition normal in each node to compute the nodal normal
            }

            noalias(NodalValues) = (CondArea/24.0) * prod(MInterp, GPValues);

            Kratos::shared_ptr< array_1d<double, 3> > pNodalValues(new array_1d<double, 3>(NodalValues) );
            pInterpValues.push_back(pNodalValues); // This is computed here because it is the constant part of RHS

            GPi += 3; // 1 Condition = 3 Gauss Points
        }

        // Iteration
        for (int k = 0; k < MaxIter; k++)
        {
            // At the begining of each iteration initialize the variable containing the assembled RHS as 0
            #pragma omp parallel for
            for (int k=0; k<static_cast<int>(mrDestinationModelPart.NumberOfNodes()); ++k)
            {
                ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin() + k;
                node_it->SetValue(MAPPER_SCALAR_PROJECTION_RHS, 0.0);
            }

            array_1d<double, 3> LocalRHS;
            array_1d<double, 3> LastSolution;
            int IV_iter = 0;

            for (ModelPart::ConditionsContainerType::iterator cond_it = mrDestinationModelPart.ConditionsBegin();
                    cond_it != mrDestinationModelPart.ConditionsEnd();
                    cond_it++)
            {
                GeometryType& rGeom = cond_it->GetGeometry();

                for(unsigned int j = 0; j < 3; j++)
                {
                    const array_1d<double,3> NormalVector = rGeom[j].GetValue(NORMAL);
                    LocalRHS[j] = 0.0;
                    LastSolution[j] = sign * inner_prod(rGeom[j].FastGetSolutionStepValue(rDestVar), NormalVector);
                }

                double CondArea = 0.0;
                mGaussPointList[3 * IV_iter]->GetArea(CondArea);

                noalias(LocalRHS) = *pInterpValues[IV_iter] - (CondArea/12.0) * prod(MCons, LastSolution);
                // We are taking advantage of 1 Condition = 3 Gauss Points to iterate the interpolation results and the Gauss Point Vector Again

                for (unsigned int j = 0; j < 3 ; j++)
                {
                    rGeom[j].GetValue(MAPPER_SCALAR_PROJECTION_RHS) += LocalRHS[j];
                }

                IV_iter++;
            }

            // Solve
            double dValNorm      = 0.0;
            double ValNorm       = 0.0;
            const unsigned int NodeNum = mrDestinationModelPart.NumberOfNodes();

            #pragma omp parallel for reduction(+ : dValNorm,ValNorm)
            for (int i=0; i<static_cast<int>(mrDestinationModelPart.NumberOfNodes()); ++i)
            {
                ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin() + i;
                const array_1d<double,3> NormalVector = node_it->GetValue(NORMAL);
                const double NodeArea = node_it->GetValue(NODAL_MAUX);
                double dVal = node_it->GetValue(MAPPER_SCALAR_PROJECTION_RHS)/NodeArea;
                node_it->FastGetSolutionStepValue(rDestVar) += sign * dVal * NormalVector;

                // Variables for convergence check
                dValNorm += dVal * dVal;
                ValNorm  += std::pow(inner_prod(node_it->FastGetSolutionStepValue(rDestVar), NormalVector), 2);
            }
            //std::cout << "ScalarToNormalVectorMap iteration: " << k+1 << std::endl;

            // Check Convergence
            const double RelativeError = (ValNorm > 10e-15) ? dValNorm / ValNorm : 0.0;
            if( (ValNorm/NodeNum < 0.00001 * std::pow(TolIter, 2.00)) || (RelativeError < std::pow(TolIter, 2.00)) )
            {
                // std::cout << "ScalarToNormalVectorMap converged in " << k + 1 << " iterations." << std::endl;
                break;
            }
            else if ((k + 1) == MaxIter)
            {
                KRATOS_WARNING("AdvancedNMPointsMapper") << "ScalarToNormalVectorMap did not converge in " << k + 1 << " iterations." << std::endl;
            }
        } // End of Iteration
    }
} // End of Map (scalar to normal vector version)

/**
 * It maps a normal vector variable to a scalar from a model part to other
 * @param rOriginVar: The original value (normal vector) of the variable
 * @param rDestVar: The variable (scalar) in the destiny modelpart
 * @param MaxIter: Maximum number of iteration allowed
 * @param TolIter: Tolerance accepted in the iteration
 * @param sign_pos: Positive or negative projection
 */
void AdvancedNMPointsMapper::NormalVectorToScalarMap(const Variable<array_1d<double,3> >& rOriginVar,
                                                     const Variable<double> & rDestVar,
                                                     const int MaxIter,
                                                     const double TolIter,
                                                     const bool sign_pos)
{
    const unsigned int dimension = mrDestinationModelPart.ConditionsBegin()->GetGeometry().WorkingSpaceDimension();

    // Define if the mapping swaps the sign of the variable values
    const double sign = (sign_pos == false) ? -1.0 : 1.0;

    // Initialize destination variable
    VariableUtils().SetToZero_ScalarVar(rDestVar, mrDestinationModelPart.Nodes());

    // Compute nodal lengths/areas (NODAL_MAUX) in both origin and destination modelparts
    ComputeNodalLengthArea();


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

        std::vector< Kratos::shared_ptr<array_1d<double,6> > > pInterpValues;

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
            array_1d<double,3> NormalVector = ZeroVector(3);
            mGaussPointList[GPi]->GetNormal(NormalVector);

            // Store the condition normal in each node to compute the nodal normal
            cond_it->GetGeometry()[0].GetValue(NORMAL) += NormalVector;
            cond_it->GetGeometry()[1].GetValue(NORMAL) += NormalVector;

            // Store in GPValues the projected value in the destiny condition Gauss Points
            // Note that currently the implementation is valid for only 2 GP
            array_1d<double, 3> TempValue = ZeroVector(3);

            for (unsigned int j = 0; j < 3; j++)
            {
                mGaussPointList[GPi]->GetProjectedValue(rOriginVar, TempValue);
                GPValues[j] = TempValue[j];

                mGaussPointList[GPi + 1]->GetProjectedValue(rOriginVar, TempValue);
                GPValues[j + 3] = TempValue[j];
            }

            const double K = CondLength/2.0;

            for (unsigned int j = 0; j < 3; j++)
            {
                NodalValues[j] = K*(MInterp(0,0)*GPValues[j] + MInterp(0,1)*GPValues[3 + j]);
                NodalValues[3 + j] = K*(MInterp(1,0)*GPValues[j] + MInterp(1,1)*GPValues[3 + j]);
            }

            Kratos::shared_ptr< array_1d<double,6> > pNodalValues(new array_1d<double,6>(NodalValues) );
            pInterpValues.push_back(pNodalValues); // This is computed here because it is the constant part of RHS

            GPi += 2; // 1 Condition = 2 Gauss Points
        }

        // Compute the nodal normals (normalised sum of the conditions normals) in the destination modelpart
        #pragma omp parallel for
        for (int i=0; i<static_cast<int>(mrDestinationModelPart.NumberOfNodes()); ++i)
        {
            ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin() + i;
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
                node_it->SetValue(MAPPER_VECTOR_PROJECTION_RHS, ZeroVector(3));
            }

            array_1d<double, 3> LocalRHS0, LocalRHS1; // Local RHS for each node
            array_1d<double, 6> LastSolution;
            int IV_iter = 0;

            for (ModelPart::ConditionIterator cond_it = mrDestinationModelPart.ConditionsBegin();
                    cond_it != mrDestinationModelPart.ConditionsEnd();
                    cond_it++)
            {
                double CondLength = 0.0;
                GeometryType& rGeom = cond_it->GetGeometry();
                mGaussPointList[2 * IV_iter]->GetArea(CondLength);

                LocalRHS0 = ZeroVector(3);
                LocalRHS1 = ZeroVector(3);

                array_1d<double, 3> NormalVector0 = rGeom[0].GetValue(NORMAL);
                array_1d<double, 3> NormalVector1 = rGeom[1].GetValue(NORMAL);

                for(unsigned int j = 0; j < 3; j++)
                {
                    LastSolution[j]     = sign * rGeom[0].FastGetSolutionStepValue(rDestVar) * NormalVector0[j];
                    LastSolution[3 + j] = sign * rGeom[1].FastGetSolutionStepValue(rDestVar) * NormalVector1[j];
                }

                const double K = CondLength/6.0;
                array_1d<double,6> CondValues = *pInterpValues[IV_iter];

                for(unsigned int j = 0; j < 3; j++)
                {
                    LocalRHS0[j] = CondValues[j] - K * (2.0*LastSolution[j] + 1.0*LastSolution[3 + j]);
                    LocalRHS1[j] = CondValues[3 + j] - K * (1.0*LastSolution[j] + 2.0*LastSolution[3 + j]);
                }
                // We are taking advantage of 1 Condition = 2 Gauss Points to iterate the interpolation results and the Gauss Point Vector Again

                rGeom[0].GetValue(MAPPER_VECTOR_PROJECTION_RHS) += LocalRHS0;
                rGeom[1].GetValue(MAPPER_VECTOR_PROJECTION_RHS) += LocalRHS1;

                IV_iter++;
            }

            // Solve
            array_1d<double,3> dVal = ZeroVector(3);
            double dValNorm      = 0.0;
            double ValNorm       = 0.0;
            const unsigned int NodeNum = mrDestinationModelPart.NumberOfNodes();

            #pragma omp parallel for reduction(+ : dValNorm,ValNorm)
            for (int i=0; i<static_cast<int>(mrDestinationModelPart.NumberOfNodes()); ++i)
            {
                ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin() + i;
                const array_1d<double,3> NormalVector = node_it->GetValue(NORMAL);
                const double NodeLength = node_it->GetValue(NODAL_MAUX);
                dVal = node_it->GetValue(MAPPER_VECTOR_PROJECTION_RHS)/NodeLength;
                node_it->FastGetSolutionStepValue(rDestVar) += sign * inner_prod(dVal, NormalVector);

                // Variables for convergence check
                for (unsigned int j = 0; j < 3; j++)
                {
                    dValNorm += dVal[j] * dVal[j];
                    ValNorm += std::pow(node_it->FastGetSolutionStepValue(rDestVar) * NormalVector[j], 2);
                }
            }
            //std::cout << "NormalVectorToScalarMap iteration: " << k+1 << std::endl;

            // Check Convergence
            const double RelativeError = (ValNorm > 10e-15) ? dValNorm / ValNorm : 0.0;
            if( (ValNorm/NodeNum < 0.00001 * std::pow(TolIter, 2.00)) || (RelativeError < std::pow(TolIter, 2.00)) )
            {
                // std::cout << "NormalVectorToScalarMap converged in " << k + 1 << " iterations." << std::endl;
                break;
            }
            else if ((k + 1) == MaxIter)
            {
                KRATOS_WARNING("AdvancedNMPointsMapper") << "NormalVectorToScalarMap did not converge in " << k + 1 << " iterations." << std::endl;
            }
        }
    }
    else
    {

        std::vector< Kratos::shared_ptr<array_1d<double,9> > > pInterpValues;

        // Here we loop both the Destination Model Part and the Gauss Point List, using the
        // fact that the GP list was created by a loop over the Dest. Model Part's conditions
        int GPi = 0;
        for (ModelPart::ConditionIterator cond_it = mrDestinationModelPart.ConditionsBegin();
                cond_it != mrDestinationModelPart.ConditionsEnd();
                cond_it++)
        {
            array_1d<double, 9> GPValues;
            array_1d<double, 9> NodalValues;

            const unsigned int ngauss = 3;
            double CondArea = 0.0;
            mGaussPointList[GPi]->GetArea(CondArea);
            const double K = CondArea/24.0;

            array_1d<double,3> NormalVector = ZeroVector(3);
            mGaussPointList[GPi]->GetNormal(NormalVector);

            for (unsigned int i = 0; i < 3; i++)
            {
                cond_it->GetGeometry()[i].GetValue(NORMAL) += NormalVector;
            }

            // Get the Gauss pts. values
            array_1d<double,3> TempValues = ZeroVector(3);

            for (unsigned int j=0; j<ngauss; ++j)
            {
                mGaussPointList[GPi + j]->GetProjectedValue(rOriginVar, TempValues);
                for (unsigned int i = 0; i < 3; i++)
                {
                    GPValues[ngauss*j + i] = TempValues[i];
                }
            }

            // Nodal values from GP projections
            for (unsigned int i = 0; i < 3; i++)
            {
                NodalValues[i]     = K*(6.0 * GPValues[i] + 1.0 * GPValues[3 + i] + 1.0 * GPValues[6 + i]);
                NodalValues[3 + i] = K*(1.0 * GPValues[i] + 6.0 * GPValues[3 + i] + 1.0 * GPValues[6 + i]);
                NodalValues[6 + i] = K*(1.0 * GPValues[i] + 1.0 * GPValues[3 + i] + 6.0 * GPValues[6 + i]);
            }

            Kratos::shared_ptr< array_1d<double,9> > pNodalValues(new array_1d<double,9>(NodalValues) );
            pInterpValues.push_back(pNodalValues); // This is computed here because it is the constant part of RHS

            GPi += 3; // 1 Condition = 3 Gauss Points
        }

        // Compute the destination interface nodal unit normals
        #pragma omp parallel for
        for (int i=0; i<static_cast<int>(mrDestinationModelPart.NumberOfNodes()); ++i)
        {
            ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin() + i;
            const array_1d<double,3> NormalVector = node_it->GetValue(NORMAL);
            node_it->GetValue(NORMAL) = NormalVector/norm_2(NormalVector);
        }

        // Iteration
        for (int k = 0; k < MaxIter; k++)
        {
            // At the begining of each iteration initialize the variable containing the assembled RHS as 0
            #pragma omp parallel for
            for (k=0; k<static_cast<int>(mrDestinationModelPart.NumberOfNodes()); ++k)
            {
                ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin() + k;
                node_it->SetValue(MAPPER_VECTOR_PROJECTION_RHS, ZeroVector(3));
            }

            array_1d<double, 3> LocalRHS0, LocalRHS1, LocalRHS2; // Local RHS for each node
            array_1d<double, 9> LastSolution;
            int IV_iter = 0;

            for (ModelPart::ConditionIterator cond_it = mrDestinationModelPart.ConditionsBegin();
                    cond_it != mrDestinationModelPart.ConditionsEnd();
                    cond_it++)
            {
                double CondArea = 0.0;
                GeometryType& rGeom = cond_it->GetGeometry();
                mGaussPointList[3 * IV_iter]->GetArea(CondArea);
                const double K = CondArea/12.0;

                LocalRHS0 = ZeroVector(3);
                LocalRHS1 = ZeroVector(3);
                LocalRHS2 = ZeroVector(3);

                for(unsigned int j = 0; j < 3; j++)
                {
                    array_1d<double,3> NormalVector = rGeom[0].GetValue(NORMAL);
                    LastSolution[j    ] = sign * rGeom[0].FastGetSolutionStepValue(rDestVar) * NormalVector[j];
                    NormalVector = rGeom[1].GetValue(NORMAL);
                    LastSolution[j + 3] = sign * rGeom[1].FastGetSolutionStepValue(rDestVar) * NormalVector[j];
                    NormalVector = rGeom[2].GetValue(NORMAL);
                    LastSolution[j + 6] = sign * rGeom[2].FastGetSolutionStepValue(rDestVar) * NormalVector[j];
                }

                array_1d<double,9> CondValues = *pInterpValues[IV_iter];

                for(unsigned int j = 0; j < 3; j++)
                {
                    LocalRHS0[j] = CondValues[j]     - K * (2.0 * LastSolution[j] + 1.0 * LastSolution[j + 3] + 1.0 * LastSolution[j + 6]);
                    LocalRHS1[j] = CondValues[j + 3] - K * (1.0 * LastSolution[j] + 2.0 * LastSolution[j + 3] + 1.0 * LastSolution[j + 6]);
                    LocalRHS2[j] = CondValues[j + 6] - K * (1.0 * LastSolution[j] + 1.0 * LastSolution[j + 3] + 2.0 * LastSolution[j + 6]);
                }
                // We are taking advantage of 1 Condition = 3 Gauss Points to iterate the interpolation results and the Gauss Point Vector Again

                rGeom[0].GetValue(MAPPER_VECTOR_PROJECTION_RHS) += LocalRHS0;
                rGeom[1].GetValue(MAPPER_VECTOR_PROJECTION_RHS) += LocalRHS1;
                rGeom[2].GetValue(MAPPER_VECTOR_PROJECTION_RHS) += LocalRHS2;

                IV_iter++;
            }

            // Solve
            array_1d<double,3> dVal = ZeroVector(3);
            double dValNorm      = 0.0;
            double ValNorm       = 0.0;
            const unsigned int NodeNum = mrDestinationModelPart.NumberOfNodes();

            #pragma omp parallel for reduction(+ : dValNorm,ValNorm)
            for (int i=0; i<static_cast<int>(mrDestinationModelPart.NumberOfNodes()); ++i)
            {
                ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin() + i;
                const array_1d<double,3> NormalVector = node_it->GetValue(NORMAL);
                const double NodeArea = node_it->GetValue(NODAL_MAUX);
                dVal = node_it->GetValue(MAPPER_VECTOR_PROJECTION_RHS)/NodeArea;
                node_it->FastGetSolutionStepValue(rDestVar) += sign * inner_prod(dVal, NormalVector);

                // Variables for convergence check
                for (unsigned int j = 0; j < 3; j++)
                {
                    dValNorm += dVal[j] * dVal[j];
                    ValNorm += std::pow(node_it->FastGetSolutionStepValue(rDestVar) * NormalVector[j], 2);
                }
            }
            //std::cout << "NormalVectorToScalarMap iteration: " << k+1 << std::endl;

            // Check Convergence
            const double RelativeError = (ValNorm > 10e-15) ? dValNorm / ValNorm : 0.0;
            if( (ValNorm/NodeNum < 0.00001 * std::pow(TolIter, 2.00)) || (RelativeError < std::pow(TolIter, 2.00)) )
            {
                // std::cout << "NormalVectorToScalarMap converged in " << k + 1 << " iterations." << std::endl;
                break;
            }
            else if ((k + 1) == MaxIter)
            {
                KRATOS_WARNING("AdvancedNMPointsMapper") << "NormalVectorToScalarMap did not converge in " << k + 1 << " iterations." << std::endl;
            }
        } // End of Iteration
    }
} // End of Map (normal vector to scalar version)

/**
 * It maps a variable (scalar) from a model part to other
 * @param rOriginVar: The original value of the variable
 * @param rDestVar: The variable in the destiny modelpart
 * @param MaxIter: Maximum number of iteration allowed
 * @param TolIter: Tolerance accepted in the iteration
 * @param sign_pos: Positive or negative projection
 */
void AdvancedNMPointsMapper::ScalarMap(const Variable<double> & rOriginVar,
                                       const Variable<double> & rDestVar,
                                       const int MaxIter,
                                       const double TolIter,
                                       const bool sign_pos)
{
    const unsigned int dimension = mrDestinationModelPart.ConditionsBegin()->GetGeometry().WorkingSpaceDimension();

    // Define if the mapping swaps the sign of the variable values
    const double sign = (sign_pos == false) ? -1.0 : 1.0;

    // Initialize destination variable
    VariableUtils().SetToZero_ScalarVar(rDestVar, mrDestinationModelPart.Nodes());

    // Compute nodal lengths/areas (NODAL_MAUX) in both origin and destination modelparts
    ComputeNodalLengthArea();

    if (dimension == 2) // 2D case
    {

        // Interpolation matrix obtention
        boost::numeric::ublas::bounded_matrix<double,2,2> MCons; // Elemental Consistent Mass Matrix = L/6 * MCons
        boost::numeric::ublas::bounded_matrix<double,2,2> MInterp;
        boost::numeric::ublas::bounded_matrix<double,2,2> aux_GPPos;
        boost::numeric::ublas::bounded_matrix<double,2,2> inv_aux_GPPos;

        MCons(0, 0) = 2.0/3.0; MCons(0, 1) = 1.0/3.0;
        MCons(1, 0) = 1.0/3.0; MCons(1, 1) = 2.0/3.0;

        aux_GPPos(0, 0) = (1.0+(-std::sqrt(1.0/3.0)))/2.0; aux_GPPos(0, 1) = (1.0-(-std::sqrt(1.0/3.0)))/2.0;
        aux_GPPos(1, 0) = (1.0+(+std::sqrt(1.0/3.0)))/2.0; aux_GPPos(1, 1) = (1.0-(+std::sqrt(1.0/3.0)))/2.0;

        const double det_aux_GPPos = (aux_GPPos(0,0)*aux_GPPos(1,1))-(aux_GPPos(0,1)*aux_GPPos(1,0));

        inv_aux_GPPos(0, 0) = aux_GPPos(1, 1)/det_aux_GPPos;  inv_aux_GPPos(0, 1) = -aux_GPPos(0, 1)/det_aux_GPPos;
        inv_aux_GPPos(1, 0) = -aux_GPPos(1, 0)/det_aux_GPPos; inv_aux_GPPos(1, 1) = aux_GPPos(0, 0)/det_aux_GPPos;

        MInterp = prod(MCons, inv_aux_GPPos); // Interpolation Matrix (NodalValues = (L/6)*MInterp*GaussValues)

        std::vector< Kratos::shared_ptr<array_1d<double,2> > > pInterpValues;

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

            mGaussPointList[GPi]->GetProjectedValue(rOriginVar, TempValue);
            GPValues[0] = TempValue;

            mGaussPointList[GPi + 1]->GetProjectedValue(rOriginVar, TempValue);
            GPValues[1] = TempValue;

            NodalValues[0] = K*(MInterp(0,0)*GPValues[0] + MInterp(0,1)*GPValues[1]);
            NodalValues[1] = K*(MInterp(1,0)*GPValues[0] + MInterp(1,1)*GPValues[1]);

            Kratos::shared_ptr< array_1d<double,2> > pNodalValues(new array_1d<double,2>(NodalValues));
            pInterpValues.push_back(pNodalValues); // This is computed here because it is the constant part of RHS

            GPi += 2; // 1 Condition = 2 Gauss Points
        }

        // Iteration
        for (int k = 0; k < MaxIter; k++)
        {
            // At the begining of each iteration initialize the variable containing the assembled RHS as 0
            #pragma omp parallel for
            for (int i=0; i<static_cast<int>(mrDestinationModelPart.NumberOfNodes()); ++i)
            {
                ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin() + i;
                node_it->SetValue(MAPPER_SCALAR_PROJECTION_RHS, 0.0);
            }

            double LocalRHS0, LocalRHS1; // Local RHS for each node
            array_1d<double, 2> LastSolution;
            int IV_iter = 0;

            for (ModelPart::ConditionIterator cond_it = mrDestinationModelPart.ConditionsBegin();
                    cond_it != mrDestinationModelPart.ConditionsEnd();
                    cond_it++)
            {
                double CondLength = 0.0;
                GeometryType& rGeom = cond_it->GetGeometry();
                mGaussPointList[2 * IV_iter]->GetArea(CondLength);
                const double K = CondLength/6.0;

                LocalRHS0 = 0.0;
                LocalRHS1 = 0.0;

                LastSolution[0] = sign * rGeom[0].FastGetSolutionStepValue(rDestVar);
                LastSolution[1] = sign * rGeom[1].FastGetSolutionStepValue(rDestVar);

                array_1d<double,2> CondValues = *pInterpValues[IV_iter];

                LocalRHS0 = CondValues[0] - K * (2.0*LastSolution[0] + 1.0*LastSolution[1]);
                LocalRHS1 = CondValues[1] - K * (1.0*LastSolution[0] + 2.0*LastSolution[1]);

                // We are taking advantage of 1 Condition = 2 Gauss Points to iterate the interpolation results and the Gauss Point Vector Again
                rGeom[0].GetValue(MAPPER_SCALAR_PROJECTION_RHS) += LocalRHS0;
                rGeom[1].GetValue(MAPPER_SCALAR_PROJECTION_RHS) += LocalRHS1;

                IV_iter++;
            }

            // Solve
            double dValNorm      = 0.0;
            double ValNorm       = 0.0;
            const unsigned int NodeNum = mrDestinationModelPart.NumberOfNodes();

            #pragma omp parallel for reduction(+ : dValNorm, ValNorm)
            for (int i=0; i<static_cast<int>(mrDestinationModelPart.NumberOfNodes()); ++i)
            {
                ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin() + i;
                const double NodeLength = node_it->GetValue(NODAL_MAUX);
                const double dVal = (node_it->GetValue(MAPPER_SCALAR_PROJECTION_RHS))/NodeLength;
                double& rDestinationValue = node_it->FastGetSolutionStepValue(rDestVar);
                rDestinationValue += sign * dVal;

                // Variables for convergence check
                dValNorm += dVal;
                ValNorm  += rDestinationValue*rDestinationValue;
            }

            // Check Convergence
            const double RelativeError = (ValNorm > 10e-15) ? dValNorm / ValNorm : 0.0;
            if( (ValNorm/NodeNum < 0.00001 * TolIter * TolIter) || RelativeError < TolIter * TolIter)
            {
                // std::cout << "ScalarMap converged in " << k + 1 << " iterations." << std::endl;
                break;
            }
            else if ( (k + 1) == MaxIter)
            {
                KRATOS_WARNING("AdvancedNMPointsMapper") << "VectorMap did not converge in " << k + 1 << " iterations." << std::endl;
            }
        } // End of Iteration
    }

    else // 3D case
    {
        // Define some variables that will be used in the iteration
        MatrixVar MCons; // Elemental Consistent Mass Matrix = Aelem/12 * MCons
        MCons(0, 0) = 2.0; MCons(0, 1) = 1.0; MCons(0, 2) = 1.0;
        MCons(1, 0) = 1.0; MCons(1, 1) = 2.0; MCons(1, 2) = 1.0;
        MCons(2, 0) = 1.0; MCons(2, 1) = 1.0; MCons(2, 2) = 2.0;

        MatrixVar MInterp; // Interpolation Matrix (NodalValues = (A/24)*MInterp*GaussValues)
        MInterp(0, 0) = 6.0; MInterp(0, 1) = 1.0; MInterp(0, 2) = 1.0;
        MInterp(1, 0) = 1.0; MInterp(1, 1) = 6.0; MInterp(1, 2) = 1.0;
        MInterp(2, 0) = 1.0; MInterp(2, 1) = 1.0; MInterp(2, 2) = 6.0;

        std::vector< Kratos::shared_ptr< array_1d<double, 3> > > pInterpValues;

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
                mGaussPointList[GPi + i]->GetProjectedValue(rOriginVar, GPValues[i]);
            }

            noalias(NodalValues) = (CondArea/24.0) * prod(MInterp, GPValues);

            Kratos::shared_ptr< array_1d<double, 3> > pNodalValues(new array_1d<double, 3>(NodalValues) );
            pInterpValues.push_back(pNodalValues); // This is computed here because it is the constant part of RHS

            GPi += 3; // 1 Condition = 3 Gauss Points
        }

        // Iteration
        for (int k = 0; k < MaxIter; k++)
        {
            // At the begining of each iteration initialize the variable containing the assembled RHS as 0
            #pragma omp parallel for
            for (int i=0; i<static_cast<int>(mrDestinationModelPart.NumberOfNodes()); ++i)
            {
                ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin() + i;
                node_it->SetValue(MAPPER_SCALAR_PROJECTION_RHS, 0.0);
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
            const unsigned int NodeNum = mrDestinationModelPart.NumberOfNodes();

            #pragma omp parallel for reduction(+ : dValNorm, ValNorm)
            for (int i=0; i<static_cast<int>(mrDestinationModelPart.NumberOfNodes()); ++i)
            {
                ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin() + i;
                const double NodeArea = node_it->GetValue(NODAL_MAUX);
                const double dVal = node_it->GetValue(MAPPER_SCALAR_PROJECTION_RHS)/NodeArea;
                double& rDestinationValue = node_it->FastGetSolutionStepValue(rDestVar);
                rDestinationValue += sign * dVal;

                // Variables for convergence check
                dValNorm += dVal * dVal;
                ValNorm  += rDestinationValue*rDestinationValue;

            }
            //std::cout << "ScalarMap iteration: " << k+1 << std::endl;

            // Check Convergence
            const double RelativeError = (ValNorm > 10e-15) ? dValNorm / ValNorm : 0.0;
            if( (ValNorm/NodeNum < 0.00001 * std::pow(TolIter, 2.00)) || (RelativeError < std::pow(TolIter, 2.00)) )
            {
                // std::cout << "ScalarMap converged in " << k + 1 << " iterations." << std::endl;
                break;
            }
            else if ((k + 1) == MaxIter)
            {
                KRATOS_WARNING("AdvancedNMPointsMapper") << "ScalarMap did not converge in " << k + 1 << " iterations." << std::endl;
            }
        } // End of Iteration
    }

} // End of Map (scalar version)

/**
 * It maps a variable (vector) from a model part to other
 * @param rOriginVar: The original value of the variable
 * @param rDestVar: The variable in the destiny modelpart
 * @param MaxIter: Maximum number of iteration allowed
 * @param TolIter: Tolerance accepted in the iteration
 * @param sign_pos: Positive or negative projection
 */
void AdvancedNMPointsMapper::VectorMap(const Variable<array_1d<double,3> >& rOriginVar,
                                       const Variable<array_1d<double,3> >& rDestVar,
                                       const int MaxIter,
                                       const double TolIter,
                                       const bool sign_pos)
{
    const unsigned int dimension = mrDestinationModelPart.ConditionsBegin()->GetGeometry().WorkingSpaceDimension();

    // Define if the mapping swaps the sign of the variable values
    const double sign = (sign_pos == false) ? -1.0 : 1.0;

    // Initialize destination variable
    VariableUtils().SetToZero_VectorVar(rDestVar, mrDestinationModelPart.Nodes());

    // Compute nodal lengths/areas (NODAL_MAUX) in both origin and destination modelparts
    ComputeNodalLengthArea();

    if (dimension == 2) // 2D case
    {
        // Interpolation matrix obtention
        boost::numeric::ublas::bounded_matrix<double,2,2> MCons; // Elemental Consistent Mass Matrix = L/2 * MCons
        boost::numeric::ublas::bounded_matrix<double,2,2> MInterp;
        boost::numeric::ublas::bounded_matrix<double,2,2> aux_GPPos;
        boost::numeric::ublas::bounded_matrix<double,2,2> inv_aux_GPPos;

        MCons(0, 0) = 2.0/3.0; MCons(0, 1) = 1.0/3.0;
        MCons(1, 0) = 1.0/3.0; MCons(1, 1) = 2.0/3.0;

        aux_GPPos(0, 0) = (1.0+(-std::sqrt(1.0/3.0)))/2.0; aux_GPPos(0, 1) = (1.0-(-std::sqrt(1.0/3.0)))/2.0;
        aux_GPPos(1, 0) = (1.0+(+std::sqrt(1.0/3.0)))/2.0; aux_GPPos(1, 1) = (1.0-(+std::sqrt(1.0/3.0)))/2.0;

        double det_aux_GPPos;
        det_aux_GPPos = (aux_GPPos(0,0)*aux_GPPos(1,1))-(aux_GPPos(0,1)*aux_GPPos(1,0));

        inv_aux_GPPos(0, 0) = aux_GPPos(1, 1)/det_aux_GPPos;  inv_aux_GPPos(0, 1) = -aux_GPPos(0, 1)/det_aux_GPPos;
        inv_aux_GPPos(1, 0) = -aux_GPPos(1, 0)/det_aux_GPPos; inv_aux_GPPos(1, 1) = aux_GPPos(0, 0)/det_aux_GPPos;

        MInterp = prod(MCons, inv_aux_GPPos); // Interpolation matrix (NodalValues = (L/6)*MInterp*GaussValues)

        std::vector< Kratos::shared_ptr<array_1d<double,6> > > pInterpValues;

        // Here we loop both the Destination Model Part and the Gauss Point List, using the
        // fact that the GP list was created by a loop over the Dest. Model Part's conditions
        int GPi = 0;
        for (ModelPart::ConditionIterator cond_it = mrDestinationModelPart.ConditionsBegin();
                cond_it != mrDestinationModelPart.ConditionsEnd();
                cond_it++)
        {

            const unsigned int ngauss = 2;
            double CondLength = 0.0;
            mGaussPointList[GPi]->GetArea(CondLength); // Gets the length of the parent condition of the two points considered
            const double K = CondLength/2.0;

            // Store in GPValues the projected value in the destiny condition Gauss Points
            // Note that currently the implementation is valid for only 2 GP
            array_1d<double, 6> GPValues = ZeroVector(6);
            array_1d<double, 6> NodalValues = ZeroVector(6);
            array_1d<double, 3> TempValues = ZeroVector(3);

            // Get the Gauss pts. values
            for (unsigned int i=0; i<ngauss; ++i)
            {
                mGaussPointList[GPi + i]->GetProjectedValue(rOriginVar, TempValues);
                for (unsigned int j = 0; j < 3; j++)
                {
                     GPValues[3*i + j] = TempValues[j];
                }
            }

            // Compute the nodal values from the projected ones using the interpolation matrix
            NodalValues[0] = K*(MInterp(0,0)*GPValues[0] + MInterp(0,1)*GPValues[3]);
            NodalValues[1] = K*(MInterp(0,0)*GPValues[1] + MInterp(0,1)*GPValues[4]);
            NodalValues[2] = 0.0;
            NodalValues[3] = K*(MInterp(1,0)*GPValues[0] + MInterp(1,1)*GPValues[3]);
            NodalValues[4] = K*(MInterp(1,0)*GPValues[1] + MInterp(1,1)*GPValues[4]);
            NodalValues[5] = 0.0;

            Kratos::shared_ptr< array_1d<double,6> > pNodalValues(new array_1d<double,6>(NodalValues) );
            pInterpValues.push_back(pNodalValues); // This is computed here because it is the constant part of RHS

            GPi += 2; // 1 Condition = 2 Gauss Points
        }

        // Iteration
        for (int k = 0; k < MaxIter; k++)
        {
            // At the begining of each iteration initialize the variable containing the assembled RHS as 0
            #pragma omp parallel for
            for (int i=0; i<static_cast<int>(mrDestinationModelPart.NumberOfNodes()); ++i)
            {
                ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin() + i;
                node_it->SetValue(MAPPER_VECTOR_PROJECTION_RHS, ZeroVector(3));
            }

            array_1d<double, 3> LocalRHS0, LocalRHS1; // Local RHS for each node
            array_1d<double, 6> LastSolution;
            int IV_iter = 0;

            for (ModelPart::ConditionIterator cond_it = mrDestinationModelPart.ConditionsBegin();
                    cond_it != mrDestinationModelPart.ConditionsEnd();
                    cond_it++)
            {
                double CondLength = 0.0;
                GeometryType& rGeom = cond_it->GetGeometry();
                mGaussPointList[2 * IV_iter]->GetArea(CondLength);
                const double K = CondLength/6.0;

                LocalRHS0 = ZeroVector(3);
                LocalRHS1 = ZeroVector(3);

                for(unsigned int j = 0; j < 3; j++)
                {
                    LastSolution[j]     = sign * rGeom[0].FastGetSolutionStepValue(rDestVar)[j];
                    LastSolution[3 + j] = sign * rGeom[1].FastGetSolutionStepValue(rDestVar)[j];
                }

                array_1d<double,6> CondValues = *pInterpValues[IV_iter];

                for(unsigned int j = 0; j < 3; j++)
                {
                    LocalRHS0[j] = CondValues[j]     - K * (2.0*LastSolution[j] + 1.0*LastSolution[j + 3]);
                    LocalRHS1[j] = CondValues[j + 3] - K * (1.0*LastSolution[j] + 2.0*LastSolution[j + 3]);
                }
                // We are taking advantage of 1 Condition = 2 Gauss Points to iterate the interpolation results and the Gauss Point Vector Again

                rGeom[0].GetValue(MAPPER_VECTOR_PROJECTION_RHS) += LocalRHS0;
                rGeom[1].GetValue(MAPPER_VECTOR_PROJECTION_RHS) += LocalRHS1;

                IV_iter++;
            }
            
            // Solve
            array_1d<double,3> dVal = ZeroVector(3);
            double dValNorm      = 0.0;
            double ValNorm       = 0.0;
            const unsigned int NodeNum = mrDestinationModelPart.NumberOfNodes();

            #pragma omp parallel for reduction(+ : dValNorm, ValNorm) private(dVal)
            for (int i=0; i<static_cast<int>(mrDestinationModelPart.NumberOfNodes()); ++i)
            {
                ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin() + i;
                const double & NodeLength = node_it->GetValue(NODAL_MAUX);
                noalias(dVal) = node_it->GetValue(MAPPER_VECTOR_PROJECTION_RHS)/NodeLength;
                array_1d<double, 3>& rDestinationValue = node_it->FastGetSolutionStepValue(rDestVar);
                rDestinationValue += sign * dVal;

                // Variables for convergence check
                for (unsigned int j = 0; j < 3; j++)
                {
                    dValNorm += dVal[j] * dVal[j];
                    ValNorm += rDestinationValue[j] * rDestinationValue[j];
                }
            }

            // Check Convergence
            const double RelativeError = (ValNorm > 10e-15) ? dValNorm / ValNorm : 0.0;
            if( (ValNorm/NodeNum < 0.00001 * TolIter * TolIter) || RelativeError < TolIter * TolIter)
            {
                // std::cout << "VectorMap converged in " << k + 1 << " iterations." << std::endl;
                break;
            }
            else if ( (k + 1) == MaxIter)
            {
                KRATOS_WARNING("AdvancedNMPointsMapper") << "VectorMap did not converge in " << k + 1 << " iterations." << std::endl;
            }
        }
        // End of Iteration

    }
    else // 3D case
    {
        std::vector< Kratos::shared_ptr<array_1d<double,9> > > pInterpValues;

        // Here we loop both the Destination Model Part and the Gauss Point List, using the
        // fact that the GP list was created by a loop over the Dest. Model Part's conditions
        int GPi = 0;
        for (ModelPart::ConditionIterator cond_it = mrDestinationModelPart.ConditionsBegin();
                cond_it != mrDestinationModelPart.ConditionsEnd();
                cond_it++)
        {
            array_1d<double, 9> GPValues;
            array_1d<double, 9> NodalValues;

            const unsigned int ngauss = 3;
            double CondArea = 0.0;
            mGaussPointList[GPi]->GetArea(CondArea);
            const double K = CondArea/24.0;

            array_1d<double,3> TempValues = ZeroVector(3);

            // Get the Gauss pts. values
            for (unsigned int i=0; i<ngauss; ++i)
            {
                mGaussPointList[GPi + i]->GetProjectedValue(rOriginVar, TempValues);
                for (unsigned int j = 0; j < 3; j++)
                {
                     GPValues[3*i + j] = TempValues[j];
                }
            }

            // Compute the nodal values from the projected ones using the interpolation matrix
            for (unsigned int i = 0; i < 3; i++)
            {
                NodalValues[i]     = K*(6.0 * GPValues[i] + 1.0 * GPValues[3 + i] + 1.0 * GPValues[6 + i]);
                NodalValues[3 + i] = K*(1.0 * GPValues[i] + 6.0 * GPValues[3 + i] + 1.0 * GPValues[6 + i]);
                NodalValues[6 + i] = K*(1.0 * GPValues[i] + 1.0 * GPValues[3 + i] + 6.0 * GPValues[6 + i]);
            }

            Kratos::shared_ptr< array_1d<double,9> > pNodalValues(new array_1d<double,9>(NodalValues) );
            pInterpValues.push_back(pNodalValues); // This is computed here because it is the constant part of RHS

            GPi += 3; // 1 Condition = 3 Gauss Points
        }

        // Iteration
        for (int k = 0; k < MaxIter; k++)
        {
            // At the begining of each iteration initialize the variable containing the assembled RHS as 0
            #pragma omp parallel for
            for (k=0; k<static_cast<int>(mrDestinationModelPart.NumberOfNodes()); ++k)
            {
                ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin() + k;
                node_it->SetValue(MAPPER_VECTOR_PROJECTION_RHS, ZeroVector(3));
            }

            array_1d<double, 3> LocalRHS0, LocalRHS1, LocalRHS2; // Local RHS for each node
            array_1d<double, 9> LastSolution;
            int IV_iter = 0;

            for (ModelPart::ConditionIterator cond_it = mrDestinationModelPart.ConditionsBegin();
                    cond_it != mrDestinationModelPart.ConditionsEnd();
                    cond_it++)
            {
                double CondArea = 0.0;
                GeometryType& rGeom = cond_it->GetGeometry();
                mGaussPointList[3 * IV_iter]->GetArea(CondArea);
                const double K = CondArea/12.0;

                LocalRHS0 = ZeroVector(3);
                LocalRHS1 = ZeroVector(3);
                LocalRHS2 = ZeroVector(3);

                for(unsigned int j = 0; j < 3; j++)
                {
                    LastSolution[j]     = sign * rGeom[0].FastGetSolutionStepValue(rDestVar)[j];
                    LastSolution[3 + j] = sign * rGeom[1].FastGetSolutionStepValue(rDestVar)[j];
                    LastSolution[6 + j] = sign * rGeom[2].FastGetSolutionStepValue(rDestVar)[j];
                }

                array_1d<double,9> CondValues = *pInterpValues[IV_iter];

                for(unsigned int j = 0; j < 3; j++)
                {
                    LocalRHS0[j] = CondValues[j]     - K * (2.0 * LastSolution[j] + 1.0 * LastSolution[j + 3] + 1.0 * LastSolution[j + 6]);
                    LocalRHS1[j] = CondValues[j + 3] - K * (1.0 * LastSolution[j] + 2.0 * LastSolution[j + 3] + 1.0 * LastSolution[j + 6]);
                    LocalRHS2[j] = CondValues[j + 6] - K * (1.0 * LastSolution[j] + 1.0 * LastSolution[j + 3] + 2.0 * LastSolution[j + 6]);
                }
                // We are taking advantage of 1 Condition = 3 Gauss Points to iterate the interpolation results and the Gauss Point Vector Again

                rGeom[0].GetValue(MAPPER_VECTOR_PROJECTION_RHS) += LocalRHS0;
                rGeom[1].GetValue(MAPPER_VECTOR_PROJECTION_RHS) += LocalRHS1;
                rGeom[2].GetValue(MAPPER_VECTOR_PROJECTION_RHS) += LocalRHS2;

                IV_iter++;
            }

            // Solve
            array_1d<double,3> dVal = ZeroVector(3);
            double dValNorm      = 0.0;
            double ValNorm       = 0.0;
            const unsigned int NodeNum = mrDestinationModelPart.NumberOfNodes();

            #pragma omp parallel for reduction(+ : dValNorm, ValNorm) private(dVal)
            for (int i=0; i<static_cast<int>(mrDestinationModelPart.NumberOfNodes()); ++i)
            {
                ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin() + i;
                const double & NodeArea = node_it->GetValue(NODAL_MAUX);
                noalias(dVal) = node_it->GetValue(MAPPER_VECTOR_PROJECTION_RHS)/NodeArea;
                array_1d<double, 3>& rDestinationValue = node_it->FastGetSolutionStepValue(rDestVar);
                rDestinationValue += sign * dVal;

                // Variables for convergence check
                for (unsigned int j = 0; j < 3; j++)
                {
                    dValNorm += dVal[j] * dVal[j];
                    ValNorm += rDestinationValue[j]*rDestinationValue[j];
                }
            }

            //std::cout << "VectorMap iteration: " << k+1 << std::endl;

            // Check Convergence
            const double RelativeError = (ValNorm > 10e-15) ? dValNorm / ValNorm : 0.0;
            if( (ValNorm/NodeNum < 0.00001 * TolIter * TolIter) || RelativeError < TolIter * TolIter)
            {
                // std::cout << "VectorMap converged in " << k + 1 << " iterations." << std::endl;
                break;
            }
            else if ( (k + 1) == MaxIter)
            {
                KRATOS_WARNING("AdvancedNMPointsMapper") << "VectorMap did not converge in " << k + 1 << " iterations." << std::endl;
            }
        }
        // End of Iteration
    }
} // End of Map (vector version)

/**
 * Test function, stores the distance between a Gauss Point and its projection
 */
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

/**
 * Auxiliar function to compute the nodal length/area of each node in both origin and destiny model parts.
 */
void AdvancedNMPointsMapper::ComputeNodalLengthArea()
{
    const unsigned int dimension = mrDestinationModelPart.ConditionsBegin()->GetGeometry().WorkingSpaceDimension();
    ModelPart::NodesContainerType& rOriginModelPartNodes = mrOriginModelPart.Nodes();
    ModelPart::NodesContainerType& rDestinationModelPartNodes = mrDestinationModelPart.Nodes();
    ModelPart::ConditionsContainerType& rOriginModelPartConditions = mrOriginModelPart.Conditions();
    ModelPart::ConditionsContainerType& rDestinationModelPartConditions = mrDestinationModelPart.Conditions();

    // NODAL_MAUX initialization
    // OriginModelPart
    #pragma omp parallel for
    for (int k = 0; k < static_cast<int>(rOriginModelPartNodes.size()); ++k)
    {
        ModelPart::NodesContainerType::iterator itNode = rOriginModelPartNodes.begin() + k;
        itNode->SetValue(NODAL_MAUX, 0.0);
    }

    // DestinationModelPart
    #pragma omp parallel for
    for (int k = 0; k < static_cast<int>(rDestinationModelPartNodes.size()); ++k)
    {
        ModelPart::NodesContainerType::iterator itNode = rDestinationModelPartNodes.begin() + k;
        itNode->SetValue(NODAL_MAUX, 0.0);
    }

    // NODAL_MAUX computation
    if (dimension == 2)
    {
        // OriginModelPart
        #pragma omp parallel for
        for (int k = 0; k < static_cast<int>(rOriginModelPartConditions.size()); ++k)
        {
            ModelPart::ConditionsContainerType::iterator itCond = rOriginModelPartConditions.begin() + k;
            GeometryType& rGeometry = itCond->GetGeometry();
            const double CondLength = rGeometry.Length();

            for (unsigned int i = 0; i < 2; i++)
            {
                double& nm = rGeometry[i].GetValue(NODAL_MAUX);
                #pragma omp atomic
                nm += 0.5 * CondLength;
            }
        }

        // DestinationModelPart
        #pragma omp parallel for
        for (int k = 0; k < static_cast<int>(rDestinationModelPartConditions.size()); ++k)
        {
            ModelPart::ConditionsContainerType::iterator itCond = rDestinationModelPartConditions.begin() + k;
            GeometryType& rGeometry = itCond->GetGeometry();
            const double CondLength = rGeometry.Length();

            for (unsigned int i = 0; i < 2; i++)
            {
                double& nm = rGeometry[i].GetValue(NODAL_MAUX);
                #pragma omp atomic
                nm += 0.5 * CondLength;
            }
        }
    }
    else
    {
        // OriginModelPart
        #pragma omp parallel for
        for (int k = 0; k < static_cast<int>(rOriginModelPartConditions.size()); ++k)
        {
            ModelPart::ConditionsContainerType::iterator itCond = rOriginModelPartConditions.begin() + k;
            GeometryType& rGeometry = itCond->GetGeometry();
            const double CondArea = rGeometry.Area();

            for (unsigned int i = 0; i < 3; i++)
            {
                double& nm = rGeometry[i].GetValue(NODAL_MAUX);
                #pragma omp atomic
                nm += (1.0/3.0) * CondArea;
            }
        }

        // DestinationModelPart
        #pragma omp parallel for
        for (int k = 0; k < static_cast<int>(rDestinationModelPartConditions.size()); ++k)
        {
            ModelPart::ConditionsContainerType::iterator itCond = rDestinationModelPartConditions.begin() + k;
            GeometryType& rGeometry = itCond->GetGeometry();
            const double CondArea = rGeometry.Area();

            for (unsigned int i = 0; i < 3; i++)
            {
                double& nm = rGeometry[i].GetValue(NODAL_MAUX);
                #pragma omp atomic
                nm += (1.0/3.0) * CondArea;
            }
        }
    }
}

} // Namespace Kratos.
