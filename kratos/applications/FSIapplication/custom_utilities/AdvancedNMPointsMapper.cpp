#include "AdvancedNMPointsMapper.h"

namespace Kratos
{
    // GaussPointItem Methods
    void GaussPointItem::Project(Condition::Pointer pOriginCond,
            array_1d<double,2>& Coords, double& Dist)
    {
        // xi,yi,zi are Nodal Coordinates, n is the destination condition's unit normal
        // and d is the distance along n from the point to its projection in the condition
        // | DestX-x3 |   | x1-x3 x2-x3 nx |   | Chi |
        // | DestY-y3 | = | y1-y3 y2-y3 ny | . | Eta |
        // | DestZ-z3 |   | z1-z3 z2-z3 nz |   |  d  |

        Matrix ChangeMatrix(3,3);
        Matrix InvChange(3,3);
        double det;

        array_1d<double,3> RHS, Res;

        RHS[0] = this->Coordinate(1) - pOriginCond->GetGeometry()[2].X();
        RHS[1] = this->Coordinate(2) - pOriginCond->GetGeometry()[2].Y();
        RHS[2] = this->Coordinate(3) - pOriginCond->GetGeometry()[2].Z();

        ChangeMatrix(0,0) = pOriginCond->GetGeometry()[0].X() - pOriginCond->GetGeometry()[2].X();
        ChangeMatrix(1,0) = pOriginCond->GetGeometry()[0].Y() - pOriginCond->GetGeometry()[2].Y();
        ChangeMatrix(2,0) = pOriginCond->GetGeometry()[0].Z() - pOriginCond->GetGeometry()[2].Z();

        ChangeMatrix(0,1) = pOriginCond->GetGeometry()[1].X() - pOriginCond->GetGeometry()[2].X();
        ChangeMatrix(1,1) = pOriginCond->GetGeometry()[1].Y() - pOriginCond->GetGeometry()[2].Y();
        ChangeMatrix(2,1) = pOriginCond->GetGeometry()[1].Z() - pOriginCond->GetGeometry()[2].Z();

        ChangeMatrix(0,2) = mNormal[0];
        ChangeMatrix(1,2) = mNormal[1];
        ChangeMatrix(2,2) = mNormal[2];

        MathUtils<double>::InvertMatrix3(ChangeMatrix,InvChange,det);
        noalias(Res) = prod(InvChange,RHS);

        Coords[0] = Res[0];
        Coords[1] = Res[1];
        // Keep distance positive, regardless of normal orientation
        Dist = (Res[2]<0)? -Res[2] : Res[2] ;

    }

    void GaussPointItem::GetProjectedValue(const Variable<double> & rOriginVar, double& Value)
    {
        if (mProjStatus == 1) // Get Interpolated value from origin condition
        {
            Value = ( (mpOriginCond.lock())->GetGeometry()[0].FastGetSolutionStepValue(rOriginVar)*mOriginCoords[0]
                    + (mpOriginCond.lock())->GetGeometry()[1].FastGetSolutionStepValue(rOriginVar)*mOriginCoords[1]
                    + (mpOriginCond.lock())->GetGeometry()[2].FastGetSolutionStepValue(rOriginVar)
                    *(1.0 - mOriginCoords[0] - mOriginCoords[1]) );
        } else if (mProjStatus == 2) // Get Value from origin node
        {
            Value = (mpOriginNode.lock())->FastGetSolutionStepValue(rOriginVar);
        } else // mProjStatus == 0: Return 0
        {
            Value = 0.0;
        }
    }

    void GaussPointItem::GetProjectedValue(const Variable<array_1d<double,3> >& rOriginVar, array_1d<double,3>& Value)
    {

        if (mProjStatus == 1) // Get Interpolated value from origin condition
        {
            for (int i=0; i<3; i++)
                Value[i] = ((mpOriginCond.lock())->GetGeometry()[0].FastGetSolutionStepValue(rOriginVar)[i]*mOriginCoords[0]
                        + (mpOriginCond.lock())->GetGeometry()[1].FastGetSolutionStepValue(rOriginVar)[i]*mOriginCoords[1]
                        + (mpOriginCond.lock())->GetGeometry()[2].FastGetSolutionStepValue(rOriginVar)[i]
                        *(1.0 - mOriginCoords[0] - mOriginCoords[1]) );
        } else if (mProjStatus == 2) // Get Value from origin node
        {
            Value = (mpOriginNode.lock())->FastGetSolutionStepValue(rOriginVar);
        } else // mProjStatus == 0: Return 0
        {
            Value[0] = Value[1] = Value[2] = 0.0;
        }
    }

    // Mapper Methods
    AdvancedNMPointsMapper::AdvancedNMPointsMapper(const ModelPart& rOriginModelPart,
            ModelPart& rDestinationModelPart):

            mrOriginModelPart(rOriginModelPart),
            mrDestinationModelPart(rDestinationModelPart),
            mBucketSize(4)
    {
        // The constructor defines the Gauss Points in the destination interface, stores them as
        // GaussPointItem instances and builds a std::vector of pointers to their position in memory
        // The Gauss Point Coordinates (Gix,Giy,Gyz) are obtained from the nodes' (A,B,C) coordinates
        // as follows:
        // | G1x G1y G1z |   | 0.6 0.2 0.2 |   | Ax Ay Az |
        // | G2x G2y G2z | = | 0.2 0.6 0.2 | . | Bx By Bz |
        // | G3x G3y G3z |   | 0.2 0.2 0.6 |   | Cx Cy Cz |

        array_1d<double,3> GPCoord; // Will store the coordinates of a condition's Gauss Points
        GPCoord[0] = GPCoord[1] = GPCoord[2] = 0.0;
        array_1d<double,3> Normal;
        double Area;

        MatrixVar Nodes, GaussPoints, GPPos;
        GPPos(0,0) = 0.6; GPPos(0,1) = 0.2; GPPos(0,2) = 0.2;
        GPPos(1,0) = 0.2; GPPos(1,1) = 0.6; GPPos(1,2) = 0.2;
        GPPos(2,0) = 0.2; GPPos(2,1) = 0.2; GPPos(2,2) = 0.6;

        for (
                ModelPart::ConditionsContainerType::iterator cond_it = rDestinationModelPart.ConditionsBegin();
                cond_it != rDestinationModelPart.ConditionsEnd();
                cond_it++)
        {
            CalcNormalAndArea(*cond_it.base(),Normal,Area);

            for(unsigned int i=0; i<3; i++)
                for(unsigned int j=0; j<3; j++)
                {
                    Nodes(i,j) = cond_it->GetGeometry()[i].Coordinate(j+1);
                }

            noalias(GaussPoints) = prod(GPPos,Nodes);

            for(unsigned int k=0; k<3; k++)
            {
                for(unsigned int l=0; l<3; l++)
                {
                    GPCoord[l] = GaussPoints(k,l);
                }
                GaussPointItem::Pointer pGP = GaussPointItem::Pointer(new GaussPointItem(GPCoord,Area,Normal) );
                mGaussPointList.push_back( pGP );
            }
        }
    }

    void AdvancedNMPointsMapper::FindNeighbours(double SearchRadiusFactor)
    {
        // Initialize some values
        const unsigned int MaxResults = 5000; // Maximum number of points to find in a single tree search
        GaussPointVector Results(MaxResults);
        std::vector<double> ResultDistances(MaxResults);
        array_1d<double,3> ZeroVect;
        ZeroVect[0] = ZeroVect[1] = ZeroVect[2] = 0.0;

        // Create a tree
        // It will use a copy of mGaussPoinList (a std::vector which contains pointers)
        // Copying the list is required because the tree will reorder it for efficiency
        GaussPointVector pGaussPoints = mGaussPointList;
        tree Tree_conds(pGaussPoints.begin(),pGaussPoints.end(),mBucketSize);

        double Radius, SearchRadius;
        double MaxRadius = 0.0;

        for(
                ModelPart::ConditionsContainerType::const_iterator cond_it = mrOriginModelPart.ConditionsBegin();
                cond_it != mrOriginModelPart.ConditionsEnd();
                cond_it++)
        {
            Point<3> Center;
            TriangleCenterAndRadius(*cond_it.base(),Center,Radius);

            // Create a fake GaussPoint with the same coordinates
            // (to ensure that the object used in the kd-tree search is of the same type as the tree contents)
            GaussPointItem CenterGP(Center.Coordinates(),0,ZeroVect);

            if (Radius > MaxRadius) MaxRadius = Radius;

            // Count Gauss Points within SearchRadius
            SearchRadius = SearchRadiusFactor*Radius;
            int Found = Tree_conds.SearchInRadius(CenterGP,SearchRadius,Results.begin(),ResultDistances.begin(),MaxResults);

            for(int i=0; i<Found;i++) SetProjectionToCond(*Results[i], *cond_it.base());
        }

        // Try to find reasonable origin values in nearby nodes for Gauss Points that couldn't be projected to a condition
        GaussPointVector ProjectionlessGP;

        for (
                GaussPointIterator gauss_it = mGaussPointList.begin();
                gauss_it != mGaussPointList.end();
                gauss_it++)
        {
            int Status = 0;
            (*gauss_it)->GetProjStatus(Status);

            if( Status != 1) ProjectionlessGP.push_back( *gauss_it );
        }

        // Try to use a nearby node as reference for points that couldn't be projected to a condition
        if ( ProjectionlessGP.size() != 0)
        {
            std::cout << "AdvancedNMPointsMapper: " << ProjectionlessGP.size()
                    << " Gauss points could not be projected to a condition.\n"
                    << "    Attempring to get a reasonable value for them from a node..." << std::endl;

            GaussPointVector NodeResults(MaxResults);
            std::vector<double> NodeResultDistances(MaxResults);

            tree kdtree_nodes(ProjectionlessGP.begin(),ProjectionlessGP.end(),mBucketSize);
            SearchRadius = SearchRadiusFactor*MaxRadius;

            for (
                    ModelPart::NodesContainerType::const_iterator node_it = mrOriginModelPart.NodesBegin();
                    node_it != mrOriginModelPart.NodesEnd();
                    node_it++)
            {
                // Our tree uses GaussPointItem objects as input to sort
                GaussPointItem NodePos(node_it->Coordinates(),0,ZeroVect);
                int Found = kdtree_nodes.SearchInRadius(NodePos,SearchRadius,NodeResults.begin(),NodeResultDistances.begin(),MaxResults);

                for (int i=0; i<Found;i++) SetProjectionToNode(*NodeResults[i],*node_it.base(),NodeResultDistances[i]);
            }
            //DistanceCheck(); // Test function

            // Count how many points without a projection remain
            unsigned int counter = 0;
            for (
                    GaussPointIterator gauss_it = mGaussPointList.begin();
                    gauss_it != mGaussPointList.end();
                    gauss_it++)
            {
                int Status = 0;
                (*gauss_it)->GetProjStatus(Status);

                if( Status == 0) counter++;
            }
            std::cout << "   ... " << counter << " Gauss Points without a reference value remain." << std::endl;
        }
    }

    void AdvancedNMPointsMapper::CalcNormalAndArea(Condition::Pointer Cond,
            array_1d<double,3> & Normal, double & Area)
    {
        array_1d<double,3> v1,v2;

        v1[0] = Cond->GetGeometry()[1].X() - Cond->GetGeometry()[0].X();
        v1[1] = Cond->GetGeometry()[1].Y() - Cond->GetGeometry()[0].Y();
        v1[2] = Cond->GetGeometry()[1].Z() - Cond->GetGeometry()[0].Z();

        v2[0] = Cond->GetGeometry()[2].X() - Cond->GetGeometry()[0].X();
        v2[1] = Cond->GetGeometry()[2].Y() - Cond->GetGeometry()[0].Y();
        v2[2] = Cond->GetGeometry()[2].Z() - Cond->GetGeometry()[0].Z();

        MathUtils<double>::CrossProduct(Normal,v1,v2);

        double NNorm = sqrt(Normal[0]*Normal[0] + Normal[1]*Normal[1]
            + Normal[2]*Normal[2]);

        for(int i=0; i<3; i++) Normal[i] = Normal[i]/NNorm;
        Area = 0.5*NNorm;
        //std::cout << Normal[0] << " " << Normal[1] << " " << Normal[2] << std::endl;
    }

    void AdvancedNMPointsMapper::TriangleCenterAndRadius(const Condition::Pointer Cond,
                Point<3>& Center, double& Radius)
    {
        Radius = 0;
        Center.Coordinate(1) = 0.33333333333333333333*(Cond->GetGeometry()[0].X() + Cond->GetGeometry()[1].X() + Cond->GetGeometry()[2].X());
        Center.Coordinate(2) = 0.33333333333333333333*(Cond->GetGeometry()[0].Y() + Cond->GetGeometry()[1].Y() + Cond->GetGeometry()[2].Y());
        Center.Coordinate(3) = 0.33333333333333333333*(Cond->GetGeometry()[0].Z() + Cond->GetGeometry()[1].Z() + Cond->GetGeometry()[2].Z());

        for(int i=0;i<3;i++)
        {
            double dx = Center.Coordinate(1) - Cond->GetGeometry()[i].X();
            double dy = Center.Coordinate(2) - Cond->GetGeometry()[i].Y();
            double dz = Center.Coordinate(3) - Cond->GetGeometry()[i].Z();

            double tmp = dx*dx + dy*dy + dz*dz;

            if(tmp > Radius) Radius = tmp;
        }
        Radius = sqrt(Radius);
    }

    void AdvancedNMPointsMapper::SetProjectionToCond(GaussPointItem& GaussPoint,
            Condition::Pointer pCandidateCond)
    {
        array_1d<double,2> LocalCoords;
        double Dist;
        GaussPoint.Project(pCandidateCond,LocalCoords,Dist);

        if (// Point belongs to Condition
                LocalCoords[0]>=0 &&
                LocalCoords[1]>=0 &&
                (1.0 - LocalCoords[0] - LocalCoords[1])>=0)
        {
            int Status = 0;
            GaussPoint.GetProjStatus(Status);

            if ( Status != 1) // No good projection found (yet)
                GaussPoint.SetProjection(pCandidateCond,LocalCoords,Dist);
            else // A good projection was found in a previous iteration
            {
                double CurrentDist = 0.0;
                GaussPoint.GetDist(CurrentDist);

                if (Dist < CurrentDist)
                    GaussPoint.SetProjection(pCandidateCond,LocalCoords,Dist);
            }
        }
    }

    void AdvancedNMPointsMapper::SetProjectionToNode(GaussPointItem& GaussPoint,
            Node<3>::Pointer pCandidateNode, double& Dist)
    {
        int Status = 0;
        GaussPoint.GetProjStatus(Status);

        if ( Status == 0)
            GaussPoint.SetProjection(pCandidateNode,Dist);
        else
        {
            double CurrentDist = 0.0;
            GaussPoint.GetDist(CurrentDist);

            if (Dist < CurrentDist)
            {
//                Node<3>::WeakPointer pNewNode(*pCandidateNode);
//                pNewNode.weak_ptr(pCandidateNode,)
                GaussPoint.SetProjection(pCandidateNode,Dist);
            }
        }
    }

    void AdvancedNMPointsMapper::ScalarMap(const Variable<double> & rOriginVar, Variable<double> & rDestVar,
            const int MaxIter, const double TolIter)
    {
        // Define some variables that will be used in the iteration
        MatrixVar MCons; // Elemental Consistent Mass Matrix = Aelem/12 * MCons
        MCons(0,0) = 2.0; MCons(0,1) = 1.0; MCons(0,2) = 1.0;
        MCons(1,0) = 1.0; MCons(1,1) = 2.0; MCons(1,2) = 1.0;
        MCons(2,0) = 1.0; MCons(2,1) = 1.0; MCons(2,2) = 2.0;

        MatrixVar MInterp; // Interpolation Matrix (NodalValues = (A/24)*MInterp*GaussValues)
        MInterp(0,0) = 6.0; MInterp(0,1) = 1.0; MInterp(0,2) = 1.0;
        MInterp(1,0) = 1.0; MInterp(1,1) = 6.0; MInterp(1,2) = 1.0;
        MInterp(2,0) = 1.0; MInterp(2,1) = 1.0; MInterp(2,2) = 6.0;

        std::vector< boost::shared_ptr< array_1d<double,3> > > pInterpValues;

        // Build (Diagonal) System Matrix and initialize results
        for ( ModelPart::NodesContainerType::iterator node_it = mrDestinationModelPart.NodesBegin();
                node_it != mrDestinationModelPart.NodesEnd();
                node_it++)
        {
            node_it->GetValue(NODAL_MAUX) = 0.0;
            node_it->FastGetSolutionStepValue(rDestVar) = 0.0;
        }

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

            for (int i=0;i<3;i++)
            {
                cond_it->GetGeometry()[i].GetValue(NODAL_MAUX) += 0.333333333333333333*CondArea;
                mGaussPointList[GPi+i]->GetProjectedValue(rOriginVar,GPValues[i]);
            }

            noalias(NodalValues) = (CondArea/24.0)*prod(MInterp,GPValues);

            boost::shared_ptr< array_1d<double,3> > pNodalValues(new array_1d<double,3>(NodalValues) );
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
                node_it->GetValue(AUX)=0.0;

            array_1d<double,3> LocalRHS;
            array_1d<double,3> LastSolution;
            int IV_iter = 0;

            for (ModelPart::ConditionsContainerType::iterator cond_it = mrDestinationModelPart.ConditionsBegin();
                    cond_it != mrDestinationModelPart.ConditionsEnd();
                    cond_it++)
            {
                for(int j=0; j<3; j++)
                {
                    LocalRHS[j] = 0.0;
                    LastSolution[j] = cond_it->GetGeometry()[j].FastGetSolutionStepValue(rDestVar);
                }

                double CondArea = 0.0;
                mGaussPointList[3*IV_iter]->GetArea(CondArea);

                noalias(LocalRHS) = *pInterpValues[IV_iter] - (CondArea/12.0)*prod(MCons,LastSolution);
                // We are taking advantage of 1 Condition = 3 Gauss Points to iterate the interpolation results and the Gauss Point Vector Again

                for (int j=0; j<3 ; j++) cond_it->GetGeometry()[j].GetValue(AUX)+=LocalRHS[j];

                IV_iter++;
            }

            // Solve
            double dValNorm = 0.0;
            double ValNorm = 0.0;
            double RelativeError = 0.0;
            unsigned int NodeNum = 0;

            for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
                    node_it != mrDestinationModelPart.NodesEnd();
                    node_it++)
            {
                const double & NodeArea = node_it->GetValue(NODAL_MAUX);
                //node_it->FastGetSolutionStepValue(NODAL_MAUX)=NodeArea; // TEST: store nodal area so GiD can paint it later
                double dVal = node_it->GetValue(AUX)/NodeArea;
                node_it->FastGetSolutionStepValue(rDestVar) += dVal;

                // Variables for convergence check
                dValNorm += dVal*dVal;
                ValNorm += pow(node_it->FastGetSolutionStepValue(rDestVar),2);
                NodeNum++;
            }
            //std::cout << "ScalarMap iteration: " << k+1 << std::endl;

            // Check Convergence
            if(ValNorm > 10e-15) RelativeError = dValNorm / ValNorm;
            if( (ValNorm/NodeNum < 0.00001*TolIter*TolIter) ||
                    (RelativeError < TolIter*TolIter) )
            {
                std::cout << "ScalarMap converged in " << k+1 << " iterations." << std::endl;
                break;
            } else if ((k+1) == MaxIter) std::cout << "WARNING: ScalarMap did not converge in " << k+1 << " iterations." << std::endl;
        } // End of Iteration
    } // End of Map (scalar version)

    void AdvancedNMPointsMapper::VectorMap(const Variable<array_1d<double,3> >& rOriginVar, Variable<array_1d<double,3> >& rDestVar,
            const int MaxIter, const double TolIter)
    {
        array_1d<double,3> ZeroVect;
        ZeroVect[0] = ZeroVect[1] = ZeroVect[2] = 0.0;

        std::vector< boost::shared_ptr<array_1d<double,9> > > pInterpValues;

        // Build (Diagonal) System Matrix and initialize results
        for ( ModelPart::NodesContainerType::iterator node_it = mrDestinationModelPart.NodesBegin();
                node_it != mrDestinationModelPart.NodesEnd();
                node_it++)
        {
            node_it->GetValue(NODAL_MAUX) = 0.0;
            node_it->FastGetSolutionStepValue(rDestVar) = ZeroVect;
        }

        // Here we loop both the Destination Model Part and the Gauss Point List, using the
        // fact that the GP list was created by a loop over the Dest. Model Part's conditions
        int GPi = 0;
        for (ModelPart::ConditionIterator cond_it = mrDestinationModelPart.ConditionsBegin();
                cond_it != mrDestinationModelPart.ConditionsEnd();
                cond_it++)
        {
            array_1d<double,9> GPValues;
            array_1d<double,9> NodalValues;

            double CondArea = 0.0;
            mGaussPointList[GPi]->GetArea(CondArea);
            const double K = CondArea/24.0;

            for (int i=0;i<3;i++) cond_it->GetGeometry()[i].GetValue(NODAL_MAUX) += 0.333333333333333333*CondArea;

            array_1d<double,3> TempValues = ZeroVect;
            mGaussPointList[GPi]->GetProjectedValue(rOriginVar,TempValues);
            for (int i=0;i<3;i++) GPValues[i] = TempValues[i];

            mGaussPointList[GPi+1]->GetProjectedValue(rOriginVar,TempValues);
            for (int i=0;i<3;i++) GPValues[3+i] = TempValues[i];

            mGaussPointList[GPi+2]->GetProjectedValue(rOriginVar,TempValues);
            for (int i=0;i<3;i++) GPValues[6+i] = TempValues[i];

            for (int i=0;i<3;i++)
            {
                NodalValues[i] = K*(6.0*GPValues[i] + 1.0*GPValues[3+i] + 1.0*GPValues[6+i]);
                NodalValues[3+i] = K*(1.0*GPValues[i] + 6.0*GPValues[3+i] + 1.0*GPValues[6+i]);
                NodalValues[6+i] = K*(1.0*GPValues[i] + 1.0*GPValues[3+i] + 6.0*GPValues[6+i]);
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
                node_it->GetValue(VAUX) = ZeroVect;

            array_1d<double,3> LocalRHS0, LocalRHS1, LocalRHS2; // Local RHS for each node
            array_1d<double,9> LastSolution;
            int IV_iter = 0;

            for (ModelPart::ConditionIterator cond_it = mrDestinationModelPart.ConditionsBegin();
                    cond_it != mrDestinationModelPart.ConditionsEnd();
                    cond_it++)
            {
                double CondArea = 0.0;
                mGaussPointList[3*IV_iter]->GetArea(CondArea);
                const double K = CondArea/12.0;

                for(int j=0; j<3; j++)
                {
                    LocalRHS0[j] = LocalRHS1[j] = LocalRHS2[j] = 0.0;
                    LastSolution[j] = cond_it->GetGeometry()[0].FastGetSolutionStepValue(rDestVar)[j];
                    LastSolution[3+j] = cond_it->GetGeometry()[1].FastGetSolutionStepValue(rDestVar)[j];
                    LastSolution[6+j] = cond_it->GetGeometry()[2].FastGetSolutionStepValue(rDestVar)[j];
                }

                array_1d<double,9> CondValues = *pInterpValues[IV_iter];

                for(int j=0;j<3;j++)
                {
                    LocalRHS0[j] = CondValues[j] - K*(2.0*LastSolution[j] + 1.0*LastSolution[j+3] + 1.0*LastSolution[j+6]);
                    LocalRHS1[j] = CondValues[j+3] - K*(1.0*LastSolution[j] + 2.0*LastSolution[j+3] + 1.0*LastSolution[j+6]);
                    LocalRHS2[j] = CondValues[j+6] - K*(1.0*LastSolution[j] + 1.0*LastSolution[j+3] + 2.0*LastSolution[j+6]);
                }
                // We are taking advantage of 1 Condition = 3 Gauss Points to iterate the interpolation results and the Gauss Point Vector Again

                cond_it->GetGeometry()[0].GetValue(VAUX)+=LocalRHS0;
                cond_it->GetGeometry()[1].GetValue(VAUX)+=LocalRHS1;
                cond_it->GetGeometry()[2].GetValue(VAUX)+=LocalRHS2;

                IV_iter++;
            }
            // Solve
            array_1d<double,3> dVal = ZeroVect;
            double dValNorm = 0.0;
            double ValNorm = 0.0;
            double RelativeError = 0.0;
            unsigned int NodeNum = 0;

            for (ModelPart::NodeIterator node_it = mrDestinationModelPart.NodesBegin();
                    node_it != mrDestinationModelPart.NodesEnd();
                    node_it++)
            {
                const double & NodeArea = node_it->GetValue(NODAL_MAUX);
                noalias(dVal) = node_it->GetValue(VAUX)/NodeArea;
                node_it->FastGetSolutionStepValue(rDestVar) += dVal;

                // Variables for convergence check
                for (int j=0; j<3; j++)
                {
                    dValNorm += dVal[j]*dVal[j];
                    ValNorm += pow(node_it->FastGetSolutionStepValue(rDestVar)[j],2);
                }
                NodeNum++;
            }
            //std::cout << "VectorMap iteration: " << k+1 << std::endl;

            // Check Convergence
            if(ValNorm > 10e-15) RelativeError = dValNorm / ValNorm;
            if( (ValNorm/NodeNum < 0.00001*TolIter*TolIter) ||
                    RelativeError < TolIter*TolIter)
            {
                std::cout << "VectorMap converged in " << k+1 << " iterations." << std::endl;
                break;
            } else if ( (k+1) == MaxIter) std::cout << "WARNING: VectorMap did not converge in " << k+1 << " iterations." << std::endl;
        } // End of Iteration
    } // End of Map (vector version)

    void AdvancedNMPointsMapper::DistanceCheck() // Test function, stores the distance between a Gauss Point and its projection
    {
        unsigned int GPiter = 0;
        for (ModelPart::ConditionsContainerType::iterator cond_it = mrDestinationModelPart.ConditionsBegin();
                cond_it != mrDestinationModelPart.ConditionsEnd();
                cond_it++)
        {
            for (int i = 0; i<3; i++)
            {
                double dist = 0.0;
                int Status = 0;
                mGaussPointList[GPiter+i]->GetProjStatus(Status);
                std::cout << GPiter+i << " " << Status << std::endl;

                if (Status != 0)
                {
                    mGaussPointList[GPiter+i]->GetDist(dist);
                    if (Status == 2) dist = -sqrt(dist);
                } else {
                    dist = -10000;
                }

                cond_it->GetGeometry()[i].FastGetSolutionStepValue(FICTITIOUS_FLUID_DENSITY) = dist;
            }
            GPiter += 3;
        }
    }
}
