//
//   Project Name:        Kratos Wind Turbine Application Utility
//   Last Modified by:    $Author: efusto $
//   Date:                $Date: 2012-05-16 11:45:23 $
//   Revision:            $Revision: 0.1 $
//
//

// System includes
#include <string>
#include <iostream>
#include <algorithm>
#include <vector>
#include <math.h>

// file printout for debug
#include <iostream>
#include <fstream>
#include <sstream>

//////////////////////////////////////////////////////////////
// Activate this #define to compile with MPI parallel extension
//
//#define WIND_TURBINE_USE_PARALLEL_EXTENSION
//
//////////////////////////////////////////////////////////////

#if defined( WIND_TURBINE_USE_PARALLEL_EXTENSION )
    #undef REAL
    #include <mpi.h>
#endif

#include "wind_turbine_rotation_utilities.h"

namespace Kratos
{

WindTurbineRegionMultiplicity::WindTurbineRegionMultiplicity()
	:mRegion(0)
	,mMultiplicity(0)
{}

WindTurbineRegionMultiplicity::WindTurbineRegionMultiplicity(int region)
{
	mRegion = region;
	mMultiplicity = 0;
}

WindTurbineRegionMultiplicity::WindTurbineRegionMultiplicity(int region, int multiplicity)
{
        mRegion = region;
        mMultiplicity = multiplicity;
}

bool WindTurbineRegionMultiplicity::operator<(WindTurbineRegionMultiplicity otherRegionMult) const
{
	return mMultiplicity < otherRegionMult.Multiplicity();
}


WindTurbineRegionMultiplicity WindTurbineRegionMultiplicity::operator++(int)
{
	WindTurbineRegionMultiplicity copy = *this;
	mMultiplicity++;
	return copy;
}


WindTurbineRotationUtilities::WindTurbineRotationUtilities(ModelPart& rAllTheModelPart)
	:mrGlobalModelPart(rAllTheModelPart)
{
	mFirstOuterInterfaceNodeOffset = 0;
        mNumberOfBoundaryFaces = 0;
        mAngleHistory.reserve(2); // Old rotation angles
        mEchoLevel = WIND_TURBINE_ECHOLEVEL_DEEPINFO;   //_ALL;

        FillNodeRegions();
        FillElementRegions();

#if defined( WIND_TURBINE_USE_PARALLEL_EXTENSION )
        Parallel_DecideRemeshingProcessor();
        Parallel_MigrateEntities(mInterfaceNodes);
        // the processors different from the remesher don't need the crown elements anymore
        if (mThisRank != mRemeshingRank)
        {
            DestroyCrownElements();
            DestroyCrownNodes();
        }

        MPI_Barrier(MPI_COMM_WORLD);
        RemoveLocalNodesWithNoElements();
        MPI_Barrier(MPI_COMM_WORLD);
#else
        // these are unused though, in this case
        mThisRank = 0;
        mRemeshingRank = 0;
        mNumberOfRanks = 1;
#endif

}

/// Percolate the whole ModelPart and put the elements in the proper region lists
void WindTurbineRotationUtilities::FillElementRegions()
{
	ModelPart::ElementsContainerType::Pointer pWholeElements = mrGlobalModelPart.pElements();

	// getting Elems geometry to initialize the FLAG_VARIABLE marker object (we must track the regions)
	ModelPart::ElementsContainerType::iterator itr = pWholeElements->begin();
        mNumberOfNodesPerElement = itr->GetGeometry().PointsNumber();
	std::vector<WindTurbineRegionMultiplicity> elemBelongings(WIND_TURBINE_REGION_NUMBER);

        if (mEchoLevel > WIND_TURBINE_ECHOLEVEL_NONE)
        {
            std::cout << "Number of Nodes per Element is " << mNumberOfNodesPerElement << std::endl;
        }

	while ( itr != pWholeElements->end() )
	{
		// initialize at each cycle...
		for (unsigned int i=0; i < WIND_TURBINE_REGION_NUMBER; i++)
                    elemBelongings[i] = WindTurbineRegionMultiplicity(i);

                if (mEchoLevel > WIND_TURBINE_ECHOLEVEL_DEEPINFO)
                {
                    std::cout << "Processing Element " << itr->Id() << std::endl;
                }
                // set region multiplicity of each region
                for( unsigned int n = 0; n < mNumberOfNodesPerElement; n++ )
		{
                    int region = (int)itr->GetGeometry()[n].FastGetSolutionStepValue(FLAG_VARIABLE);
                    elemBelongings[region]++;

                    if (mEchoLevel > WIND_TURBINE_ECHOLEVEL_DEEPINFO)
                    {
                        std::cout << "      Region is " << region << std::endl;
                        std::cout << "      Node " << itr->GetGeometry()[n].Id() << ": FLAG_VARIABLE " << elemBelongings[region].Region() << std::endl;
                    }

		}

		// sorting the vector for following usefulness
		std::sort(elemBelongings.begin(), elemBelongings.end());
		std::reverse(elemBelongings.begin(), elemBelongings.end());

		unsigned int destRegion = WIND_TURBINE_UNKNOWN_REGION;
		unsigned int oppositeVertex = WIND_TURBINE_UNKNOWN_REGION;
                bool warn = false;
                destRegion = DecideElementRegion(mNumberOfNodesPerElement, elemBelongings, oppositeVertex, warn);

                // pushing the element in the proper container member
                if (warn)
                {
                    std::cout << "Undefined belongings (" << destRegion << ")for element " << itr->Id() << std::endl << ", with regions in the multiplicity vector: " << std::endl;
                    for (unsigned int i = 0; i < mNumberOfNodesPerElement; i++)
                        std::cout << elemBelongings[i].Region() << ", " << std::endl;
                    exit(-1);
                }
                else
                {
                    if (mEchoLevel > WIND_TURBINE_ECHOLEVEL_INFO)
                    {
                        std::cout << "Assigning Element " << itr->Id() << " to region " << destRegion;
                        if (oppositeVertex)
                            std::cout << " with vertex in " << oppositeVertex << std::endl;
                    }
                }

                switch (destRegion)
		{
                    case WIND_TURBINE_INNER_REGION:
                        mInnerElems.push_back(*(itr.base()));
                        break;
                    case WIND_TURBINE_INNER_INTERF_REGION:
                        mInnerInterfElems.push_back(*(itr.base()));
                        break;
                    case WIND_TURBINE_OUTER_INTERF_REGION:
                        mOuterInterfElems.push_back(*(itr.base()));
                        break;
                    case WIND_TURBINE_OUTER_REGION:
                        mOuterElems.push_back(*(itr.base()));
                        break;
                    case WIND_TURBINE_CROWN_REGION:
                        mCrownElems.push_back(*(itr.base()));
                        break;
                    default:
                        if (mEchoLevel > WIND_TURBINE_ECHOLEVEL_WARNING)
                        {
                            std::cout << "Unidentified region!!" << std::endl;
                        }
                        break;
		}

                // maybe there are also useful edges to be retrieved
                if (oppositeVertex != WIND_TURBINE_UNKNOWN_REGION)
                {
                    Geometry< Node<3> > elemGeom = itr->GetGeometry();
                    mConstrainedBoundaryNodeAuxIndices.reserve(mNumberOfNodesPerElement - 1);

                    if (mNumberOfNodesPerElement == 3)   //Triangle
                    {
                        // this means the opposite vertex lays in WIND_TURBINE_{INNER/OUTER}_INTERF_REGION,
                        // finding it I've also found the counterclockwise oriented basis

                        //loop on all the nodes in the geometry and get the neighbour elements of each node
                        for (unsigned int n=0; n < mNumberOfNodesPerElement; n++)
                        {
                            if ((unsigned int)(elemGeom[n].GetSolutionStepValue(FLAG_VARIABLE)) == oppositeVertex)
                            {
                                switch (n)
                                {
                                    case 0:
                                        mConstrainedBoundaryNodeAuxIndices.push_back( elemGeom[1].GetValue(AUX_ID) );
                                        mConstrainedBoundaryNodeAuxIndices.push_back( elemGeom[2].GetValue(AUX_ID) );
                                        break;
                                    case 1:
                                        mConstrainedBoundaryNodeAuxIndices.push_back( elemGeom[2].GetValue(AUX_ID) );
                                        mConstrainedBoundaryNodeAuxIndices.push_back( elemGeom[0].GetValue(AUX_ID) );
                                        break;
                                    case 2:
                                        mConstrainedBoundaryNodeAuxIndices.push_back( elemGeom[0].GetValue(AUX_ID) );
                                        mConstrainedBoundaryNodeAuxIndices.push_back( elemGeom[1].GetValue(AUX_ID) );
                                        break;
                                }
                            }
                        }
                    }
                    else if (mNumberOfNodesPerElement == 4)  //Thetraedra
                    {
                        //loop on all the nodes in the geometry and get the neighbour elements of each
                        for (unsigned int n=0; n < mNumberOfNodesPerElement; n++)
                        {
                            if ((unsigned int)(elemGeom[n].GetSolutionStepValue(FLAG_VARIABLE)) == oppositeVertex)
                            {
                                switch (n)
                                {
                                    // vertex : base
                                    // 0 ------ 1 2 3
                                    // 1 ------ 0 3 2
                                    // 2 ------ 0 1 3
                                    // 3 ------ 0 2 1
                                    case 0:
                                        mConstrainedBoundaryNodeAuxIndices.push_back( elemGeom[1].GetValue(AUX_ID) );
                                        mConstrainedBoundaryNodeAuxIndices.push_back( elemGeom[2].GetValue(AUX_ID) );
                                        mConstrainedBoundaryNodeAuxIndices.push_back( elemGeom[3].GetValue(AUX_ID) );
                                        break;
                                    case 1:
                                        mConstrainedBoundaryNodeAuxIndices.push_back( elemGeom[0].GetValue(AUX_ID) );
                                        mConstrainedBoundaryNodeAuxIndices.push_back( elemGeom[3].GetValue(AUX_ID) );
                                        mConstrainedBoundaryNodeAuxIndices.push_back( elemGeom[2].GetValue(AUX_ID) );
                                        break;
                                    case 2:
                                        mConstrainedBoundaryNodeAuxIndices.push_back( elemGeom[0].GetValue(AUX_ID) );
                                        mConstrainedBoundaryNodeAuxIndices.push_back( elemGeom[1].GetValue(AUX_ID) );
                                        mConstrainedBoundaryNodeAuxIndices.push_back( elemGeom[3].GetValue(AUX_ID) );
                                        break;
                                    case 3:
                                        mConstrainedBoundaryNodeAuxIndices.push_back( elemGeom[0].GetValue(AUX_ID) );
                                        mConstrainedBoundaryNodeAuxIndices.push_back( elemGeom[2].GetValue(AUX_ID) );
                                        mConstrainedBoundaryNodeAuxIndices.push_back( elemGeom[1].GetValue(AUX_ID) );
                                        break;
                                }
                                mNumberOfBoundaryFaces++;
                            }
                        }
                    }
                }
		// move on the iterator
		itr++;
	}


        if (mEchoLevel > WIND_TURBINE_ECHOLEVEL_DEEPINFO)
        {
            std::cout << std::endl << "--- Assignments----------------------" << std::endl;
            std::cout << "Element regions filled." << std::endl;
            std::cout << "Inner Interface: " << mInnerInterfElems.size() << " elems" << std::endl;
            std::cout << "Outer Interface: " << mOuterInterfElems.size() << " elems" << std::endl << std::endl;

            std::cout << "Inner Region: " << mInnerElems.size() << " elems" << std::endl;
            std::cout << "Outer Region: " << mOuterElems.size() << " elems" << std::endl;
            std::cout << "Crown Region: " << mCrownElems.size()  << " elems" << std::endl << std::endl;

            std::cout << "Constrained Components: " << mNumberOfBoundaryFaces << ". Constrained nodes :" << mConstrainedBoundaryNodeAuxIndices.size() << std::endl;
        }
}

void WindTurbineRotationUtilities::FillNodeRegions()
{
        int interfNode_AUX_ID = 1;
	ModelPart::NodesContainerType outerInterfaceNodes;

	for ( ModelPart::NodesContainerType::iterator itr = mrGlobalModelPart.NodesBegin();
	itr != mrGlobalModelPart.NodesEnd();
	itr++)
	{
            switch((int)(itr->GetSolutionStepValue(FLAG_VARIABLE)))
            {
                case WIND_TURBINE_INNER_REGION:
                    mInnerNodes.push_back(*(itr.base()));
                    break;

                case WIND_TURBINE_INTERNAL_CYL_BASE:
                case WIND_TURBINE_INNER_INTERF_REGION:
                    itr->GetValue(AUX_ID) = interfNode_AUX_ID;
                    mInterfaceNodes.push_back(*(itr.base()));
                    interfNode_AUX_ID++;
                    break;

                case WIND_TURBINE_EXTERNAL_CYL_BASE:
                case WIND_TURBINE_OUTER_INTERF_REGION:
                    itr->GetValue(AUX_ID) = interfNode_AUX_ID;
                    outerInterfaceNodes.push_back(*(itr.base()));
                    interfNode_AUX_ID++;
                    break;

                default:
                    break;
            }
	}

	mFirstOuterInterfaceNodeOffset = mInterfaceNodes.size();
        mInterfaceNodes.reserve( outerInterfaceNodes.size() );
	for ( ModelPart::NodesContainerType::iterator itr = outerInterfaceNodes.begin();
        itr != outerInterfaceNodes.end();
        itr++)
        {
		mInterfaceNodes.push_back(*(itr.base()));
	}
}

// rotation is made around the origin O(0,0)
void WindTurbineRotationUtilities::DoRotationAndRemesh(const int& dimensions, const double& rotAngle, const double& timeStep)
{
#if defined ( WIND_TURBINE_USE_PARALLEL_EXTENSION )
    // in this case, only the remesher's subdomain can contain crown elements
    if ( mThisRank == mRemeshingRank )
    {
#endif
	// erasing elements in the crown region
	DestroyCrownElements();
#if defined ( WIND_TURBINE_USE_PARALLEL_EXTENSION )
    }
#endif
        // erasing also the nodes in the crown region (this just for 3D, but could be also for 2D)
        if (dimensions == 3)
            DestroyCrownNodes();

#if defined ( WIND_TURBINE_USE_PARALLEL_EXTENSION )
    MPI_Barrier(MPI_COMM_WORLD);
#endif
	// rotating inner nodes and consequently the elements identified by them
        RotateEntities((double)WIND_TURBINE_INNER_REGION, rotAngle, timeStep);
        RotateEntities((double)WIND_TURBINE_INNER_INTERF_REGION, rotAngle, timeStep);

#if defined ( WIND_TURBINE_USE_PARALLEL_EXTENSION )
    // in this case, only the remeshing processor can perform the regeneration
    if ( mThisRank == mRemeshingRank )
    {
#endif
        // remeshing elements in the crown region
        if (dimensions == 3)
            RegenerateCrownElements3D();
        else if (dimensions == 2)
            RegenerateCrownElements2D();
        else
            std::cout << "Geometry not supported." << std::endl;
#if defined ( WIND_TURBINE_USE_PARALLEL_EXTENSION )
    }

    MPI_Barrier(MPI_COMM_WORLD);
#endif

}

void WindTurbineRotationUtilities::RotateEntities(const double& regionFlag, const double& rotAngle, const double& timeStep)
{
	ModelPart::NodesContainerType* pNodes = &mInnerNodes;
	int itrOffset = 0;
	switch ((int)regionFlag)
	{
		case WIND_TURBINE_INNER_REGION:
			pNodes = &mInnerNodes;
			break;
		case WIND_TURBINE_INNER_INTERF_REGION:
			pNodes = &mInterfaceNodes;
                        itrOffset = mInterfaceNodes.size() - mFirstOuterInterfaceNodeOffset;  // mind that in the 3D case, this offset includes also the two cylindrical bases
			break;
		default:
			std::cout << "WARNING: unhandled region " << regionFlag << "... won't rotate anything." << std::endl;
                        return;
			break;
	}

        double cosRotAngle = cos(rotAngle);
        double sinRotAngle = sin(rotAngle);

        // set mesh velocity
        double Omega = CalculateRotationVelocity(rotAngle);

        // now i am iterating in a container already fully stored, this means i can change any content properties
        // operating in-iterator-step (without offsetting from .ElementBegin())...
        for (ModelPart::NodesContainerType::iterator itr = pNodes->begin(); itr != (pNodes->end() - itrOffset); itr++)
        {
                double X = itr->X();
                double X0 = itr->X0();
                double Y = itr->Y();
                double Y0 = itr->Y0();

                // calculating x and y offset due to rotation 
                double Xoffset = X*cosRotAngle - Y*sinRotAngle - X0;
                double Yoffset = X*sinRotAngle + Y*cosRotAngle - Y0;

                // setting new coordinates through their offset at the actual step of the history
                itr->GetSolutionStepValue(DISPLACEMENT_X, 0) = Xoffset;
                itr->GetSolutionStepValue(DISPLACEMENT_Y, 0) = Yoffset;

                X = X0 + Xoffset;
                Y = Y0 + Yoffset;
                itr->X() = X;
                itr->Y() = Y;

                array_1d<double,3> MeshVelocity(3,0.0);
                double Radius = sqrt( X*X + Y*Y );
                if ( Radius > 1.0e-12 )
                {
                    double Vtang = Omega * Radius;
                    MeshVelocity[0] = - Vtang * Y / Radius; //sinRotAngle;
                    MeshVelocity[1] = Vtang * X / Radius; //cosRotAngle;
                }
                itr->FastGetSolutionStepValue(MESH_VELOCITY) = MeshVelocity;
                if ( itr->IsFixed(VELOCITY_X) || itr->IsFixed(VELOCITY_Y) || itr->IsFixed(VELOCITY_Z) ) itr->FastGetSolutionStepValue(VELOCITY) = MeshVelocity;
        }
}


void WindTurbineRotationUtilities::DestroyCrownElements()
{
	ModelPart::MeshType& rMesh = mrGlobalModelPart.GetMesh();

	for (ModelPart::ElementsContainerType::iterator itr = mCrownElems.begin(); itr != mCrownElems.end(); itr++)
	{
		// removing the element from the model part
		rMesh.RemoveElement( itr->Id() );
	}

	// now cleaning the crownRegion reference list
#if defined(WIND_TURBINE_USE_PARALLEL_EXTENSION)
        mLastKratosGlobalElementId -= mCrownElems.size();
#endif
	mCrownElems.clear();
}

void WindTurbineRotationUtilities::DestroyCrownNodes()
{
        ModelPart::MeshType& rMesh = mrGlobalModelPart.GetMesh();

        for ( ModelPart::NodesContainerType::iterator itr = mrGlobalModelPart.NodesBegin();
        itr != mrGlobalModelPart.NodesEnd();
        itr++)
        {
                // removing the node from the model part
                if ((int)(itr->GetSolutionStepValue(FLAG_VARIABLE)) == WIND_TURBINE_UNKNOWN_REGION   // <- because the new nodes by Tetgen do not have flag associated to them!!
                 || (int)(itr->GetSolutionStepValue(FLAG_VARIABLE)) == WIND_TURBINE_CROWN_REGION)
                    rMesh.RemoveNode( itr->Id() );
        }
}

void WindTurbineRotationUtilities::RegenerateCrownElements2D()
{
	// ****** preparing food for Trigen
        char trigenOptsNormal[] = "PpcYYQj";
        char trigenOptsVerbose[] = "PpcYYVVj";
        char* trigenOpts = trigenOptsNormal;

        if (mEchoLevel > WIND_TURBINE_ECHOLEVEL_NONE)
        {
            trigenOpts = trigenOptsVerbose;    // setting verbosity to the stars...
        }

	// initializing i/o containers
	struct triangulateio inData;
        struct triangulateio outData;
        struct triangulateio vorOutData;

        InitTriangulationDataStructure(inData);
        InitTriangulationDataStructure(outData);
        InitTriangulationDataStructure(vorOutData);

        inData.numberofpoints = mInterfaceNodes.size();
        inData.pointlist = (REAL*) malloc(inData.numberofpoints * 2 * sizeof(REAL));
        inData.pointmarkerlist = (int*) malloc(inData.numberofpoints * sizeof(int));

        inData.holelist = new REAL[2];
        *(inData.holelist) = (REAL) 0.0;
        *(inData.holelist + 1) = (REAL) 0.0;
        inData.numberofholes = 1;

        inData.numberofsegments = mInterfaceNodes.size();   // constrained segments equal the number of the inner plus outer boundary points
        inData.segmentlist = (int*) malloc(inData.numberofsegments * 2 * sizeof(int));
        inData.segmentmarkerlist = (int*) malloc(inData.numberofsegments * sizeof(int));

        ModelPart::NodesContainerType::iterator boundaryNodesItr = mInterfaceNodes.begin();
        for (unsigned int idx = 0; idx < mInterfaceNodes.size(); idx++)
        {
            int auxId = boundaryNodesItr->GetValue(AUX_ID);
            int base = (auxId - 1)*2;
            inData.pointlist[base] = boundaryNodesItr->X();
            inData.pointlist[base + 1] = boundaryNodesItr->Y();

            if (mEchoLevel > WIND_TURBINE_ECHOLEVEL_INFO)
            {
                std::cout << "from interface nodes, the number "<< idx << " has Kratos Id " << boundaryNodesItr->Id() << ", coords. (" << boundaryNodesItr->X() << "," << boundaryNodesItr->Y() << "," << boundaryNodesItr->Z() << "), and AUX_ID " << boundaryNodesItr->GetValue(AUX_ID) << std::endl;
            }

            boundaryNodesItr++;            
        }

        // now parsing the segment counterclockwised scheme and filling the segment list,
        // finding the index of the source container mInterfaceNodes.
        // (NOTE: not using the '-z' switch, this means dereferencing in triangleio starts from 1)

        unsigned int sgmntListIdx = 0;
        for (unsigned int auxIdxPos = 0;
             auxIdxPos < mConstrainedBoundaryNodeAuxIndices.size();
             )
        {
            inData.segmentmarkerlist[sgmntListIdx++] = 1;
            inData.segmentlist[auxIdxPos] = mConstrainedBoundaryNodeAuxIndices.at(auxIdxPos);
            auxIdxPos++;
            inData.segmentlist[auxIdxPos] = mConstrainedBoundaryNodeAuxIndices.at(auxIdxPos);
            auxIdxPos++;
        }

	// ****** FEEDING TRIGEN
        triangulate(trigenOpts, &inData, &outData, &vorOutData);

	// ****** CREATING THE NEW KRATOS ELEMENTS AND UPDATING THE MODEL PART
        unsigned int newElemsNumber = outData.numberoftriangles;
        Properties::Pointer properties = mrGlobalModelPart.GetMesh().pGetProperties(1);
        unsigned int lastElemId = (mrGlobalModelPart.ElementsEnd()-1)->Id();
#if defined(WIND_TURBINE_USE_PARALLEL_EXTENSION)
        lastElemId = mLastKratosGlobalElementId;
        std::cout << "RegenerateCrownElements2D(): first assigned id for this cycle is " << lastElemId+1 << std::endl;
#endif

	mrGlobalModelPart.Elements().reserve(newElemsNumber);	// reserving new space to avoid reallocation while iterating 

        boundaryNodesItr = mInterfaceNodes.begin(); // restart from container base element
        for(unsigned int i = 0; i < newElemsNumber; i++)
        {
            int id = lastElemId + i + 1;
            unsigned int base = i * 3;

            Condition::NodesArrayType geom;
            geom.reserve(3);

            for (unsigned int point = base; point < base+3; point++)
            {
                for (ModelPart::NodesContainerType::iterator refItr = boundaryNodesItr; refItr != mInterfaceNodes.end(); refItr++)
                {
                    if ( refItr->GetValue(AUX_ID) == outData.trianglelist[point] )
                        geom.push_back( *(refItr.base()) );
                }
            }

	    Element::Pointer pElem = mInnerInterfElems.begin()->Create(id, geom, properties);
            (mrGlobalModelPart.Elements()).push_back(pElem);

            mCrownElems.push_back(pElem);  // update crown region element reference list
        }

#if defined(WIND_TURBINE_USE_PARALLEL_EXTENSION)
            std::cout << "RegenerateCrownElements2D(): for next cycle the first global id will be " << mLastKratosGlobalElementId +1 << std::endl;
            mLastKratosGlobalElementId += mCrownElems.size();
#endif

        // how can i relabel all? ..could be enough by sorting with Kratos Container Sort()?
        mrGlobalModelPart.Elements().Sort();

	CleanTriangulationDataStructure(vorOutData);
	CleanTriangulationDataStructure(inData);
	CleanTriangulationDataStructure(outData);
}

void WindTurbineRotationUtilities::RegenerateCrownElements3D()
{
        char tetgenOptsNormal[] = "pYYQ";
        char tetgenOptsVerbose[] = "pYYVV";
        char* tetgenOpts = tetgenOptsNormal;

        if (mEchoLevel > WIND_TURBINE_ECHOLEVEL_NONE)
        {
            tetgenOpts = tetgenOptsVerbose;    // setting verbosity to the stars...
        }

        tetgenio inData, outData;
        inData.initialize();
        outData.initialize();
        inData.firstnumber = 1;   // specifying that internal reference indices start from 1, not 0

        inData.numberofpoints = mInterfaceNodes.size();
        inData.pointlist = new REAL[inData.numberofpoints * 3];   // a base is a triangle
        inData.pointmarkerlist = new int[inData.numberofpoints];

        inData.holelist = new REAL[3];
        *(inData.holelist) = (REAL)0.0;
        *(inData.holelist + 1) = (REAL)0.0;
        *(inData.holelist + 2) = (REAL)0.0; //(-0.75); TODO: clean this owful hard-coded setting
        inData.numberofholes = 1;

        inData.numberoffacets = mNumberOfBoundaryFaces;   // constrained facets
        inData.facetlist = new tetgenio::facet[inData.numberoffacets];
        inData.facetmarkerlist = new int[inData.numberoffacets];

        ModelPart::NodesContainerType::iterator boundaryNodesItr = mInterfaceNodes.begin();
        for (unsigned int idx = 0; idx < mInterfaceNodes.size(); idx++)
        {
            int auxId = boundaryNodesItr->GetValue(AUX_ID);
            int base = (auxId - 1)*3;
            inData.pointlist[base] = boundaryNodesItr->X();
            inData.pointlist[base + 1] = boundaryNodesItr->Y();
            inData.pointlist[base + 2] = boundaryNodesItr->Z();

            inData.pointmarkerlist[auxId-1] = 1;

            if (mEchoLevel > WIND_TURBINE_ECHOLEVEL_INFO)
            {
                std::cout << "setting Kratos point " << boundaryNodesItr->Id()
                          << " coord. (" << inData.pointlist[base] << ", " << inData.pointlist[base+1] << ", " << inData.pointlist[base+2]
                          << ")" << " as reference idx " << auxId << std::endl;
            }
            boundaryNodesItr++;
        }

        // Preparing base faces for Tetgen:
        // parsing the facet counterclockwised scheme and filling the facet list by allocating polygons,
        // finding the index of the source container mInterfaceNodes.

        unsigned int facetListIdx = 0;
        for (unsigned int auxIdxPos = 0;
             auxIdxPos < mConstrainedBoundaryNodeAuxIndices.size();
             )
        {
            tetgenio::polygon polyg;
            polyg.numberofvertices = 3;
            polyg.vertexlist = new int[3];

            tetgenio::facet face;
            face.numberofholes = 0;
            face.holelist = (REAL*)NULL;
            face.numberofpolygons = 1;
            face.polygonlist = new tetgenio::polygon[1];
            face.polygonlist[0] = polyg;

            face.polygonlist[0].vertexlist[0] = mConstrainedBoundaryNodeAuxIndices.at(auxIdxPos++);
            face.polygonlist[0].vertexlist[1] = mConstrainedBoundaryNodeAuxIndices.at(auxIdxPos++);
            face.polygonlist[0].vertexlist[2] = mConstrainedBoundaryNodeAuxIndices.at(auxIdxPos++);

            inData.facetlist[facetListIdx] = face;
            inData.facetmarkerlist[facetListIdx] = 1;
            facetListIdx++;
        }

        // FEEDING TETGEN!
        tetrahedralize(tetgenOpts, &inData, &outData);

        // ****** CREATING THE NEW KRATOS ELEMENTS AND UPDATING THE MODEL PART
        unsigned int newElemsNumber = outData.numberoftetrahedra;

        Properties::Pointer properties = mrGlobalModelPart.GetMesh().pGetProperties(1);
        unsigned int lastElemId = (mrGlobalModelPart.ElementsEnd()-1)->Id();
#if defined(WIND_TURBINE_USE_PARALLEL_EXTENSION)
        lastElemId = mLastKratosGlobalElementId;
        std::cout << "RegenerateCrownElements3D(): first assigned id for this cycle is " << lastElemId+1 << std::endl;
#endif

        if (mEchoLevel > WIND_TURBINE_ECHOLEVEL_INFO)
        {
            std::cout << "Reserving space for " << newElemsNumber << " tetrahedra." << std::endl;
        }
        mrGlobalModelPart.Elements().reserve(newElemsNumber);	// reserving new space to avoid reallocation while iterating

        boundaryNodesItr = mInterfaceNodes.begin(); // restart from container base element
        for(unsigned int i = 0; i < newElemsNumber; i++)
        {
            int id = lastElemId + i + 1;
            unsigned int base = i * 4;

            Condition::NodesArrayType geom;
            geom.reserve(4);

            for (unsigned int point = base; point < base+4; point++)
            {
                for (ModelPart::NodesContainerType::iterator refItr = boundaryNodesItr; refItr != mInterfaceNodes.end(); refItr++)
                {
                    if ( refItr->GetValue(AUX_ID) == outData.tetrahedronlist[point] )
                        geom.push_back( *(refItr.base()) );
                }
            }

            Element::Pointer pElem = mInnerInterfElems.begin()->Create(id, geom, properties);
            (mrGlobalModelPart.Elements()).push_back(pElem);

            mCrownElems.push_back(pElem);  // update crown region element reference list

        }
#if defined(WIND_TURBINE_USE_PARALLEL_EXTENSION)
            std::cout << "RegenerateCrownElements3D(): for next cycle the first global id will be " << mLastKratosGlobalElementId << std::endl;
            mLastKratosGlobalElementId += mCrownElems.size();
#endif

        // resorting elements
        mrGlobalModelPart.Elements().Sort();

        if (mEchoLevel > WIND_TURBINE_ECHOLEVEL_INFO)
        {
            std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> END OF GENERATION" << std::endl;
        }
        // free some memory!!
        //inData.deinitialize();
        //outData.deinitialize();
}


unsigned int WindTurbineRotationUtilities::DecideElementRegion(const unsigned int& nodesNumberPerElement, std::vector<WindTurbineRegionMultiplicity>& regionMultiplicity, unsigned int& edgeOppositeVertex, bool& warning) const
{
	unsigned int destRegion = WIND_TURBINE_UNKNOWN_REGION;
	edgeOppositeVertex = WIND_TURBINE_UNKNOWN_REGION;
        warning = false;

	if (nodesNumberPerElement == 3)
	{
		switch (regionMultiplicity[0].Multiplicity())
                {
                	case 3:
                        	destRegion = regionMultiplicity[0].Region();
                                break;

                        case 2:
                                switch (regionMultiplicity[0].Region())
                                {
                                	case WIND_TURBINE_OUTER_REGION:
                                        	destRegion = WIND_TURBINE_OUTER_INTERF_REGION;
                                                break;
                                        case WIND_TURBINE_INNER_REGION:
                                                destRegion = WIND_TURBINE_INNER_INTERF_REGION;
                                                break;
                                        case WIND_TURBINE_INNER_INTERF_REGION:
                                                if (regionMultiplicity[1].Region() == WIND_TURBINE_INNER_REGION)
                                                	destRegion = WIND_TURBINE_INNER_INTERF_REGION;
                                                else if (regionMultiplicity[1].Region() == WIND_TURBINE_CROWN_REGION || regionMultiplicity[1].Region() == WIND_TURBINE_OUTER_INTERF_REGION)
						{
                                                	destRegion = WIND_TURBINE_CROWN_REGION;
							edgeOppositeVertex = WIND_TURBINE_OUTER_INTERF_REGION;
						}
                                                else
                                                {
                                                        std::cout << "Warning: not handled (switched in case WIND_TURBINE_INNER_INTERF_REGION)" << std::endl;
                                                }
                                                break;
                                        case WIND_TURBINE_OUTER_INTERF_REGION:
                                        	if (regionMultiplicity[1].Region() == WIND_TURBINE_OUTER_REGION)
                                                	destRegion = WIND_TURBINE_OUTER_INTERF_REGION;
                                                else if (regionMultiplicity[1].Region() == WIND_TURBINE_CROWN_REGION || regionMultiplicity[1].Region() == WIND_TURBINE_INNER_INTERF_REGION)
						{
                                                	destRegion = WIND_TURBINE_CROWN_REGION;
							edgeOppositeVertex = WIND_TURBINE_INNER_INTERF_REGION;
						}
                                               	else
                                                        std::cout << "Warning: not handled (switched in case WIND_TURBINE_OUTER_INTERF_REGION)" << std::endl;
                                                break;
					default:
                                                std::cout << "Warning: not handled (switched in the case of \"2\" Multiplicity)" << std::endl;
						break;
                                }
                                break;

			case 1:  // particular case (to be more roboust and avoid fault by GiD user THIS CASE SHOULD BE IMPROVED)
                        	destRegion = WIND_TURBINE_CROWN_REGION;
                                break;

                        default:
                                std::cout << "There's something nasty in your FLAG_VARIABLE assignment in pre-processing phase..." << std::endl;
                        	break;
                }
	}
        else if (nodesNumberPerElement == 4)
        {
            switch (regionMultiplicity[0].Multiplicity())
            {
                    case 4:
                        destRegion = regionMultiplicity[0].Region();
                        break;

                    case 3:   // this case is useful to determine the crown-boundary triangular facets
                        switch (regionMultiplicity[0].Region())
                        {
                            case WIND_TURBINE_INNER_REGION:
                                destRegion = WIND_TURBINE_INNER_INTERF_REGION;
                                break;

                            case WIND_TURBINE_CROWN_REGION:
                                destRegion = WIND_TURBINE_CROWN_REGION;
                                break;

                            case WIND_TURBINE_OUTER_REGION:
                                destRegion = WIND_TURBINE_OUTER_INTERF_REGION;
                                break;

                            case WIND_TURBINE_INNER_INTERF_REGION:
                                if ( regionMultiplicity[1].Region() == WIND_TURBINE_INNER_REGION
                                     || regionMultiplicity[1].Region() == WIND_TURBINE_INTERNAL_CYL_BASE)
                                    destRegion = WIND_TURBINE_INNER_INTERF_REGION;
                                else if (regionMultiplicity[1].Region() == WIND_TURBINE_OUTER_INTERF_REGION   // very likely case
                                     || regionMultiplicity[1].Region() == WIND_TURBINE_EXTERNAL_CYL_BASE
                                     || regionMultiplicity[1].Region() == WIND_TURBINE_CROWN_REGION )   // very likely case
                                {
                                    destRegion = WIND_TURBINE_CROWN_REGION;
                                    edgeOppositeVertex = regionMultiplicity[1].Region();
                                }
                                else
                                {
                                    std::cout << "Warning: not handled (switched in case WIND_TURBINE_INNER_INTERF_REGION)!" << std::endl;
                                    warning = true;
                                }
                                break;

                            case WIND_TURBINE_OUTER_INTERF_REGION:
                                if ( regionMultiplicity[1].Region() == WIND_TURBINE_OUTER_REGION )
                                    destRegion = WIND_TURBINE_OUTER_INTERF_REGION;
                                else if (regionMultiplicity[1].Region() == WIND_TURBINE_INNER_INTERF_REGION
                                     || regionMultiplicity[1].Region() == WIND_TURBINE_INTERNAL_CYL_BASE
                                     || regionMultiplicity[1].Region() == WIND_TURBINE_CROWN_REGION
                                  /* || regionMultiplicity[1].Region() == WIND_TURBINE_EXTERNAL_CYL_BASE*/)   // in this case two faces will be extracted from the element
                                {
                                    destRegion = WIND_TURBINE_CROWN_REGION;
                                    edgeOppositeVertex = regionMultiplicity[1].Region();
                                }
                                else
                                {
                                    std::cout << "Warning: not handled (switched in case WIND_TURBINE_OUTER_INTERF_REGION)!" << std::endl;
                                    warning = true;
                                }
                                break;

                            case WIND_TURBINE_EXTERNAL_CYL_BASE:
                                if ( regionMultiplicity[1].Region() == WIND_TURBINE_OUTER_REGION )
                                    destRegion = WIND_TURBINE_OUTER_INTERF_REGION;
                                else if (regionMultiplicity[1].Region() == WIND_TURBINE_INNER_INTERF_REGION
                                     || regionMultiplicity[1].Region() == WIND_TURBINE_INTERNAL_CYL_BASE
                                     || regionMultiplicity[1].Region() == WIND_TURBINE_CROWN_REGION
                                     || regionMultiplicity[1].Region() == WIND_TURBINE_OUTER_INTERF_REGION) // <- this case brings two faces to be extracted!!!
                                {                                                                                // NEIGHBOUR ELEMENTS TO BE INVOKED:
                                    destRegion = WIND_TURBINE_CROWN_REGION;                                      // the one with only 1 node in the OUTER_REGION
                                    edgeOppositeVertex = regionMultiplicity[1].Region();                         // and matching 3 nodes of the present element
                                                                                                                 // is the one from which to extract the second face correctly
                                }
                                else
                                {
                                    std::cout << "Warning: not handled (switched in case WIND_TURBINE_EXTERNAL_CYL_BASE)!" << std::endl;
                                    warning = true;
                                }
                                break;

                            case WIND_TURBINE_INTERNAL_CYL_BASE:
                                if ( regionMultiplicity[1].Region() == WIND_TURBINE_INNER_REGION )
                                    destRegion = WIND_TURBINE_INNER_INTERF_REGION;
                                else if ( regionMultiplicity[1].Region() == WIND_TURBINE_OUTER_INTERF_REGION
                                     || regionMultiplicity[1].Region() == WIND_TURBINE_EXTERNAL_CYL_BASE
                                     || regionMultiplicity[1].Region() == WIND_TURBINE_CROWN_REGION)
                                {
                                    destRegion = WIND_TURBINE_CROWN_REGION;
                                    edgeOppositeVertex = regionMultiplicity[1].Region();
                                }
                                else
                                {
                                    std::cout << "Warning: not handled (switched in case WIND_TURBINE_INTERNAL_CYL_BASE)!" << std::endl;
                                    warning = true;
                                }
                                break;
                        }
                        break;

                    case 2:   // here other limit configurations are handled, useful mainly to identify crown region elements to be destroyed
                            switch (regionMultiplicity[0].Region())
                            {
                                case WIND_TURBINE_OUTER_REGION:   // unuseful
                                    destRegion = WIND_TURBINE_OUTER_INTERF_REGION;
                                    break;

                                case WIND_TURBINE_CROWN_REGION:
                                    destRegion = WIND_TURBINE_CROWN_REGION;
                                    break;

                                case WIND_TURBINE_INNER_REGION:   // unuseful
                                    destRegion = WIND_TURBINE_INNER_INTERF_REGION;
                                    break;

                                case WIND_TURBINE_INNER_INTERF_REGION:
                                    switch (regionMultiplicity[1].Multiplicity())
                                    {
                                        case 2:
                                            if (regionMultiplicity[1].Region() == WIND_TURBINE_INNER_REGION)
                                                destRegion = WIND_TURBINE_INNER_INTERF_REGION;
                                            else if (regionMultiplicity[1].Region() == WIND_TURBINE_OUTER_INTERF_REGION
                                                     || regionMultiplicity[1].Region() == WIND_TURBINE_EXTERNAL_CYL_BASE
                                                     || regionMultiplicity[1].Region() == WIND_TURBINE_CROWN_REGION)
                                                destRegion = WIND_TURBINE_CROWN_REGION;
                                            else
                                            {
                                                std::cout << "Warning: not handled (part 1 WIND_TURBINE_OUTER_INTERF_REGION)!" << std::endl;
                                                warning = true;
                                            }
                                            break;

                                        case 1:
                                            if ((regionMultiplicity[1].Region() == WIND_TURBINE_INNER_REGION && regionMultiplicity[2].Region() == WIND_TURBINE_INTERNAL_CYL_BASE)
                                                || (regionMultiplicity[2].Region() == WIND_TURBINE_INNER_REGION && regionMultiplicity[1].Region() == WIND_TURBINE_INTERNAL_CYL_BASE))
                                                destRegion = WIND_TURBINE_INNER_INTERF_REGION;
                                            else if ((regionMultiplicity[1].Region() == WIND_TURBINE_EXTERNAL_CYL_BASE && regionMultiplicity[2].Region() == WIND_TURBINE_INTERNAL_CYL_BASE)
                                                || (regionMultiplicity[2].Region() == WIND_TURBINE_EXTERNAL_CYL_BASE && regionMultiplicity[1].Region() == WIND_TURBINE_INTERNAL_CYL_BASE))
                                            {
                                                destRegion = WIND_TURBINE_CROWN_REGION;
                                                edgeOppositeVertex = WIND_TURBINE_EXTERNAL_CYL_BASE;
                                            }
                                            else if ((regionMultiplicity[1].Region() == WIND_TURBINE_OUTER_INTERF_REGION && regionMultiplicity[2].Region() == WIND_TURBINE_INTERNAL_CYL_BASE)
                                                || (regionMultiplicity[2].Region() == WIND_TURBINE_OUTER_INTERF_REGION && regionMultiplicity[1].Region() == WIND_TURBINE_INTERNAL_CYL_BASE))
                                            {
                                                destRegion = WIND_TURBINE_CROWN_REGION;
                                                edgeOppositeVertex = WIND_TURBINE_OUTER_INTERF_REGION;
                                            }
                                            else if ((regionMultiplicity[1].Region() == WIND_TURBINE_INTERNAL_CYL_BASE && regionMultiplicity[2].Region() == WIND_TURBINE_CROWN_REGION)
                                                || (regionMultiplicity[2].Region() == WIND_TURBINE_INTERNAL_CYL_BASE && regionMultiplicity[1].Region() == WIND_TURBINE_CROWN_REGION))
                                            {
                                                destRegion = WIND_TURBINE_CROWN_REGION;
                                                edgeOppositeVertex = WIND_TURBINE_CROWN_REGION;
                                            }
                                            else if ((regionMultiplicity[1].Region() == WIND_TURBINE_OUTER_INTERF_REGION && regionMultiplicity[2].Region() == WIND_TURBINE_CROWN_REGION)
                                                || (regionMultiplicity[2].Region() == WIND_TURBINE_OUTER_INTERF_REGION && regionMultiplicity[1].Region() == WIND_TURBINE_CROWN_REGION)
                                                || (regionMultiplicity[1].Region() == WIND_TURBINE_EXTERNAL_CYL_BASE && regionMultiplicity[2].Region() == WIND_TURBINE_CROWN_REGION)
                                                || (regionMultiplicity[2].Region() == WIND_TURBINE_EXTERNAL_CYL_BASE && regionMultiplicity[1].Region() == WIND_TURBINE_CROWN_REGION)
                                                || (regionMultiplicity[1].Region() == WIND_TURBINE_EXTERNAL_CYL_BASE && regionMultiplicity[2].Region() == WIND_TURBINE_OUTER_INTERF_REGION)
                                                || (regionMultiplicity[2].Region() == WIND_TURBINE_EXTERNAL_CYL_BASE && regionMultiplicity[1].Region() == WIND_TURBINE_OUTER_INTERF_REGION))
                                            {
                                                destRegion = WIND_TURBINE_CROWN_REGION;
                                            }
                                            else
                                            {
                                                std::cout << "Warning: not handled (part 2 WIND_TURBINE_INNER_INTERF_REGION)!" << std::endl;
                                                warning = true;
                                            }
                                            break;

                                        default:
                                            std::cout << "Something wrong in node region assignment in pre-processing!" << std::endl;
                                            warning = true;
                                            break;
                                    }
                                    break;

                                case WIND_TURBINE_OUTER_INTERF_REGION:
                                    switch (regionMultiplicity[1].Multiplicity())
                                    {
                                        case 2:
                                            if (regionMultiplicity[1].Region() == WIND_TURBINE_OUTER_REGION)
                                                destRegion = WIND_TURBINE_OUTER_INTERF_REGION;
                                            else if (regionMultiplicity[1].Region() == WIND_TURBINE_INNER_INTERF_REGION
                                                     || regionMultiplicity[1].Region() == WIND_TURBINE_INTERNAL_CYL_BASE)
                                                destRegion = WIND_TURBINE_CROWN_REGION;
                                            else
                                                std::cout << "Warning: (part 1 WIND_TURBINE_OUTER_INTERF_REGION)!" << std::endl;
                                            break;

                                        case 1:
                                            if ((regionMultiplicity[1].Region() == WIND_TURBINE_INNER_INTERF_REGION && regionMultiplicity[2].Region() == WIND_TURBINE_EXTERNAL_CYL_BASE)
                                                || (regionMultiplicity[2].Region() == WIND_TURBINE_INNER_INTERF_REGION && regionMultiplicity[1].Region() == WIND_TURBINE_EXTERNAL_CYL_BASE))
                                            {
                                                destRegion = WIND_TURBINE_CROWN_REGION;
                                                edgeOppositeVertex = WIND_TURBINE_INNER_INTERF_REGION;
                                            }
                                            else if ((regionMultiplicity[1].Region() == WIND_TURBINE_EXTERNAL_CYL_BASE && regionMultiplicity[2].Region() == WIND_TURBINE_OUTER_REGION)
                                                || (regionMultiplicity[2].Region() == WIND_TURBINE_EXTERNAL_CYL_BASE && regionMultiplicity[1].Region() == WIND_TURBINE_OUTER_REGION))
                                            {
                                                destRegion = WIND_TURBINE_OUTER_INTERF_REGION;
                                            }
                                            else if ((regionMultiplicity[1].Region() == WIND_TURBINE_INNER_INTERF_REGION && regionMultiplicity[2].Region() == WIND_TURBINE_CROWN_REGION)
                                                || (regionMultiplicity[2].Region() == WIND_TURBINE_INNER_INTERF_REGION && regionMultiplicity[1].Region() == WIND_TURBINE_CROWN_REGION))
                                            {
                                                destRegion = WIND_TURBINE_CROWN_REGION;
                                            }
                                            else if ((regionMultiplicity[1].Region() == WIND_TURBINE_CROWN_REGION && regionMultiplicity[2].Region() == WIND_TURBINE_EXTERNAL_CYL_BASE)
                                                || (regionMultiplicity[2].Region() == WIND_TURBINE_CROWN_REGION && regionMultiplicity[1].Region() == WIND_TURBINE_EXTERNAL_CYL_BASE))
                                            {
                                                edgeOppositeVertex = WIND_TURBINE_CROWN_REGION;
                                                destRegion = WIND_TURBINE_CROWN_REGION;
                                            }
                                            else
                                            {
                                                std::cout << "Warning: (part 2 WIND_TURBINE_OUTER_INTERF_REGION)!" << std::endl;
                                                warning = true;
                                            }
                                            break;

                                        default:
                                            std::cout << "Something wrong in node region assignment in pre-processing!" << std::endl;
                                            warning = true;
                                            break;
                                    }
                                    break;

                                case WIND_TURBINE_EXTERNAL_CYL_BASE:
                                    switch (regionMultiplicity[1].Multiplicity())
                                    {
                                        case 2:
                                            if (regionMultiplicity[1].Region() == WIND_TURBINE_OUTER_REGION)
                                                destRegion = WIND_TURBINE_OUTER_INTERF_REGION;
                                            else if (regionMultiplicity[1].Region() == WIND_TURBINE_CROWN_REGION
                                                     || regionMultiplicity[1].Region() == WIND_TURBINE_INTERNAL_CYL_BASE
                                                     || regionMultiplicity[1].Region() == WIND_TURBINE_INNER_INTERF_REGION
                                                     /*|| regionMultiplicity[1].Region() == WIND_TURBINE_OUTER_INTERF_REGION*/)
                                                destRegion = WIND_TURBINE_CROWN_REGION;
                                            else
                                            {
                                                std::cout << "Warning: WIND_TURBINE_EXTERNAL_CYL_BASE (part 1)!" << std::endl;
                                                warning = true;
                                            }
                                            break;

                                        case 1:
                                            if ((regionMultiplicity[1].Region() == WIND_TURBINE_OUTER_INTERF_REGION && regionMultiplicity[2].Region() == WIND_TURBINE_INNER_INTERF_REGION)
                                                || (regionMultiplicity[2].Region() == WIND_TURBINE_OUTER_INTERF_REGION && regionMultiplicity[1].Region() == WIND_TURBINE_INNER_INTERF_REGION))
                                            {
                                                destRegion = WIND_TURBINE_CROWN_REGION;
                                                edgeOppositeVertex = WIND_TURBINE_INNER_INTERF_REGION;
                                            }
                                            else if ((regionMultiplicity[1].Region() == WIND_TURBINE_OUTER_INTERF_REGION && regionMultiplicity[2].Region() == WIND_TURBINE_CROWN_REGION)
                                                || (regionMultiplicity[2].Region() == WIND_TURBINE_OUTER_INTERF_REGION && regionMultiplicity[1].Region() == WIND_TURBINE_CROWN_REGION))
                                            {
                                                     destRegion = WIND_TURBINE_CROWN_REGION;
                                                     edgeOppositeVertex = WIND_TURBINE_CROWN_REGION;
                                            }
                                            else if ((regionMultiplicity[1].Region() == WIND_TURBINE_INTERNAL_CYL_BASE && regionMultiplicity[2].Region() == WIND_TURBINE_INNER_INTERF_REGION)
                                                || (regionMultiplicity[2].Region() == WIND_TURBINE_INTERNAL_CYL_BASE && regionMultiplicity[1].Region() == WIND_TURBINE_INNER_INTERF_REGION)
                                                || (regionMultiplicity[1].Region() == WIND_TURBINE_INTERNAL_CYL_BASE && regionMultiplicity[2].Region() == WIND_TURBINE_CROWN_REGION)
                                                || (regionMultiplicity[2].Region() == WIND_TURBINE_INTERNAL_CYL_BASE && regionMultiplicity[1].Region() == WIND_TURBINE_CROWN_REGION)
                                                || (regionMultiplicity[1].Region() == WIND_TURBINE_INNER_INTERF_REGION && regionMultiplicity[2].Region() == WIND_TURBINE_CROWN_REGION)
                                                || (regionMultiplicity[2].Region() == WIND_TURBINE_INNER_INTERF_REGION && regionMultiplicity[1].Region() == WIND_TURBINE_CROWN_REGION))
                                            {
                                                     destRegion = WIND_TURBINE_CROWN_REGION;
                                            }
                                            else if ((regionMultiplicity[1].Region() == WIND_TURBINE_OUTER_INTERF_REGION && regionMultiplicity[2].Region() == WIND_TURBINE_OUTER_REGION)
                                                || (regionMultiplicity[2].Region() == WIND_TURBINE_OUTER_INTERF_REGION && regionMultiplicity[1].Region() == WIND_TURBINE_OUTER_REGION))
                                               destRegion = WIND_TURBINE_OUTER_INTERF_REGION;
                                            else
                                            {
                                                std::cout << "Warning: WIND_TURBINE_EXTERNAL_CYL_BASE (part 2)!" << std::endl;
                                                warning = true;
                                            }
                                            break;
                                    }
                                    break;

                                case WIND_TURBINE_INTERNAL_CYL_BASE:
                                    switch (regionMultiplicity[1].Multiplicity())
                                    {
                                        case 2:
                                            if (regionMultiplicity[1].Region() == WIND_TURBINE_INNER_REGION
                                              || regionMultiplicity[1].Region() == WIND_TURBINE_INNER_INTERF_REGION)
                                                destRegion = WIND_TURBINE_INNER_INTERF_REGION;
                                            else if (regionMultiplicity[1].Region() == WIND_TURBINE_EXTERNAL_CYL_BASE
                                              || regionMultiplicity[1].Region() == WIND_TURBINE_CROWN_REGION
                                              || regionMultiplicity[1].Region() == WIND_TURBINE_OUTER_INTERF_REGION)
                                                destRegion = WIND_TURBINE_CROWN_REGION;
                                            else
                                            {
                                                std::cout << "Warning: WIND_TURBINE_INTERNAL_CYL_BASE (part 1)!" << std::endl;
                                                warning = true;
                                            }
                                            break;

                                        case 1:
                                            if ((regionMultiplicity[1].Region() == WIND_TURBINE_INNER_REGION && regionMultiplicity[2].Region() == WIND_TURBINE_INNER_INTERF_REGION)
                                                || (regionMultiplicity[2].Region() == WIND_TURBINE_INNER_REGION && regionMultiplicity[1].Region() == WIND_TURBINE_INNER_INTERF_REGION))
                                                destRegion = WIND_TURBINE_INNER_INTERF_REGION;
                                            else if ((regionMultiplicity[1].Region() == WIND_TURBINE_CROWN_REGION && regionMultiplicity[2].Region() == WIND_TURBINE_EXTERNAL_CYL_BASE)
                                                || (regionMultiplicity[2].Region() == WIND_TURBINE_CROWN_REGION && regionMultiplicity[1].Region() == WIND_TURBINE_EXTERNAL_CYL_BASE))
                                            {
                                                destRegion = WIND_TURBINE_CROWN_REGION;
                                            }
                                            else if ((regionMultiplicity[1].Region() == WIND_TURBINE_CROWN_REGION && regionMultiplicity[2].Region() == WIND_TURBINE_INNER_INTERF_REGION)
                                                || (regionMultiplicity[2].Region() == WIND_TURBINE_CROWN_REGION && regionMultiplicity[1].Region() == WIND_TURBINE_INNER_INTERF_REGION))
                                            {
                                                edgeOppositeVertex = WIND_TURBINE_CROWN_REGION;
                                                destRegion = WIND_TURBINE_CROWN_REGION;
                                            }
                                            else if ((regionMultiplicity[1].Region() == WIND_TURBINE_EXTERNAL_CYL_BASE && regionMultiplicity[2].Region() == WIND_TURBINE_INNER_INTERF_REGION)
                                                || (regionMultiplicity[2].Region() == WIND_TURBINE_EXTERNAL_CYL_BASE && regionMultiplicity[1].Region() == WIND_TURBINE_INNER_INTERF_REGION))
                                            {
                                                 destRegion = WIND_TURBINE_CROWN_REGION;
                                                 edgeOppositeVertex = WIND_TURBINE_EXTERNAL_CYL_BASE;
                                            }
                                            else if ((regionMultiplicity[1].Region() == WIND_TURBINE_OUTER_INTERF_REGION && regionMultiplicity[2].Region() == WIND_TURBINE_INNER_INTERF_REGION)
                                                || (regionMultiplicity[2].Region() == WIND_TURBINE_OUTER_INTERF_REGION && regionMultiplicity[1].Region() == WIND_TURBINE_INNER_INTERF_REGION))
                                            {
                                                destRegion = WIND_TURBINE_CROWN_REGION;
                                                edgeOppositeVertex = WIND_TURBINE_OUTER_INTERF_REGION;
                                            }
                                            else
                                            {
                                                std::cout << "Warning: WIND_TURBINE_INTERNAL_CYL_BASE (part 2)!" << std::endl;
                                                warning = true;
                                            }
                                            break;
                                    }
                                    break;

                                default:
                                        std::cout << "Warning: not handled (switched in the case of \"2\" Multiplicity)" << std::endl;
                                        warning = true;
                                        break;
                            }
                            break;

                    case 1:  // particular case (to be more roboust and avoid fault by GiD user THIS CASE SHOULD BE IMPROVED)
                            destRegion = WIND_TURBINE_CROWN_REGION;
                            break;

                    default:
                            std::cout << "There's something nasty in your FLAG_VARIABLE assignment in pre-processing phase..." << std::endl;
                            warning = true;
                            break;
            }
        }
        else
        	std::cout << "Warning: Element type not handled..." << std::endl;

        return destRegion;
}

void WindTurbineRotationUtilities::InitTriangulationDataStructure( triangulateio& tr )
{
        tr.pointlist                  = (REAL*) NULL;
        tr.pointattributelist         = (REAL*) NULL;
        tr.pointmarkerlist            = (int*) NULL;
        tr.numberofpoints             = 0;
        tr.numberofpointattributes    = 0;
        tr.trianglelist               = (int*) NULL;
        tr.triangleattributelist      = (REAL*) NULL;
        tr.trianglearealist           = (REAL*) NULL;
        tr.neighborlist               = (int*) NULL;
        tr.numberoftriangles          = 0;
        tr.numberofcorners            = 3;
        tr.numberoftriangleattributes = 0;
        tr.segmentlist                = (int*) NULL;
        tr.segmentmarkerlist          = (int*) NULL;
        tr.numberofsegments           = 0;
        tr.holelist                   = (REAL*) NULL;
        tr.numberofholes              = 0;
        tr.regionlist                 = (REAL*) NULL;
        tr.numberofregions            = 0;
        tr.edgelist                   = (int*) NULL;
        tr.edgemarkerlist             = (int*) NULL;
        tr.normlist                   = (REAL*) NULL;
        tr.numberofedges              = 0;
}

void WindTurbineRotationUtilities::CleanTriangulationDataStructure( triangulateio& tr )
{
        if(tr.pointlist != NULL) free(tr.pointlist );
        if(tr.pointattributelist != NULL) free(tr.pointattributelist );
        if(tr.pointmarkerlist != NULL) free(tr.pointmarkerlist   );
        if(tr.trianglelist != NULL) free(tr.trianglelist  );
        if(tr.triangleattributelist != NULL) free(tr.triangleattributelist );
        if(tr.trianglearealist != NULL) free(tr.trianglearealist );
        if(tr.neighborlist != NULL) free(tr.neighborlist   );
        if(tr.segmentlist != NULL) free(tr.segmentlist    );
        if(tr.segmentmarkerlist != NULL) free(tr.segmentmarkerlist   );
        if(tr.holelist != NULL) delete[] tr.holelist;
        if(tr.regionlist != NULL) free(tr.regionlist  );
        if(tr.edgelist != NULL) free(tr.edgelist   );
        if(tr.edgemarkerlist != NULL) free(tr.edgemarkerlist   );
        if(tr.normlist != NULL) free(tr.normlist  );
}

double WindTurbineRotationUtilities::CalculateRotationVelocity(double NewRotAngle)
{
    static int Step = 0;
    double RotVelocity;
    double Dt = mrGlobalModelPart.GetProcessInfo().GetValue(DELTA_TIME);

    switch (Step)
    {
    case 0:
        Step++;
        mAngleHistory.push_back(NewRotAngle);
        mAngleHistory.push_back(0.0);
        RotVelocity = NewRotAngle / Dt;
        break;
    case 1:
        Step++;
        mAngleHistory[0] += NewRotAngle;
        mAngleHistory[1] += NewRotAngle;
        RotVelocity = NewRotAngle / Dt;
        break;
    default:
        double OldDt = mrGlobalModelPart.GetProcessInfo().GetPreviousTimeStepInfo(1)[DELTA_TIME];
        double Rho = OldDt / Dt;
        double TimeCoeff = 1.0 / (Dt * Rho * Rho + Dt * Rho);

        RotVelocity = TimeCoeff * (Rho * Rho + 2.0 * Rho) * (NewRotAngle + mAngleHistory[0]);
        RotVelocity -= TimeCoeff * (Rho * Rho + 2.0 * Rho + 1.0) * mAngleHistory[0];
        RotVelocity += TimeCoeff * mAngleHistory[1];

        mAngleHistory[0] += NewRotAngle;
        mAngleHistory[1] += NewRotAngle;
        break;
    }
    return RotVelocity;
}


// this function is for debug purposes
void WindTurbineRotationUtilities::DoExtractFaceNodes(ModelPart& auxModelPart, const int& domainSize)
{
    if (mThisRank == mRemeshingRank)
    {
        //construct a new auxiliary model part
        auxModelPart.SetBufferSize(mrGlobalModelPart.GetBufferSize());
        //mspalart_model_part.Nodes() = mr_model_part.Nodes();
        auxModelPart.SetProcessInfo(mrGlobalModelPart.pGetProcessInfo());
        auxModelPart.SetProperties(mrGlobalModelPart.pProperties());

        std::string ConditionName;
        if (domainSize == 2)
            ConditionName = std::string("Condition2D");
        else
            ConditionName = std::string("Condition3D");

        const Condition& rReferenceCondition = KratosComponents<Condition>::Get(ConditionName);

        for (ModelPart::NodesContainerType::ptr_iterator itr = mInterfaceNodes.ptr_begin();
             itr != mInterfaceNodes.ptr_end();
             itr++)
        {
            auxModelPart.Nodes().push_back(*itr);
        }

        //generating the conditions
        unsigned int condIndex = 0;

        Properties::Pointer properties = auxModelPart.GetMesh().pGetProperties(1);
        ModelPart::NodesContainerType::iterator boundaryNodesItr = mInterfaceNodes.begin();

        for (unsigned int idx = 0;
             idx < mConstrainedBoundaryNodeAuxIndices.size();
            )
        {
            //generate a face condition
            condIndex++;
            Condition::NodesArrayType temp;
            temp.reserve(domainSize);

            for (int pback=0; pback < domainSize; pback++)
            {
                int auxIdx = mConstrainedBoundaryNodeAuxIndices.at(idx++);
                for (ModelPart::NodesContainerType::iterator refItr = boundaryNodesItr; refItr != mInterfaceNodes.end(); refItr++)
                {
                    if ( refItr->GetValue(AUX_ID) == auxIdx )
                        temp.push_back( *(refItr.base()) );
                }
            }

            Condition::Pointer pCond = rReferenceCondition.Create(condIndex, temp, properties);
            auxModelPart.Conditions().push_back( pCond );
        }

        auxModelPart.Elements() = mCrownElems;
    }
    // else
      //auxModelPart.Clear();
}


/// parallel functions
#if defined( WIND_TURBINE_USE_PARALLEL_EXTENSION )

void WindTurbineRotationUtilities::Parallel_DecideRemeshingProcessor()
{
    MPI_Comm_rank(MPI_COMM_WORLD, &mThisRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mNumberOfRanks);

    if (mEchoLevel > WIND_TURBINE_ECHOLEVEL_INFO)
    {
        std::cout << "I am the ID " << mThisRank << " of " << mNumberOfRanks << " processors " <<  std::endl;

        std::cout << "   rank " << mThisRank << "(out of " << mNumberOfRanks << "): " << mInnerNodes.size() << " inner nodes." << std::endl;
        std::cout << "   rank " << mThisRank << "(out of " << mNumberOfRanks << "): " << mCrownElems.size() << " crown elems." << std::endl;
        std::cout << "   rank " << mThisRank << "(out of " << mNumberOfRanks << "): " << mInterfaceNodes.size() << " interface nodes." << std::endl;
        std::cout << "   rank " << mThisRank << "(out of " << mNumberOfRanks << "): " << mConstrainedBoundaryNodeAuxIndices.size() << " nodes aux indices." << std::endl << std::endl;
    }

    int sendData = mInterfaceNodes.size();
    int recvData[mNumberOfRanks];
    MPI_Allgather(&sendData, 1, MPI_INT, recvData, 1, MPI_INT, MPI_COMM_WORLD);
    std::cout << "   rank " << mThisRank << "(out of " << mNumberOfRanks << "): sendData => " << sendData << std::endl;

    if (mEchoLevel > WIND_TURBINE_ECHOLEVEL_INFO)
    {
        std::cout << "   rank " << mThisRank << "(out of " << mNumberOfRanks << "): recvData[] = { ";
        for (int r = 0; r < mNumberOfRanks; r++)
            std::cout << recvData[r] << " ";
        std::cout << "}" << std::endl;
    }

    // here I must be sure that if domains exist, with the same number of boundary nodes, the same remeshing processor ID is picked by each processor
    std::vector<WindTurbineRegionMultiplicity> ordRecvData;
    for (int r = 0; r < mNumberOfRanks; r++)
    {
        WindTurbineRegionMultiplicity procMult(r, recvData[r]);
        ordRecvData.push_back( procMult );
    }

    std::sort(ordRecvData.begin(), ordRecvData.end());
    WindTurbineRegionMultiplicity remesherMult = ordRecvData.at( mNumberOfRanks - 1 );
    for (int p = mNumberOfRanks - 2; p > -1; p--)
    {
        if ( remesherMult.Multiplicity() == ordRecvData.at(p).Multiplicity() && remesherMult.Region() > ordRecvData.at(p).Region())
            remesherMult = ordRecvData.at(p);
    }

    mRemeshingRank = remesherMult.Region();

    if (mEchoLevel > WIND_TURBINE_ECHOLEVEL_NONE)
    {
        std::cout << "   rank " << mThisRank << "(out of " << mNumberOfRanks << "): the target remeshing processor is: " << mRemeshingRank << std::endl;
    }

    // now remarking the interface nodes AUX_IDs
    int startAUX_ID = 1;
    for (int rankTurn = 0; rankTurn < mThisRank; rankTurn++)
    {
        startAUX_ID += recvData[rankTurn];
    }

    std::vector<bool> alreadyMarked(mConstrainedBoundaryNodeAuxIndices.size());
    for (unsigned int markId=0; markId < alreadyMarked.size(); markId++) alreadyMarked[markId] = false;

    for (ModelPart::NodesContainerType::iterator itr = mInterfaceNodes.begin();
         itr != mInterfaceNodes.end();
         itr++)
    {
        for (unsigned int n=0; n < mConstrainedBoundaryNodeAuxIndices.size(); n++)
        {
            if (itr->GetValue(AUX_ID) == mConstrainedBoundaryNodeAuxIndices.at(n) && !alreadyMarked[n])
            {
                mConstrainedBoundaryNodeAuxIndices.at(n) = startAUX_ID;
                alreadyMarked[n] = true;
            }
        }
        itr->GetValue(AUX_ID) = startAUX_ID++;
    }

}

template <class EntitiesContainer>
void WindTurbineRotationUtilities::Parallel_MigrateEntities(const EntitiesContainer& container)
{
    ///////////////////////////////////// Serializing node rough data
    std::string sendData;         // buffer of serialized entities for this rank
    int depth = 0;                // sendData char buffer depth, included the '\0' terminator character
    int* recvRankDepths = NULL;   // number of serialized chars received at the remeshing rank, from each process!!!

    if (mEchoLevel > WIND_TURBINE_ECHOLEVEL_DEEPINFO)
    {
        ////////////////////////////////////////
        // FILE PRINTOUT FOR DEBUGGING
        std::string s;
        std::stringstream out;
        out << mThisRank;
        s = out.str();
        const std::string outFilename = "/tmp/SENDING_nodeinfo_proc_" + s + ".log";
        InspectNodeContainerAndLogToFile(mInterfaceNodes, outFilename);
        // END OF FILE PRINTOUT FOR DEBUGGING
        ////////////////////////////////////////
    }

    Parallel_SerializerSave(container, sendData);
    depth = sendData.size() + 1;
    if (mThisRank == mRemeshingRank)
        recvRankDepths = new int[mNumberOfRanks];

    // make the remesher receive the buffer length from all the other ranks
    MPI_Gather(&depth, 1, MPI_INT, recvRankDepths, 1, MPI_INT, mRemeshingRank, MPI_COMM_WORLD);

    char *recvData = NULL;         // the receiving recipient (this is influent only for the remeshing process)
    int  *recvRankOffsets = NULL;  // offsets from receiving recipient, for each process!!!

    if ( mThisRank == mRemeshingRank )
    {
        recvRankOffsets  = new int[mNumberOfRanks];

        // let's allocate enough space to receive everything from the other processes
        int totalDepth = 0;
        recvRankOffsets[0] = 0;
        for (int rankTurn = 0; rankTurn < mNumberOfRanks; rankTurn++)
        {
            std::cout << ">>> REMESHING PROC >>> allocating " << recvRankDepths[rankTurn] << " chars for rank " << rankTurn << std::endl;
            if (rankTurn) recvRankOffsets[rankTurn] = recvRankOffsets[rankTurn-1] + recvRankDepths[rankTurn-1];
            totalDepth += recvRankDepths[rankTurn];
        }
        recvData = new char[totalDepth];
    }

    char *pData = (char*)sendData.c_str(); // i must pass a non-const pointer to char
    MPI_Gatherv(pData, depth, MPI_CHAR, recvData, recvRankDepths, recvRankOffsets, MPI_CHAR, mRemeshingRank, MPI_COMM_WORLD);

    ///////////////////////////////////////////////
    // before rebuilding and fusing the received structures
    // i also have to gather the number of entities used for reenumeration purpose through unique global model part IDs
    int recvUniqueQuantities[mNumberOfRanks];   // contains the number of elements (of every rank)
                                                // then the inner/outer interface node container offsets (of every rank)
                                                // and finally the depths of the constrained nodes AUX_ID array (of every rank)
    // Determining the first Unique Global Kratos element ID, to be used later by the remesher:
    // gathering the number of elements of every rank
    Parallel_FindLastKratosGlobalElementId();

    // gathering the offset in the inner/outer interface node array, of every rank
    depth = mFirstOuterInterfaceNodeOffset;
    MPI_Gather(&depth, 1, MPI_INT, recvUniqueQuantities, 1, MPI_INT, mRemeshingRank, MPI_COMM_WORLD);
    if (mThisRank == mRemeshingRank)
    {
        for (int rankTurn = 0; rankTurn < mNumberOfRanks; rankTurn++)
        {
            mRankInnerInterfaceOffsets.push_back(recvUniqueQuantities[rankTurn]);
        }
    }

    // gathering the number of the constrained edges AUX_ID (and finally their container) identifying the faces for remeshing, of every rank
    depth = mConstrainedBoundaryNodeAuxIndices.size();
    MPI_Gather(&depth, 1, MPI_INT, recvUniqueQuantities, 1, MPI_INT, mRemeshingRank, MPI_COMM_WORLD);
    int *recvConstrainedNodes = NULL;
    int recvConstrainedNodesOffsets[mNumberOfRanks];
    if (mThisRank == mRemeshingRank)
    {
        // let's allocate enough space to receive everything from the other processes
        int totalDepth = 0;
        recvConstrainedNodesOffsets[0] = 0;
        for (int rankTurn = 0; rankTurn < mNumberOfRanks; rankTurn++)
        {
            if (rankTurn) recvConstrainedNodesOffsets[rankTurn] = recvConstrainedNodesOffsets[rankTurn-1] + recvUniqueQuantities[rankTurn-1];
            totalDepth += recvUniqueQuantities[rankTurn];
        }
        recvConstrainedNodes = new int[totalDepth];
    }
    int *sourceIntegers = (int*)(mConstrainedBoundaryNodeAuxIndices.data());  //warning! std::vector could not be contiguous...;
    MPI_Gatherv(sourceIntegers, depth, MPI_INT, recvConstrainedNodes, recvUniqueQuantities, recvConstrainedNodesOffsets, MPI_INT, mRemeshingRank, MPI_COMM_WORLD);

    // rebuild structures from the series of chars
    if ( mThisRank == mRemeshingRank )
    {
        // extracting base structures of real node objects
        ModelPart::NodesContainerType* receivedInterfaceNodeContainers[mNumberOfRanks];
        for (int rankTurn = 0; rankTurn < mNumberOfRanks; rankTurn++)
        {
            ModelPart::NodesContainerType* nodesContainer = new ModelPart::NodesContainerType();
            ModelPart::NodesContainerType::iterator nItr;
            Parallel_SerializerLoad(*nodesContainer, &recvData[recvRankOffsets[rankTurn]], recvRankDepths[rankTurn]-1);
            std::cout << "remeshing rank: from rank " << rankTurn << " I received " << nodesContainer->size() << " nodes" << std::endl;
            nItr = nodesContainer->begin();
            if (nodesContainer->size())
                std::cout << "   with first/last Kratos IDs(AUX_IDs): " << nItr->Id() << "("<< nItr->GetValue(AUX_ID) << ")/" << (--(nodesContainer->end()))->Id() << "(" <<  (--(nodesContainer->end()))->GetValue(AUX_ID) << ")"<< std::endl;
            receivedInterfaceNodeContainers[rankTurn] = nodesContainer;
        }

        if (mEchoLevel > WIND_TURBINE_ECHOLEVEL_DEEPINFO)
        {
            ////////////////////////////////////////
            // FILE PRINTOUT FOR DEBUGGING
            for (int rankTurn = 0; rankTurn < mNumberOfRanks; rankTurn++)
            {
                std::string s;
                std::stringstream out;
                out << rankTurn;
                s = out.str();
                std::string outFilename = "/tmp/recv_nodeinfo_proc_" + s + ".log";
                InspectNodeContainerAndLogToFile(*(receivedInterfaceNodeContainers[rankTurn]), outFilename);
            }
            // END OF FILE PRINTOUT FOR DEBUGGING
            ////////////////////////////////////////
        }

        // callapsing constrained boundary node AUX_IDs series in one list
        std::vector<int> duplicatePurgeList;
        mConstrainedBoundaryNodeAuxIndices.clear(); // this will be the target list
        for (int rankTurn = 0; rankTurn < mNumberOfRanks; rankTurn++)
        {
            int numberOfRankConstraints = recvUniqueQuantities[rankTurn];
            int rankStartPos = recvConstrainedNodesOffsets[rankTurn];
            std::cout << "Parsing constrained node chain of rank " << rankTurn << "! size: " << numberOfRankConstraints << ", offset: " << rankStartPos << std::endl;
            // START OF DEBUGGING CHAIN PRINT OUT
            for (int j=rankStartPos; j < rankStartPos + numberOfRankConstraints; j+=(mNumberOfNodesPerElement-1))
            {
                switch (mNumberOfNodesPerElement)
                {
                    case 4:
                        std::cout << "(" << recvConstrainedNodes[j] << ", " << recvConstrainedNodes[j+1] << ", " << recvConstrainedNodes[j+2] << ")-";
                        break;

                    case 3:
                        std::cout << "(" << recvConstrainedNodes[j] << ", " << recvConstrainedNodes[j+1] << ")-";
                        break;
                }
            }
            std::cout << std::endl;
            // END OF DEBUGGING CHAIN PRINT OUT

            for (int i = rankStartPos; i < rankStartPos + numberOfRankConstraints; i++)
            {
                bool matched = false;
                std::cout << "     Managing the original AUX_ID " << recvConstrainedNodes[i] << ", cycling rank " << rankTurn << ". Does it match through these rank interface nodes?" << std::endl;
                for ( ModelPart::NodesContainerType::iterator rankNodesItr = receivedInterfaceNodeContainers[rankTurn]->begin();
                      rankNodesItr != receivedInterfaceNodeContainers[rankTurn]->end();
                      rankNodesItr++
                     )
                {
                    int auxId = rankNodesItr->GetValue(AUX_ID);
                    if (recvConstrainedNodes[i] == auxId)
                    {
                        matched = true;
                        int owner = rankNodesItr->GetSolutionStepValue(PARTITION_INDEX);
                        std::cout << "         Found a match (" << auxId << ") with owner (" << owner << ")";
                        if (owner != rankTurn) // this is a ghost node! I must pick the right owner node AUX_ID
                        {
                            for (ModelPart::NodesContainerType::iterator ownerItr = receivedInterfaceNodeContainers[owner]->begin();
                                  ownerItr != receivedInterfaceNodeContainers[owner]->end();
                                  ownerItr++
                                 )
                            {
                                if (rankNodesItr->Id() == ownerItr->Id())
                                {
                                    auxId = ownerItr->GetValue(AUX_ID);
                                    std::cout << " with AUX_ID " << auxId << std::endl;
                                    if (rankTurn == mRemeshingRank)
                                    {
                                        std::cout << ">>> changing the AUX_ID of this ghost node for the remeshing rank (" << rankTurn << ")";
                                        std::cout << "    was " << rankNodesItr->GetValue(AUX_ID);
                                        // changing the AUX_ID of this ghost node for the remeshing rank
                                        ModelPart::NodesContainerType::iterator iNode = mrGlobalModelPart.Nodes().find(rankNodesItr->Id());
                                        if (iNode != mrGlobalModelPart.Nodes().end())
                                            mrGlobalModelPart.Nodes().find(rankNodesItr->Id())->GetValue(AUX_ID) = auxId;
                                        else
                                        {
                                            std::cout << "NOT POSSIBLE! SOMETHING REALLY REALLY NASTY HAS COME!" << std::endl;
                                        }
                                        std::cout << "    is " << mrGlobalModelPart.Nodes().find(rankNodesItr->Id())->GetValue(AUX_ID);
                                        // appending the AUX_ID in the
                                        duplicatePurgeList.push_back(auxId);
                                        std::cout << ">> [rank " << rankTurn << "] adding AUX_ID " << auxId << " of Kratos Id = "<< rankNodesItr->Id() << ", whose owner is " << owner << ", in the duplicateList" << std::endl;
                                    }
                                    break;
                                }
                            }
                        }
                        mConstrainedBoundaryNodeAuxIndices.push_back(auxId);
                        std::cout << std::endl << "[node AUX_ID " << auxId << " added to the Global Constrained series]." << std::endl;
                        break;
                    }

                }
                if (!matched)
                    std::cout << "Warning: AUX_ID = " << recvConstrainedNodes[i] << " didn't match with any entry in rank " << rankTurn << "!" << std::endl;
            }
        }


        // now i must reforge the inner/outer interface nodes container with only the proprietary nodes
        // received from every cluster rank and, accordingly recalculate the global offset pivot
        ModelPart::NodesContainerType outerNodesRefsContainer;
        ModelPart::NodesContainerType innerNodesRefsContainer;

        outerNodesRefsContainer.reserve(mInterfaceNodes.size() - mFirstOuterInterfaceNodeOffset);
        int itemCount = 0;

        std::cout << "" << std::endl;
        for (ModelPart::NodesContainerType::iterator itr = mInterfaceNodes.begin();
             itr != mInterfaceNodes.end();
             itr++)
        {
            std::cout << "APPPENDING old remeshing processor node with Kratos Id " << itr->Id() << ", AUX_ID = " << itr->GetValue(AUX_ID) << ", coords. ("<< itr->X() << ", " << itr->Y() << ", " << itr->Z() << ")";

            if (itemCount < mFirstOuterInterfaceNodeOffset)
            {
                std::cout << " as INNER interface node" << std::endl;
                innerNodesRefsContainer.push_back(*(itr.base()));
            }
            else
            {
                std::cout << " as OUTER interface node" << std::endl;
                outerNodesRefsContainer.push_back(*(itr.base()));
            }

            itemCount++;
        }

        mInterfaceNodes.clear();
        mInterfaceNodes = innerNodesRefsContainer;

        for (int rankTurn = 0; rankTurn < mNumberOfRanks; rankTurn++)
        {
            int itemCount = 0;
            for (ModelPart::NodesContainerType::iterator itr = receivedInterfaceNodeContainers[rankTurn]->begin();
                 itr != receivedInterfaceNodeContainers[rankTurn]->end();
                 itr++)
            {
                if ( itr->GetSolutionStepValue(PARTITION_INDEX) == rankTurn )
                {
                    if (rankTurn != mRemeshingRank) // the remeshing rank already keeps its node references
                    {
                        int auxId = itr->GetValue(AUX_ID);
                        bool isPresent = false;
                        for (unsigned int findIdx = 0; findIdx < duplicatePurgeList.size(); findIdx++)
                        {
                            if (duplicatePurgeList.at(findIdx) == auxId)
                            {
                                isPresent = true;
                                std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!! Found AUX_ID " << auxId << " in the duplicateList" << std::endl;
                                break;
                            }
                        }
                        if (!isPresent)
                        {
                            // creating the new node
                            ModelPart::NodeType::Pointer pNewNode = *(itr.base());  // getting the shared pointer inside the node container

                            std::cout << "Appending new node with older Kratos ID " << itr->Id() << " (AUX_ID = " << pNewNode->GetValue(AUX_ID) << "), coords. ("<< pNewNode->X() << ", " << pNewNode->Y() << ", " << pNewNode->Z() << ")";
                            if ( itemCount < mRankInnerInterfaceOffsets.at(rankTurn) )
                            {
                                std::cout << " as INNER interface node" << std::endl;
                                mInterfaceNodes.push_back(pNewNode);
                            }
                            else
                            {
                                std::cout << " as OUTER interface node" << std::endl;
                                outerNodesRefsContainer.push_back(pNewNode);
                            }
                            // i must attach the new nodes to the model part
//                            mrGlobalModelPart.AddNode(pNewNode);          // trying substituting this .AddNode() with a .push_back() (and then a .Sort() suddently out of this cycle)
                            mrGlobalModelPart.Nodes().push_back(pNewNode);
                        }
                    }
                }
                itemCount++;
            }
        }

        mrGlobalModelPart.Nodes().Sort();

        // appending the outer interface nodes and update the list pivot offset
        mFirstOuterInterfaceNodeOffset = mInterfaceNodes.size();
        mInterfaceNodes.reserve( outerNodesRefsContainer.size() );
        for ( ModelPart::NodesContainerType::iterator itr = outerNodesRefsContainer.begin();
        itr != outerNodesRefsContainer.end();
        itr++)
        {
            mInterfaceNodes.push_back(*(itr.base()));
        }

        // final re-enumeration of nodes and setting of the new ownership
        int startAUX_ID = 1;
        std::vector<bool> alreadyMarked(mConstrainedBoundaryNodeAuxIndices.size());
        for (unsigned int markId=0; markId < alreadyMarked.size(); markId++) alreadyMarked[markId] = false;

        for (ModelPart::NodesContainerType::iterator itr = mInterfaceNodes.begin();
             itr != mInterfaceNodes.end();
             itr++)
        {
            for (unsigned int n=0; n < mConstrainedBoundaryNodeAuxIndices.size(); n++)
            {
                if (itr->GetValue(AUX_ID) == mConstrainedBoundaryNodeAuxIndices.at(n) && !alreadyMarked[n])
                {
                    mConstrainedBoundaryNodeAuxIndices.at(n) = startAUX_ID;
                    alreadyMarked[n] = true;
                }
            }
            itr->GetValue(AUX_ID) = startAUX_ID++;
            itr->GetSolutionStepValue(PARTITION_INDEX) = mThisRank; //marking them with the REMESHING_RANK as the new owner!
        }

        // then i must regenerate the crown elements based on the nodes arrived from the other ranks
        // and assign them to the remeshing model part (I redo it in a whole using the remesher itself)
        mLastKratosGlobalElementId += mCrownElems.size();
        DestroyCrownElements();
        DestroyCrownNodes();
        switch (mNumberOfNodesPerElement)
        {
            case 3:
                RegenerateCrownElements2D();
                break;

            case 4:
                mNumberOfBoundaryFaces = mConstrainedBoundaryNodeAuxIndices.size()/3;
                RegenerateCrownElements3D();
                break;
        }

    }
    else
    {
        // setting the new ownership
        for (ModelPart::NodesContainerType::iterator itr = mInterfaceNodes.begin();
             itr != mInterfaceNodes.end();
             itr++)
        {
            itr->GetSolutionStepValue(PARTITION_INDEX) = mRemeshingRank; //marking them with the REMESHING_RANK as the new owner!
        }
    }

    delete[] recvRankDepths;
    delete[] recvData;
    delete[] recvRankOffsets;
    delete[] recvConstrainedNodes;
}

// Determining the first Unique Global Kratos element ID, to be used later by the remesher
void WindTurbineRotationUtilities::Parallel_FindLastKratosGlobalElementId()
{
    // gathering the number of crown elements of every rank
    int depth = mrGlobalModelPart.Elements().size();
    int numberOfElementsPerRank[mNumberOfRanks];
    MPI_Gather(&depth, 1, MPI_INT, numberOfElementsPerRank, 1, MPI_INT, mRemeshingRank, MPI_COMM_WORLD);
    mLastKratosGlobalElementId = 0;
    if (mThisRank == mRemeshingRank)
    {
        for (int rankTurn = 0; rankTurn < mNumberOfRanks; rankTurn++)
        {
            mLastKratosGlobalElementId += numberOfElementsPerRank[rankTurn];
            // debugging couts
            std:: cout << "°°° rank " << rankTurn << "'s elements number " << numberOfElementsPerRank[rankTurn] << std::endl;
            // end of debugging couts
        }
    }
}

void WindTurbineRotationUtilities::InspectNodeContainerAndLogToFile(ModelPart::NodesContainerType& rNodeContainer, const std::string& outFilename)
{
    std::ofstream logFile;
    const char* ptrFilenameStr = outFilename.c_str();
    logFile.open(ptrFilenameStr, std::ios::out | std::ios::app);

    for ( ModelPart::NodesContainerType::iterator rankItr = rNodeContainer.begin();
          rankItr != rNodeContainer.end();
          rankItr++
         )
    {
        logFile << "Id: " << rankItr->Id() << "\tAUX_ID: " << rankItr->GetValue(AUX_ID) << std::endl;
        logFile << "     dofs: ";
        Node < 3 > ::DofsContainerType& referenceDofs = rankItr->GetDofs();
        for (Node < 3 > ::DofsContainerType::iterator iii = referenceDofs.begin(); iii != referenceDofs.end(); iii++)
        {
            std::cout << "\t\t" << iii->GetVariable() << "(key: " << iii->GetVariable().Key() << "), " << std::endl;
        }
        logFile << std::endl << "           historical database: ";

        // intepolating the data
        unsigned int bufferSize = rankItr->GetBufferSize();
        for (unsigned int step = 0; step < bufferSize; step++)
        {
            double* stepData = rankItr->SolutionStepData().Data(step);
            for (unsigned int j = 0; j < mrGlobalModelPart.GetNodalSolutionStepDataSize(); j++)
            {
                std::cout << stepData[j] << ", ";
            }
        }

        logFile << std::endl << "           non-historical database: " << rankItr->Data() << std::endl;
        // copying non-historical database data and position into the new node


    }
    logFile.flush();
    logFile.close();
}


// function suggested by Riccardo, to mark for deletion possibly the orphan nodes
void WindTurbineRotationUtilities::RemoveLocalNodesWithNoElements()
{
    for (ModelPart::NodesContainerType::iterator itr = mrGlobalModelPart.NodesBegin();
        itr != mrGlobalModelPart.NodesEnd();
        itr++)
    {
        itr->Set(TO_ERASE, true);
    }

    // marking nodes belonging to elements and conditions with TO_ERASE = false
    for (ModelPart::ElementsContainerType::iterator itr = mrGlobalModelPart.ElementsBegin();
        itr != mrGlobalModelPart.ElementsEnd();
        itr++)
    {
        for( unsigned int n = 0; n < mNumberOfNodesPerElement; n++ )
        {
            itr->GetGeometry()[n].Set(TO_ERASE, false);
        }
    }

    for (ModelPart::ConditionsContainerType::iterator itr = mrGlobalModelPart.ConditionsBegin();
        itr != mrGlobalModelPart.ConditionsEnd();
        itr++)
    {
        for( unsigned int n = 0; n < mNumberOfNodesPerElement-1; n++ )
        {
            itr->GetGeometry()[n].Set(TO_ERASE, false);
        }
    }

    int orphanCounter = 0;
    for (ModelPart::NodesContainerType::iterator itr = mrGlobalModelPart.NodesBegin();
        itr != mrGlobalModelPart.NodesEnd();
        itr++)
    {
        if (itr->Is(TO_ERASE) == true)
            orphanCounter++;
    }

    std::cout << orphanCounter << " orphan nodes have been marked for deletion";
    if (mNumberOfRanks > 1)
        std::cout << " in rank " << mThisRank << std::endl;
}

// next two Load/Save functions make use of the Kratos Serializer (which transforms any Kratos object in a series of chars).
template <class EntitiesContainer>
void WindTurbineRotationUtilities::Parallel_SerializerSave(EntitiesContainer& rInputEntities, std::string& buffer)
{
    Kratos::Serializer entitySerializer;
    std::stringstream *serializer_buffer;


    entitySerializer.save("entities", rInputEntities);

    serializer_buffer = (std::stringstream *)entitySerializer.pGetBuffer();
    buffer = std::string(serializer_buffer->str());
}

template <class EntitiesContainer>
void WindTurbineRotationUtilities::Parallel_SerializerLoad(EntitiesContainer& rOutputEntities, std::string& buffer)
{
    Kratos::Serializer entitySerializer;
    std::stringstream *serializer_buffer;

    serializer_buffer = (std::stringstream *)entitySerializer.pGetBuffer();
    serializer_buffer->write((char*)(buffer.c_str()), buffer.size());

    entitySerializer.load("entities", rOutputEntities);
}

template <class EntitiesContainer>
void WindTurbineRotationUtilities::Parallel_SerializerLoad(EntitiesContainer& rOutputEntities, const char* buffer, const int& bufferSize)
{
    Kratos::Serializer entitySerializer;
    std::stringstream *serializer_buffer;

    serializer_buffer = (std::stringstream *)entitySerializer.pGetBuffer();
    serializer_buffer->write(buffer, bufferSize);

    entitySerializer.load("entities", rOutputEntities);
}


#endif  // WIND_TURBINE_USE_PARALLEL_EXTENSION

}
