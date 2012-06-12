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

        FillNodeRegions();
        FillElementRegions();

        mEchoLevel = WIND_TURBINE_ECHOLEVEL_NONE;   //_ALL;
}

/// Percolate the whole ModelPart and put the elements in the proper region lists
void WindTurbineRotationUtilities::FillElementRegions()
{
	ModelPart::ElementsContainerType::Pointer pWholeElements = mrGlobalModelPart.pElements();

	// getting Elems geometry to initialize the FLAG_VARIABLE marker object (we must track the regions)
	ModelPart::ElementsContainerType::iterator itr = pWholeElements->begin();
	unsigned int nodesNumber = itr->GetGeometry().PointsNumber();
	std::vector<WindTurbineRegionMultiplicity> elemBelongings(WIND_TURBINE_REGION_NUMBER);

        if (mEchoLevel > WIND_TURBINE_ECHOLEVEL_NONE)
        {
            std::cout << "Number of Nodes per Element is " << nodesNumber << std::endl;
        }

	while ( itr != pWholeElements->end() )
	{
		// initialize at each cycle...
		for (unsigned int i=0; i < WIND_TURBINE_REGION_NUMBER; i++)
                    elemBelongings[i] = WindTurbineRegionMultiplicity(i);

                if (mEchoLevel > WIND_TURBINE_ECHOLEVEL_NONE)
                {
                    std::cout << "Processing Element " << itr->Id() << std::endl;
                }
                // set region multiplicity of each region
		for( unsigned int n = 0; n < nodesNumber; n++ )
		{
                    int region = (int)itr->GetGeometry()[n].FastGetSolutionStepValue(FLAG_VARIABLE);
                    elemBelongings[region]++;

                    if (mEchoLevel > WIND_TURBINE_ECHOLEVEL_NONE)
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
                destRegion = DecideElementRegion(nodesNumber, elemBelongings, oppositeVertex, warn);

                if (mEchoLevel > WIND_TURBINE_ECHOLEVEL_INFO)
                {
                    // pushing the element in the proper container member
                    if (warn)
                        std::cout << "Undefined belongings (" << destRegion << ")for element " << itr->Id() << std::endl;
                    else
                    {
                        std::cout << "Assigning Element " << itr->Id() << " to region " << destRegion;
                        if (oppositeVertex)
                            std::cout << " with vertex in " << oppositeVertex;
                    }
                    std::cout << std::endl;
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
                    mConstrainedBoundaryNodeAuxIndices.reserve(nodesNumber - 1);

                    if (nodesNumber == 3)   //Triangle
                    {
                        // this means the opposite vertex lays in WIND_TURBINE_{INNER/OUTER}_INTERF_REGION,
                        // finding it I've also found the counterclockwise oriented basis

                        //loop on all the nodes in the geometry and get the neighbour elements of each node
                        for (unsigned int n=0; n < nodesNumber; n++)
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
                    else if (nodesNumber == 4)  //Thetraedra
                    {
                        //loop on all the nodes in the geometry and get the neighbour elements of each
                        for (unsigned int n=0; n < nodesNumber; n++)
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

//        if (mEchoLevel > WIND_TURBINE_ECHOLEVEL_NONE)
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
	// erasing elements in the crown region
	DestroyCrownElements();

        // erasing also the nodes in the crown region (this just for 3D, but could be also for 2D)
        if (dimensions == 3)
            DestroyCrownNodes();

	// rotating inner nodes and consequently the elements identified by them
        RotateEntities((double)WIND_TURBINE_INNER_REGION, rotAngle, timeStep);
        RotateEntities((double)WIND_TURBINE_INNER_INTERF_REGION, rotAngle, timeStep);

        // remeshing elements in the crown region
        if (dimensions == 3)
            RegenerateCrownElements3D();
        else if (dimensions == 2)
            RegenerateCrownElements2D();
        else
            std::cout << "Geometry not supported." << std::endl;
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
            inData.segmentlist[auxIdxPos] = mConstrainedBoundaryNodeAuxIndices[auxIdxPos];
            auxIdxPos++;
            inData.segmentlist[auxIdxPos] = mConstrainedBoundaryNodeAuxIndices[auxIdxPos];
            auxIdxPos++;
        }

	// ****** FEEDING TRIGEN
        triangulate(trigenOpts, &inData, &outData, &vorOutData);

	// ****** CREATING THE NEW KRATOS ELEMENTS AND UPDATING THE MODEL PART
        unsigned int newElemsNumber = outData.numberoftriangles;
        Properties::Pointer properties = mrGlobalModelPart.GetMesh().pGetProperties(1);
        unsigned int lastElemId = (mrGlobalModelPart.ElementsEnd()-1)->Id();

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
        *(inData.holelist + 2) = (REAL)(-0.75);
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

            face.polygonlist[0].vertexlist[0] = mConstrainedBoundaryNodeAuxIndices[auxIdxPos++];
            face.polygonlist[0].vertexlist[1] = mConstrainedBoundaryNodeAuxIndices[auxIdxPos++];
            face.polygonlist[0].vertexlist[2] = mConstrainedBoundaryNodeAuxIndices[auxIdxPos++];

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
                                                std::cout << "Warning: not handled (switched in the case of \"2\" Multiplicity" << std::endl;
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
        if(tr.holelist != NULL) delete tr.holelist;
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
                int auxIdx = mConstrainedBoundaryNodeAuxIndices[idx++];
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

}
