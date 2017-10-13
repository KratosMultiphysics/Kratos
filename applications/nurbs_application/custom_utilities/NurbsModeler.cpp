//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.3 $
//
// 



// System includes


// External includes 


// Project includes
#include "custom_utilities/NurbsModeler.h"
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d.h"
#include "geometries/nurbs_2d.h"
#include "nurbs_testcase_application.h"
#include "includes/variables.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"


namespace Kratos
{


/**
 * @brief NurbsModeler::ReadModelPart
 * This function reads in a NURBS-Geometry using an existing ModelPart and stores
 * the NURBS-Object as a member variable
 *
 * @param NurbsModelPart: Contains all Control Points
 * @param ConnectivitiesStart: Has to provided the first Control Point which belongs to the NURBS-Patch
 * @param ConnectivitiesEnd: Has to provided the last Control Point which belongs to the NURBS-Patch
 * @param PolynomialDegreeXi: Polynomial Degree in U-/Xi-direction
 * @param PolynomialDegreeEta: Polynomial Degree in V-/Eta-direction
 * @param KnotsXi: Knots in U-/Xi-direction
 * @param KnotsEta: Knots in V-/Eta-direction
 * @param NumberOfCpsU: Number of Control Points in U-/Xi-direction
 * @param NumberOfCpsV: Number of Control Points in V-/Eta-direction
 * @param Weights: Weights of the Control Points
 */
void NurbsModeler::ReadModelPart(ModelPart &NurbsModelPart,
                                 int dimension,
                                 int ConnectivitiesStart,
                                 int ConnectivitiesEnd,
                                 int PolynomialDegreeXi,
                                 int PolynomialDegreeEta,
                                 Vector KnotsXi,
                                 Vector KnotsEta,
                                 int NumberOfCpsU,
                                 int NumberOfCpsV,
                                 Vector Weights)
{

     PointerVector< Node<3> > PatchControlPoints;

     //All CPs are in one .mdpa file, therefore the connectivities (_Start and _End) define the correct CPs corressponding
     //to one patch
      for (int i=ConnectivitiesStart; i< ConnectivitiesEnd+1; i++)
      {
               PatchControlPoints.push_back(NurbsModelPart.pGetNode(i));
      }

NurbsPatchGeometry< Node<3> >::Pointer pNurbsSurface;
if (dimension ==2)
{
    pNurbsSurface = NurbsPatchGeometry< Node<3> >::Pointer( new NurbsPatchGeometry2D< Node<3> >(PatchControlPoints,
                                                                                                                                       Weights,
                                                                                                                                       KnotsXi,
                                                                                                                                       KnotsEta,
                                                                                                                                       PolynomialDegreeXi,
                                                                                                                                       PolynomialDegreeEta,
                                                                                                                                       NumberOfCpsU,
                                                                                                                                       NumberOfCpsV));
}

else if (dimension ==3)
{
    pNurbsSurface = NurbsPatchGeometry< Node<3> >::Pointer( new NurbsPatchGeometry3D< Node<3> >(PatchControlPoints,
                                                                                                                                       Weights,
                                                                                                                                       KnotsXi,
                                                                                                                                       KnotsEta,
                                                                                                                                       PolynomialDegreeXi,
                                                                                                                                       PolynomialDegreeEta,
                                                                                                                                       NumberOfCpsU,
                                                                                                                                       NumberOfCpsV));
}

else
{
    std::cout <<"Error in NurbsModeler::ReadModelPart() : This implementation only support dimension 2 or 3."<<std::endl;
                return;
}

mNurbsPatchGeometry.push_back( pNurbsSurface );

     std::vector< Geometry< Node<3> >::Pointer > Geometries;
     pNurbsSurface->DefineGeometries(Geometries);

    // Add a property
    Properties::Pointer Prop = Properties::Pointer(new Properties(mNurbsPatchGeometry.size()-1) );
    Prop->SetValue(CONDUCTIVITY,1);
    NurbsModelPart.AddProperties(Prop);


    // Add elements
    for (int i=0; i<pNurbsSurface->GeometryNumber(); i++)
    {
        // Iterate over geometries, create an element with each, add to model part
        Element::Pointer NurbsElement = Element::Pointer( new NurbsPoisson2D(NurbsModelPart.NumberOfElements()+1, Geometries[i],Prop));
        NurbsModelPart.AddElement( NurbsElement );

    }


} //NurbsModeler::ReadModelPart




/**
 * @brief NurbsModeler::ClosestPoint
 * Calculates the closest distance of selected Points of a fluid ModelPart (maybe any ModelPart)
 * to a NURBS-Model Part. Therefore the distance will be calculated only for Points which are "close" to
 * the NURBS-Geometry. This closeness will be checked by a BIN-search, which takes for each NURBS-Element half
 * of the maximum distance from 2 Control Points as search criteria. The search itself is based on the
 * Control Points of the NURBS-Model Part and the Nodes of the other ModelPart (FluidMesh)
 * For all Nodes which will be not detected by the BIN-search the highest possible value for a double will be
 * assigned as Distance.
 * This algorithm calculates the Distance to the Fluid-Nodes on element level, so the distance to some Nodes
 * will be calculated to different elements. At the end only the smallest distance will be safed as nodal values
 * in the fluid ModelPart.
 *
 * @param FluidMesh: ModelPart which provides a mesh to which nodes the distance to the NURBS-Geometry shall
 * be calculat * note: Previously the function ReadModelPart has to be called to initialize the NURBS Geometry
 *       (stored in mNurbsPatchGeometry).ed
 * @param NurbsModelPart: NURBS-Geometry to which the distance will be calculated
 *
 * note: Previously the function ReadModelPart has to be called to initialize the NURBS Geometry
 *       (stored in mNurbsPatchGeometry).
 */
void NurbsModeler::ClosestPoint(ModelPart &FluidMesh, ModelPart & NurbsModelPart, double Max_search_radius )
{


 unsigned int NumberOfFluidNodes = FluidMesh.Nodes().size();
 for(
     ModelPart::NodesContainerType::iterator inode = FluidMesh.NodesBegin();
     inode != FluidMesh.NodesEnd();
     inode++)
 {
     inode->GetSolutionStepValue(DISTANCE,0) = Max_search_radius;
 }

 typename StaticBins::SizeType CellSize = 20;
 typename StaticBins::Pointer mp_search_structure = typename StaticBins::Pointer(new StaticBins(FluidMesh.Nodes().ptr_begin(),FluidMesh.Nodes().ptr_end(),CellSize) );

 const unsigned int  max_results = NumberOfFluidNodes;
 PointVector FoundPoints = PointVector(max_results);
 DistanceVector PointDistances = DistanceVector(max_results);

int elecount(0);
for(
    ModelPart::ElementsContainerType::iterator iel = NurbsModelPart.ElementsBegin();
    iel != NurbsModelPart.ElementsEnd();
    iel++
    ) //Loop over Elements
{
    elecount++;
    std::cout<<elecount<< ".Element"<<std::endl;
    std::clock_t start;
    double duration(0.0);
    start = std::clock();


double  x_max(std::numeric_limits<double>::min()),
        y_max(std::numeric_limits<double>::min()),
        z_max(std::numeric_limits<double>::min()),
        x_min(std::numeric_limits<double>::max()),
        y_min(std::numeric_limits<double>::max()),
        z_min(std::numeric_limits<double>::max());

Node<3> work_point;

for (unsigned int i = 0; i < iel->GetGeometry().PointsNumber(); i++) //Loop over all Control Points in the Element Geometry
{
    ModelPart::NodeType& rNode = iel->GetGeometry()[i];

    work_point.X() += rNode.X();
    work_point.Y() += rNode.Y();
    work_point.Z() += rNode.Z();

    //Min and Max Values for Coordinates of the Control Points
    if(rNode.X() < x_min) x_min = rNode.X();
    if(rNode.Y() < y_min) y_min = rNode.Y();
    if(rNode.Z() < z_min) z_min = rNode.Z();
    if(rNode.X() > x_max) x_max = rNode.X();
    if(rNode.Y() > y_max) y_max = rNode.Y();
    if(rNode.Z() > z_max) z_max = rNode.Z();
}


//Center Point of the Control Polygon
work_point.X() = work_point.X()/iel->GetGeometry().PointsNumber();
work_point.Y() = work_point.Y()/iel->GetGeometry().PointsNumber();
work_point.Z() = work_point.Z()/iel->GetGeometry().PointsNumber();

//Search Radius defined as max. distance between two Control Points
double rSearchRadius = sqrt(    (x_max-x_min)*(x_max-x_min)+
                                (y_max-y_min)*(y_max-y_min)+
                                (z_max-z_min)*(z_max-z_min));



int NumberOfNodesInRadius;

if (Max_search_radius < rSearchRadius){
NumberOfNodesInRadius = mp_search_structure->SearchInRadius(work_point, Max_search_radius, FoundPoints.begin(), PointDistances.begin(), max_results);
}
else
{
 NumberOfNodesInRadius = mp_search_structure->SearchInRadius(work_point, rSearchRadius, FoundPoints.begin(), PointDistances.begin(), max_results);
}
//KRATOS_WATCH(work_point);
//KRATOS_WATCH(*iel);
//KRATOS_WATCH(rSearchRadius);

for( PointIterator it_found = FoundPoints.begin(); it_found != FoundPoints.begin() + NumberOfNodesInRadius; it_found++)
{
//    KRATOS_WATCH(**it_found);
    // calculate real distance
    //Initialize Values for computation of exact distance
    Vector GlobalCoordinates(3), ShapeFunctions;
    GlobalCoordinates = ZeroVector(3);
    array_1d<double,3> DistanceCoordinates;
    double distance(0);
    //Calculate closest local Coordinates to Candidates
    iel->GetGeometry().PointLocalCoordinates(DistanceCoordinates,*(*it_found));
    //Calculate shape functions values at the closest possible local coordinates
    iel->GetGeometry().ShapeFunctionsValues(ShapeFunctions, DistanceCoordinates);

    //Calculate Global Coordinates of the closest Point on the NURBS-Surface
    for (unsigned int k=0; //Loop over Control Points in ElementGeometry
         k< iel->GetGeometry().Points().size();
         k++)
    {
        GlobalCoordinates[0] +=  iel->GetGeometry()[k].X() * ShapeFunctions[k];
        GlobalCoordinates[1] +=  iel->GetGeometry()[k].Y() * ShapeFunctions[k];
        GlobalCoordinates[2] +=  iel->GetGeometry()[k].Z() * ShapeFunctions[k];
    }

    //Distance of Candidate to NURBS-Element
    distance = sqrt( (GlobalCoordinates[0] - (*it_found)->X()) * (GlobalCoordinates[0] - (*it_found)->X())
                    +(GlobalCoordinates[1] - (*it_found)->Y()) * (GlobalCoordinates[1] - (*it_found)->Y())
                    +(GlobalCoordinates[2] - (*it_found)->Z()) * (GlobalCoordinates[2] - (*it_found)->Z()));


    //Check sign of Distance
        //Compute base vectors
        Vector gXi, gEta, Normal;
        Normal.resize(3);
        Normal = ZeroVector(3);
        gXi.resize(3);
        gXi = ZeroVector(3);
        gEta.resize(3);
        gEta = ZeroVector(3);

        Matrix Derivatives;
        iel->GetGeometry().ShapeFunctionsLocalGradients(Derivatives,DistanceCoordinates);
        gXi[0] = Derivatives(0,0);
        gXi[1] = Derivatives(0,1);
        gXi[2] = Derivatives(0,2);
        gEta[0] = Derivatives(1,0);
        gEta[1] = Derivatives(1,1);
        gEta[2] = Derivatives(1,2);

        Normal = MathUtils<double>::UnitCrossProduct(gXi,gEta);

        //KRATOS_WATCH(Normal);

        //Check Scalar product of Normal and Distance Vector.
        Vector DistanceVector;
        DistanceVector.resize(3);
        DistanceVector[0] = (*it_found)->X() - GlobalCoordinates[0];
        DistanceVector[1] = (*it_found)->Y() - GlobalCoordinates[1];
        DistanceVector[2] = (*it_found)->Z() - GlobalCoordinates[2];
        double ScalarProduct =  DistanceVector[0] * Normal[0] +
                                DistanceVector[1] * Normal[1] +
                                DistanceVector[2] * Normal[2];

        if (ScalarProduct < 0)
        {distance = distance * (-1);}



    //Check if computed distance is smaller to a maybe previous calculated distance
    if(MathUtils<double>::Abs(distance) < MathUtils<double>::Abs((*it_found)->GetSolutionStepValue(DISTANCE,0)))
    {(*it_found)->GetSolutionStepValue(DISTANCE,0) = distance;


}


        }

duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

std::cout<<elecount<<".Element Closest Point search needed "<< duration<<" seconds"<<std::endl;
    }//Loop over Elements


//Print all Distances
//double dis;
//for(
//    ModelPart::NodesContainerType::iterator inode = FluidMesh.NodesBegin();
//    inode != FluidMesh.NodesEnd();
//    inode++)
//{
//    dis = inode->GetSolutionStepValue(DISTANCE,0);
//   KRATOS_WATCH(*inode);
//   KRATOS_WATCH(dis);
//}


}//Closest Point




/**
 * @brief NurbsModeler::InterpolateDesignVariables
 * Function which is used for post-processing of NURBS-Temperature-Computations.
 * Here the NURBS-local-coordinates of the TriangleModelPart will be used to calculated
 * the correct Temperature-Values at the Triangle-Nodes. Therefore the normal post-process
 * of GiD can be used to illustrate the Results of the NURBS-Temperature-Computation.
 * @param NurbsModelPart: The NURBS-ModelPart on which the computation already has been performed
 * @param TriangleModelPart: TriangelModelPart which contains the local coordinates of the corressponding
 * NURBS-geometry at the Triangle Nodes.
 *
 * note: Previously the function ReadModelPart has to be called to initialize the NURBS Geometry
 *       (stored in mNurbsPatchGeometry).
 */
void NurbsModeler::InterpolateDesignVariables(ModelPart &NurbsModelPart, ModelPart &TriangleModelPart)
{


    double Xi(0.0), Eta(0.0);
    int ElementId(0);
    Element::Pointer NurbsElement;
    double TemperatureAtTriangleNode(0.0);
    Vector NurbsValues;


    for (ModelPart::NodesContainerType::iterator inode = TriangleModelPart.NodesBegin();
         inode != TriangleModelPart.NodesEnd();
         inode++)
        {
        array_1d<double,3>coords;

        coords =inode->FastGetSolutionStepValue(NURBS_COORDINATES,0);



        ElementId = mNurbsPatchGeometry[0]->FindGeometryId(Xi,Eta);
        NurbsElement = NurbsModelPart.pGetElement(ElementId,0);
        Geometry< Node<3> > & NurbsGeometry = NurbsElement->GetGeometry();

        NurbsGeometry.ShapeFunctionsValues(NurbsValues,coords);
        TemperatureAtTriangleNode = 0;
        for(unsigned int i=0; i<NurbsValues.size();i++)
        {
            TemperatureAtTriangleNode += NurbsValues[i] * NurbsGeometry[i].FastGetSolutionStepValue(TEMPERATURE,0);
        }
        inode->GetSolutionStepValue(TEMPERATURE,0) = TemperatureAtTriangleNode;
        }




}




NurbsModeler::NurbsModeler()
{
    NurbsPatchGeometry2D< Node<3> > mNurbsPatchGeometry();

}

NurbsModeler::~NurbsModeler()
{}

}  // namespace Kratos.


