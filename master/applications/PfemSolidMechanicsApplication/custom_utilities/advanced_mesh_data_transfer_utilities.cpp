//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

// Project includes
#include "includes/kratos_flags.h"
#include "custom_utilities/mesh_data_transfer_utilities.hpp"
#include "pfem_solid_mechanics_application_variables.h"

// functions to compute the "advanced"  mesh data transfer
// follows the scheme to compute the super mesh and from the mapping interpolate in a Least squares.


namespace Kratos 
{


   // crec que la solució és

   void MeshDataTransferUtilities::PerformAnAdvancedTransfer( ElementsContainerType& rTemporal_elements, ModelPart& rModelPart, const int& rLeastSquareInterpolation)
   {
      ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();


      std::cout << " [ ADVANCED MESH INTERPOLATION: Beggin " << std::endl;
      struct triangulateio Mesh1; 
      std::vector< std::vector<int> > EdgeNeigh1;
      ConvertAnElementContainerToTriangulate( rTemporal_elements, Mesh1, EdgeNeigh1);
      std::cout << "   Mesh1Written: " << Mesh1.numberofpoints << " nodes "<<  std::endl;

      struct triangulateio Mesh2;
      std::vector< std::vector<int> > EdgeNeigh2;
      ElementsContainerType& rSecondMesh = rModelPart.Elements();
      ConvertAnElementContainerToTriangulate( rSecondMesh, Mesh2, EdgeNeigh2);
      std::cout << "   Mesh2Written: " << Mesh2.numberofpoints << " nodes "<<  std::endl;

      struct triangulateio SuperMesh;
      std::vector<int> EdgeMarkers1, EdgeMarkers2;
      bool fail = GenerateTheSuperMesh( Mesh1, Mesh2, SuperMesh,EdgeNeigh1, EdgeNeigh2, EdgeMarkers1, EdgeMarkers2);
      
      // it is not a problem that the superMesh fail since the information has been interpolated in the previous function.
      if (fail)
      {
         for (unsigned int i = 0; i < 100; i++)
         {
            std::cout << "   FAILED TO GENERATE THE SUPER MESH " << std::endl;
            std::cout << "   ADVANCED MESH INTERPOLATION FAIL]" << std::endl;
         }
         return;

      }
      std::cout << "   Computed the SuperMesh: " << SuperMesh.numberofpoints << " nodes "  << std::endl;

      std::vector<int> Mapping1, Mapping2;
      GenerateTheMapping( Mesh1, SuperMesh, EdgeMarkers1, Mapping1);
      std::cout << "   Generated the first Mapping " << std::endl;
      GenerateTheMapping( Mesh2, SuperMesh, EdgeMarkers2, Mapping2);
      std::cout << "   Generated the second Mapping" << std::endl;


      // DEFINITION OF THE IMPORTANT VECTORS...
      int from, to;
      std::vector<double> Areas(Mesh2.numberoftriangles, 0);
      std::vector<double> MatrixOfDeterminants( Mesh2.numberoftriangles, 0);
      std::vector<double> MatrixOfPlastic( Mesh2.numberoftriangles, 0);
      std::vector<Matrix> MatrixOfStress;
      for ( int elem = 0; elem < Mesh2.numberoftriangles; elem++)
         MatrixOfStress.push_back( ZeroMatrix(3,3));

      ///// VERY IMPORTANT: Complete cauchy stress tensor is the complete kirchoff stress!!!!!!!!!!
      ElementsContainerType::iterator firstelemMesh1 = rTemporal_elements.begin();
      ElementsContainerType::iterator firstelemMesh2 = rSecondMesh.begin();
   
      Variable<Matrix> VariableTensor = KIRCHHOFF_STRESS_TENSOR; // the one that we Get
      Variable<Vector> VariableVector = ELASTIC_LEFT_CAUCHY_FROM_KIRCHHOFF_STRESS; // the one that we Set
      //Variable<double> PlasticVariable = PLASTIC_STRESS_LIKE;
      bool DetF0Interpolation = false;
      if ( rLeastSquareInterpolation == 1)
      {
         // es el base case
      }
      else if ( rLeastSquareInterpolation == 2)
      {
         DetF0Interpolation = true;
      }
      else if ( rLeastSquareInterpolation == 3 )
      {
         VariableTensor = ELASTIC_LEFT_CAUCHY_GREEN_TENSOR;
         VariableVector = ELASTIC_LEFT_CAUCHY_GREEN_VECTOR;
      }
      else if ( rLeastSquareInterpolation == 4)
      {
         VariableTensor = ELASTIC_LEFT_CAUCHY_GREEN_TENSOR;
         VariableVector = ELASTIC_LEFT_CAUCHY_GREEN_VECTOR;
         DetF0Interpolation = true;
      }

      // definitionOfStupidQuantities
      double Area;
      std::vector<Matrix> Result1, Result2;
      std::vector<double> ResultJ; 

      std::cout << "   Previous to compute the mapping. Case: " << rLeastSquareInterpolation << std::endl;
      int nSMelem = SuperMesh.numberoftriangles;
      int exception1 = 0;
      int exception2 = 0;
      for ( int Selem = 0; Selem < nSMelem; Selem++)
      {
         from = Mapping1[Selem]-1;
         to = Mapping2[Selem]-1;

         if ( from < 0) {
            exception1 +=1;
            break;
         }
         if ( to < 0) {
            exception2 +=1;
            break;
         }

         Area = ComputeElementArea( SuperMesh, Selem);
         Areas[to] += Area;

         (firstelemMesh1 + from)->GetValueOnIntegrationPoints( VariableTensor, Result2, CurrentProcessInfo);
         MatrixOfStress[to] += Result2[0]* Area;

         (firstelemMesh1 + from)->GetValueOnIntegrationPoints( DETERMINANT_F, ResultJ, CurrentProcessInfo);
         MatrixOfDeterminants[to] += ResultJ[0]*Area;

         //(firstelemMesh1 + from)->GetValueOnIntegrationPoints( PlasticVariable, ResultJ, CurrentProcessInfo);
         MatrixOfPlastic[to] += ResultJ[0]*Area;
      }
      std::cout << "   End with graving information: " << exception1 << " without from and " << exception2 << " without to " << std::endl;

      ResultJ.resize(1);
      for (int elem = 0; elem < Mesh2.numberoftriangles; elem++)
      {
         if (Areas[elem] < 1e-5)
         {
            break;
         }

         std::vector<Matrix> Result(1);

         Result[0] = (MatrixOfStress[elem] / Areas[elem]);
         // ahora una manera un poco rara de preparar la información para pasar-la
         std::vector<Vector> ThisInput;
         Vector V = ZeroVector(6);
         Matrix SM = Result[0];
         for (int i = 0; i < 3; ++i)
            V(i) = SM(i,i); 
         V(3) = SM(0,1);
         V(4) = SM(0,2);
         V(5) = SM(1,2);
         ThisInput.push_back(V);

         (firstelemMesh2 + elem)->SetValuesOnIntegrationPoints( VariableVector, ThisInput, CurrentProcessInfo);

         if ( DetF0Interpolation) {
            ResultJ[0] = (MatrixOfDeterminants[elem]) / Areas[elem];
            (firstelemMesh2 + elem)->SetValuesOnIntegrationPoints( DETERMINANT_F, ResultJ, CurrentProcessInfo);
         }

         ResultJ[0] = MatrixOfPlastic[elem] / Areas[elem];
         //(firstelemMesh2 + elem)->SetValuesOnIntegrationPoints( PlasticVariable, ResultJ, CurrentProcessInfo);
      }

      ClearTrianglesList( SuperMesh );
      ClearTrianglesList( Mesh1 );
      ClearTrianglesList( Mesh2 );
      std::cout << " end with AdvancedInterpolation ];; " << std::endl; 


   }

   double MeshDataTransferUtilities::ComputeElementArea( const struct triangulateio& rMesh, const int& elem)
   {
      // First I have to get the nodes and the coordinates....
      Matrix Coords = ZeroMatrix(3,3);
      int nod;
      for (unsigned int i = 0; i < 3; i++)
      {
         nod = rMesh.trianglelist[elem*3 +i] - 1;
         Coords(i,0) = rMesh.pointlist[2*nod]; // not sure
         Coords(i,1) = rMesh.pointlist[2*nod+1]; // (nod-1) ?
         Coords(i,2) = 1.0;

      }

      double Area = MathUtils<double>::Det(Coords) /2.0; 
      return Area;

   }

   
   void MeshDataTransferUtilities::GenerateTheMapping( struct triangulateio& rMesh, struct triangulateio& rSuperMesh, std::vector<int>& rEdgeMarkers, std::vector<int>& rMapping)
   {

      int BoundaryTag = 6666;
      int NoInfoTag = 9999;


      rMapping.clear();
      rMapping.resize( rSuperMesh.numberoftriangles, -1);

      int nLeft = rSuperMesh.numberoftriangles;
      int nZeroCase = 0;
      int nFirstCase = 0;
      int nSecondCase = 0;
      int nThirdCase = 0;
      int nLastCall = 0;

      for ( int elem = 0; elem != rSuperMesh.numberoftriangles; elem++ )
      {
         bool found = false;
         // 1. Find the edges....
         std::vector<int> ThisElem(3);
         std::vector<std::vector<int> > PossibleCandidates;
         ThisElem[0] = rSuperMesh.trianglelist[elem*3];
         ThisElem[1] = rSuperMesh.trianglelist[elem*3 + 1];
         ThisElem[2] = rSuperMesh.trianglelist[elem*3 + 2];
         // 2.
         for (int edg = 0; edg != rSuperMesh.numberofedges; edg++) {
            if( EdgeInElement( ThisElem, rSuperMesh.edgelist[edg*2], rSuperMesh.edgelist[edg*2+1])){
               std::vector<int> EdgeCandidates;
               UncodeUniqueFunction( rEdgeMarkers[edg], EdgeCandidates);
               PossibleCandidates.push_back( EdgeCandidates);
            }
            if (PossibleCandidates.size() == 3) {
               //break;
            }
         }

         if ( (PossibleCandidates.size() != 3)  )
         {
            std::cout << " I HAVE CLEARLY A PROBLEM "<< PossibleCandidates.size() << std::endl;
            std::cout << "        CLEARLY ELEM " << ThisElem[0] << " " << ThisElem[1] << " "<< ThisElem[2] << std::endl;
            int thisind = 0;
            for ( int p = 0; p != rSuperMesh.numberofedges  ; p++) {
               if ( EdgeInElement(ThisElem, rSuperMesh.edgelist[p*2], rSuperMesh.edgelist[p*2+1] ) ) {
                  std::cout << "        CLEARLY EDGE " << rSuperMesh.edgelist[2*p] << " " << rSuperMesh.edgelist[2*p+1] << std::endl;
                  std::cout << " the problem might be " << PossibleCandidates[thisind].size() << std::endl;
                  std::cout << "    possible candidates " << PossibleCandidates[thisind][0] << " and " << PossibleCandidates[thisind][1] << std::endl;
                  thisind += 1;
               }
            }
         }

         for (unsigned int p = 0; p != PossibleCandidates.size(); ++p)
         {
            if ( PossibleCandidates[p][1] == BoundaryTag) // the other is good informationi.... a repassar
            {

               found = true;
               rMapping[elem] = PossibleCandidates[p][0];
               nLeft -= 1;
               nZeroCase += 1;
               break;
            }
         }
         if ( ! found)
         {
            std::vector<int> CandidatesV;
            for (unsigned int p = 0; p != PossibleCandidates.size(); p++)
            {
               for (int q = 0; q < 2; q++) {
                  if ( (PossibleCandidates[p][q] != NoInfoTag ) && (PossibleCandidates[p][q] != BoundaryTag) )
                     CandidatesV.push_back( PossibleCandidates[p][q] );
               }
            }
            if ( CandidatesV.size() > 0) {  // first case, all edges were.
               std::sort(CandidatesV.begin(), CandidatesV.end() );
               for (unsigned int p = 0; p != CandidatesV.size()-2; p++)
               {
                  if ( (CandidatesV[p] == CandidatesV[p+1] ) && ( CandidatesV[p] == CandidatesV[p+2] ) )
                  {
                     found = true;
                     rMapping[elem] = CandidatesV[p];
                     nLeft -= 1;
                     nFirstCase += 1;
                     break;
                  }

               }
               if  (!found) { // two edges were

                  for (unsigned int p = 0; p != CandidatesV.size()-1; p++)
                  {
                     if ( (CandidatesV[p] == CandidatesV[p+1] ) )
                     {
                        found = true;
                        rMapping[elem] = CandidatesV[p];
                        nLeft -= 1;
                        nSecondCase += 1;
                        break;
                     }

                  }
               }
               if ( (!found) && (CandidatesV.size() == 2) ) {
                  int Chosen;
                  found = ChooseBetweenTwoCandidates( rSuperMesh, elem, rMesh, CandidatesV , Chosen);
                  //found = ChooseBetweenTwoCandidates( rSuperMesh, elem, rMesh, CandidatesV, Chosen);
                  if (found) {
                     rMapping[elem] = Chosen;
                     nLeft -= 1;
                     nThirdCase += 1;
                     found = true;
                  }
               }
            }
            else {
               rMapping[elem] = -3;
            }
         } // end Found  ( looking for 2 repeated info)
      } // end superMesh elem Loop


      int ThisLoopCleaned = 1;
      int NumberOfLoop = 0;

      while ( ThisLoopCleaned != 0 ) 
      {
         ThisLoopCleaned = 0;
         for (int elem = 0; elem != rSuperMesh.numberoftriangles; elem++) 
         {

            if ( rMapping[elem] < -2)
            {
               //std::cout << " elem " << elem << " rMap " << rMapping[elem] << std::endl;
               int ThisNeigh;
               for (int i = 0; i < 3; ++i)
               {
                  ThisNeigh = rSuperMesh.neighborlist[3*(elem)+i]; // same error
                  if (ThisNeigh > 0)
                  {
                     if ( rMapping[ThisNeigh-1] > 0)
                     {
                        rMapping[elem] = rMapping[ThisNeigh-1];
                        nLeft -= 1;
                        nLastCall +=1;
                        ThisLoopCleaned += 1;
                        break;
                     }
                  }
               }
               //std::cout << " elem " << elem << " rMap " << rMapping[elem] << std::endl;
            }
         }
         std::cout << "      FinallLoop: number " << NumberOfLoop << " Found " << ThisLoopCleaned << " we have left " << nLeft << std::endl;
         NumberOfLoop += 1;

      }

      std::cout << "     nZeroCase " << nZeroCase << " firstCase " << nFirstCase << " secondCase " << nSecondCase  << " thirdCase " << nThirdCase << " lastCall " << nLastCall << std::endl;
   }


   bool MeshDataTransferUtilities::ChooseBetweenTwoCandidates( const struct triangulateio& rSuperMesh, const int& rElem, const struct triangulateio& rMesh, const std::vector<int>& rCandidates, int& rChosen)
   {
      bool found = true;

      int firstCandidate = rCandidates[0];
      int secondCandidate = rCandidates[1];

      std::vector<int> firstElem;
      std::vector<int> secondElem;

      for (int i = 0; i < 3; i++)
      {
         firstElem.push_back( rMesh.trianglelist[ (firstCandidate-1)*3 + i] );
         secondElem.push_back( rMesh.trianglelist[ (secondCandidate-1)*3 + i] );
      }


      std::vector<int> ThisLine;

      for (int i = 0; i < 3; ++i)
      {
         int Nod1 = firstElem[i];
         for (int j = 0 ; j < 3; ++j)
         {
            int Nod2 = secondElem[j];
            if ( Nod1 == Nod2)
               ThisLine.push_back(Nod1);
         }
      }

      if ( ThisLine.size() != 2) {
         std::cout << " TWO MEGA ERROR " << std::endl;
         return false;
      }

      std::vector<double> FirstPoint(2);
      FirstPoint[0] = rMesh.pointlist[ 2*(ThisLine[0]-1)];
      FirstPoint[1] = rMesh.pointlist[ 2*(ThisLine[0]-1) + 1];
      std::vector<double> SecondPoint(2);
      SecondPoint[0] = rMesh.pointlist[ 2*(ThisLine[1]-1)];
      SecondPoint[1] = rMesh.pointlist[ 2*(ThisLine[1]-1) + 1];


      std::vector<double> Line(2);
      Line[0] = SecondPoint[0] - FirstPoint[0];
      Line[1] = SecondPoint[1] - FirstPoint[1];


      std::vector<double> InsidePoint(2, 0.0);
      for ( int i = 0; i < 3; i++) 
      {
         for (int j = 0; j < 2; j++)
         {
            InsidePoint[j] += rSuperMesh.pointlist[ 2*(rSuperMesh.trianglelist[3*(rElem)+i]-1) +j ] / 3.0; // SI?? SEGURO??
         }
      }

      std::vector<double> InsideTo1(2);
      InsideTo1[0] = InsidePoint[0] - FirstPoint[0];
      InsideTo1[1] = InsidePoint[1] - FirstPoint[1];

      std::vector<double> Normal(2);
      Normal[0] = -Line[1];
      Normal[1] = Line[0];

      double product = 0;
      product = Normal[0]*InsideTo1[0] + Normal[1]*InsideTo1[1];

      if ( product > 0)
      {
         rChosen = rCandidates[0];
      }
      else
      {
         rChosen = rCandidates[1];
      }



      return found;

   }


   int MeshDataTransferUtilities::ComputeUniqueFunction( const int& rA, const int& rB)
   {
      int A, B;
      A = rA; B = rB;
      if ( rA < 1 )
         A = 6666;
      if (rB < 1)
         B = 6666;

      int Res = 0;
      Res = A + (A + B - 2) * ( A + B - 1) / 2;

      return Res;

   }

   void MeshDataTransferUtilities::UncodeUniqueFunction( const int& rCoded, std::vector<int>& rDecoded)
   {

      double Coded = double(rCoded);
      rDecoded.clear();
      rDecoded.resize(2);
      if (Coded < 1)
      {
         rDecoded[0] = 9999;
         rDecoded[1] = 9999;
         return;
      }

      double L, M;

      L = 0.5 + sqrt( 2*Coded - 1);
      L = floor(L);
      M = Coded - (L-1)*(L)*0.5;

      rDecoded[0] = int(M);
      rDecoded[1] = int( 1 + L - M);
      std::sort( rDecoded.begin(), rDecoded.end() );

   }

   
   bool MeshDataTransferUtilities::EdgeInElement( std::vector<int> & rElement, int& rN1, int& rN2)
   {
      int NumberOfCoin = 0;
      bool IsInside = false;

      for ( int i = 0; i < 3; i++)
      {
         if (rElement[i] == rN1) {
            NumberOfCoin += 1;
         }
         else if (rElement[i] == rN2) {
            NumberOfCoin += 1;
         }


      }

      if ( NumberOfCoin == 2)
         IsInside = true;

      return IsInside;

   }

   bool MeshDataTransferUtilities::GenerateTheSuperMesh(struct triangulateio& rMesh1, struct triangulateio& rMesh2, struct triangulateio& rSuperMesh, 
         std::vector< std::vector<int> > & rEdgNeigh1, std::vector< std::vector<int> > & rEdgNeigh2,
         std::vector<int>& rEdgeMarkers1, std::vector<int>& rEdgeMarkers2)
   {

      bool fail;

      double Tolerance = 2.0e-6;


      // 1. Generate A List of Unique Nodes
      int size =  rMesh1.numberofpoints + rMesh2.numberofpoints;
      Matrix Nodes = ZeroMatrix(size,size);
      for (int i = 0; i != rMesh1.numberofpoints; i++) {
         Nodes(i,0) = rMesh1.pointlist[2*i];
         Nodes(i,1) = rMesh1.pointlist[2*i+1];
      }


      Vector ThisNode = ZeroVector(2);

      std::vector<int> NodeDictionary( rMesh2.numberofpoints, -1);
      int nUniqueNodes = rMesh1.numberofpoints;
      bool repeated;

      for ( int i = 0; i != rMesh2.numberofpoints; i++) {
         repeated = false;
         ThisNode(0) = rMesh2.pointlist[i*2];
         ThisNode(1) = rMesh2.pointlist[i*2+1];
         // Look if the node is repeated or new. I have to think how i will save the information.
         double distance;
         for ( int prev = 0; prev != rMesh1.numberofpoints; prev++) {
            distance = sqrt( pow(ThisNode(0)-Nodes(prev,0), 2) + pow(ThisNode(1)-Nodes(prev, 1), 2));
            if (distance < Tolerance) {
               repeated = true;
               NodeDictionary[i] = prev+1;
               break;
            }
         }
         if ( ! repeated)
         {
            Nodes( nUniqueNodes, 0) = ThisNode(0);
            Nodes( nUniqueNodes, 1) = ThisNode(1);
            NodeDictionary[i] = nUniqueNodes +1 ;
            nUniqueNodes += 1;
         }

      }

      // 2. Generate A List of Unique Edges ...
      //size = rMesh1.numberofsegments + rMesh2.numberofsegments;
      //Matrix Edges = ZeroMatrix(size,size);
      std::vector<std::vector<int > > Edges;

      int nEdges = rMesh1.numberofsegments;
      std::vector< int > ThisEdge(2);
      for ( int ed = 0; ed != rMesh1.numberofsegments; ed++)
      {
         ThisEdge[0] = rMesh1.segmentlist[ed*2];
         ThisEdge[1] = rMesh1.segmentlist[ed*2 + 1];
         std::sort (ThisEdge.begin(), ThisEdge.end() );
         Edges.push_back(ThisEdge);
      }

      // ahora la segunda malla, que es donde viene lo...
      //Vector EdgeDictionary = ZeroVector( rMesh2.numberofsegments);
      std::vector<int> EdgeDictionary (rMesh2.numberofsegments, -1);
      for ( int ed = 0; ed != rMesh2.numberofsegments; ed++)
      {
         ThisEdge[0] = rMesh2.segmentlist[ed*2];
         ThisEdge[1] = rMesh2.segmentlist[ed*2 + 1];
         // ConvertToTheSuperMesh names
         ThisEdge[0] = NodeDictionary[ ThisEdge[0]-1]; // ARA TINC LA EDGE AMB NOTACIÓ Malla SUPER
         ThisEdge[1] = NodeDictionary[ ThisEdge[1]-1];
         std::sort( ThisEdge.begin(), ThisEdge.end() );
         repeated = false;
         // look of the edge is repeated...
         for ( int pEd = 0; pEd != rMesh1.numberofsegments; pEd++)
         {
            if (  (ThisEdge[0] == Edges[pEd][0] ) && ( ThisEdge[1] ==Edges[pEd][ 1] ) ) 
            {
               repeated = true;
               EdgeDictionary[ed] = pEd + 1;
               break;
            }
         }
         if ( ! repeated) {
            Edges.push_back(ThisEdge);
            EdgeDictionary[ed] = nEdges + 1;
            nEdges += 1;

         }


      }


      // NOW i think i set all the unique edges. It only need to write now the Markers.......


      // THINK ABOUT THE MARKERS
      int ThisNumber = ComputeUniqueFunction(9999, 9999);
      std::vector<int> Markers1 ( nEdges, ThisNumber); //the value of no info
      for ( int ed = 0; ed != rMesh1.numberofsegments; ed++)
      {
         Markers1[ed] = ComputeUniqueFunction( rEdgNeigh1[ed][0], rEdgNeigh1[ed][1]);
      }
      std::vector<int> Markers2 ( nEdges, ThisNumber);
      std::vector<std::vector<int> > EdgeList;
      int  destination;
      for ( int ed = 0; ed != rMesh2.numberofsegments; ed++)
      {
         destination = EdgeDictionary[ed] - 1;
         Markers2[destination] = ComputeUniqueFunction(rEdgNeigh2[ed][0], rEdgNeigh2[ed][1]);
      }

      struct triangulateio MeshAux; 
      ClearTrianglesList( MeshAux );

      // ARA HE D'Escriure la super Meshhhhh
      MeshAux.numberofpoints = nUniqueNodes;
      MeshAux.pointlist = new REAL[MeshAux.numberofpoints*2];
      for ( int node = 0; node < nUniqueNodes; node++)
      {
         MeshAux.pointlist[node*2] = Nodes(node,0);
         MeshAux.pointlist[node*2+1] = Nodes(node,1);
      }

      MeshAux.numberofsegments = nEdges;
      MeshAux.segmentlist = new int [MeshAux.numberofsegments*2];
      for ( int edg = 0; edg < nEdges; edg++)
      {
         MeshAux.segmentlist[edg*2] = Edges[edg][0];
         MeshAux.segmentlist[edg*2+1] = Edges[edg][1];
      }

      fail = MeshWithMarkers( MeshAux, rSuperMesh, Markers1, rEdgeMarkers1); 
      if (fail)
         return fail;
      fail = MeshWithMarkers( MeshAux, rSuperMesh, Markers2, rEdgeMarkers2);
      return fail;

   }

   bool MeshDataTransferUtilities::MeshWithMarkers(struct triangulateio& rMesh, struct triangulateio& rSuperMesh, std::vector<int>& rMarkers, std::vector<int>& rEdgeMarkers)
   {

      bool failed = false;
      rMesh.segmentmarkerlist = new int [rMesh.numberofsegments];

      for ( int i = 0; i != rMesh.numberofsegments; i++)
      {
         rMesh.segmentmarkerlist[i] = rMarkers[i];
      }

      char meshingOptions[255];
      struct triangulateio VorOut;
      ClearTrianglesList(rSuperMesh);
      ClearTrianglesList(VorOut);

      strcpy ( meshingOptions, "penQ");
      try {
         triangulate (meshingOptions,&rMesh,&rSuperMesh,&VorOut);
      }
      catch (int64_t error_Code) {
         if (error_Code != 0) {
            std::cout << "       In tHE CATCH. Correct. I THink We are here and falied" << std::endl;
            failed = true;
            return failed;
         }
      }


      rEdgeMarkers.clear();
      rEdgeMarkers.resize( rSuperMesh.numberofedges);
      for (int i = 0; i != rSuperMesh.numberofedges; ++i) {
         rEdgeMarkers[i] = rSuperMesh.edgemarkerlist[i];
      }
      rMesh.segmentmarkerlist          = (int*) NULL;

      return failed;

   }

   // function copied from somewhere
   void MeshDataTransferUtilities::ClearTrianglesList(struct triangulateio& tr)
   {
      KRATOS_TRY

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
      tr.numberofcorners            = 3; //for three node triangles
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

      KRATOS_CATCH(" ")
   }

   
   void MeshDataTransferUtilities::ConvertAnElementContainerToTriangulate( ElementsContainerType& rTemporal_elements, struct triangulateio& MeshAux, std::vector<std::vector<int> >& rEdNeigh)
   {
      // 1. Get The Number Of Nodes;
      unsigned int nMin = 100, nMax =  0 ;
      for (ElementsContainerType::iterator elem = rTemporal_elements.begin() ; elem != rTemporal_elements.end(); elem++)
      {
         PointsArrayType& rVertices = elem->GetGeometry().Points();
         for ( int node = 0; node < 3; node++)
         { 
            if (rVertices[node].Id() > nMax)
               nMax = rVertices[node].Id();
            if (rVertices[node].Id() < nMin )
               nMin = rVertices[node].Id();
         }
      }

      // 2. Write the Nodes And Elements
      //struct triangulateio MeshAux;
      ClearTrianglesList( MeshAux );
      MeshAux.numberoftriangles = rTemporal_elements.size();
      MeshAux.numberofpoints = nMax;
      MeshAux.pointlist = new REAL[MeshAux.numberofpoints*2];
      MeshAux.trianglelist = new int [MeshAux.numberoftriangles*3];
      int nElem = 0;
      for (ElementsContainerType::iterator elem = rTemporal_elements.begin(); elem != rTemporal_elements.end(); elem++)
      {
         PointsArrayType& rVertices = elem->GetGeometry().Points();
         for (int node = 0; node < 3; node++)
         {
            int ThisNode = rVertices[node].Id();
            MeshAux.trianglelist[nElem*3 + node] = ThisNode;
            MeshAux.pointlist[(ThisNode-1)*2] = rVertices[node].X();
            MeshAux.pointlist[(ThisNode-1)*2+1] = rVertices[node].Y();

         }
         nElem += 1;
      }  // POTSER EL PROBLEMA ÉS AQUEST, PERÒ NO CREC


      // 3. Write the Edges (requiered) and Neighbours (why?)
      std::vector<std::vector<int > > Edges;
      rEdNeigh.clear();

      int nRepeated;

      int nEdges = 0;
      for (int elem = 0; elem != MeshAux.numberoftriangles ; elem++)
      {
         for (int ed = 0; ed < 3; ed++)
         {
            std::vector<int> ThisEdge;
            if ( ed < 2) {
               ThisEdge.push_back( MeshAux.trianglelist[elem*3 + ed] );
               ThisEdge.push_back( MeshAux.trianglelist[elem*3 + ed + 1] );
            }
            else {
               ThisEdge.push_back( MeshAux.trianglelist[elem*3 + ed]);
               ThisEdge.push_back( MeshAux.trianglelist[elem*3 ] );
            }
            std::sort(ThisEdge.begin(), ThisEdge.end() );
            bool repeated = false;
            for (int doneEdges = 0; doneEdges < nEdges; doneEdges++)
            {
               if ( ThisEdge[0] == Edges[doneEdges][0] ) {
                  if ( ThisEdge[1] == Edges[doneEdges][1]) {
                     repeated = true;
                     nRepeated = doneEdges;
                     break;
                  }
               }
            }
            if ( ! repeated) {
               Edges.push_back(ThisEdge);
               std::vector<int> ThisNeigh(2);
               ThisNeigh[0] = elem+1;
               ThisNeigh[1] = -1;
               rEdNeigh.push_back(ThisNeigh);
               nEdges += 1;
            }
            else {
               rEdNeigh[nRepeated][ 1] = elem+1;
            }
         }
      }

      // WRITE THE SEGMENTS
      MeshAux.numberofsegments = nEdges;
      MeshAux.segmentmarkerlist = new int [MeshAux.numberofsegments];
      MeshAux.segmentlist = new int[MeshAux.numberofsegments*2];

      for ( int  Edge = 0; Edge < nEdges; Edge++)
      {
         MeshAux.segmentlist[2*Edge]   = Edges[Edge][ 0];
         MeshAux.segmentlist[2*Edge+1] = Edges[Edge][ 1];
      }


      return;

   }


   // function to view things in the screen.
   void MeshDataTransferUtilities::WriteAMesh( struct triangulateio& rMesh)
   {

      std::cout << " WRITTING A MESH. Nodes " << std::endl;
      for (int i = 0; i != rMesh.numberofpoints; i++)
         std::cout << i << " POINT " << rMesh.pointlist[2*i] << " , " << rMesh.pointlist[2*i+1] << std::endl;

      std::cout << " WRITTING A MESH. Elems " << std::endl;
      for (int i = 0; i != rMesh.numberoftriangles; i++) {
         std::cout << i << " ELEM " << rMesh.trianglelist[3*i] << " , " <<  rMesh.trianglelist[3*i+1] << " , " <<
            rMesh.trianglelist[3*i +2] << std::endl;
      }

      std::cout << " WRITTING A MESH. Eges " << std::endl;
      for (int i = 0; i != rMesh.numberofedges; i++)
      {
         std::cout << i  << " SEGMENT " << rMesh.edgelist[2*i] << " , " << rMesh.edgelist[2*i+1] << std::endl;
      }


   }
} // end namespace Kratos
