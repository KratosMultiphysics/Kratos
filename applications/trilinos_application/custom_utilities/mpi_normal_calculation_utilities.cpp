#include "mpi.h"
#include "custom_utilities/mpi_normal_calculation_utilities.h"
#include "utilities/math_utils.h"
#include "includes/deprecated_variables.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////
// public functions
///////////////////////////////////////////////////////////////////////////////


MPINormalCalculationUtils::MPINormalCalculationUtils()
{

}

///////////////////////////////////////////////////////////////////////////////

MPINormalCalculationUtils::~MPINormalCalculationUtils()
{

}

///////////////////////////////////////////////////////////////////////////////

int MPINormalCalculationUtils::Check(ModelPart& rModelPart)
{
    KRATOS_TRY;

    if(NORMAL.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"MPINormalCalculationUtils ERROR: NORMAL Key is 0. Kratos variables were not correctly registered.","");
    if ( !rModelPart.GetNodalSolutionStepVariablesList().Has(NORMAL) )
        KRATOS_THROW_ERROR(std::invalid_argument,"MPINormalCalculationUtils ERROR: ModelPart does not contain NORMAL as nodal solution step data.","");


    if(PARTITION_INDEX.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"MPINormalCalculationUtils ERROR: PARTITION_INDEX Key is 0. Kratos variables were not correctly registered.","");
    if ( !rModelPart.GetNodalSolutionStepVariablesList().Has(PARTITION_INDEX) )
        KRATOS_THROW_ERROR(std::invalid_argument,"MPINormalCalculationUtils ERROR: ModelPart does not contain PARTITION_INDEX as nodal solution step data.","");

    if(AUX_INDEX.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"MPINormalCalculationUtils ERROR: AUX_INDEX Key is 0. Kratos variables were not correctly registered.","");
    if ( !rModelPart.GetNodalSolutionStepVariablesList().Has(AUX_INDEX) )
        KRATOS_THROW_ERROR(std::invalid_argument,"MPINormalCalculationUtils ERROR: ModelPart does not contain AUX_INDEX as nodal solution step data.","");

    if(NODAL_PAUX.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"MPINormalCalculationUtils ERROR: NODAL_PAUX Key is 0. Kratos variables were not correctly registered.","");
    if ( !rModelPart.GetNodalSolutionStepVariablesList().Has(NODAL_PAUX) )
        KRATOS_THROW_ERROR(std::invalid_argument,"MPINormalCalculationUtils ERROR: ModelPart does not contain NODAL_PAUX as nodal solution step data.","");

    // If we reached this point without throwing an error, return 0
    return 0;

    KRATOS_CATCH("");
}


///////////////////////////////////////////////////////////////////////////////

void MPINormalCalculationUtils::OrientFaces(ModelPart& rModelPart,
                                            bool OutwardsPositive)
{
    KRATOS_TRY;

    // Initialize normals as zero
    array_1d<double,3> Zero(3,0.0);
    for(ModelPart::NodeIterator itNode =  rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); itNode++)
    {
        noalias(itNode->FastGetSolutionStepValue(NORMAL)) = Zero;
    }

    // Main loop for elements
    unsigned int ElemSwitchCount = 0;

    for (ModelPart::ElementIterator itElem = rModelPart.ElementsBegin(); itElem != rModelPart.ElementsEnd(); itElem++)
    {
        Geometry< Node<3> >& rGeom = itElem->GetGeometry();
        GeometryData::KratosGeometryType GeoType = rGeom.GetGeometryType();

        if (GeoType == GeometryData::Kratos_Tetrahedra3D4  || GeoType == GeometryData::Kratos_Triangle2D3)
        {
            bool Switched = this->OrientElement(rGeom);
            if (Switched)
                ElemSwitchCount++;

            this->NormalContribution(rGeom);
        }
    }

    // Generate output message, throw error if necessary
    std::stringstream OutMsg;
    if (ElemSwitchCount > 0)
    {
        OutMsg << "Mesh orientation check found " << ElemSwitchCount << " inverted elements." << std::endl;
    }
    else
    {
        OutMsg << "No inverted elements found" << std::endl;
    }


    rModelPart.GetCommunicator().AssembleCurrentData(NORMAL);

    // Main loop for faces
    unsigned int CondSwitchCount = 0;

    for (ModelPart::ConditionIterator itCond = rModelPart.ConditionsBegin(); itCond != rModelPart.ConditionsEnd(); itCond++)
    {
        Geometry< Node<3> >& rGeom = itCond->GetGeometry();
        GeometryData::KratosGeometryType GeoType = rGeom.GetGeometryType();
        array_1d<double,3> FaceNormal(3,0.0);

        if ( GeoType == GeometryData::Kratos_Triangle3D3 )
            FaceNormal3D(FaceNormal,rGeom);
        else if ( GeoType == GeometryData::Kratos_Line2D2 )
            FaceNormal2D(FaceNormal,rGeom);

        const unsigned int NumNodes = rGeom.PointsNumber();
        unsigned int MismatchCount = 0;
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            const array_1d<double,3>& rNormal = rGeom[i].FastGetSolutionStepValue(NORMAL);
            double Dot = FaceNormal[0]*rNormal[0] + FaceNormal[1]*rNormal[1] + FaceNormal[2]*rNormal[2];

            // Note: after orienting elements, nodal normals will always point outwards, therefore
            // (Dot < 0.0) == True => face normal points inside
            if ( (Dot < 0.0) == OutwardsPositive )
                MismatchCount++;
        }

        // Re-orient if the face normal is not aligned with node normals
        if (MismatchCount == NumNodes)
        {
            rGeom(0).swap(rGeom(1));
            CondSwitchCount++;
        }
    }

    if (CondSwitchCount > 0)
    {
        OutMsg << "Mesh orientation check found " << CondSwitchCount << " inverted conditions." << std::endl;
    }
    else
    {
        OutMsg << "No inverted conditions found" << std::endl;
    }

    KRATOS_CATCH("");

}

///////////////////////////////////////////////////////////////////////////////

void MPINormalCalculationUtils::CalculateOnSimplex(ModelPart& rModelPart,
                                                   int Dimension,
                                                   const Variable<double>& rVariable,
                                                   const double rAlpha)
{
    KRATOS_TRY;

    if (Dimension != 3)
        KRATOS_THROW_ERROR(std::invalid_argument,"MPINormalCalculationUtils::CalculateOnSimplex only implemented for 3D, but got Dimension ==  ",Dimension);

    /*
     * Nodes that lie on marked surfaces will be identified with NODAL_PAUX > 0
     * (NODAL_PAUX stores the number of faces on the surface that contain the node)
     * If that is the case, the node will also have a value of AUX_INDEX, which
     * can be used to find its data in the arrays pNormals and pActiveNeigh
     */
    int MaxNeighs = 0;
    int NumNodes = 0;
    this->IdentifyFaces(rModelPart,rVariable,MaxNeighs,NumNodes);

    //
    // Calculate condition normals and store contributions to nodes
    //

    double *pNormals = new double[3*MaxNeighs*NumNodes];
    int *pActiveNeigh = new int[NumNodes];

    this->InitializeNormalData(rModelPart,rVariable,MaxNeighs,pNormals,pActiveNeigh,NumNodes);

    //
    // Mark edges and corners on the mesh
    //

    this->DetectEdges(rModelPart,Dimension,rAlpha,pNormals,pActiveNeigh,NumNodes,MaxNeighs);

    //
    // Deallocate normal data, we no longer need it
    //

    this->FreeNormalData(pNormals,pActiveNeigh);

    //
    // Update nodal normals (using only contributions from marked faces)
    //

    this->UpdateNodeNormals(rModelPart,Dimension,rVariable);

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////

std::string MPINormalCalculationUtils::Info() const
{
    return "MPINormalCalculationUtils";
}


void MPINormalCalculationUtils::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void MPINormalCalculationUtils::PrintData(std::ostream& rOStream) const
{
    rOStream << this->Info();
}


///////////////////////////////////////////////////////////////////////////////
// private functions
///////////////////////////////////////////////////////////////////////////////

void MPINormalCalculationUtils::IdentifyFaces(ModelPart& rModelPart,
                                              const Variable<double>& rVariable,
                                              int& MaxNeigh,
                                              int& NumNodes)
{
    for (ModelPart::NodeIterator iNode = rModelPart.NodesBegin(); iNode != rModelPart.NodesEnd(); iNode++)
        iNode->FastGetSolutionStepValue(NODAL_PAUX) = 0.0;

    for (ModelPart::ConditionIterator iCond = rModelPart.ConditionsBegin(); iCond != rModelPart.ConditionsEnd(); iCond++)
    {
        if ( iCond->GetValue(rVariable) != 0.0 ) // Only use marked conditions to calculate nodal normals
        {
            Geometry< Node<3> >& rGeom = iCond->GetGeometry();
            for (Geometry< Node<3> >::iterator iNode = rGeom.begin(); iNode != rGeom.end(); iNode++)
                (iNode->FastGetSolutionStepValue(NODAL_PAUX) )++;
        }
    }

    rModelPart.GetCommunicator().AssembleCurrentData(NODAL_PAUX);

    // Determine the maximum amount of neighbour faces on a node
    MaxNeigh = 0;
    NumNodes = 0; // Assing node cound as an auxiliary index to identify nodal data
    for (ModelPart::NodeIterator iNode = rModelPart.NodesBegin(); iNode != rModelPart.NodesEnd(); iNode++)
    {
        int NNeigh = static_cast<int>(iNode->FastGetSolutionStepValue(NODAL_PAUX));
        MaxNeigh = NNeigh > MaxNeigh ? NNeigh : MaxNeigh;

        if (NNeigh > 0)
            iNode->FastGetSolutionStepValue(AUX_INDEX) = NumNodes++;
    }

    // Get the maximum from all processes
    rModelPart.GetCommunicator().MaxAll(MaxNeigh);
}


///////////////////////////////////////////////////////////////////////////////


void MPINormalCalculationUtils::InitializeNormalData(ModelPart& rModelPart,
                                                     const Variable<double>& rVariable,
                                                     int MaxNeigh,
                                                     double* pNormals,
                                                     int* pActiveNeigh,
                                                     const int NumNodes)
{
    //
    // Calculate local normal data
    // pNormals: an array containing the normals of all faces around a given node
    // pActiveNeigh: an array containing the number of faces around a given node
    //

    for (int i = 0; i < NumNodes; i++)
        pActiveNeigh[i] = 0;

    // Fill local information
    array_1d<double,3> Normal(3,0.0);
    for (ModelPart::ConditionIterator iCond = rModelPart.ConditionsBegin(); iCond != rModelPart.ConditionsEnd(); iCond++)
    {
        Geometry< Node<3> >& rGeom = iCond->GetGeometry();
        this->FaceNormal3D(Normal,rGeom);
        Normal *= 0.5; // Triangle area is 1/2 of the cross product of its sides
        iCond->SetValue(NORMAL,Normal);

        if ( iCond->GetValue(rVariable) != 0.0 )
        {
            for (Geometry< Node<3> >::iterator iNode = rGeom.begin(); iNode != rGeom.end(); iNode++)
            {
                const int NodeIndex = static_cast<const int>(iNode->FastGetSolutionStepValue(AUX_INDEX));
                //if (NodeIndex == 0) std::cout << rModelPart.GetCommunicator().MyPID() << " found AUX_INDEX 0: node " << iNode->Id() << std::endl;

                int Offset = 3*( MaxNeigh*NodeIndex + pActiveNeigh[NodeIndex] );
                pNormals[Offset]   = Normal[0];
                pNormals[Offset+1] = Normal[1];
                pNormals[Offset+2] = Normal[2];
                pActiveNeigh[NodeIndex]++;
            }
        }
    }

    //
    // Communicate non-local data to owner node
    // I use the communicator's communication schedule to send lists of normals to the node's owner
    //

    Communicator& rComm = rModelPart.GetCommunicator();
    int NumColors = rComm.GetNumberOfColors();

    for (int step = 0; step < NumColors; step++)
    {
        int Destination = rComm.NeighbourIndices()[step];

        if (Destination >= 0) // If Destination == -1 this process skips this communication step
        {
            Communicator::NodesContainerType& rLocalNodes = rComm.LocalMesh(step).Nodes();
            Communicator::NodesContainerType& rGhostNodes = rComm.GhostMesh(step).Nodes();

            //std::cout << rComm.MyPID() << " (talking to " << Destination << ") " << rLocalNodes.size() << " local nodes, " << rGhostNodes.size() << " ghost nodes." << std::endl;

            // Compute length of arrays to send and receive
            int NodeSendSize = 0;
            int DataSendSize = 0;
            for (Communicator::NodesContainerType::iterator iNode = rGhostNodes.begin(); iNode != rGhostNodes.end(); iNode++)
                if (iNode->FastGetSolutionStepValue(NODAL_PAUX) > 0.0)
                {
                    NodeSendSize++;
                    const int NodeIndex = static_cast<const int>( iNode->FastGetSolutionStepValue(AUX_INDEX) );
                    DataSendSize += 3*pActiveNeigh[ NodeIndex ];
                }

            int NodeRecvSize = 0;

            for (Communicator::NodesContainerType::iterator iNode = rLocalNodes.begin(); iNode != rLocalNodes.end(); iNode++)
                if (iNode->FastGetSolutionStepValue(NODAL_PAUX) > 0.0)
                    NodeRecvSize++;

            // Allocate arrays to send/recv, fill send data
            // The last entry in SizeSendBuff is DataSendSize (so the receiving process knows how many doubles to expect)
            int* SizesSendBuff = new int[NodeSendSize+1];
            double* DataSendBuff = new double[DataSendSize];

            SizesSendBuff[NodeSendSize] = DataSendSize;


            int ii = 0;
            int jj = 0;
            for (Communicator::NodesContainerType::iterator iNode = rGhostNodes.begin(); iNode != rGhostNodes.end(); iNode++)
                if (iNode->FastGetSolutionStepValue(NODAL_PAUX) > 0.0)
                {
                    int i = static_cast<int>( iNode->FastGetSolutionStepValue(AUX_INDEX) );
                    int offset = 3*MaxNeigh*i;
                    int n = pActiveNeigh[i];

                    SizesSendBuff[ii++] = n;

                    for (int j = 0; j < 3*n; j++)
                        DataSendBuff[jj++] = pNormals[offset+j];
                }


            // Send/Receive normal data

            MPI_Status Status;

            //std::cout << rComm.MyPID() << " sending " << NodeSendSize+1 << " ints to " << Destination << " (Expecting " << NodeRecvSize+1 << " ints in return)." << std::endl;

            int* SizesRecvBuff = new int[NodeRecvSize+1];

            MPI_Sendrecv(SizesSendBuff, NodeSendSize+1, MPI_INT, Destination, step,
                         SizesRecvBuff, NodeRecvSize+1, MPI_INT, Destination, step,
                         MPI_COMM_WORLD, &Status);

            int DataRecvSize = SizesRecvBuff[NodeRecvSize];
            double* DataRecvBuff = new double[ DataRecvSize ];

            //std::cout << rComm.MyPID() << " sending " << DataSendSize << " doubles to " << Destination << " (Expecting " << DataRecvSize << " doubles in return)." << std::endl;

            MPI_Sendrecv(DataSendBuff, DataSendSize, MPI_DOUBLE, Destination, step,
                         DataRecvBuff, DataRecvSize, MPI_DOUBLE, Destination, step,
                         MPI_COMM_WORLD, &Status);

            //std::cout << rComm.MyPID() << " finished communication with " << Destination << std::endl;

            // Merge received data to local array

            ii = 0;
            jj = 0;
            for (Communicator::NodesContainerType::iterator iNode = rLocalNodes.begin(); iNode != rLocalNodes.end(); iNode++)
                if (iNode->FastGetSolutionStepValue(NODAL_PAUX) > 0.0)
                {
                    int i = static_cast<int>( iNode->FastGetSolutionStepValue(AUX_INDEX) );

                    // Append normal data
                    int NumRecvNormals = SizesRecvBuff[ii++];
                    int offset = 3*( MaxNeigh*i + pActiveNeigh[i] );

                    for (int j = 0; j < 3*NumRecvNormals; j++)
                        pNormals[offset+j] = DataRecvBuff[jj++];


                    // Add number of extra normals
                    pActiveNeigh[i]  += NumRecvNormals;
                }

            // Delete auxiliary arrays
            delete [] SizesSendBuff;
            delete [] SizesRecvBuff;
            delete [] DataSendBuff;
            delete [] DataRecvBuff;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void MPINormalCalculationUtils::FreeNormalData(double* pNormals,
                                               int* pActiveNeigh)
{
    delete [] pNormals;
    delete [] pActiveNeigh;

    pNormals = NULL;
    pActiveNeigh = NULL;
}

///////////////////////////////////////////////////////////////////////////////

void MPINormalCalculationUtils::DetectEdges(ModelPart &rModelPart,
                                            unsigned int Dimension,
                                            double MaxAngle,
                                            double const *pNormals,
                                            int const *pActiveNeigh,
                                            const int NumNodes,
                                            const int MaxNeighs)
{
    const double to_radians = 3.14159265358979323846/180.;
    const double MaxCos = std::cos(to_radians*MaxAngle); // If the cosine between two normals is larger than this, they are considered to belong to different surfaces
    const double NodeFactor = 1.0/Dimension; // Total area is distributed between nodes with this ratio

    std::vector< array_1d<double,3> > ActiveNormals;
    ActiveNormals.reserve(MaxNeighs);

    for (ModelPart::NodeIterator iNode = rModelPart.GetCommunicator().LocalMesh().NodesBegin();
         iNode != rModelPart.GetCommunicator().LocalMesh().NodesEnd(); iNode++)
    {
        if ( iNode->FastGetSolutionStepValue(NODAL_PAUX) > 0.0)
        {
            int AuxID = static_cast<int>( iNode->FastGetSolutionStepValue(AUX_INDEX));
            int NumNormals = pActiveNeigh[AuxID];
            int Offset = 3*MaxNeighs*AuxID;

            ActiveNormals.resize(0);
            array_1d<double,3> AuxN(3,0.0);

            for (int n = 0; n < NumNormals; n++)
            {
                AuxN[0] = pNormals[Offset + 3*n];
                AuxN[1] = pNormals[Offset + 3*n + 1];
                AuxN[2] = pNormals[Offset + 3*n + 2];

                if (ActiveNormals.size() > 0)
                {
                    double norm = std::sqrt( AuxN[0]*AuxN[0] + AuxN[1]*AuxN[1] + AuxN[2]*AuxN[2] );

                    bool Added = false;

                    for (unsigned int m = 0; m < ActiveNormals.size(); m++)
                    {
                        array_1d<double,3>& ThisNormal = ActiveNormals[m];
                        double norm_i = std::sqrt( ThisNormal[0]*ThisNormal[0] + ThisNormal[1]*ThisNormal[1] + ThisNormal[2]*ThisNormal[2] );

                        double Cos = AuxN[0]*ThisNormal[0] + AuxN[1]*ThisNormal[1] + AuxN[2]*ThisNormal[2];
                        Cos /= norm*norm_i;

                        if (Cos > MaxCos) // We consider that the normals belong to the same surface
                        {
                            ThisNormal += AuxN * NodeFactor;
                            Added = true; // I think at this point I should add a break too, but the non-mpi version does exactly this (JC).
                        }
                    }

                    if (!Added)
                        ActiveNormals.push_back(AuxN*NodeFactor);

                }
                else
                {
                    ActiveNormals.push_back(AuxN*NodeFactor);
                }

            }

            //
            // Set flags to detect edges and corners
            //

            switch ( ActiveNormals.size() )
            {
            case 0:

                iNode->FastGetSolutionStepValue(IS_SLIP) = 0.0;
                break;
            case 1:
                iNode->FastGetSolutionStepValue(IS_SLIP) = 10.0;
                break;
            case 2:
                iNode->FastGetSolutionStepValue(IS_SLIP) = 20.0;
                break;
            case 3:
                iNode->FastGetSolutionStepValue(IS_SLIP) = 30.0;
                break;
            default:
                iNode->FastGetSolutionStepValue(IS_SLIP) = 30.0;
                break;
            }
        }
    }

    //
    // Communicate results to other processes
    //

    rModelPart.GetCommunicator().SynchronizeVariable(IS_SLIP);
}

///////////////////////////////////////////////////////////////////////////////

bool MPINormalCalculationUtils::OrientElement(Geometry<Node<3> > &rGeom)
{
    const unsigned int PointIndex = 0;
    const GeometryData::IntegrationMethod Method = GeometryData::GI_GAUSS_1;

    // Re-orient the element if needed
    double DetJ = rGeom.DeterminantOfJacobian(PointIndex,Method);
    if (DetJ < 0.0)
    {
        // swap two nodes to change orientation
        rGeom(0).swap(rGeom(1));
        return true;
    }
    else
        return false;
}

///////////////////////////////////////////////////////////////////////////////

void MPINormalCalculationUtils::NormalContribution(Geometry<Node<3> > &rGeom)
{
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int Dim = rGeom.WorkingSpaceDimension();

    const unsigned int PointIndex = 0;
    const GeometryData::IntegrationMethod Method = GeometryData::GI_GAUSS_1;
    double DetJ = rGeom.DeterminantOfJacobian(PointIndex,Method);

    Geometry< Node<3> >::ShapeFunctionsGradientsType DN_DX;
    rGeom.ShapeFunctionsIntegrationPointsGradients(DN_DX,Method);
    Matrix& rDN_DX = DN_DX[0];

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        array_1d<double,3>& rNormal = rGeom[i].FastGetSolutionStepValue(NORMAL);
        for (unsigned int d = 0; d < Dim; d++)
            rNormal[d] += DetJ*rDN_DX(i,d);
    }
}

///////////////////////////////////////////////////////////////////////////////

void MPINormalCalculationUtils::FaceNormal2D(array_1d<double,3> &An,
                                             Geometry<Node<3> > &rGeometry)
{
    An[0] =   rGeometry[1].Y() - rGeometry[0].Y();
    An[1] = - (rGeometry[1].X() - rGeometry[0].X());
    An[2] =    0.00;

}

///////////////////////////////////////////////////////////////////////////////

void MPINormalCalculationUtils::FaceNormal3D(array_1d<double,3> &An,
                                             Geometry<Node<3> > &rGeometry)
{

    array_1d<double,3> v1,v2;
    v1[0] = rGeometry[1].X() - rGeometry[0].X();
    v1[1] = rGeometry[1].Y() - rGeometry[0].Y();
    v1[2] = rGeometry[1].Z() - rGeometry[0].Z();

    v2[0] = rGeometry[2].X() - rGeometry[0].X();
    v2[1] = rGeometry[2].Y() - rGeometry[0].Y();
    v2[2] = rGeometry[2].Z() - rGeometry[0].Z();

    MathUtils<double>::CrossProduct(An,v1,v2);
}

///////////////////////////////////////////////////////////////////////////////


void MPINormalCalculationUtils::UpdateNodeNormals(ModelPart &rModelPart,
                                                  const unsigned int Dimension,
                                                  const Variable<double> &rVariable)
{
    //
    // Also calculating the nodal area because the serial algorithm also does it
    // NODAL_PAUX (which I use to count the number of faces that contain the node) is overwritten.
    // CHECK IF IT IS ACTUALLY USED
    //
    const double NodeFactor = 1.0/Dimension;

    const array_1d<double,3> Zero(3,0.0);
    for (ModelPart::NodeIterator iNode = rModelPart.NodesBegin(); iNode != rModelPart.NodesEnd(); iNode++)
    {
        iNode->FastGetSolutionStepValue(NODAL_PAUX) = 0.0;
        iNode->FastGetSolutionStepValue(NORMAL) = Zero;
    }

    array_1d<double,3> Normal(3,0.0);
    for (ModelPart::ConditionIterator iCond = rModelPart.ConditionsBegin(); iCond != rModelPart.ConditionsEnd(); iCond++)
    {
        if ( iCond->GetValue(rVariable) != 0.0 ) // Only use marked conditions to calculate nodal normals
        {
            Geometry< Node<3> >& rGeom = iCond->GetGeometry();
            if (Dimension == 2)
            {
                this->FaceNormal2D(Normal,rGeom);
                Normal *= 2.0*NodeFactor;
            }
            else
            {
                this->FaceNormal3D(Normal,rGeom);
                Normal *= 0.5*NodeFactor; // Triangle area is 1/2 of the cross product of its sides
            }

            double NodalArea = std::sqrt(Normal[0]*Normal[0] + Normal[1]*Normal[1] + Normal[2]*Normal[2]);

            for (Geometry< Node<3> >::iterator iNode = rGeom.begin(); iNode != rGeom.end(); iNode++)
            {
                iNode->FastGetSolutionStepValue(NORMAL) += Normal;
                iNode->FastGetSolutionStepValue(NODAL_PAUX) += NodalArea;
            }
        }

    }

    rModelPart.GetCommunicator().AssembleCurrentData(NORMAL);
    rModelPart.GetCommunicator().AssembleCurrentData(NODAL_PAUX);

}

///////////////////////////////////////////////////////////////////////////////
}
