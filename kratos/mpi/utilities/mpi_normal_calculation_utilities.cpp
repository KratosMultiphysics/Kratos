//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

#include "mpi/utilities/mpi_normal_calculation_utilities.h"
#include "includes/checks.h"
#include "includes/data_communicator.h"
#include "includes/deprecated_variables.h"
#include "utilities/math_utils.h"

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

    const auto& r_node = *(rModelPart.NodesBegin());
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NORMAL, r_node);
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PARTITION_INDEX, r_node);
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(AUX_INDEX, r_node);
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NODAL_PAUX, r_node);

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
    const array_1d<double,3> Zero(3,0.0);
    for(ModelPart::NodeIterator itNode =  rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); itNode++)
    {
        noalias(itNode->FastGetSolutionStepValue(NORMAL)) = Zero;
    }

    // Main loop for elements
    unsigned int ElemSwitchCount = 0;

    for (ModelPart::ElementIterator itElem = rModelPart.ElementsBegin(); itElem != rModelPart.ElementsEnd(); itElem++)
    {
        Geometry< Node >& rGeom = itElem->GetGeometry();
        GeometryData::KratosGeometryType GeoType = rGeom.GetGeometryType();

        if (GeoType == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4 || GeoType == GeometryData::KratosGeometryType::Kratos_Triangle2D3)
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
        Geometry< Node >& rGeom = itCond->GetGeometry();
        GeometryData::KratosGeometryType GeoType = rGeom.GetGeometryType();
        array_1d<double,3> FaceNormal(3,0.0);

        if (GeoType == GeometryData::KratosGeometryType::Kratos_Triangle3D3)
            FaceNormal3D(FaceNormal,rGeom);
        else if (GeoType == GeometryData::KratosGeometryType::Kratos_Line2D2)
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

    KRATOS_ERROR_IF_NOT(Dimension == 3)
    << "MPINormalCalculationUtils::CalculateOnSimplex only implemented for 3D, but got Dimension ==  "
    << Dimension << "." << std::endl;

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

    std::vector<double> normals(3*MaxNeighs*NumNodes);
    std::vector<int> active_neigh(NumNodes);

    this->InitializeNormalData(rModelPart,rVariable,MaxNeighs,normals,active_neigh,NumNodes);

    //
    // Mark edges and corners on the mesh
    //

    this->DetectEdges(rModelPart,Dimension,rAlpha,normals,active_neigh,NumNodes,MaxNeighs);

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
            Geometry< Node >& rGeom = iCond->GetGeometry();
            for (Geometry< Node >::iterator iNode = rGeom.begin(); iNode != rGeom.end(); iNode++)
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
    MaxNeigh = rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(MaxNeigh);
}


///////////////////////////////////////////////////////////////////////////////


void MPINormalCalculationUtils::InitializeNormalData(ModelPart& rModelPart,
                                                     const Variable<double>& rVariable,
                                                     int MaxNeigh,
                                                     std::vector<double>& rNormals,
                                                     std::vector<int>& rActiveNeigh,
                                                     const int NumNodes)
{
    //
    // Calculate local normal data
    // rNormals: an array containing the normals of all faces around a given node
    // rActiveNeigh: an array containing the number of faces around a given node
    //

    for (int i = 0; i < NumNodes; i++)
        rActiveNeigh[i] = 0;

    // Fill local information
    array_1d<double,3> normal(3,0.0);
    for (ModelPart::ConditionIterator i_condition = rModelPart.ConditionsBegin(); i_condition != rModelPart.ConditionsEnd(); i_condition++)
    {
        Geometry< Node >& r_geometry = i_condition->GetGeometry();
        this->FaceNormal3D(normal,r_geometry);
        normal *= 0.5; // Triangle area is 1/2 of the cross product of its sides
        i_condition->SetValue(NORMAL,normal);

        if ( i_condition->GetValue(rVariable) != 0.0 )
        {
            for (Geometry< Node >::iterator i_node = r_geometry.begin(); i_node != r_geometry.end(); i_node++)
            {
                const int node_index = static_cast<int>(i_node->FastGetSolutionStepValue(AUX_INDEX));
                //if (node_index == 0) std::cout << rModelPart.GetCommunicator().MyPID() << " found AUX_INDEX 0: node " << i_node->Id() << std::endl;

                int offset = 3*( MaxNeigh*node_index + rActiveNeigh[node_index] );
                rNormals[offset]   = normal[0];
                rNormals[offset+1] = normal[1];
                rNormals[offset+2] = normal[2];
                rActiveNeigh[node_index]++;
            }
        }
    }

    //
    // Communicate non-local data to owner node
    // I use the communicator's communication schedule to send lists of normals to the node's owner
    //

    Communicator& r_comm = rModelPart.GetCommunicator();
    int NumColors = r_comm.GetNumberOfColors();

    for (int step = 0; step < NumColors; step++)
    {
        const int destination = r_comm.NeighbourIndices()[step];

        if (destination >= 0) // If destination == -1 this process skips this communication step
        {
            Communicator::NodesContainerType& r_local_nodes = r_comm.LocalMesh(step).Nodes();
            Communicator::NodesContainerType& r_ghost_nodes = r_comm.GhostMesh(step).Nodes();

            //std::cout << r_comm.MyPID() << " (talking to " << destination << ") " << r_local_nodes.size() << " local nodes, " << r_ghost_nodes.size() << " ghost nodes." << std::endl;

            // Compute length of arrays to send and receive
            int node_send_size = 0;
            int data_send_size = 0;
            for (Communicator::NodesContainerType::iterator i_node = r_ghost_nodes.begin(); i_node != r_ghost_nodes.end(); i_node++)
                if (i_node->FastGetSolutionStepValue(NODAL_PAUX) > 0.0)
                {
                    node_send_size++;
                    const int node_index = static_cast<int>( i_node->FastGetSolutionStepValue(AUX_INDEX) );
                    data_send_size += 3*rActiveNeigh[ node_index ];
                }

            int node_recv_size = 0;

            for (Communicator::NodesContainerType::iterator i_node = r_local_nodes.begin(); i_node != r_local_nodes.end(); i_node++)
                if (i_node->FastGetSolutionStepValue(NODAL_PAUX) > 0.0)
                    node_recv_size++;

            // Allocate arrays to send/recv, fill send data
            // The last entry in SizeSendBuff is data_send_size (so the receiving process knows how many doubles to expect)
            std::vector<int> sizes_send_buffer(node_send_size+1);
            std::vector<double> data_send_buffer(data_send_size);

            sizes_send_buffer[node_send_size] = data_send_size;

            int ii = 0;
            int jj = 0;
            for (Communicator::NodesContainerType::iterator i_node = r_ghost_nodes.begin(); i_node != r_ghost_nodes.end(); i_node++)
                if (i_node->FastGetSolutionStepValue(NODAL_PAUX) > 0.0)
                {
                    int i = static_cast<int>( i_node->FastGetSolutionStepValue(AUX_INDEX) );
                    int offset = 3*MaxNeigh*i;
                    int n = rActiveNeigh[i];

                    sizes_send_buffer[ii++] = n;

                    for (int j = 0; j < 3*n; j++)
                        data_send_buffer[jj++] = rNormals[offset+j];
                }


            // Send/Receive normal data

            const DataCommunicator& r_data_comm = r_comm.GetDataCommunicator();

            std::vector<int> sizes_recv_buffer(node_recv_size+1);

            //std::cout << r_comm.MyPID() << " sending " << node_send_size+1 << " ints to " << destination << " (Expecting " << node_recv_size+1 << " ints in return)." << std::endl;

            r_data_comm.SendRecv(sizes_send_buffer, destination, step, sizes_recv_buffer, destination, step);

            int data_recv_size = sizes_recv_buffer[node_recv_size];
            std::vector<double> data_recv_buffer(data_recv_size);

            //std::cout << r_comm.MyPID() << " sending " << data_send_size << " doubles to " << destination << " (Expecting " << data_recv_size << " doubles in return)." << std::endl;

            r_data_comm.SendRecv(data_send_buffer, destination, step, data_recv_buffer, destination, step);

            //std::cout << r_comm.MyPID() << " finished communication with " << destination << std::endl;

            // Merge received data to local array

            ii = 0;
            jj = 0;
            for (Communicator::NodesContainerType::iterator i_node = r_local_nodes.begin(); i_node != r_local_nodes.end(); i_node++)
                if (i_node->FastGetSolutionStepValue(NODAL_PAUX) > 0.0)
                {
                    int i = static_cast<int>( i_node->FastGetSolutionStepValue(AUX_INDEX) );

                    // Append normal data
                    int num_recv_normals = sizes_recv_buffer[ii++];
                    int offset = 3*( MaxNeigh*i + rActiveNeigh[i] );

                    for (int j = 0; j < 3*num_recv_normals; j++)
                        rNormals[offset+j] = data_recv_buffer[jj++];


                    // Add number of extra normals
                    rActiveNeigh[i]  += num_recv_normals;
                }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void MPINormalCalculationUtils::DetectEdges(ModelPart &rModelPart,
                                            unsigned int Dimension,
                                            double MaxAngle,
                                            const std::vector<double>& rNormals,
                                            const std::vector<int>& rActiveNeigh,
                                            const int NumNodes,
                                            const int MaxNeighs)
{
    constexpr double to_radians = 3.14159265358979323846/180.;
    const double max_cosine = std::cos(to_radians*MaxAngle); // If the cosine between two normals is larger than this, they are considered to belong to different surfaces
    const double node_factor = 1.0/Dimension; // Total area is distributed between nodes with this ratio

    std::vector< array_1d<double,3> > active_normals;
    active_normals.reserve(MaxNeighs);

    for (ModelPart::NodeIterator i_node = rModelPart.GetCommunicator().LocalMesh().NodesBegin();
         i_node != rModelPart.GetCommunicator().LocalMesh().NodesEnd(); i_node++)
    {
        if ( i_node->FastGetSolutionStepValue(NODAL_PAUX) > 0.0)
        {
            int AuxID = static_cast<int>( i_node->FastGetSolutionStepValue(AUX_INDEX));
            int NumNormals = rActiveNeigh[AuxID];
            int Offset = 3*MaxNeighs*AuxID;

            active_normals.resize(0);
            array_1d<double,3> aux_normal(3,0.0);

            for (int n = 0; n < NumNormals; n++)
            {
                aux_normal[0] = rNormals[Offset + 3*n];
                aux_normal[1] = rNormals[Offset + 3*n + 1];
                aux_normal[2] = rNormals[Offset + 3*n + 2];

                if (active_normals.size() > 0)
                {
                    double norm = std::sqrt( aux_normal[0]*aux_normal[0] + aux_normal[1]*aux_normal[1] + aux_normal[2]*aux_normal[2] );

                    bool added = false;

                    for (unsigned int m = 0; m < active_normals.size(); m++)
                    {
                        array_1d<double,3>& this_normal = active_normals[m];
                        double norm_i = std::sqrt( this_normal[0]*this_normal[0] + this_normal[1]*this_normal[1] + this_normal[2]*this_normal[2] );

                        double cosine = aux_normal[0]*this_normal[0] + aux_normal[1]*this_normal[1] + aux_normal[2]*this_normal[2];
                        cosine /= norm*norm_i;

                        if (cosine > max_cosine) // We consider that the normals belong to the same surface
                        {
                            this_normal += aux_normal * node_factor;
                            added = true; // I think at this point I should add a break too, but the non-mpi version does exactly this (JC).
                        }
                    }

                    if (!added)
                        active_normals.push_back(aux_normal*node_factor);

                }
                else
                {
                    active_normals.push_back(aux_normal*node_factor);
                }

            }

            //
            // Set flags to detect edges and corners
            //

            switch ( active_normals.size() )
            {
            case 0:

                i_node->FastGetSolutionStepValue(IS_SLIP) = 0.0;
                break;
            case 1:
                i_node->FastGetSolutionStepValue(IS_SLIP) = 10.0;
                break;
            case 2:
                i_node->FastGetSolutionStepValue(IS_SLIP) = 20.0;
                break;
            case 3:
                i_node->FastGetSolutionStepValue(IS_SLIP) = 30.0;
                break;
            default:
                i_node->FastGetSolutionStepValue(IS_SLIP) = 30.0;
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

bool MPINormalCalculationUtils::OrientElement(Geometry<Node > &rGeom)
{
    const unsigned int PointIndex = 0;
    const GeometryData::IntegrationMethod Method = GeometryData::IntegrationMethod::GI_GAUSS_1;

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

void MPINormalCalculationUtils::NormalContribution(Geometry<Node > &rGeom)
{
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int Dim = rGeom.WorkingSpaceDimension();

    const unsigned int PointIndex = 0;
    const GeometryData::IntegrationMethod Method = GeometryData::IntegrationMethod::GI_GAUSS_1;
    double DetJ = rGeom.DeterminantOfJacobian(PointIndex,Method);

    Geometry< Node >::ShapeFunctionsGradientsType DN_DX;
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
                                             Geometry<Node > &rGeometry)
{
    An[0] =   rGeometry[1].Y() - rGeometry[0].Y();
    An[1] = - (rGeometry[1].X() - rGeometry[0].X());
    An[2] =    0.00;

}

///////////////////////////////////////////////////////////////////////////////

void MPINormalCalculationUtils::FaceNormal3D(array_1d<double,3> &An,
                                             Geometry<Node > &rGeometry)
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
            Geometry< Node >& rGeom = iCond->GetGeometry();
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

            for (Geometry< Node >::iterator iNode = rGeom.begin(); iNode != rGeom.end(); iNode++)
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
