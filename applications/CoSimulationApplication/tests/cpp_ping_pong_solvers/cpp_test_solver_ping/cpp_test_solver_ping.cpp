// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//					 license: CoSimulationApplication/license.txt
//
//  Main authors:    Aditya Ghantasala
//


//NOTE: This can use any implementation of IO which communicates with IO of CoSimulationApplication's PingPongIO.py
//      Here the header-only version of co_sim_EMPIRE_API.h is used as it is already available.
#include "custom_io/co_sim_EMPIRE_API.h" 
#include <assert.h>
#include <stdio.h>

int main(int argc, char **argv) {
    EMPIRE_API_Connect("ping.xml");

    std::cout<<"Solver 1 : "<<"------------ Defining the mesh"<<std::endl;
    int numNodes = 4;
    int numElems = 1;
    double nodes[] = {
                        0,0,0,
                        1,0,0,
                        1,1,0,
                        0,1,0
                        };
    int nodeIDs[] = {1,2,3,4};
    int elems[] = {1,2,3,4};
    int numNodesPerElem[] = {4};
    std::cout<<"Solver 1 : "<< "------------ Setting the mesh"<<std::endl;
    EMPIRE_API_sendMesh("mesh_ping", numNodes, numElems, nodes, nodeIDs, numNodesPerElem, elems);

    double toSend = 111;
    std::cout<<"Solver PING : "<<"Sending ... \n";
    EMPIRE_API_sendDataField("ping_send_data", 1, &toSend);
    std::cout<<"Solver PING : "<<"Sent: "<< toSend<<std::endl;

    double toReceive = -1;
    std::cout<<"Solver PING : "<< "Receiving ... \n";
    EMPIRE_API_recvDataField("ping_recv_data", 1, &toReceive);
    std::cout<<"Solver PING : "<<"Received: "<<toReceive<<std::endl;

    return (0);
}

