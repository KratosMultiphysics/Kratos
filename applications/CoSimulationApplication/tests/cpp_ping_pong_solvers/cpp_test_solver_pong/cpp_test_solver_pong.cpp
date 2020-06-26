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
    EMPIRE_API_Connect("pong.xml");

    std::cout<<"Solver 2 : "<< "------------ Defining the mesh"<<std::endl;
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
    std::cout<<"Solver PONG : "<< "------------ Setting the mesh"<<std::endl;
    EMPIRE_API_sendMesh("mesh_pong", numNodes, numElems, nodes, nodeIDs, numNodesPerElem, elems);    

    std::cout<<"Solver PONG : "<<"Receiving ... \n";
    double toReceive = -1;
    EMPIRE_API_recvDataField("pong_recv_data", 1, &toReceive);
    std::cout<<"Solver PONG : "<<"Received: "<<toReceive<<std::endl;

    double toSend = 222;
    std::cout<<"Solver PONG : "<<"Sending ... \n";
    EMPIRE_API_sendDataField("pong_send_data", 1, &toSend);
    std::cout<<"Solver PONG : "<<"Sent: "<<toSend<<std::endl;

    return (0);
}
