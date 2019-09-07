#include "custom_io/co_sim_EMPIRE_API.h"
#include <assert.h>
#include <stdio.h>

int main(int argc, char **argv) {
    EMPIRE_API_Connect("ping.xml");

    // Setting mesh in solver 1
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
    EMPIRE_API_sendMesh("mesh1_cpp", numNodes, numElems, nodes, nodeIDs, numNodesPerElem, elems);

    double toSend = 111;
    std::cout<<"Solver 1 : "<<"Sending ... \n";
    EMPIRE_API_sendDataField("ping_send_data", 1, &toSend);
    std::cout<<"Solver 1 : "<<"Sent: "<< toSend<<std::endl;

    double toReceive = -1;
    std::cout<<"Solver 1 : "<< "Receiving ... \n";
    EMPIRE_API_recvDataField("ping_recv_data", 1, &toReceive);
    std::cout<<"Solver 1 : "<<"Received: "<<toReceive<<std::endl;

    return (0);
}

