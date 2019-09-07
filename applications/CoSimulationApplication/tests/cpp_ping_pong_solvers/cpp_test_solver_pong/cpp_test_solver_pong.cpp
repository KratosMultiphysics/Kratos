#include "custom_io/co_sim_EMPIRE_API.h"
#include <assert.h>
#include <stdio.h>

int main(int argc, char **argv) {
    EMPIRE_API_Connect("pong.xml");


    // Setting mesh in solver 2
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
    std::cout<<"Solver 2 : "<< "------------ Setting the mesh"<<std::endl;
    EMPIRE_API_sendMesh("mesh2_cpp", numNodes, numElems, nodes, nodeIDs, numNodesPerElem, elems);    

    std::cout<<"Solver 2 : "<<"Receiving ... \n";
    double toReceive = -1;
    EMPIRE_API_recvDataField("pong_recv_data", 1, &toReceive);
    std::cout<<"Solver 2 : "<<"Received: "<<toReceive<<std::endl;

    double toSend = 222;
    std::cout<<"Solver 2 : "<<"Sending ... \n";
    EMPIRE_API_sendDataField("pong_send_data", 1, &toSend);
    std::cout<<"Solver 2 : "<<"Sent: "<<toSend<<std::endl;

    return (0);
}
