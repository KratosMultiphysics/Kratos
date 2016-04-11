#include <iostream>
#include "eig3.h"

void post(double A[3][3],double V[3][3], double d[3])
{
    std::cout << "Eigenvalues:";
    for(int i = 0; i < 3; ++i)
    {
        std::cout << " " << d[i];
    }
    std::cout << std::endl;
    
    std::cout << "Eigenvectors:" << std::endl;
    for(int i = 0; i < 3; ++i)
    {
        std::cout << "Eigenvector " << i << ":";
        for(int j = 0; j < 3; ++j)
        {
            std::cout << " " << V[j][i];
        }
        std::cout << std::endl;
    }
    std::cout << "----------------------------------------------------" << std::endl;
}

int main(void)
{
    using namespace Kratos;

    double A[3][3];
    double V[3][3];
    double d[3];
    
    A[0][0] = 1.0;
    A[0][1] = 0.0;
    A[0][2] = 0.0;
    
    A[1][0] = 0.0;
    A[1][1] = 1.0;
    A[1][2] = 0.0;
    
    A[2][0] = 0.0;
    A[2][1] = 0.0;
    A[2][2] = 1.0;
    
    eigen_decomposition(A, V, d);
    post(A, V, d);
    
    A[0][0] = 1.0000;
    A[0][1] = 0.9794;
    A[0][2] = 0.6967;
    
    A[1][0] = 0.9794;
    A[1][1] = 1.0000;
    A[1][2] = 0.7742;
    
    A[2][0] = 0.6967;
    A[2][1] = 0.7742;
    A[2][2] = 1.0000;
    
    eigen_decomposition(A, V, d);
    post(A, V, d);
    
    A[0][0] = 178;
    A[0][1] = 206;
    A[0][2] = 234;
    
    A[1][0] = 206;
    A[1][1] = 238;
    A[1][2] = 270;
    
    A[2][0] = 234;
    A[2][1] = 270;
    A[2][2] = 306;
    
    eigen_decomposition(A, V, d);
    post(A, V, d);
    
    return 0;
}
