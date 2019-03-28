#ifndef EIGEN_QR_UTILITY_HPP
#define EIGEN_QR_UTILITY_HPP

// System includes
#include <vector>

// External includes
#include <Eigen/QR>

// Project includes
#include "includes/define.h"
#include "custom_utilities/ublas_wrapper.h"
    

/*void Eigentest()
{
    std::cout<<"Eigen Utility called"<<std::endl;
}*/

namespace Kratos
{
class EigenQrUtility
{

public:


KRATOS_CLASS_POINTER_DEFINITION(EigenQrUtility);

EigenQrUtility(unsigned int system_size, const std::size_t n_sampling_points)
{
    std::cout<<"Eigen constructor called\n";

    std::cout<<"row size "<<system_size<<std::endl;
    std::cout<<"column size "<<n_sampling_points<<std::endl;
    Eigen::MatrixXd B(system_size, n_sampling_points);
     std::cout<<B.size()<<std::endl;

     //Eigen::Map<>
}

~EigenQrUtility(){}

private:

protected:

};

} // end class 
 // end namespace
#endif