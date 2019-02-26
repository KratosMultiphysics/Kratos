// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//

#ifndef OUTPUT_UTILITY_H
#define OUTPUT_UTILITY_H

// System includes

// External includes

// Project includes
#include "includes/define.h"


namespace Kratos
{

/** \brief OutputUtility
*
* 
* 
* matching the data type of the response variable (array_1d<double, 3>, Matrix etc.) 
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) OutputUtility
{
public:

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;


    static void OutputOnTerminal(const std::string output_name, const std::vector<std::vector<array_1d<double, 3>>>& output_vector)
    {
        std::cout << output_name << ":  " << std::endl;
        //std::cout << std::setw(6);
        for(IndexType i = 0; i < output_vector.size(); ++i)
        {
            std::cout << "       ";
            for(IndexType j = 0; j < output_vector[0].size(); ++j)
            {
                std::cout << "    |";
                for(IndexType dir_it = 0; dir_it < 3; ++dir_it)
                {
                    if(dir_it == 2)
                        std::cout << std::setw(8) << output_vector[i][j][dir_it] << "|    ";
                    else
                        std::cout << std::setw(8) << output_vector[i][j][dir_it] << "  :  ";
                }             
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;          
    }    
    
    static void OutputOnTerminal(const std::string output_name, const std::vector<array_1d<double, 3>>& output_vector)
    {
        std::cout << output_name << ":  " << std::endl;
        for(IndexType i = 0; i < output_vector.size(); ++i)
        {
            std::cout << "           |";
            for(IndexType dir_it = 0; dir_it < 3; ++dir_it)
            {
                if(dir_it == 2)
                    std::cout << std::setw(8) << output_vector[i][dir_it] << "|";
                else
                    std::cout << std::setw(8) << output_vector[i][dir_it] << "  :  ";
            } 
            std::cout << std::endl;
        } 
        std::cout << std::endl;         
    }

    static void OutputOnTerminal(const std::string output_name, const Vector& output_vector)
    {
        std::cout << std::setw(18) << output_name << ":   |";
        for(IndexType i = 0; i < output_vector.size(); ++i)
        {
            if (i == output_vector.size()-1)
                std::cout << std::setw(8) << output_vector[i] << "|";
            else
                std::cout << std::setw(8) << output_vector[i] << " : ";
        } 
        std::cout << std::endl;
        std::cout << std::endl;
    }
    

    static void OutputOnTerminal(const std::string output_name, const std::vector<Matrix>& output_vector)
    {
        std::cout << output_name << std::endl; 
        for(IndexType k = 0; k < output_vector.size(); ++k)
        {   
            for(IndexType i = 0; i < output_vector[k].size2(); ++i)
            {
                std::cout << "    |";
                for(IndexType j = 0; j < output_vector[k].size1(); ++j)
                {   
                    if (j == output_vector[k].size1()-1)
                        std::cout << std::setw(8) << output_vector[k](j, i) << "|";
                    else
                        std::cout << std::setw(8) << output_vector[k](j, i) << " : ";
                }
                std::cout << std::endl; 
            }                
            std::cout << "    |                                  |" << std::endl;
            std::cout << "    |                                  |" <<std::endl;
        }
    } 
    
    static void OutputOnTerminal(const std::string output_name, const Matrix& output_vector)
    {
        std::cout << output_name << std::endl; 
        for(IndexType i = 0; i < output_vector.size2(); ++i)
        {
            std::cout << "    |";
            for(IndexType j = 0; j < output_vector.size1(); ++j)
            {   
                if (j == output_vector.size1()-1)
                    std::cout << std::setw(8) << output_vector(j, i) << "|";
                else
                    std::cout << std::setw(8) << output_vector(j, i) << " : ";
            }
            std::cout << std::endl; 
        }
        std::cout << std::endl;
        std::cout << std::endl;
    }


private:

    
    
}; // class OutputUtil


} /* namespace Kratos.*/

#endif /* OUTPUT_UTILITY_H defined */