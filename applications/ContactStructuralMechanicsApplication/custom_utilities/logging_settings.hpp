// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Mohamed Khalil
// 
#ifndef LOGGING_SETTINGS_HPP_
#define LOGGING_SETTINGS_HPP_

// only tested on mac console
// for linux, prefix \033 should be changed to \e
// colors
#define RESET      "\033[0m"
#define BLACK      "\033[30m"
#define RED        "\033[31m"
#define GREEN      "\033[32m"
#define YELLOW     "\033[33m"
#define BLUE       "\033[34m"
#define MAGENTA    "\033[35m"
#define CYAN       "\033[36m"
#define WHITE      "\033[37m"
#define DK_GREY    "\033[90m"
#define LT_GREY    "\033[37m"
#define LT_RED     "\033[91m"
#define LT_GREEN   "\033[92m"
#define LT_YELLOW  "\033[93m"
#define LT_BLUE    "\033[94m"
#define LT_MAGENTA "\033[95m"
#define LT_CYAN    "\033[96m"

// formatting
#define BOLD       "\033[1m"
#define DIM        "\033[2m"
#define UNDERLINE  "\033[4m"

// output level
enum LoggingLevel { none, info_warnings, system_matrices, debug };

// Logging macros
#define RESET_LOG_SETTINGS \
    std::cout << RESET << std::endl; 

#define LOG_GENERAL(color, message, variable)\
    std::cout << color << message << variable;            \
    RESET_LOG_SETTINGS

#define LOG_DEBUG(variable)\
    std::cout << BOLD << YELLOW << #variable << " : " << LT_YELLOW << variable;            \
    RESET_LOG_SETTINGS
                                                        
#define LOG_INFO(variable)\
    std::cout << BOLD << GREEN  << #variable << " : " << LT_GREEN  << variable;            \
    RESET_LOG_SETTINGS

#define LOG_WARNING(variable)\
    std::cout << BOLD << RED    << #variable << " : " << LT_RED    << variable;            \
    RESET_LOG_SETTINGS
                                                        
#define LOG_DEBUG_HEADER(header, variable)\
    std::cout << BOLD << UNDERLINE << YELLOW << #header;                                   \
    RESET_LOG_SETTINGS                                                                          \
    std::cout << BOLD <<              YELLOW << #variable << " : " << LT_YELLOW << variable\
    RESET_LOG_SETTINGS

#define LOG_INFO_HEADER(header, variable)\
    std::cout << BOLD << UNDERLINE << GREEN  << #header;                                   \
    RESET_LOG_SETTINGS                                                                          \
    std::cout << BOLD <<              GREEN  << #variable << " : " <<  LT_GREEN << variable\
    RESET_LOG_SETTINGS
                                             
#define LOG_WARNING_HEADER(header, variable)\
    std::cout << BOLD << UNDERLINE << RED    << #header;                                   \
    RESET_LOG_SETTINGS                                                                          \
    std::cout << BOLD <<              RED    << #variable << " : " <<  LT_RED   << variable\
    RESET_LOG_SETTINGS
                     
////////////////////////////////////////////////////////////////////////////////////
//////////////////////////    MESSAGES OUTPUT SETTINGS    //////////////////////////
////////////////////////////////////////////////////////////////////////////////////
    
#define DEBUG_MSG( message )\
    std::cout << BOLD << GREEN << UNDERLINE << "DEBUG" << RESET << std::endl;\
    std::cout << BOLD << GREEN << message << std::endl;\
    RESET_LOG_SETTINGS

#define INFO_MSG( message )\
    std::cout << BOLD << DK_GREY << UNDERLINE << "INFO : " << RESET << std::endl;\
    std::cout << BOLD << DK_GREY << message << std::endl;\
    RESET_LOG_SETTINGS

#define WARNING_MSG( message )\
    std::cout << BOLD << YELLOW << UNDERLINE << "WARNING" << RESET << std::endl;\
    std::cout << BOLD << YELLOW << message << std::endl;\
    RESET_LOG_SETTINGS
    
#define ERROR_MSG( message, function )\
    std::cout << BOLD << RED << UNDERLINE << "ERROR from " << function << RESET << std::endl;\
    std::cout << BOLD << RED << message << std::endl;\
    RESET_LOG_SETTINGS
    

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////    TENSOR OUTPUT SETTINGS    //////////////////////////
//////////////////////////////////////////////////////////////////////////////////

/*
 * DEFINE LOG SETTINGS
 */
#define TENSOR_LOG_SETTINGS \
    std::cout << CYAN;\
    std::cout.precision(2);\

/*
 * DEFINE FUNCTIONALITIES
 */
#define LOG_MATRIX_PRETTY( matrix )\
    TENSOR_LOG_SETTINGS \
    std::cout << #matrix << " : [ " << matrix.size1() << " x " << matrix.size2() << " ] :" << std::endl;\
    for ( unsigned int i = 0; i < matrix.size1( ); ++i )\
    {\
        for ( unsigned int j = 0; j < matrix.size2( ); ++j )\
            std::cout << std::scientific << matrix(i,j) << "\t";\
        std::cout << std::endl;\
    }\
    RESET_LOG_SETTINGS
    
#define LOG_VECTOR_PRETTY( vector )\
    TENSOR_LOG_SETTINGS \
    std::cout << #vector << " : " << std::scientific << vector;\
    RESET_LOG_SETTINGS

#define LOG_VECTOR3( array )\
    TENSOR_LOG_SETTINGS \
    std::cout << #array << " : " << std::scientific << array[0] << ", " << std::scientific << array[1] << ", " << std::scientific << array[2];\
    RESET_LOG_SETTINGS
    
#define LOG_VECTOR2( array )\
    TENSOR_LOG_SETTINGS \
    std::cout << #array << " : " << std::scientific << array[0] << ", " << std::scientific << array[1];\
    RESET_LOG_SETTINGS

#define LOG_SCALAR( scalar )\
    TENSOR_LOG_SETTINGS \
    std::cout << #scalar << " : " << std::scientific << scalar;\
    RESET_LOG_SETTINGS


/////////////////////////////////////////////////////////////////////////////////////
//////////////////////////    CONDITION OUTPUT SETTINGS    //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

/*
 * DEFINE LOG SETTINGS
 */
#define CONDITION_LOG_SETTINGS \
    std::cout << BLUE << "Condition : " << std::endl; \
    std::cout << LT_BLUE;\

/*
 * DEFINE FUNCTIONALITIES
 */
#define LOG_CONDITION_HEADER( master, slave ) \
    CONDITION_LOG_SETTINGS \
    std::cout << "|_ Master : "; \
    for( unsigned int i = 0; i < master.PointsNumber( ) - 1; ++i ) \
        std::cout << master[i].Id( ) << ", "; \
    std::cout << master[master.PointsNumber( ) - 1].Id( ) << "\n";\
    std::cout << "|_ Slave  : "; \
    for( unsigned int i = 0; i < master.PointsNumber( ) - 1; ++i ) \
        std::cout << slave[i].Id( ) << ", "; \
    std::cout << slave[slave.PointsNumber( ) - 1].Id( ) << "\n";\
    RESET_LOG_SETTINGS


#endif /* LOGGING_SETTINGS_HPP_ */


