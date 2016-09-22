/*
 * logging_colors.hpp
 *
 *  Created on: Aug 4, 2016
 *      Author: mohamedkhalil
 */

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
#define LOG_GENERAL(color, message, variable)\
    std::cout << color << message << variable;            \
    LOG_RESET("")

#define LOG_RESET(suffix)\
    std::cout << RESET << suffix << std::endl; 

#define LOG_DEBUG(variable)\
    std::cout << BOLD << YELLOW << #variable << " : " << LT_YELLOW << variable;            \
    LOG_RESET("")
                                                        
#define LOG_INFO(variable)\
    std::cout << BOLD << GREEN  << #variable << " : " << LT_GREEN  << variable;            \
    LOG_RESET("")

#define LOG_WARNING(variable)\
    std::cout << BOLD << RED    << #variable << " : " << LT_RED    << variable;            \
    LOG_RESET("")
                                                        
#define LOG_DEBUG_HEADER(header, variable)\
    std::cout << BOLD << UNDERLINE << YELLOW << #header;                                   \
    LOG_RESET("")                                                                          \
    std::cout << BOLD <<              YELLOW << #variable << " : " << LT_YELLOW << variable\
    LOG_RESET("")

#define LOG_INFO_HEADER(header, variable)\
    std::cout << BOLD << UNDERLINE << GREEN  << #header;                                   \
    LOG_RESET("")                                                                          \
    std::cout << BOLD <<              GREEN  << #variable << " : " <<  LT_GREEN << variable\
    LOG_RESET("")
                                             
#define LOG_WARNING_HEADER(header, variable)\
    std::cout << BOLD << UNDERLINE << RED    << #header;                                   \
    LOG_RESET("")                                                                          \
    std::cout << BOLD <<              RED    << #variable << " : " <<  LT_RED   << variable\
    LOG_RESET("")
                                                        
#endif /* LOGGING_SETTINGS_HPP_ */


