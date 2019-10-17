// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//                   license: CoSimulationApplication/license.txt
//
//  Main authors:    Philipp Bucher
//

#ifndef KRATOS_CO_SIM_API_INTERNALS_H_INCLUDED
#define KRATOS_CO_SIM_API_INTERNALS_H_INCLUDED

// System includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>

// Project includes
#include "co_sim_io_define.h"

namespace CoSim {
namespace Internals {

#define CS_LOG std::cout << "[CoSimIO] "
#define CS_LOG_IF(condition) if(condition) CS_LOG

typedef std::unordered_map<std::string, std::string> SettingsType;

inline void AddMissingSettings(const SettingsType& rDefaultSettings, SettingsType& rSettings)
{
    for (const auto& r_setting : rDefaultSettings) {
        if (rSettings.count(r_setting.first) == 0) {
            rSettings[r_setting.first] = r_setting.second;
        }
    }
}

inline SettingsType ReadSettingsFile(const std::string& rSettingsFileName)
{
    std::ifstream settings_file(rSettingsFileName);

    if (!settings_file.good()) {
        std::cout << "Input file \"" << rSettingsFileName << "\" could not be read, using default configuration" << std::endl;
        return SettingsType();
    }

    SettingsType settings;
    std::string current_line;
    while (std::getline(settings_file, current_line)) {
        // TODO implement this
    }
    return settings;
}

} // namespace Internals
} // namespace CoSim

#endif /* KRATOS_CO_SIM_API_INTERNALS_H_INCLUDED */
