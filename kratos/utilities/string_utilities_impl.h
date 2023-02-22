//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

#pragma once

// --- Core Includes ---
#include "utilities/string_utilities.h"


namespace Kratos {


template <class TOutputIterator>
void PlaceholderPattern::Glob(TOutputIterator it) const
{
    KRATOS_TRY

    // Create a path from the regex string omitting the leading '^' and trailing '$\0'
    PathType pattern(this->GetRegexString().substr(1, this->GetRegexString().size()-2));

    auto it_pattern_part = pattern.begin();
    std::vector<PathType> paths;

    // Decide where to begin globbing
    if (pattern.is_absolute()) { // the pattern is absolute => begin globbing at the filesystem root
        paths.emplace_back(pattern.root_path());
        ++it_pattern_part;
    } else { // the pattern is relative => begin globbing at the current working directory
        paths.emplace_back(std::filesystem::current_path());
    }

    // Compare the pattern parts to the globbed files'/directories' parts
    for ( ; it_pattern_part!=pattern.end(); ++it_pattern_part) {
        if (paths.empty()) break;

        std::regex pattern_part_regex(it_pattern_part->string());
        std::vector<PathType> tmp_paths;

        for (const auto& r_path : paths) {
            if (std::filesystem::is_directory(r_path)) {
                for (auto item : std::filesystem::directory_iterator(r_path)) {
                    if (std::regex_match(item.path().filename().string(), pattern_part_regex)) {
                        tmp_paths.emplace_back(item);
                    }
                } // for r_item in r_path.glob(*)
            } // if r_path.is_directory()
        } // for r_path in paths

        paths = tmp_paths;
    } // for pattern_part in pattern.parts

    // Output
    for (auto& r_path : paths) {
        *it++ = std::move(r_path);
    }

    KRATOS_CATCH("");
}


} // namespace Kratos
