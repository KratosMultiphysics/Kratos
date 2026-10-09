//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <format>
#include <string>
#include <string_view>

// External includes

// Project includes

namespace Kratos::ColorUtilities {

/// ANSI escape code constants (empty on Windows)
#if !defined(_WIN32)
    inline constexpr std::string_view kReset         = "\x1B[0m";
    inline constexpr std::string_view kRed           = "\x1B[31m";
    inline constexpr std::string_view kGreen         = "\x1B[32m";
    inline constexpr std::string_view kYellow        = "\x1B[33m";
    inline constexpr std::string_view kBlue          = "\x1B[34m";
    inline constexpr std::string_view kMagenta       = "\x1B[35m";
    inline constexpr std::string_view kCyan          = "\x1B[36m";
    inline constexpr std::string_view kWhite         = "\x1B[37m";
    inline constexpr std::string_view kBold          = "\x1B[1m";
    inline constexpr std::string_view kFaint         = "\x1B[2m";
    inline constexpr std::string_view kItalic        = "\x1B[3m";
    inline constexpr std::string_view kUnderline     = "\x1B[4m";
    inline constexpr std::string_view kInverse       = "\x1B[7m";
    inline constexpr std::string_view kStrikethrough = "\x1B[9m";
#else
    inline constexpr std::string_view kReset         = "";
    inline constexpr std::string_view kRed           = "";
    inline constexpr std::string_view kGreen         = "";
    inline constexpr std::string_view kYellow        = "";
    inline constexpr std::string_view kBlue          = "";
    inline constexpr std::string_view kMagenta       = "";
    inline constexpr std::string_view kCyan          = "";
    inline constexpr std::string_view kWhite         = "";
    inline constexpr std::string_view kBold          = "";
    inline constexpr std::string_view kFaint         = "";
    inline constexpr std::string_view kItalic        = "";
    inline constexpr std::string_view kUnderline     = "";
    inline constexpr std::string_view kInverse       = "";
    inline constexpr std::string_view kStrikethrough = "";
#endif

/// Foreground color formatting
[[nodiscard]] inline std::string Red(std::string_view text)     { return std::format("{}{}{}", kRed, text, kReset); }
[[nodiscard]] inline std::string Green(std::string_view text)   { return std::format("{}{}{}", kGreen, text, kReset); }
[[nodiscard]] inline std::string Yellow(std::string_view text)  { return std::format("{}{}{}", kYellow, text, kReset); }
[[nodiscard]] inline std::string Blue(std::string_view text)    { return std::format("{}{}{}", kBlue, text, kReset); }
[[nodiscard]] inline std::string Magenta(std::string_view text) { return std::format("{}{}{}", kMagenta, text, kReset); }
[[nodiscard]] inline std::string Cyan(std::string_view text)    { return std::format("{}{}{}", kCyan, text, kReset); }
[[nodiscard]] inline std::string White(std::string_view text)   { return std::format("{}{}{}", kWhite, text, kReset); }

/// Font style formatting
[[nodiscard]] inline std::string Bold(std::string_view text)          { return std::format("{}{}{}", kBold, text, kReset); }
[[nodiscard]] inline std::string Faint(std::string_view text)         { return std::format("{}{}{}", kFaint, text, kReset); }
[[nodiscard]] inline std::string Italic(std::string_view text)        { return std::format("{}{}{}", kItalic, text, kReset); }
[[nodiscard]] inline std::string Underline(std::string_view text)     { return std::format("{}{}{}", kUnderline, text, kReset); }
[[nodiscard]] inline std::string Inverse(std::string_view text)       { return std::format("{}{}{}", kInverse, text, kReset); }
[[nodiscard]] inline std::string Strikethrough(std::string_view text) { return std::format("{}{}{}", kStrikethrough, text, kReset); }

} // namespace Kratos::ColorUtilities

// Backward-compatible macros — prefer using Kratos::ColorUtilities directly
#define RST   Kratos::ColorUtilities::kReset
#define KRED  Kratos::ColorUtilities::kRed
#define KGRN  Kratos::ColorUtilities::kGreen
#define KYEL  Kratos::ColorUtilities::kYellow
#define KBLU  Kratos::ColorUtilities::kBlue
#define KMAG  Kratos::ColorUtilities::kMagenta
#define KCYN  Kratos::ColorUtilities::kCyan
#define KWHT  Kratos::ColorUtilities::kWhite

#define FRED(x)     Kratos::ColorUtilities::Red(x)
#define FGRN(x)     Kratos::ColorUtilities::Green(x)
#define FYEL(x)     Kratos::ColorUtilities::Yellow(x)
#define FBLU(x)     Kratos::ColorUtilities::Blue(x)
#define FMAG(x)     Kratos::ColorUtilities::Magenta(x)
#define FCYN(x)     Kratos::ColorUtilities::Cyan(x)
#define FWHT(x)     Kratos::ColorUtilities::White(x)

#define BOLDFONT(x)  Kratos::ColorUtilities::Bold(x)
#define FAINTFONT(x) Kratos::ColorUtilities::Faint(x)
#define ITAFONT(x)   Kratos::ColorUtilities::Italic(x)
#define UNDL(x)      Kratos::ColorUtilities::Underline(x)
#define INVFONT(x)   Kratos::ColorUtilities::Inverse(x)
#define STRFONT(x)   Kratos::ColorUtilities::Strikethrough(x)