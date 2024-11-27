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
//                   The implementations come from CppBenchmark, with MIT license, see https://github.com/chronoxor/CppBenchmark/
//

// System includes
#include <sstream>
#include <codecvt>
#include <cstring>
#include <locale>
#include <regex>
#include <fstream>
#include <set>
// #include <thread>

// External includes

// Project includes
#include "input_output/logger.h"
#include "includes/system_information.h"

// OS specific includes
#if defined(KRATOS_COMPILED_IN_OS)
#include <mach/mach.h>
#include <mach/mach_time.h>
#include <sys/sysctl.h>
#include <math.h>
#include <pthread.h>
#elif defined(KRATOS_COMPILED_IN_LINUX)
#include <sys/stat.h>
#include <sys/utsname.h>
#include <pthread.h>
#include <unistd.h>
#elif defined(KRATOS_COMPILED_IN_WINDOWS)
#include <windows.h>
#include <winternl.h>
#include <memory>
#define STATUS_SUCCESS 0x00000000
#endif

namespace Kratos
{
std::string SystemInformation::OSVersion()
{
#if defined(KRATOS_COMPILED_IN_OS)
    char result[1024];
    std::size_t size = sizeof(result);
    if (sysctlbyname("kern.osrelease", result, &size, nullptr, 0) == 0) {
        return result;
    }

    return "<apple>";
#elif defined(__CYGWIN__)
    struct utsname name;
    if (uname(&name) == 0) {
        std::string result(name.sysname);
        result.append(" ");
        result.append(name.release);
        result.append(" ");
        result.append(name.version);
        return result;
    }

    return "<cygwin>";
#elif defined(linux) || defined(__linux) || defined(__linux__)
    static std::regex pattern("DISTRIB_DESCRIPTION=\"(.*)\"");

    std::string line;
    std::ifstream stream("/etc/lsb-release");
    while (getline(stream, line)) {
        std::smatch matches;
        if (std::regex_match(line, matches, pattern)) {
            return matches[1];
        }
    }

    return "<linux>";
#elif defined(KRATOS_COMPILED_IN_WINDOWS)
    static NTSTATUS(__stdcall *RtlGetVersion)(OUT PRTL_OSVERSIONINFOEXW lpVersionInformation) = (NTSTATUS(__stdcall*)(PRTL_OSVERSIONINFOEXW))GetProcAddress(GetModuleHandle("ntdll.dll"), "RtlGetVersion");
    static void(__stdcall *GetNativeSystemInfo)(OUT LPSYSTEM_INFO lpSystemInfo) = (void(__stdcall*)(LPSYSTEM_INFO))GetProcAddress(GetModuleHandle("kernel32.dll"), "GetNativeSystemInfo");
    static BOOL(__stdcall *GetProductInfo)(IN DWORD dwOSMajorVersion, IN DWORD dwOSMinorVersion, IN DWORD dwSpMajorVersion, IN DWORD dwSpMinorVersion, OUT PDWORD pdwReturnedProductType) = (BOOL(__stdcall*)(DWORD, DWORD, DWORD, DWORD, PDWORD))GetProcAddress(GetModuleHandle("kernel32.dll"), "GetProductInfo");

    OSVERSIONINFOEXW osvi;
    ZeroMemory(&osvi, sizeof(OSVERSIONINFOEXW));
    osvi.dwOSVersionInfoSize = sizeof(OSVERSIONINFOEXW);

    if (RtlGetVersion != nullptr) {
        NTSTATUS ntRtlGetVersionStatus = RtlGetVersion(&osvi);
        if (ntRtlGetVersionStatus != STATUS_SUCCESS) {
            return "<windows>";
        }
    } else {
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable: 4996) // C4996: 'function': was declared deprecated
#endif
        BOOL bOsVersionInfoEx = GetVersionExW((OSVERSIONINFOW*)&osvi);
        if (bOsVersionInfoEx == 0) {
            return "<windows>";
        }
#if defined(_MSC_VER)
#pragma warning(pop)
#endif
    }

    SYSTEM_INFO si;
    ZeroMemory(&si, sizeof(SYSTEM_INFO));

    if (GetNativeSystemInfo != nullptr) {
        GetNativeSystemInfo(&si);
    } else {
        GetSystemInfo(&si);
    }

    if ((osvi.dwPlatformId != VER_PLATFORM_WIN32_NT) || (osvi.dwMajorVersion <= 4)) {
        return "<windows>";
    }

    std::stringstream os;

    // Windows version
    os << "Microsoft ";
    if (osvi.dwMajorVersion >= 6) {
        if (osvi.dwMajorVersion == 10) {
            if (osvi.dwMinorVersion == 0) {
                if (osvi.wProductType == VER_NT_WORKSTATION) {
                    if (osvi.dwBuildNumber >= 22000) {
                        os << "Windows 11 ";
                    } else {
                        os << "Windows 10 ";
                    }
                } else {
                    if (osvi.dwBuildNumber >= 20348) {
                        os << "Windows Server 2022 ";
                    } else if (osvi.dwBuildNumber >= 17763) {
                        os << "Windows Server 2019 ";
                    } else {
                        os << "Windows Server 2016 ";
                    }
                }
            }
        } else if (osvi.dwMajorVersion == 6) {
            if (osvi.dwMinorVersion == 3) {
                if (osvi.wProductType == VER_NT_WORKSTATION) {
                    os << "Windows 8.1 ";
                } else {
                    os << "Windows Server 2012 R2 ";
                }
            } else if (osvi.dwMinorVersion == 2) {
                if (osvi.wProductType == VER_NT_WORKSTATION) {
                    os << "Windows 8 ";
                } else {
                    os << "Windows Server 2012 ";
                }
            } else if (osvi.dwMinorVersion == 1) {
                if (osvi.wProductType == VER_NT_WORKSTATION) {
                    os << "Windows 7 ";
                } else {
                    os << "Windows Server 2008 R2 ";
                }
            } else if (osvi.dwMinorVersion == 0) {
                if (osvi.wProductType == VER_NT_WORKSTATION) {
                    os << "Windows Vista ";
                } else {
                    os << "Windows Server 2008 ";
                }
            }
        }

        DWORD dwType;
        if ((GetProductInfo != nullptr) && GetProductInfo(osvi.dwMajorVersion, osvi.dwMinorVersion, 0, 0, &dwType)) {
            switch (dwType) {
                case PRODUCT_ULTIMATE:
                    os << "Ultimate Edition";
                    break;
                case PRODUCT_PROFESSIONAL:
                    os << "Professional";
                    break;
                case PRODUCT_HOME_PREMIUM:
                    os << "Home Premium Edition";
                    break;
                case PRODUCT_HOME_BASIC:
                    os << "Home Basic Edition";
                    break;
                case PRODUCT_ENTERPRISE:
                    os << "Enterprise Edition";
                    break;
                case PRODUCT_BUSINESS:
                    os << "Business Edition";
                    break;
                case PRODUCT_STARTER:
                    os << "Starter Edition";
                    break;
                case PRODUCT_CLUSTER_SERVER:
                    os << "Cluster Server Edition";
                    break;
                case PRODUCT_DATACENTER_SERVER:
                    os << "Datacenter Edition";
                    break;
                case PRODUCT_DATACENTER_SERVER_CORE:
                    os << "Datacenter Edition (core installation)";
                    break;
                case PRODUCT_ENTERPRISE_SERVER:
                    os << "Enterprise Edition";
                    break;
                case PRODUCT_ENTERPRISE_SERVER_CORE:
                    os << "Enterprise Edition (core installation)";
                    break;
                case PRODUCT_ENTERPRISE_SERVER_IA64:
                    os << "Enterprise Edition for Itanium-based Systems";
                    break;
                case PRODUCT_SMALLBUSINESS_SERVER:
                    os << "Small Business Server";
                    break;
                case PRODUCT_SMALLBUSINESS_SERVER_PREMIUM:
                    os << "Small Business Server Premium Edition";
                    break;
                case PRODUCT_STANDARD_SERVER:
                    os << "Standard Edition";
                    break;
                case PRODUCT_STANDARD_SERVER_CORE:
                    os << "Standard Edition (core installation)";
                    break;
                case PRODUCT_WEB_SERVER:
                    os << "Web Server Edition";
                    break;
            }
        }
    } else if ((osvi.dwMajorVersion == 5) && (osvi.dwMinorVersion == 2)) {
        if (GetSystemMetrics(SM_SERVERR2)) {
            os << "Windows Server 2003 R2, ";
        } else if (osvi.wSuiteMask & VER_SUITE_STORAGE_SERVER) {
            os << "Windows Storage Server 2003";
        } else if (osvi.wSuiteMask & VER_SUITE_WH_SERVER) {
            os << "Windows Home Server";
        } else if ((osvi.wProductType == VER_NT_WORKSTATION) && (si.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_AMD64)) {
            os << "Windows XP Professional x64 Edition";
        } else {
            os << "Windows Server 2003, ";
        }
        if (osvi.wProductType != VER_NT_WORKSTATION) {
            if (si.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_IA64) {
                if (osvi.wSuiteMask & VER_SUITE_DATACENTER) {
                    os << "Datacenter Edition for Itanium-based Systems";
                } else if (osvi.wSuiteMask & VER_SUITE_ENTERPRISE) {
                    os << "Enterprise Edition for Itanium-based Systems";
                }
            } else if (si.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_AMD64) {
                if (osvi.wSuiteMask & VER_SUITE_DATACENTER) {
                    os << "Datacenter x64 Edition";
                } else if (osvi.wSuiteMask & VER_SUITE_ENTERPRISE) {
                    os << "Enterprise x64 Edition";
                } else {
                    os << "Standard x64 Edition";
                }
            } else {
                if (osvi.wSuiteMask & VER_SUITE_COMPUTE_SERVER) {
                    os << "Compute Cluster Edition";
                } else if (osvi.wSuiteMask & VER_SUITE_DATACENTER) {
                    os << "Datacenter Edition";
                } else if (osvi.wSuiteMask & VER_SUITE_ENTERPRISE) {
                    os << "Enterprise Edition";
                } else if (osvi.wSuiteMask & VER_SUITE_BLADE) {
                    os << "Web Edition";
                } else {
                    os << "Standard Edition";
                }
            }
        }
    } else if ((osvi.dwMajorVersion == 5) && (osvi.dwMinorVersion == 1)) {
        os << "Windows XP ";
        if (osvi.wSuiteMask & VER_SUITE_PERSONAL) {
            os << "Home Edition";
        } else {
            os << "Professional";
        }
    } else if ((osvi.dwMajorVersion == 5) && (osvi.dwMinorVersion == 0)) {
        os << "Windows 2000 ";
        if (osvi.wProductType == VER_NT_WORKSTATION) {
            os << "Professional";
        } else {
            if (osvi.wSuiteMask & VER_SUITE_DATACENTER) {
                os << "Datacenter Server";
            } else if (osvi.wSuiteMask & VER_SUITE_ENTERPRISE) {
                os << "Advanced Server";
            } else {
                os << "Server";
            }
        }
    }

    // Windows Service Pack version
    std::wstring sp_version(osvi.szCSDVersion);
    std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> convert;
    if (std::wcslen(osvi.szCSDVersion) > 0) {
        os << " " << convert.to_bytes(sp_version.data(), sp_version.data() + sp_version.size());
    }

    // Windows build
    os << " (build " << osvi.dwBuildNumber << ")";

    // Windows architecture
    if (osvi.dwMajorVersion >= 6) {
        if (si.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_INTEL) {
            os << ", 32-bit";
        } else if (si.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_AMD64) {
            os << ", 64-bit";
        } else if (si.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_IA64) {
            os << ", Intel Itanium";
        } else if (si.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_ARM) {
            os << ", ARM";
#if !defined(__MINGW32__) && !defined(__MINGW64__)
        } else if (si.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_ARM64) {
            os << ", ARM64";
        }
#endif
    }

    return os.str();
#else
    #error Unsupported platform
#endif
}

/***********************************************************************************/
/***********************************************************************************/

std::string SystemInformation::CPUArchitecture()
{
#if defined(KRATOS_COMPILED_IN_OS)
    char result[1024];
    std::size_t size = sizeof(result);
    if (sysctlbyname("machdep.cpu.brand_string", result, &size, nullptr, 0) == 0) {
        return result;
    }

    return "<unknown>";
#elif defined(KRATOS_COMPILED_IN_LINUX)
    static std::regex pattern("model name(.*): (.*)");

    std::string line;
    std::ifstream stream("/proc/cpuinfo");
    while (getline(stream, line)) {
        std::smatch matches;
        if (std::regex_match(line, matches, pattern)) {
            return matches[2];
        }
    }

    return "<unknown>";
#elif defined(KRATOS_COMPILED_IN_WINDOWS)
    HKEY hKeyProcessor;
    LONG lError = RegOpenKeyExA(HKEY_LOCAL_MACHINE, "HARDWARE\\DESCRIPTION\\System\\CentralProcessor\\0", 0, KEY_READ, &hKeyProcessor);
    if (lError != ERROR_SUCCESS)
        return "<unknown>";

    // Smart resource cleaner pattern
    auto clearer = [](HKEY hKey) { RegCloseKey(hKey); };
    auto key = std::unique_ptr<std::remove_pointer<HKEY>::type, decltype(clearer)>(hKeyProcessor, clearer);

    CHAR pBuffer[_MAX_PATH] = { 0 };
    DWORD dwBufferSize = sizeof(pBuffer);
    lError = RegQueryValueExA(key.get(), "ProcessorNameString", nullptr, nullptr, (LPBYTE)pBuffer, &dwBufferSize);
    if (lError != ERROR_SUCCESS) {
        return "<unknown>";
    }

    return std::string(pBuffer);
#else
    #error Unsupported platform
#endif
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t SystemInformation::CPULogicalCores()
{
    const int cpu_logical_cores = CPUTotalCores().first;
    if (cpu_logical_cores > 0) {
        return static_cast<std::size_t>(cpu_logical_cores);
    } else {
        KRATOS_WARNING("SystemInformation") << "The number of logical cores is not available" << std::endl;
        return 0;
    }
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t SystemInformation::CPUPhysicalCores()
{
    const int cpu_physical_cores = CPUTotalCores().second;
    if (cpu_physical_cores > 0) {
        return static_cast<std::size_t>(cpu_physical_cores);
    } else {
        KRATOS_WARNING("SystemInformation") << "The number of physical cores is not available" << std::endl;
        return 0;
    }
}

/***********************************************************************************/
/***********************************************************************************/

#if defined(KRATOS_COMPILED_IN_WINDOWS)
namespace Internals
{
// Helper function to count set bits in the processor mask
DWORD CountSetBits(ULONG_PTR pBitMask)
{
    DWORD dwLeftShift = sizeof(ULONG_PTR) * 8 - 1;
    DWORD dwBitSetCount = 0;
    ULONG_PTR pBitTest = (ULONG_PTR)1 << dwLeftShift;

    for (DWORD i = 0; i <= dwLeftShift; ++i) {
        dwBitSetCount += ((pBitMask & pBitTest) ? 1 : 0);
        pBitTest /= 2;
    }

    return dwBitSetCount;
}
} // namespace Internals
#endif

std::pair<int, int> SystemInformation::CPUTotalCores()
{
    // This code is generic, cane be used to ant platform
    // int logical = std::thread::hardware_concurrency();
    
#if defined(KRATOS_COMPILED_IN_OS)
    int logical = 0;
    std::size_t logical_size = sizeof(logical);
    if (sysctlbyname("hw.logicalcpu", &logical, &logical_size, nullptr, 0) != 0) {
        logical = -1;
    }

    int physical = 0;
    std::size_t physical_size = sizeof(physical);
    if (sysctlbyname("hw.physicalcpu", &physical, &physical_size, nullptr, 0) != 0) {
        physical = -1;
    }

    return std::make_pair(logical, physical);
#elif defined(KRATOS_COMPILED_IN_LINUX)
    std::set<int> unique_cores;
    std::ifstream cpuinfo("/proc/cpuinfo");
    std::string line;
    while (std::getline(cpuinfo, line)) {
        if (line.find("core id") != std::string::npos) {
            const int core_id = std::stoi(line.substr(line.find(":") + 1));
            unique_cores.insert(core_id);
        }
    }
    const int physical = unique_cores.size();
    const long logical = sysconf(_SC_NPROCESSORS_ONLN);
    return std::make_pair(logical, physical);
#elif defined(KRATOS_COMPILED_IN_WINDOWS)
    BOOL allocated = FALSE;
    PSYSTEM_LOGICAL_PROCESSOR_INFORMATION pBuffer = nullptr;
    DWORD dwLength = 0;

    while (!allocated) {
        BOOL bResult = GetLogicalProcessorInformation(pBuffer, &dwLength);
        if (bResult == FALSE) {
            if (GetLastError() == ERROR_INSUFFICIENT_BUFFER) {
                if (pBuffer != nullptr) {
                    std::free(pBuffer);
                }
                pBuffer = (PSYSTEM_LOGICAL_PROCESSOR_INFORMATION)std::malloc(dwLength);
                if (pBuffer == nullptr) {
                    return std::make_pair(-1, -1);
                }
            } else {
                return std::make_pair(-1, -1);
            }
        } else {
            allocated = TRUE;
        }
    }

    std::pair<int, int> result(0, 0);
    PSYSTEM_LOGICAL_PROCESSOR_INFORMATION pCurrent = pBuffer;
    DWORD dwOffset = 0;

    while (dwOffset + sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION) <= dwLength) {
        switch (pCurrent->Relationship) {
            case RelationProcessorCore:
                result.first += Internals::CountSetBits(pCurrent->ProcessorMask);
                result.second += 1;
                break;
            case RelationNumaNode:
            case RelationCache:
            case RelationProcessorPackage:
                break;
            default:
                return std::make_pair(-1, -1);
        }
        dwOffset += sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION);
        pCurrent++;
    }

    std::free(pBuffer);

    return result;
#else
    #error Unsupported platform
#endif
}

/***********************************************************************************/
/***********************************************************************************/

int64_t SystemInformation::CPUClockSpeed()
{
#if defined(KRATOS_COMPILED_IN_OS)
    uint64_t frequency = 0;
    size_t size = sizeof(frequency);
    if (sysctlbyname("hw.cpufrequency", &frequency, &size, nullptr, 0) == 0)
        return frequency;

    return -1;
#elif defined(KRATOS_COMPILED_IN_LINUX)
    static std::regex pattern("cpu MHz(.*): (.*)");

    std::string line;
    std::ifstream stream("/proc/cpuinfo");
    while (getline(stream, line)) {
        std::smatch matches;
        if (std::regex_match(line, matches, pattern)) {
            return (int64_t)(atof(matches[2].str().c_str()) * 1000000);
        }
    }

    return -1;
#elif defined(KRATOS_COMPILED_IN_WINDOWS)
    HKEY hKeyProcessor;
    long lError = RegOpenKeyExA(HKEY_LOCAL_MACHINE, "HARDWARE\\DESCRIPTION\\System\\CentralProcessor\\0", 0, KEY_READ, &hKeyProcessor);
    if (lError != ERROR_SUCCESS) {
        return -1;
    }

    // Smart resource cleaner pattern
    auto clearer = [](HKEY hKey) { RegCloseKey(hKey); };
    auto key = std::unique_ptr<std::remove_pointer<HKEY>::type, decltype(clearer)>(hKeyProcessor, clearer);

    DWORD dwMHz = 0;
    DWORD dwBufferSize = sizeof(DWORD);
    lError = RegQueryValueExA(key.get(), "~MHz", nullptr, nullptr, (LPBYTE)&dwMHz, &dwBufferSize);
    if (lError != ERROR_SUCCESS)
        return -1;

    return dwMHz * 1000000;
#else
    #error Unsupported platform
#endif
}

/***********************************************************************************/
/***********************************************************************************/

bool SystemInformation::CPUHyperThreading()
{
    const std::pair<int, int> cores = CPUTotalCores();
    return (cores.first != cores.second);
}

/***********************************************************************************/
/***********************************************************************************/

int64_t SystemInformation::RamTotal()
{
#if defined(KRATOS_COMPILED_IN_OS)
    int64_t memsize = 0;
    std::size_t size = sizeof(memsize);
    if (sysctlbyname("hw.memsize", &memsize, &size, nullptr, 0) == 0) {
        return memsize;
    }

    return -1;
#elif defined(KRATOS_COMPILED_IN_LINUX)
    const int64_t pages = sysconf(_SC_PHYS_PAGES);
    const int64_t page_size = sysconf(_SC_PAGESIZE);
    if ((pages > 0) && (page_size > 0)) {
        return pages * page_size;
    }

    return -1;
#elif defined(KRATOS_COMPILED_IN_WINDOWS)
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return status.ullTotalPhys;
#else
    #error Unsupported platform
#endif
}

/***********************************************************************************/
/***********************************************************************************/

int64_t SystemInformation::RamFree()
{
#if defined(KRATOS_COMPILED_IN_OS)
    mach_port_t host_port = mach_host_self();
    if (host_port == MACH_PORT_NULL) {
        return -1;
    }

    vm_size_t page_size = 0;
    host_page_size(host_port, &page_size);

    vm_statistics_data_t vmstat;
    mach_msg_type_number_t count = HOST_VM_INFO_COUNT;
    kern_return_t kernReturn = host_statistics(host_port, HOST_VM_INFO, (host_info_t)&vmstat, &count);
    if (kernReturn != KERN_SUCCESS) {
        return -1;
    }

    [[maybe_unused]] int64_t used_mem = (vmstat.active_count + vmstat.inactive_count + vmstat.wire_count) * page_size;
    int64_t free_mem = vmstat.free_count * page_size;
    return free_mem;
#elif defined(KRATOS_COMPILED_IN_LINUX)
    const int64_t pages = sysconf(_SC_AVPHYS_PAGES);
    const int64_t page_size = sysconf(_SC_PAGESIZE);
    if ((pages > 0) && (page_size > 0)) {
        return pages * page_size;
    }

    return -1;
#elif defined(KRATOS_COMPILED_IN_WINDOWS)
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return status.ullAvailPhys;
#else
    #error Unsupported platform
#endif
}

/***********************************************************************************/
/***********************************************************************************/

std::string SystemInformation::GenerateClockSpeed(const int64_t Hertz)
{
    // Define the stream
    std::ostringstream stream;

    // Absolute value
    std::size_t abs_hertz = static_cast<std::size_t>(Hertz);

    // The logic
    if (abs_hertz >= 1000000000) {
        const std::size_t gigahertz = abs_hertz / 1000000000;
        const std::size_t megahertz = (abs_hertz % 1000000000) / 1000000;
        stream << gigahertz << '.' << ((megahertz < 100) ? "0" : "") << ((megahertz < 10) ? "0" : "") << megahertz << " GHz";
    } else if (abs_hertz >= 1000000) {
        const std::size_t megahertz = abs_hertz / 1000000;
        const std::size_t kilohertz = (abs_hertz % 1000000) / 1000;
        stream << megahertz << '.' << ((kilohertz < 100) ? "0" : "") << ((kilohertz < 10) ? "0" : "") << kilohertz << " MHz";
    } else if (abs_hertz >= 1000) {
        const std::size_t kilohertz = abs_hertz / 1000;
        const std::size_t hertz_1000 = abs_hertz % 1000;
        stream << kilohertz << '.' << ((hertz_1000 < 100) ? "0" : "") << ((hertz_1000 < 10) ? "0" : "") << hertz_1000 << " kHz";
    } else {
        stream << abs_hertz << " Hz";
    }

    return stream.str();
}

/***********************************************************************************/
/***********************************************************************************/

std::string SystemInformation::GenerateDataSize(const int64_t Bytes)
{
    // Define the stream
    std::ostringstream stream;

    // Absolute value
    std::size_t abs_bytes = static_cast<std::size_t>(Bytes);

    // The logic
    if (abs_bytes >= (1024ll * 1024ll * 1024ll * 1024ll)) {
        const std::size_t tb = abs_bytes / (1024ll * 1024ll * 1024ll * 1024ll);
        const std::size_t gb = (abs_bytes % (1024ll * 1024ll * 1024ll * 1024ll)) / (1024 * 1024 * 1024);
        stream << tb << '.' << ((gb < 100) ? "0" : "") << ((gb < 10) ? "0" : "") << gb << " TiB";
    } else if (abs_bytes >= (1024 * 1024 * 1024)) {
        const std::size_t gb = abs_bytes / (1024 * 1024 * 1024);
        const std::size_t mb = (abs_bytes % (1024 * 1024 * 1024)) / (1024 * 1024);
        stream << gb << '.' << ((mb < 100) ? "0" : "") << ((mb < 10) ? "0" : "") << mb << " GiB";
    } else if (abs_bytes >= (1024 * 1024)) {
        const std::size_t mb = abs_bytes / (1024 * 1024);
        const std::size_t kb = (abs_bytes % (1024 * 1024)) / 1024;
        stream << mb << '.' << ((kb < 100) ? "0" : "") << ((kb < 10) ? "0" : "") << kb << " MiB";
    } else if (abs_bytes >= 1024) {
        const std::size_t kb = abs_bytes / 1024;
        const std::size_t bytes_1024 = abs_bytes % 1024;
        stream << kb << '.' << ((bytes_1024 < 100) ? "0" : "") << ((bytes_1024 < 10) ? "0" : "") << bytes_1024 << " KiB";
    } else {
        stream << abs_bytes << " bytes";
    }

    return stream.str();
}

} // namespace Kratos