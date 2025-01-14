#!/usr/bin/env bash

scriptName="$(basename ${BASH_SOURCE[0]})"
scriptDir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
projectName="$(basename $scriptDir)"

print_help() {
    echo "$scriptName - Configure, build, and install $projectName."
    echo "Usage: $scriptName [OPTION [ARGUMENT]]"
    echo "-h                : print this help and exit."
    echo "-p                : package after building."
    echo "-t build-type     : build type [Debug, Release, RelWithDebInfo] (default: Release)."
    echo "-j job-count      : number of build jobs to launch in parallel (Default: use as many threads as the machine supports)."
    echo "-b build-dir      : build directory."
    echo "-i install-dir    : install directory (local directory by default)."
    echo "-o [opts]         : options/arguments to pass on to CMake. Semicolon (;) delimited, or defined repeatedly."
}

# Default arguments
package=0
buildType="Release"
buildDir="$scriptDir/build"
installDir="$scriptDir/install"
cmakeArguments=""
jobCount=""

while getopts ":h p t: j: b: i: o:" arg; do
    case "$arg" in
        h)  # Print help and exit without doing anything
            print_help
            exit 0
            ;;
        p)  # Package after building
            package=1
            ;;
        t)  # Set build type
            buildType="$OPTARG"
            if ! [[ "${buildType}" = "Debug"            \
                 || "${buildType}" = "RelWithDebInfo"   \
                 || "${buildType}" = "Release" ]]; then
                 print_help
                 echo "Invalid build type: $buildType"
                 exit 1
            fi
            ;;
        j)  # Set the number of build jobs
            if [[ "$OPTARG" =~ ^[1-9][0-9]*$ ]]; then
                jobCount="$OPTARG"
            else
                echo "Error: invalid number of jobs requested: '$OPTARG'"
                exit 1
            fi
            ;;
        b)  # Set build directory
            buildDir="$OPTARG"
            ;;
        i)  # Set install directory
            installDir="$OPTARG"
            ;;
        o)  # Append CMake arguments
            cmakeArguments="$cmakeArguments;$OPTARG"
            ;;
        \?) # Unrecognized argument
            print_help
            echo "Error: unrecognized argument -$OPTARG"
            exit 1
            ;;
    esac
done

case "$(uname -s)" in
    Linux*)
        export cxx=g++
        ;;
    Darwin*)
        if [ -z "$jobCount" ]; then
            jobCount=$(sysctl -n machdep.cpu.thread_count)
        fi

        if [ -z "$cxx" ]; then
            # Set clang from homebrew
            if ! command -v brew &> /dev/null; then
                echo "Error: $scriptName requires Homebrew"
                exit 1
            fi

            if ! brew list llvm >/dev/null 2>&1; then
                echo "Error: missing dependency: llvm"
                echo "Consider running 'brew install llvm'"
                exit 1
            fi

            toolchainRoot="$(brew --prefix llvm)"
            toolchainBin="${toolchainRoot}/bin"
            toolchainLib="${toolchainRoot}/lib"
            toolchainInclude="${toolchainRoot}/include"
            export cc="$toolchainBin/clang"
            export cxx="$toolchainBin/clang++"
        fi
        ;;
    \?)
        echo "Error: unsupported OS $(uname -s)"
        exit 1
esac

# Create or clear the build directory
if ! [ -d "$buildDir" ]; then
    mkdir -p "$buildDir"
else
    rm -f "$buildDir/cmake_install.cmake"
    rm -f "$buildDir/CMakeCache.txt"
fi

# Generate with CMake
if ! cmake                                                  \
    "-H$scriptDir"                                          \
    "-B$buildDir"                                           \
    "-DCMAKE_INSTALL_PREFIX:STRING=$installDir"             \
    "-DCMAKE_BUILD_TYPE:STRING=$buildType"                  \
    "-DCMAKE_CXX_COMPILER:STRING=$cxx"                      \
    "-DCMAKE_COLOR_DIAGNOSTICS:BOOL=ON"                     \
    "-DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON"               \
    "-DMCGS_BUILD_TESTS:BOOL=ON"                            \
    $(echo $cmakeArguments | tr '\;' '\n')                  \
    ; then
    exit 1
fi

# Build and install
if ! cmake --build "$buildDir" --config "$buildType" --target install -j$jobCount; then
    exit 1
fi

# Package
if [ $package -eq 1 ]; then
    cd "$buildDir"
    echo make package
    make package
    make package_source
fi

# Success
exit 0
