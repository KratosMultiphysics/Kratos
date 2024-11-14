#!/bin/bash
# authors: Máté Kelemen
# Run this script with the -h flag for more info.

# Make sure the working directory is where the script is,
# which should be <kratos_root>/documents.
script_dir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
script_name="$(basename ${BASH_SOURCE[0]})"
kratos_root_dir="$(dirname "${script_dir}")"
cd "$script_dir"

# Function for printing usage information.
print_help() {
    echo "$script_name - Configure and generate doxygen documentation for KratosMultiphysics."
    echo "Usage: $script_name [-h] [-a <path-to-application>] [-A] [-c] [-C]"
    echo "-h                            : print this help and exit."
    echo "-a <path-to-app>              : add an application to generate documentation for (pass repeatedly for multiple apps)."
    echo "-A                            : generate docs for all Kratos applications."
    echo "-e <path-to-external-app>     : add an external application to generate docs for (pass repeatedly for multiple apps)."
    echo "-c                            : add all compiled applications to the list to generate docs for."
    echo "-C                            : clear existing docs."
}

# Parse command line arguments.
applications="" ##< semicolon separated list of app name - app parent dir pairs (components separated by commas).
clean=0

while getopts ":h a: A e: c C" arg; do
    case "$arg" in
        h)  # Print help and exit without doing anything else.
            print_help
            exit 0
            ;;
        a)  # Add an application from the repo apps directory.
            applications="$applications $OPTARG,$kratos_root_dir/applications"
            ;;
        A)  # Add all Kratos applications.
            for app_path in $kratos_root_dir/applications/*Application; do
                app_name="$(basename $app_path)"
                applications="$applications $app_name,$kratos_root_dir/applications"
            done
            ;;
        e)  # Add an application from an external directory. The passed argument
            # must be an app name - app parent directory pair that are separated
            # by a comma.
            app_name=$(echo $OPTARG | cut -d "," -f 1)
            app_dir=$(realpath $(echo $OPTARG | cut -d "," -f 2) )

            if ! [ -d "$app_dir" ]; then
                >&2 echo "$app_dir is not a directory."
                exit 1
            fi

            if ! [ -d "$app_dir/$app_name" ]; then
                >&2 echo "no application named $app_name in $app_dir."
            fi

            applications="$applications $app_name,$app_dir"
            ;;
        c)  # Add all compiled applications.
            if python3 -c "
import KratosMultiphysics
with open('compiled_kratos_applications.csv', 'w') as file:
    for app_name in KratosMultiphysics.python_registry.GetListOfAvailableApplications():
        file.write(f'{app_name};')
            " &> /dev/null; then
                extra=$(cat compiled_kratos_applications.csv | sed "s=;=,$kratos_root_dir/applications =g")
                applications="$applications $extra"
                rm compiled_kratos_applications.csv
            fi
            ;;
        C)  # Set the clear flag.
            clean=1
            ;;
        \?) # Unrecognized argument.
            >&2 echo "Error: unrecognized argument: -$OPTARG"
            >&2 print_help
            exit 1
    esac
done

# Clear existing output if requested, then exit.
if [ $clean -ne 0 ]; then
    if [ -d "$kratos_root_dir/documents/html" ]; then
        rm -rf "$kratos_root_dir/documents/html"
    fi

    if [ -f "$kratos_root_dir/documents/KratosCore.tag" ]; then
        rm "$kratos_root_dir/documents/KratosCore.tag"
    fi

    for pair in $applications; do
        application_name=$(echo $pair | cut -d "," -f 1)
        application_dir=$(echo $pair | cut -d "," -f 2)
        app_doc_dir="$application_dir/$application_name/documents"

        if [ -d "$app_doc_dir/html" ]; then
            rm -rf "$app_doc_dir/html"
        fi

        if [ -f "$app_doc_dir/$application_name.tag" ]; then
            rm "$app_doc_dir/$application_name.tag"
        fi
    done
    exit 0
fi

# First pass:
# Generate the doxygen tag files for core and all applications.
# These will be used by every application to link relevant references
# from their dependencies (most commonly from core).
cd "$kratos_root_dir/documents"
if ! { ( cat doxyfile ; echo "GENERATE_HTML=NO" ) | doxygen - &>tag.log & } ; then
    >&2 echo "Error: generating tag file for core failed."
    exit 1
fi
cd "$kratos_root_dir"

for pair in $applications; do
    application_name=$(echo $pair | cut -d "," -f 1)
    application_dir=$(echo $pair | cut -d "," -f 2)
    if [ -d "$application_dir/$application_name/documents" ]; then
        cd "$application_dir/$application_name/documents"
        if [ -f doxyfile ]; then
            # Generate the tag file.
            if ! { ( cat doxyfile ; echo "GENERATE_HTML=NO" ) | doxygen - &>tag.log & } ; then
                >&2 echo "Error: generating tag file for $application_name failed."
                exit 1
            fi
        else
            # Error if no doxyfile is defined for the application.
            >&2 echo "Error: expecting a doxyfile for $application_name in $(pwd), but found none."
            exit 1
        fi
    else
        # Error if the application has no directory for documentation.
        >&2 echo "Error: $application_name at $application_dir lacks documentation!"
        exit 1
    fi
    cd "$kratos_root_dir"
done

wait
echo "Status: finished generating tag files."

# Second pass:
# Generate HTML documentation.

# Collect all apps that docs get generated for and create an index page for them.
cd "$kratos_root_dir/documents"

echo "# List of Applications" > index.md
echo "- [KratosCore](../README.md)" >> index.md

application_tagfiles=
for pair in $applications; do
    application_name=$(echo $pair | cut -d "," -f 1)
    application_dir=$(echo $pair | cut -d "," -f 2)
    echo "- [$application_name](../applications/$application_name/README.md)" >> index.md
    application_tagfiles="${application_tagfiles}../applications/$application_name/documents/$application_name.tag=../../applications/$application_name/documents/html "
done

# Run doxygen for core
if ! { ( cat doxyfile ; echo "TAGFILES=$application_tagfiles" ; echo "INPUT+=index.md" ) | doxygen - &>doxygen.log & } ; then
    >&2 echo "Error: generating HTML docs for core failed."
    exit 1
fi
cd "$kratos_root_dir"

for pair in $applications; do
    application_name=$(echo $pair | cut -d "," -f 1)
    application_dir=$(echo $pair | cut -d "," -f 2)
    if [ -d "$kratos_root_dir/applications/$application_name/documents" ]; then
        cd "$kratos_root_dir/applications/$application_name/documents"
        if [ -f doxyfile ]; then
            # Run doxygen.
            if ! { doxygen doxyfile &>doxygen.log & } ; then
                >&2 echo "Error: generating HTML docs for $application_name failed."
                exit 1
            fi
        fi
    fi
    cd "$kratos_root_dir"
done

wait
echo "Status: finished generating documentation."
