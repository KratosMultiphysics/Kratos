# This function automatically configures a given application to build its benchmarks
macro(kratos_add_benchmarks)
    set(options USE_MPI USE_CUSTOM_MAIN)
    set(oneValueArgs TARGET WORKING_DIRECTORY)
    set(multiValueArgs SOURCES)

    cmake_parse_arguments(KRATOS_ADD_BENCHMARK "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(KRATOS_ADD_BENCHMARK_SOURCES)
        include(GoogleBenchmark)
    endif(KRATOS_ADD_BENCHMARK_SOURCES)

endmacro(kratos_add_benchmarks)