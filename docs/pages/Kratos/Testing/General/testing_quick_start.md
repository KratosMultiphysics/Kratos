---
title: Testing Quick Start
keywords: 
tags: [Testing Test Cpp_Tests Python_Tests Quick Start]
sidebar: kratos_testing
summary: 
---

Kratos includes multiple layers of integrated testing to ensure functional correctness and maintain code stability, even after modifications by collaborators.

This section provides a concise overview of the test structure, key features, and usage guidelines.

## Running Tests

This section outlines the steps to execute various Kratos tests and describes the available options for configuring test runs.

- [Cpp Tests](#cpp-tests)
- [Python Tests](#python-tests)

### CPP Tests

There are several ways to run the different cpp level tests in Kratos:

- [Direct Execution](#direct-execution)
- [Python Launcher](#python-launcher)

#### Direct Execution

The simplest way to run a test is by executing it directly from the corresponding binary file. You can locate the test binaries in the following directory:
`kratos/bin/[Release|Debug|FullDebug]/test`.

The available binaries may vary depending on the applications and components you have compiled.

For example, compiling the `FluidDynamicsApplication`, `MappingApplication`, and `TrilinosApplication` with MPI support enabled will generate the following binaries:

```console
$ ls -la
-rwxr-xr-x 1 user group 60834304 Jun 13 09:56 KratosCoreTest
-rwxr-xr-x 1 user group 10199936 Jun 13 09:56 KratosFluidDynamicsCoreTest
-rwxr-xr-x 1 user group  4140904 Jun 13 09:57 KratosMappingCoreTest
-rwxr-xr-x 1 user group  2633736 Jun 13 09:57 KratosMappingMPICoreTest
-rwxr-xr-x 1 user group  2344424 Jun 13 09:57 KratosMeshMovingCoreTest
-rwxr-xr-x 1 user group  8239136 Jun 13 09:56 KratosMPICoreTest
-rwxr-xr-x 1 user group  4928160 Jun 13 09:56 KratosTrilinosCoreTest
```

To run serial tests (those without MPI in the name), use the following command:

```console
./KratosCoreTest
```

For MPI-based tests, execute them by invoking the binary through the MPI wrapper:

```console
mpiexec -np [N] ./KratosMPICoreTest
```

These commands will produce an output similar to the following:

```console
[==========] Running 1495 tests from 10 test suites.
[----------] Global test environment set-up.
[----------] 936 tests from KratosCoreFastSuite
[ RUN      ] KratosCoreFastSuite.ConstitutiveLawHasMethods
[       OK ] KratosCoreFastSuite.ConstitutiveLawHasMethods (1 ms)
...
[----------] Global test environment tear-down
[==========] 1495 tests from 10 test suites ran. (379719 ms total)
[  PASSED  ] 1488 tests.
```

Several options are available to give you greater control over which tests within a suite are executed, what information is displayed, and to customize the test run according to your needs. For a full list of options, refer to the  [Running Test Programs: Advanced Options](https://google.github.io/googletest/advanced.html#running-test-programs-advanced-options) documentation.

Below is a summary of commonly used options in Kratos:

- **--gtest_list_tests**: Displays the list of available tests in the binary.
- **--gtest_filter=\***: Runs only a specified subset of tests.
- **--gtest_repeat=N**: Repeats the selected tests multiple times (useful for detecting race conditions).
- **--gtest_shuffle**: Randomizes the execution order of tests. When combined with `--gtest_repeat`, each repetition runs in a different order.
- **--gtest_brief=N**: Adjusts the verbosity of the output.

For more details, we also provide several [examples]() demonstrating the usage of these options.

#### Debugging

While historically it has been hard for the Kratos tests to be run alongisde debugging tools due to its reliance on python launcher, this is no longer the case anymore with GTest as they provide a direct access to the C++ layer of Kratos. 

#### Python Launcher

You may also run the tests through the python launcher located in 

### Python Tests

