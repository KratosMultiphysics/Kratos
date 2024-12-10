---
title: Testing Quick Start
keywords: 
tags: [Testing Test Cpp_Tests Python_Tests Quick Start]
sidebar: kratos_testing
summary: 
---

## How to run tests

This section provides insight on how to run the different Kratos tests and some of the options available. 

- [Cpp Tests](#cpp-tests)
- [Python Tests](#python-tests)

### CPP Tests

There are several ways to run the different cpp level tests in Kratos:

- [Direct Execution](#direct-execution)
- [Python Launcher](#python-launcher)

#### Direct Execution

The most straigh forward way to run a test is execute it directly from the binary file containing it. You will be able to find the test binaries in `kratos/bin/[Release|Debug|FullDebug]/test`.

You may see different binaries based on the applications and components that you have compiled. 

For exampe having compiled the `FluidDynamicsApplication`, the `MappingApplication` and `TrilinosApplication` with `MPI` support enabled, will generate the following binaries:

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

You may execute the serial tests (no `MPI` in the name with the following command):

```console
./KratosCoreTest
```

For MPI based tests, you may execute them invoking the binary through the mpi wrapper:

```console
mpiexec -np [N] ./KratosMPICoreTest
```

This commands will result in an output like this:

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

There are serveral options that can be used to have more control over which tests inside a suite are run, what information is shown, and overall tailor the run to your needs. You can find a complete list of options in [Running Test Programs: Advanced Options](https://google.github.io/googletest/advanced.html#running-test-programs-advanced-options)

As a resume, a brief list of some commonly used ones for Kratos are:

- **--gtest_list_tests**: Prints the list of available tests in a binary.
- **--gtest_filter=\***: Runs only a subset of the specified tests.
- **--gtest_repeat=N**: Repeats the selected tests multiple times (usefull to detect race conditions).
- **--gtest_shuffle**: Randomizes the running order of tests. While combined with `--gtest_repeat` will select a different order for every repetition.
- **--gtest_brief=N**: Controls the verbosity of the output. 

We also provide different [examples]() of this different options if you are interested.

#### Debugging

While historically it has been hard for the Kratos tests to be run alongisde debugging tools due to its reliance on python launcher, this is no longer the case anymore with GTest as they provide a direct access to the C++ layer of Kratos. 

#### Python Launcher

You may also run the tests through the python launcher located in 

### Python Tests

