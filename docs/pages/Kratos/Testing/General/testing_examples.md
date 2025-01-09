---
title: Testing Examples
keywords: 
tags: [Testing Test Cpp_Tests Python_Tests Examples]
sidebar: kratos_testing
summary: 
---

# How to run tests

This section will provided some common use cases for Kratos tests.

Unless stated otherwise, we asume that we running `KratosCoreTest` or `KratosMPICoreTest`. We will consider the `Array1DTest` as an example test

Note that even if not present in this page, you can mix and combine many of the different flags shown in the examples.

## Exectuion

### Running a suite of tests

```console
./KratosCoreTests
```

### Running a single test from a suite

```console
./KratosCoreTests --gtest_filter="KratosCoreFastSuite.Array1DTest"
```

```console
Kratos::Testing::GTestMain::InitializeTesting
Note: Google Test filter = KratosCoreFastSuite.Array1DTest
[==========] Running 1 test from 1 test suite.
[----------] Global test environment set-up.
[----------] 1 test from KratosCoreFastSuite
[ RUN      ] KratosCoreFastSuite.Array1DTest
[       OK ] KratosCoreFastSuite.Array1DTest (1 ms)
[----------] 1 test from KratosCoreFastSuite (1 ms total)

[----------] Global test environment tear-down
[==========] 1 test from 1 test suite ran. (2 ms total)
[  PASSED  ] 1 test.
```

### Running all tests matching a patter from a suite

In this example we make use of a regex `KratosCoreFastSuite.*` in the filter in order to run all the tests from the `KratosCoreFastSuite`
isnide the KratosCoreTests binary.

```console
./KratosCoreTests --gtest_filter="KratosCoreFastSuite.*"
```

After executuon the command we will se something like:

```console
[       OK ] KratosCoreFastSuite.XmlOStreamWriterWriteDataElementBinaryMixed (1 ms)
[----------] 936 tests from KratosCoreFastSuite (378175 ms total)
[----------] Global test environment tear-down
[==========] 936 tests from 1 test suites ran. (379719 ms total)
```

### Running a single tests from a suite multiple times

In this example we make use of the ... to run a single test multiple times.

We take the change to encourage you to make small tests to make the best use of options like this.

Be warned that negative values will repeat the test endlessly

```console
./KratosCoreTests --gtest_repeat=10 --gtest_filter="KratosCoreFastSuite.Array1DTest"
```

If negative values are combined with `--gtest_fail_fast` which causes the test to stop at the first failure it is very usefull to catch race conditions that only happens a given a small number of times of the total amount of executions. 

```console
./KratosCoreTests --gtest_repeat=-1 --gtest_fail_fast --gtest_filter="KratosCoreFastSuite.[YourRaceConditionTest]"
```

### Randomizing the execution of a subset of tests from a suite

Historically Kratos tests have had the problem of sharing a common Initialized kernel which lead to situation in which the test execution order could result in problems. While this should be minimized in the current framework were every tests is run in its own initializes kernel, GTest still provides a mechanism to randomize the execution order in case hidden dependencies causing problems are still present in the code.

In this example we show how to randomize the execution order of the tests inside a binary.

```console
./KratosCoreTests --gtest_filter="KratosCoreFastSuite.*"
```

### Get the list of available tests from a suite

```console
./KratosCoreTests --gtest_list_tests
```

## Debugging

### Running a single test from a suite with GDB

```console
./KratosCoreTests --gtest_catch_exceptions=0 --gtest_filter="KratosCoreFastSuite.Array1DTest"
```

### Running multiple tests with GDB stoping at the first failure

```console
gdb --args ./KratosCoreTests --gtest_catch_exceptions=0 --gtest_break_on_failure --gtest_filter="KratosCoreFastSuite.*"
```

### Running a single test from a suite with valgrind

```console
valgrind ./KratosCoreTests --gtest_catch_exceptions=0 --gtest_break_on_failure --gtest_filter="KratosCoreFastSuite.Array1DTest"
```

### Running a single test from a suite with valgrind (to discover memory leaks)

```console
valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./KratosCoreTests --gtest_catch_exceptions=0 --gtest_break_on_failure --gtest_filter="KratosCoreFastSuite.*"
```

## MPI

### Running multiple test with MPI with 4 procs

```console
OMP_NUM_THREADS=1 mpirun -np 4 ./KratosMPICoreTests --gtest_filter="KratosMPICoreFastSuite.*"
```

