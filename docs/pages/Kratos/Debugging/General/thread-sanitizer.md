---
title: Thread Sanitizer
keywords: 
tags: [Thread Sanitizer Debug]
sidebar: kratos_debugging
summary: 
---

The aim of this tutorial is to briefly describe how to compile Kratos with ThreadSanitizer (or TSan) and run a test example from a python script. If someone is interested in the details of ThreadSanitizer and data race bugs, please read the manual [Thread Sanitizer Cpp Manual](https://github.com/google/sanitizers/wiki/).

Let us notice that ThreadSanitizer is supported on Linux x86_64 and is part of clang 3.2 and gcc 4.8. This tutorial has been tested in Ubuntu 16.04 with gcc 6.3.0, compiling Kratos in release mode.

First, modify the CMAKE_CXX_FLAGS in your configure.sh with `-fsanitize=thread` and `-ltsan` flags:

```console
 -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 -std=c++11 -fsanitize=thread -ltsan"     \
```
{: data-lang="configure.sh"}

Once the compilation finishes, test examples must be run by '''pre-loading the ThreadSanitizer library''' of your compiler. For instance, in our case it was:

```console
 LD_PRELOAD=/path/to/gcc-6.3.0/lib64/libtsan.so python3 MainKratos.py
```
{: data-lang="Terminal Command"}

Visit [Thread Sanitizer Report Format](https://github.com/google/sanitizers/wiki/ThreadSanitizerReportFormat) for explanation of reports format.