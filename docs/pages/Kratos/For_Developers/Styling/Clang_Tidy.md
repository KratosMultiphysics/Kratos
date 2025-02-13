---
title: Clang Tidy to automatically correct code
keywords: 
tags: [How-to-use-Clang-Tidy-to-automatically-correct-code.md]
sidebar: kratos_for_developers
summary: 
---

# Install 

*Clang-Tidy* is a tool developed and maintained by the *Clang/LLVM* community. The official documentation can be found at [http://clang.llvm.org/extra/clang-tidy/](http://clang.llvm.org/extra/clang-tidy/). When running Linux, *clang-tidy* is usually easy to get via your distribution’s package manager. On Ubuntu Linux:

```console
sudo apt-get install clang-tidy clang-tools
```
 
 Additionally you can download directly in the [project page](http://releases.llvm.org/download.html).
 
 # Introduction
 
 A typical invocation of the command-line tool looks like this:
 
 ```console
 clang-tidy test.cpp -- -Imy_project/include -DMY_DEFINES ...
 ```
 
 Executing it like this, the tool will print a bunch of warnings and notes (if applicable), in exactly the same way Clang/GCC provide diagnostics, too.
 
## Available checkers

 The avalaible checkers avaible can be found with the command `list-checks`. To correct the code we can use the command `-fix` (not all the avalaible checkers have the `-fix` function). For example for `modernize` and `performance`.
 
### Modernize

To modernize code we have the following options:

 ```console
 clang-tidy --list-checks -checks='*' | grep "modernize"
    modernize-avoid-bind
    modernize-avoid-c-arrays
    modernize-concat-nested-namespaces
    modernize-deprecated-headers
    modernize-deprecated-ios-base-aliases
    modernize-loop-convert
    modernize-make-shared
    modernize-make-unique
    modernize-pass-by-value
    modernize-raw-string-literal
    modernize-redundant-void-arg
    modernize-replace-auto-ptr
    modernize-replace-disallow-copy-and-assign-macro
    modernize-replace-random-shuffle
    modernize-return-braced-init-list
    modernize-shrink-to-fit
    modernize-unary-static-assert
    modernize-use-auto
    modernize-use-bool-literals
    modernize-use-default-member-init
    modernize-use-emplace
    modernize-use-equals-default
    modernize-use-equals-delete
    modernize-use-nodiscard
    modernize-use-noexcept
    modernize-use-nullptr
    modernize-use-override
    modernize-use-trailing-return-type
    modernize-use-transparent-functors
    modernize-use-uncaught-exceptions
    modernize-use-using
 ```
 
For example to add the missing overrides in the following code:
 
 ```c
 struct Base {
    virtual void reimplementMe(int a) {}
};
struct Derived : public Base  {
    virtual void reimplementMe(int a) {}
};
 ```
 
 We use the next command:
 
 ```console
 clang-tidy -checks='modernize-use-override' -fix test.cpp -- -std=c++11
 ```
 
### Performance

For performance checks we have the following options:

 ```console
clang-tidy --list-checks -checks='*' | grep "performance"                                                                                                                                                                   
    clang-analyzer-optin.performance.GCDAntipattern
    clang-analyzer-optin.performance.Padding
    performance-faster-string-find
    performance-for-range-copy
    performance-implicit-conversion-in-loop
    performance-inefficient-algorithm
    performance-inefficient-string-concatenation
    performance-inefficient-vector-operation
    performance-move-const-arg
    performance-move-constructor-init
    performance-no-automatic-move
    performance-no-int-to-ptr
    performance-noexcept-move-constructor
    performance-trivially-destructible
    performance-type-promotion-in-math-fn
    performance-unnecessary-copy-initialization
    performance-unnecessary-value-param
 ```

 # Using to correct a whole project
 
 The previous example will work just with a very simple example contained in one file. To correct the whole project we will need to create a `json` file containing all the file in the project. For that we add the following line to ou configure.sh:
 
 ```console
 -DCMAKE_EXPORT_COMPILE_COMMANDS=ON                                                       \
 ```
 
 This will create a file named `compile_commands.json` that we will use with the following python script from the [LLVM project](https://github.com/llvm/llvm-project/blob/main/clang-tools-extra/clang-tidy/tool/run-clang-tidy.py)).
 
 Once we have the script and the json file we can check and fix the whole project by the following way:
 
 ```console
 run-clang-tidy.py -header-filter='.*' -checks='-*,modernize-use-override' -fix
 ```

You can run simmultaneously all the possible modernize commands using the following [shell script](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Resources_files/Clang-tidy%20modernize/modernize.sh).

 ```console
sh modernize.sh
 ```
 
