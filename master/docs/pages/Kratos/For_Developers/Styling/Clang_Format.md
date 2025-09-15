---
title: How to configure clang‐format
keywords: 
tags: [How-to-configure-clang‐format.md]
sidebar: kratos_for_developers
summary: 
---

# Overview
_clang-format_ is an utility to format source code in several languages according to predefined settings.
Settings for auto-formatting of C++ source files in Kratos are summarized below, and should be saved in _.clang-format_ file in the Kratos root directory. These settings are chosen to closely match the formatting style used in Kratos.

# Configuring Editors
## Visual Studio Code
Install the C/C++ extension and customize settings.json to include:

    "C_Cpp.clang_format_path": "/path/to/<clang-format-executable>",
    "C_Cpp.clang_format_style": "file",
    "C_Cpp.clang_format_fallbackStyle": "none"


## Emacs
Add the line

    (load "/path/to/clang-format.el")


to ~/.emacs.d/init.el. Format a source file with `M-x clang-format-region`.
## CLion
- Add the .clang-format file to the Kratos root directory as explained above
- Go to File->Settings->Tools->External Tools and click on the plus sign. A window should pop up. Choose a name, for example "clang-format"
- For the Tool settings tab I'm using this configuration:
    - Program: clang-format (you should use the name of your executable here)
    - Parameters: --style=file -i $FileName$
    - Working directory: $FileDir$

Now, with your file open, you can go to Tools->External tools and run the config above. It basically calls clang-format and does inplace formatting using the style define in the first .clang-format file found in a parent directory.

# .clang-format file
(tested with clang-format version 3.5, 3.8)

    Language: Cpp
    AccessModifierOffset: -4
    AlignEscapedNewlinesLeft: true
    AlignTrailingComments: true
    AllowAllParametersOfDeclarationOnNextLine: true
    AllowShortBlocksOnASingleLine: false
    AllowShortFunctionsOnASingleLine: None
    AllowShortIfStatementsOnASingleLine: false
    AllowShortLoopsOnASingleLine: false
    AlwaysBreakBeforeMultilineStrings: true
    AlwaysBreakTemplateDeclarations: true
    BinPackParameters: false
    BreakBeforeBinaryOperators: false
    BreakBeforeBraces: Stroustrup
    BreakBeforeTernaryOperators: true
    BreakConstructorInitializersBeforeComma: false
    ColumnLimit: 80
    CommentPragmas: ''
    ConstructorInitializerAllOnOneLineOrOnePerLine: true
    ConstructorInitializerIndentWidth: 4
    ContinuationIndentWidth: 4
    DerivePointerAlignment: false
    DisableFormat: false
    ExperimentalAutoDetectBinPacking: false
    IndentCaseLabels: false
    IndentWidth: 4
    IndentWrappedFunctionNames: false
    IndentFunctionDeclarationAfterType: false
    MaxEmptyLinesToKeep: 1
    KeepEmptyLinesAtTheStartOfBlocks: false
    NamespaceIndentation: None
    PenaltyBreakBeforeFirstCallParameter: 1
    PenaltyBreakComment: 300
    PenaltyBreakString: 1000
    PenaltyBreakFirstLessLess: 120
    PenaltyExcessCharacter: 1
    PenaltyReturnTypeOnItsOwnLine: 1000
    PointerAlignment: Left
    SpaceBeforeAssignmentOperators: true
    SpaceBeforeParens: ControlStatements
    SpaceInEmptyParentheses: false
    SpacesBeforeTrailingComments: 1
    SpacesInAngles: false
    SpacesInCStyleCastParentheses: false
    SpacesInContainerLiterals: true
    SpacesInParentheses: false
    Cpp11BracedListStyle: true
    Standard: Cpp11
    TabWidth: 4
    UseTab: Never