#!/bin/bash
#
# Script to check include/test code for common pybind11 code style errors.
#
# This script currently checks for
#
# 1. use of tabs instead of spaces
# 2. MSDOS-style CRLF endings
# 3. trailing spaces
# 4. missing space between keyword and parenthesis, e.g.: for(, if(, while(
# 5. Missing space between right parenthesis and brace, e.g. 'for (...){'
# 6. opening brace on its own line. It should always be on the same line as the
#    if/while/for/do statment.
# 7. Leftover markers denoting incomplete implementations/tasks
#
# Invoke as: tools/check-style.sh
#

errors=0
IFS=$'\n'
found=
# The mt=41 sets a red background for matched tabs:
GREP_COLORS='mt=41' GREP_COLOR='41' grep $'\t' include/ src/ tests/*.{h,cpp} docs/*.rst -rn --color=always |
while read f; do
    if [ -z "$found" ]; then
        echo -e '\033[31m\033[01mError: found tabs instead of spaces in the following files:\033[0m'
        found=1
        errors=1
    fi

    echo "    $f"
done

found=
# The mt=41 sets a red background for matched MS-DOS CRLF characters
GREP_COLORS='mt=41' GREP_COLOR='41' grep -IUlr $'\r' include/ src/ tests/*.{h,cpp} docs/*.rst --color=always |
while read f; do
    if [ -z "$found" ]; then
        echo -e '\033[31m\033[01mError: found CRLF characters in the following files:\033[0m'
        found=1
        errors=1
    fi

    echo "    $f"
done

found=
# The mt=41 sets a red background for matched trailing spaces
GREP_COLORS='mt=41' GREP_COLOR='41' grep '[[:blank:]]\+$' include/ src/ tests/*.{h,cpp} docs/*.rst -rn --color=always |
while read f; do
    if [ -z "$found" ]; then
        echo -e '\033[31m\033[01mError: found trailing spaces in the following files:\033[0m'
        found=1
        errors=1
    fi

    echo "    $f"
done

found=
grep '\<\(if\|for\|while\|catch\)(\|){' include/ src/ tests/*.{h,cpp} -rn --color=always |
while read f; do
    if [ -z "$found" ]; then
        echo -e '\033[31m\033[01mError: found the following coding style problems:\033[0m'
        found=1
        errors=1
    fi

    echo "    $f"
done

found=
GREP_COLORS='mt=41' GREP_COLOR='41' grep '^\s*{\s*$' include/ src/ tests/*.{h,cpp} -rn --color=always |
while read f; do
    if [ -z "$found" ]; then
        echo -e '\033[31m\033[01mError: braces should occur on the same line as the if/while/.. statement. Found issues in the following files: \033[0m'
        found=1
        errors=1
    fi

    echo "    $f"
done

found=
GREP_COLORS='mt=41' GREP_COLOR='41' grep '\<\(TODO\|XXX\)' include/ src/ tests/*.{h,cpp} -rn --color=always |
while read f; do
    if [ -z "$found" ]; then
        echo -e '\033[31m\033[01mError: Incomplete implementation markers in code. Found issues in the following files: \033[0m'
        found=1
        errors=1
    fi

    echo "    $f"
done

exit $errors
