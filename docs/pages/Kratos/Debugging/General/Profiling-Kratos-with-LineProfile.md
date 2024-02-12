---
title: Profiling Kratos with LineProfile
keywords: 
tags: [Profiling-Kratos-with-LineProfile.md]
sidebar: kratos_debugging
summary: 
---

## Profiling Kratos with LineProfile

This page details the steps to follow in order to profile Kratos with a python line profiler. Please be aware that this is an intrusive process, and you will have to modify some of your code.

Obtaining the profiler:
You can find the python line profiler an all the relevant information here: https://github.com/pyutils/line_profiler.

If you prefer only to install it just:
```console
python -m pip install line_profiler
```

### Preparing your code.
The python line profiler needs to know beforehand which functions are you interested in debugger. In that sense, it is helpful that you already have a slight idea on what sections of your code are interested on.

In order to mark such parts, you will have to use the `@profile` decorator. A example is presented below:

```python
@profile
def my_function():
    var = calculate()
    ...  etc  ...
```

Unfortunately, this process will make your code break if you try to run it without the profiler. We are working on integrating a mechanism that solved this inside Kratos, but in the meantime we recommend you to use a launcher script to run your code. This will allow to ignore the decorators without removing them from the code if you intend to run the code normally.

Such launcher script have to look like this:

```python
import os
import sys
import builtins
import line_profiler

prof = line_profiler.LineProfiler()

builtins.__dict__['profile'] = prof

exec(open(sys.argv[1]).read())

with open("profile_"+str(os.getpid())+".prof", "w+") as prof_file:
    prof.print_stats(prof_file)

print("Profile output written in:", "profile_"+str(os.getpid())+".prof")
```

We do recommend to launch your case with this launcher as well if you intend to run with MPI, as it will automatically generate a different report for every process. An example of an invocation call will look like this:

Serial
```console
python profiler.py MainKratos.py
```

Distributed
```console
mpirun -np [N] --output-filename profile python profiler.py MainKratos.py
```
## Analyzing output

Parser.py
```python
import os
import re
import sys


def ReadValueAs(type, item):
    return type(item) if item.strip() != "" else ""

def ReadBlock(my_file, json_object):
    total_time = my_file.readline()             # Total time or end of file

    if total_time == "" or total_time == "\n":
        return True
    else:
        profiled_file = my_file.readline().rstrip().split(" ")[1]
        profiled_function = my_file.readline().split(" ")[1]

        my_file.readline()                   # Space
        my_file.readline()                   # Headings
        my_file.readline()                   # ========

        stop = False

        function_lines = {}

        if profiled_file not in json_object:
            json_object[profiled_file] = {}

        while not stop:
            profile_line = my_file.readline()

            if profile_line == "" or profile_line == "\n":
                stop = True
            else:
                json_line = {"Line":None, "Hits":None, "Time":None, "Per Hit":None, "% Time":None, "Line Content":None}

                # This may need some tunning. I am not sure how to correctly detect the separators here.
                column_size = 48
                line_items = re.split(r'\s+', profile_line[0:column_size].strip())

                # Has info
                if len(line_items) == 5:
                    json_line["Line"]           = ReadValueAs(int,   line_items[0])
                    json_line["Hits"]           = ReadValueAs(int,   line_items[1])
                    json_line["Time"]           = ReadValueAs(float, line_items[2])
                    json_line["Per Hit"]        = ReadValueAs(float, line_items[3])
                    json_line["% Time"]         = ReadValueAs(float, line_items[4])

                    json_line["Line Content"]   = profile_line[column_size:-1].rstrip()

                    function_lines[json_line["Line"]] = json_line

        json_object[profiled_file][profiled_function] = function_lines

        return False

def ParseFile(filename):
    json_object = {}

    with open(filename, "r") as parse_file:
        timer_unit = parse_file.readline()          # Timer
        parse_file.readline()                       # Space

        is_end_block = False

        while not is_end_block:
            is_end_block = ReadBlock(parse_file, json_object)    # Block Content

    return json_object

def ComputeTotalTime(trace):
    total_time = 0.0

    for filename in trace:
        for function in trace[filename]:
            for line in trace[filename][function]:
                if trace[filename][function][line]["Time"] != '':
                    total_time += trace[filename][function][line]["Time"]

    return total_time
```

Analyzer.py
```python
import parser

profile_outputs = {
    2: parser.ParseFile("profile-2"),
    4: parser.ParseFile("profile-4"),
    8: parser.ParseFile("profile-8"),
    16: parser.ParseFile("profile-16"),
    32: parser.ParseFile("profile-32"),
    64: parser.ParseFile("profile-64"),
    128: parser.ParseFile("profile-128")
}

def CalculateTotalSpeedUp(input_files, ref):
    speedups = {}
    for proc in input_files:
        speedups[proc] = parser.ComputeTotalTime(input_files[proc])

    for proc in input_files:
        if(proc != ref):
            speedups[proc] = speedups[ref]/speedups[proc]
    speedups[ref] = 1

    print(speedups)

def CalculateSpeedUpByLine(input_files, ref):
    speedups = {}
    for proc in input_files:
        trace = input_files[proc]
        for filename in trace:
            for function in trace[filename]:
                for line in trace[filename][function]:
                    if trace[filename][function][line]["Time"] != '' and trace[filename][function][line]["Time"] != 0.0:
                        trace[filename][function][line]["SpeedUp"] = input_files[ref][filename][function][line]["Time"] / input_files[proc][filename][function][line]["Time"]

def ExtractExpensiveLines(input_files, max_core_num, max_line_count, check_criteria):
    most_relevant_entry_lines = []

    trace = input_files[max_core_num]
    for filename in trace:
        for function in trace[filename]:
            for line in trace[filename][function]:
                if trace[filename][function][line][check_criteria] != '':
                    most_relevant_entry_lines.append({"Line": trace[filename][function][line], "Function":function, "File":filename})

    most_relevant_entry_lines.sort(key=lambda x: x["Line"][check_criteria], reverse=True)

    return most_relevant_entry_lines[0:max_line_count]


def GenerateFieldTable(input_files, function_lines, field, min_procs, max_procs, scale_to_min_procs):
    for e_line in function_lines:

        file_name = e_line["File"]
        function_name = e_line["Function"]
        function_line = e_line["Line"]

        buffer = str(file_name.split("/")[-1])+":"+str(function_name)+":"+str(function_line["Line"])+":"+str(function_line["Line Content"])

        for proc in input_files:
            trace = input_files[proc]
            value = trace[file_name][function_name][function_line["Line"]][field]

            if field == "SpeedUp" and scale_to_min_procs:
                value = value * min_procs

            buffer += "|" + str(value)

        buffer += "|" + str(input_files[max_procs][file_name][function_name][function_line["Line"]]["Time"]/1e6)
        buffer += "|" + str(input_files[max_procs][file_name][function_name][function_line["Line"]]["Hits"])

        print(buffer)

min_procs = 1
max_procs = 4

CalculateSpeedUpByLine(profile_outputs, min_procs)

most_relevant_entry_lines = ExtractExpensiveLines(profile_outputs, max_procs, 20, "Time")

GenerateFieldTable(profile_outputs, most_relevant_entry_lines, "SpeedUp", min_procs, max_procs, True)
CalculateTotalSpeedUp(profile_outputs, min_procs)
```
