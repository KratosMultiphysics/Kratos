---
title: Logger
keywords: 
tags: [How To Use Logger]
sidebar: kratos_for_developers
summary: describes the usage of the Kratos Logger
---

## Structure

The logging system in Kratos has 3 important parts:
- `LoggerMessage` a data class storing the message with some attributes like label, category, severity
- `Logger` is a singleton object in charge of gathering all messages produced in the code and passing them to the outputs
- `LoggerOutput` takes message and write it to the output (or file). This is the extension point of the logger and one can create its output to write only some messages (filtering by category, severity, label, etc.) with custom format. 

## Sending a Message in C++
The simplest way to send a message is to use predefined macros:
```cpp
KRATOS_INFO("Some label") << "Some message with value: " << 3.14 << " with more message";
```
This example sends a message with `INFO` severity and `STATUS` category. There are macros for each severity level:
- `KRATOS_WARNING` for reporting a warning without interrupting the simulation
- `KRATOS_INFO` is the standard level in which minimum level of information should be provided
- `KRATOS_DETAIL` is to make more detailed output. All above macros should be used with low frequencies and should be avoided in fine grain parts like elements and conditions. 
- `KRATOS_TRACE` which will be enabled only in debug mode and will be ignored completely in release modes. This is the most verbose level for debugging only.
- `KRATOS_CHECK_POINT`is an special macro which to be used for quality control and regression tests and is enabled only when `KRATOS_ENABLE_CHECK_POINT` is defined. Providing outputs with `CHECKING` category to be used in filtering output for regression tests.

Apart from severity each message has a category which can be:
- `STATUS` to be used to indicate status of solution like passing some point, if is converged, some steps started or finished, etc.
- `CRITICAL` is to indicate an exceptional nature of the message used for errors and warnings.
- `STATISTICS` indicates that the message has statistical information like number of iterations, solver residual, mesh quality and so on.
- `PROFILING` to be used for timing (not implemented yet)
- `CHECKING` is for regression tests and quality control messages.


One can change the category by passing it to the message: 
```cpp
KRATOS_INFO("Number of Iterations") << number_of_iterations << LoggerMessage::Category::STATISTICS;
```

Note that each logger output can filter the messages by their category and severity and show only the ones it wants. For example an output associated with statistical file can only print the messages with `STATISTICS` as category. 

### List of macros
The logger offers several utility macros that allow more fine grain control over when a message needs to be printed. The list of the current available macros is the following:

#### Default
The standard output, prints the label alongside the message

**Variants**:
  - `KRATOS_INFO("Some Label")`
  - `KRATOS_WARNING("Some Label")`
  - `KRATOS_DETAIL("Some Label")`
  - `KRATOS_TRACE("Some Label")`
  - `KRATOS_CHECK_POINT("Some Label")`

**Example**:
```cpp
KRATOS_INFO("Example") << "Hello World." << std::endl;
```
Output:
```console
Hello World.
```

#### IF
Prints the message only if condition is evaluated to true

**Variants**:
  - `KRATOS_INFO_IF("Some Label", condition)`
  - `KRATOS_WARNING_IF("Some Label", condition)`

**Example**:
```cpp
for(std::size i = 0; i < 4; i++) {
    KRATOS_INFO_IF("Example", i % 2) << "Iter: " << i << std::endl;
}
```
Output:
```console
Iter: 1
Iter: 3
```

#### ONCE
Prints the message only once

**Variants**:
  - `KRATOS_INFO_ONCE("Some Label")`
  - `KRATOS_WARNING_ONCE("Some Label")`

**Example**:
```cpp
for(std::size i = 0; i < 4; i++) {
    KRATOS_INFO_ONCE("Example") << "Iter: " << i << std::endl;
}
```
Output:
```console
Iter: 0
```

#### FIRST_N
Prints the message the first N times

**Variants**:
  - `KRATOS_INFO_FIRST_N("Some Label", count)`
  - `KRATOS_WARNING_FIRST_N("Some Label", count)`

**Example**:
```cpp
for(std::size i = 0; i < 4; i++) {
    KRATOS_FIRST_N("Example", 2) << "Iter: " << i << std::endl;
}
```
Output:
```console
Iter: 0
Iter: 1
```

### Sending a Message in Python
Logger is also has a limited version available. In order to use the logger you must import it with `KratosMultihpysics` module.

```python
from KratosMultiphysics import Logger
# or
from KratosMultiphysics import *
```

It has three different print functions you can use: __Print__, __PrintInfo__ or __PrintWarning__.
All functions have the same syntax and print the message in the default output channel, which is the stdout, in other words your terminal. For __PrintInfo__ and __PrintWarning__ the first argument will always be the label, while for __Print__ you will have to set it using the `label` named argument. For a message splitted in several arguments, spaces are added automatically between the parts of the message when printing.

```python
Logger.PrintInfo("Message Label", "This", "is", "a message") 
Logger.PrintWarning("Message Label", "This", "is", "a message")
Logger.Print("This", "is", "a", "message", label="Message Label") # Prints a message
```

You can specify the severity and category of the message by using the _severity_ and _category_ named arguments, similarly to the print function from python.
Both `PrintInfo` and `PrintWarning` have they severity set to __INFO__ and __WARNING__ and we recommend not to change them. `Print` function has its channel set to __INFO__ by default but this function is intended for you to change its severity and category.

`Category` and `Severity` levels are found in their respective namespaces inside the `Logger` module:

```python
Logger.PrintWarning("Message Label Warning", "This", "is", "a", "message") # Prints a message in the warning channel
Logger.Print("This", "is", "a", "message", severity=Logger.Severity.WARNING) # Does the same
```

Finally, you can also interact with the default output and change its `Category` and `Severity` but is not possible to change the output itself at the moment:

```python
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
Logger.PrintInfo("Message Label Info", "You won't see this")
Logger.PrintWarning("Message Label Warning", "But you will see this")
Logger.Print("And this", severity=Logger.severity.WARNING, label="Message Label Custom")
```

Result:
```console
Message Label Warning: But you will see this
Message Label Custom: And this
```

### Writing Log Files
Log Files can easily be created by adding a file logger to the Logging system. An arbitrary number of file loggers can be created and added. Each file logger can be configured individually. The following code shows how a file logger for warnings can be added in python: 

```python
import KratosMultiphysics as KM

file_logger = KM.FileLoggerOutput("KratosWarning.log")
file_logger.SetSeverity(KM.Logger.Severity.WARNING)
KM.Logger.AddOutput(file_logger)
```

## Writing a Logger Output
Creating a custom output is relatively easy by deriving a new output class from `LoggerOutput `class. The following code shows a sample implementation for an output which reports only the warnings in a csv format:

```cpp
class WarningCSVOutput : public LoggerOutput{
public:
  WarningCSVOutput(std::ostream& rOutputStream) : LoggerOutput(rOutputStream){}

  void WriteMessage(LoggerMessage const& TheMessage) override {
    auto message_severity = TheMessage.GetSeverity();
    if (message_severity == LoggerMessage::Severity::WARNING){
	(*this) << TheMessage.GetLabel() << " , " << TheMessage.GetMessage() << std::endl;
    }
  }
};
``` 
To add your output to the Logging system one should use the `Logger::AddOutput` method:

```cpp
WarningCSVOutput warning_csv_output(output_file);
Logger::AddOutput(warning_csv_output);
```

For another example see also the `FileLoggerOutput` class. 

