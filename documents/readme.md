# Documentation Structure

The documentation of `KratosMultiphysics` mainly consists of two parts:
- annotated source code
- standalone markdown docs

This directory is dedicated to the configuration and generation of the technical documentation from annotated sources using `doxygen`.

The `doxygen` documentation is broken up into the mostly self-contained units per application and core. References to other applications and core are handled via [tag files](https://www.doxygen.nl/manual/external.html) (essentially indexes of other documentations). To handle the partitioned generation of docs, run the [`mkdocs.sh`](mkdocs.sh) script that handles tag files and the assembly of the HTML documentation.

# Custom Applications

Should you wish to create `doxygen` documentation for your own application, run `mkdocs.sh`, create a separate configuration file (`doxyfile`) for your app and reference the tag file of core(`KratosCore.tag`) as well as any application yours depends on. Finally, run `doxygen` for your app.

# Dependencies
- `doxygen` for generating HTML documentation from annotated source code.
- `graphviz` for [`dot`](https://www.howtoinstall.co/en/ubuntu/xenial/graphviz) functionalities, which is required for graphs and UML diagrams.
