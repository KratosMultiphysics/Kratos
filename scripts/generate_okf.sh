rm -rf ../okf_py
rm -rf ../okf_cpp
okf-rs generate --lsp ../bin/Release/KratosMultiphysics/ --output ../okf_py/knowledge
okf-rs generate --lsp .. --output ../okf_cpp/knowledge
okf-rs validate ../okf_py/knowledge
okf-rs validate ../okf_cpp/knowledge
