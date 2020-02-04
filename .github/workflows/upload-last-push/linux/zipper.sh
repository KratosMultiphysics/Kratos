cd /zip
python3 -m pip install pyinstaller
python3 -m pip install scipy matplotlib

pyinstaller runkratos.py --add-data /__w/Kratos/Kratos/bin/Release/KratosMultiphysics:./KratosMultiphysics --add-binary /__w/Kratos/Kratos/bin/Release/libs/*:. -p .
