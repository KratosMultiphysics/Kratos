rm *.inp *.pyc *.rpy ../Base.inp
abaqus cae noGUI=MakeInp.py
cp -r Base.inp ../Base.inp
