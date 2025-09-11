 #!/bin/bash
 
source /VW/python_venv/bin/activate

cd /VW/Wheel_Sources
python3 -m pip install --no-index --find-links=./ -r requirements.txt

cd /data
$1