echo "running with 0 level of refinement"

time python cpu_edgebased_lev0.py 0 > cpu_lev0.txt

echo "running with 1 level of refinement"
time python cpu_edgebased_lev0.py 1 > cpu_lev1.txt

echo "running with 2 levels of refinement"
time python cpu_edgebased_lev0.py 2 > cpu_lev2.txt
