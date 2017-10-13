echo "running with 0 level of refinement"

time python gpu_edgebased_lev0.py 0 > lev0.txt

echo "running with 1 level of refinement"
time python gpu_edgebased_lev0.py 1 > lev1.txt

echo "running with 2 levels of refinement"
time python gpu_edgebased_lev0.py 2 > lev2.txt
