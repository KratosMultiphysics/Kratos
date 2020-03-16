rm -f AbaqusHosts.txt
NPROC=$(nproc --all)
NAME=$(hostname | sed -e 's/.ugent.be//')

for i in $(seq 1 $NPROC)
do
echo $NAME >> AbaqusHosts.txt
done
