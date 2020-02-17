cd CSM
ls -1 | egrep -v '(makeHostFile.sh|AbaqusHosts.txt)' | xargs -I files rm -r "files"

