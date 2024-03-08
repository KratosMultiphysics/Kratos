$RELEASE_BRANCH="Release-9.4.6"

$HOST_SWAP="C:/data_swap_host"
$GUEST_SWAP="/data_swap_guest"

docker run --name release_build -v ${HOST_SWAP}:${GUEST_SWAP} --entrypoint="/bin/bash" kratosmultiphysics/kratos-wheelbuilder-linux-mpi /workspace/scripts/start.sh ${RELEASE_BRANCH}
docker rm release_build