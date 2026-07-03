$RELEASE_BRANCH="Release-10.4.3"

$HOST_SWAP="C:/data_swap_host"
$GUEST_SWAP="/data_swap_guest"

docker run --name release_build -v ${HOST_SWAP}:${GUEST_SWAP} --entrypoint="/bin/bash" kratosmultiphysics/kratos-wheelbuilder-linux-2-28 /workspace/scripts/start.sh ${RELEASE_BRANCH}
docker rm release_build