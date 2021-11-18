$RELEASE_BRANCH="wheels/updating_triangle"

$HOST_SWAP="C:/data_swap_host"
$GUEST_SWAP="/data_swap_guest"

docker run --name release_build -v ${HOST_SWAP}:${GUEST_SWAP} --entrypoint="/bin/bash" kratos/wheels /workspace/scripts/start.sh ${RELEASE_BRANCH}
docker rm release_build