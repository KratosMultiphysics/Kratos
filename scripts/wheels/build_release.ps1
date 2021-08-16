RELEASE_BRANCH=Release-9.0

HOST_SWAP=C:/data_swap_host
GUEST_SWAP=/data_swap_guest

docker run --name release_build -v ${HOST_SWAP}:${GUEST_SWAP} -d kratos/wheels ${RELEASE_BRANCH}
docker rm release_build