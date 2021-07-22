RELEASE_BRANCH=Release-9.0

HOST_SWAP=C:/data_swap_host

docker run --name release_build -v ${HOST_SWAP}:/data_swap_guest -d kratos/wheels ${RELEASE_BRANCH}
docker rm release_build