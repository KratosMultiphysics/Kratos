$HOST_SWAP="C:/data_swap_host"
$GUEST_SWAP="C:/data_swap_guest"

docker run -it --name test_release -v ${HOST_SWAP}:${GUEST_SWAP} --cpus=24 --memory=64g --entrypoint="powershell" kratosmultiphysics/kratos-wheelbuilder-windows