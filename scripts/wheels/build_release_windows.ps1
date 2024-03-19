$RELEASE_BRANCH="Release-9.4.6"

$HOST_SWAP="C:/data_swap_host"
$GUEST_SWAP="C:/data_swap_guest"

docker run --name release_build -v ${HOST_SWAP}:${GUEST_SWAP} --cpus=24 --memory=120g  --dns 1.1.1.1 --entrypoint="powershell.exe" kratosmultiphysics/kratos-wheelbuilder-windows "c:\\scripts\\start.ps1" ${RELEASE_BRANCH}
docker rm release_build