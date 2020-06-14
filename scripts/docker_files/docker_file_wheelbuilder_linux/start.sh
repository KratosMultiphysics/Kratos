#!/bin/bash

set -e

branch=${branch:-master}
pythons=${pythons:-35,36,37,38}
cpus=${cpus:-4}
repository=${repository:-https://github.com/KratosMultiphysics/Kratos.git}
cotire=${cotire:-OFF}

while [ $# -gt 0 ]; do

  if [[ $1 == *"--"* ]]; then
    param="${1/--/}"
    declare $param="$2"
  fi

  shift
done

cd /workspace/kratos
git clone --depth 1 --single-branch -b "$branch" "$repository"

cd /workspace/kratos/Kratos/scripts/wheels/linux/
chmod +x build.sh
./build.sh "$cpus" "$pythons" "$cotire"
