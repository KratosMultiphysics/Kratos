#!/bin/bash

if [ ! -d "../_site/Page_files" ]; then
    if [ -d ../_site ]; then
        rm -rf ../_site
    fi

    git clone --depth 1 --filter=blob:none --sparse https://github.com/KratosMultiphysics/Documentation.git ../_site
    cd ../_site
    git sparse-checkout init --cone
    git sparse-checkout set Page_files

    rm ../_site/*
fi

cp -r ../pages ../_site/Page_files
cp ../_config.yml ../_site/Page_files

cd ../_site/Page_files
jekyll build