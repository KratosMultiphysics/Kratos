---
title: Documentation testing and updating
keywords: 
tags: [Testing.md]
sidebar: documentation_guide
summary: 
---


This section explains how the written content can be tested in your local PC.

## Dependencies

Followings needs to be installed in your system to locally test the changes you made in the documentation.

1. jekyll
2. python

## Testing locally

The whole kratos documentation can be hosted in your local machine so that you can visualize the changes you made by using any browser of your choice.

1. goto the "docs" folder
2. python process_pages.py -t local
3. goto "scripts" folder in "docs".
4. execute "build_site.sh"
5. Then navidate to "docs/_site/Page_files"
6. jekyll serve

The last command will host the webpage under the address "http://127.0.0.1:4000". So now you can visualize the documentation with your local changes by browsing the above address using any browser.

The generated html files for the whole web page can be found in "docs/_site/Page_files/_site" folder.
