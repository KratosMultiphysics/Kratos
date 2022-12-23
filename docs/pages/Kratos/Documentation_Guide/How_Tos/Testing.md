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

The whole kratos documentation can be hosted in your local machine so that you acn visualize the changes you made by using any browser of your choice.

1. goto the "docs" folder
2. python process_pages.py -t local
3. jekyll serve

The last command will host the webpage under the address "http://127.0.0.1:4000". So now you can visualize the documentation with your local changes by browsing the above address using any browser.

## Submitting content to repository

Please make sure you only have changes to the folder "docs/pages". All the other folders/files will be changed accordingly based on the content provided in the "docs/pages" folder. Then run the following command (without `-t local` option)

```bash
    python process_pages.py
```

This will generate the files which needs to be submitted to the github repository so the changes made can be visualized automatically.