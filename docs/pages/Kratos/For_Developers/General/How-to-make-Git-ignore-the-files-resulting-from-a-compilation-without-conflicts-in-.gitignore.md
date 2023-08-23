---
title: How to make Git ignore the files resulting from a compilation without conflicts in .gitignore
keywords: 
tags: [How-to-make-Git-ignore-the-files-resulting-from-a-compilation-without-conflicts-in-.gitignore.md]
sidebar: kratos_for_developers
summary: 
---

Go to your _Kratos_ folder:

```sh
cd kratos
``` 

And edit the file _.git/info/exclude_ (NOTE: You can use your favourite editor, not necessarily Emacs):

```sh
emacs .git/info/exclude
``` 

Add the following lines inside the file: 

```sh
KratosMultiphysics/
applications/python_scripts/ 
libpython2.7.so.1.0
runkratos
libs/
``` 

Note that some files will change according to your Python version. Feel free to add to this list every folder that you want to ignore.

Close the text editor and run the following command: 

```sh
git update-index
``` 

The folders and files that you added to the file .git/info/exclude should not appear when typing:

```sh
git status
```  