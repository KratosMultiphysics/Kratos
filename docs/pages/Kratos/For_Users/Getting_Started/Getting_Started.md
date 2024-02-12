---
title: Getting Started
keywords: 
tags: [Python Script Tutorial Getting Started]
sidebar: kratos_for_users
summary: 
---

In order to write a python script with *Kratos* the first step would be to *"Obtain"* the `KratosMultiphysics` module and then *"add"* it as an external module to *Python*. 

There are sveral ways to Obtain the module:

- Download the latest release of *Kratos* from Kratos' [release page](https://github.com/KratosMultiphysics/Kratos/releases)
- Installing from PIP 'KratosMultiphysics-all'
- Getting the *Kratos* and its **GUI** from *GiD* internet retrieve.
- Compile the source from the GitHub repo.

Remember that if you compile Kratos you will have to add the module path to the `PYTHONPATH` environment variable. 

## Linux
In Linux you can set the python path to the Kratos root directory using `export` command in terminal:

```console
$ export PYTHONPATH=/path/to/Kratos
```
{: data-lang="Console"}

### Example:
For a release version extracted your home directory

```console
$ export PYTHONPATH=$HOME/Kratos
```
{: data-lang="Console"}

In case of having a GiD installed in `GiDx64/GiD14.0.0` subfolder of home: 
```console
$ export PYTHONPATH=$HOME/GiDx64/GiD14.0.0/problemtypes/kratos.gid/exec/Kratos
```
{: data-lang="Console"}

## Windows

To create the path on *Windows 10* and having installed *Kratos* through *GiD* or either downloaded manually, follow these steps:

- In the search bar type `system` and select `System (Control Panel)`
- In the left menu, click the `Advanced system settings` link.
- Click Environment Variables Button at the bottom. In the section System Variables, find the `PYTHONPATH` and `PYTHONHOME` environment variables and select them. Edit both of them. If the `PYTHONHOME` or `PYTHONPATH` environment variable do not exist, click New and create them.

   - If you have installed *Kratos* using *GiD*:
Go to the *GiD* folder. For example, it could be in `C:\Program Files\GiD`
Inside the *GiD* folder, you should follow `\problemtypes\kratos.gid\exec\Kratos`. There are many subfolders here and can be confusing. The correct one is just when you can see the *KratosMultiphysics* folder.
Copy the current path and paste it as variable name.

   - If you have downloaded yourself or you have compiled *Kratos* in your own computer, copy the *Kratos* path, for example `C:\Kratos`

- Click `Accept`.

**Next** [Hello Kratos](Hello_World)