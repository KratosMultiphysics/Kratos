---
title: Custom Entries
keywords: 
tags: [adding custom entries, adding external links, adding kratos examples]
sidebar: documentation_guide
summary: 
---


```custom_entries``` is not a mandatory entry in ```menu_info.json``` file. This is used order the menu entries when automating creation of the navigation bars. If this entry is present, then the entries in ```custom_entries``` will be added first to the navigation bar (files and folder names are processed in the given order). What ever folder or MarkDown file not found in this ```custom_entries``` list will be sorted separately in alphabatical order and then added after the final menu entry given by ```custom_entries```. This allows content creator to have ordered menu entries if required.

There are two types of ```custom_entries``` which are allowed to have.
* String entry
* Dictionary entry

## String entry

String entries represents linking the menu bar entry with a folder or a MarkDown content page. The path of the folder or MarkDown page file relative to the corresponding ```menu_info.json``` file should be provided.

## Dictionary entry

There are two types of dictionary entries.
* ```external```
* ```kratos_example```

### External
The external type of entry allows having navigation bar entries which links with an external link. These entries must have ```type = external```, ```title``` which corresponds to entry name, ```url``` which corresponds to the external link. Please start with ```https://``` or ```http://``` for ```url``` entry.

### Kratos example
The external type of entry allows having navigation bar entries which links with a kratos example found in [https://github.com/KratosMultiphysics/Examples](https://github.com/KratosMultiphysics/Examples). These entries must have ```type = kratos_example```, ```raw_url``` which corresponds to the raw url of the MarkDown file for the specific example, ```source_url``` which corresponds to the external folder for the specific example found in the github and ```file_name``` which is used to save downloaded content from ```raw_url```. The ```raw_url``` must start with ```https://raw.githubusercontent.com/```.

## Example

```json
    "custom_entries": [
        "General",
        "How_to_get_Kratos",
        "Getting_started",
        {
            "type": "external",
            "title": "Google",
            "url": "https://www.google.com"
        },
        {
            "type": "kratos_example",
            "raw_url": "https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/shape_optimization/use_cases/02_Strain_Energy_Minimization_3D_Shell/README.md",
            "source_url": "https://github.com/KratosMultiphysics/Examples/tree/master/shape_optimization/use_cases/02_Strain_Energy_Minimization_3D_Shell",
            "file_name": "Strain_Energy_Minimization_3D_Shell.md"
        }
    ]
```