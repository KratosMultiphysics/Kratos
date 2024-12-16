---
title: Application Documentation
keywords: 
tags: [creating application documentation]
sidebar: documentation_guide
summary: 
---

Each application in *Kratos* can have its own documentation with its own side navigation bar.

## How to create application documentation guide

1. Create a folder with the application name in **docs/pages/Applications/\<YOUR_APPLICATION_NAME\>**
2. Create a file named ```menu_info.json``` in that folder.
3. Now create the folder structure and MarkDown file structure as preferred. [Documentation Pages](Documentation_Pages.html)
4. Add the following content to the ```menu_info.json``` file.

```json
{
    "side_bar_name": "shape_optimization_application",
    "landing_page": "General/Multi_objective_optimization.md",
    "additional_menu_options": {
        "title": "sidebar",
        "product": "Shape Optimization Application"
    },
    "custom_entries": [
        "General",
        "Technologies",
        "Examples"
    ]
}
```

Please refer to [Side Navigation Bar](../Modify_Navigation_Bars/Modify_side_nav.html) for more information on how to modify the side navigation bar menu entries.

## How to add external links to navigation bars

Please refer to [External entries](../Modify_Navigation_Bars/Custom_entries.html#external)

## How to add kratos examples

Please refer to [Kratos Example entries](../Modify_Navigation_Bars/Custom_entries.html#kratos-example)