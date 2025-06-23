---
title: Top Navigation Bar
keywords: 
sidebar: documentation_guide
summary: 
tags: [modify top navigation bar]
---

This will eloborate on how a user/developer can create/modify top navigation bar entries.

## Creating navigation bar entries

Creation of the top navigation bar entries are automated by the **python** file **process_pages.py**. When creating the top navigation bar, it recursively iterates through the *docs/pages* folder maximum up to two sub-levels and creates top navigation bar entries. This automation can be controlled by editing the *docs/pages/menu_info.json*.

The section ```additional_menu_options``` has ```title``` entry which referes to the type of the navigation bar. This should not be changed in the *docs/pages/menu_info.json*.

A content creator is allowed to not to have the entry ```custom_entries```. In that case, the folders found in recursive search will be first sorted alphabatically and then will be put to the navigation bar. If the content creator wants to have custom ordering, then he/she can add the order of the folder entries in the ```custom_entries```. Each of the leaf folders in the recusive search for maximum upto two levels should have a ```menu_info.json``` file. As an example, all the following folders needs to have their own ```menu_info.json```.
* *docs/pages/Kratos/For_Users*
* *docs/pages/Kratos/For_Developers*
* *docs/pages/Kratos/Debugging*
* *docs/pages/Kratos/Documentation_Guide*

This is because each folder will be a leaf entry in the top navigation bar and they need to provide landing content page information when navigated to that entry via top navigation bar. [Modify Side Navigation Bar](Modify_side_nav.html)

The automation process will replace "_" with " " when creating menu entries from folder or MarkDown file names.

## Custom entries

Please refer to [Custom entries](Custom_entries.html)

## Example

```json
{
    "additional_menu_options": {
        "title": "Topnav dropdowns"
    },
    "custom_entries": [
        "Kratos",
        "Applications"
    ]
}
```






