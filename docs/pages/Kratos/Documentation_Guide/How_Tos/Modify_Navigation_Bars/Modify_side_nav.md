---
title: Side Navigation Bar
keywords: 
tags: [modify side navigation bar]
sidebar: documentation_guide
summary: 
---

This elaborates how to modify side navigation bar for a specific content page.

This will eloborate on how a user/developer can create/modify top navigation bar entries.

## Creating navigation side bar entries

Creation of the side navigation bar entries aCustom_entries.mdre automated by the **python** file **process_pages.py**. When creating the side navigation bar, it recursively iterates through the *docs/pages/*/*/* folder maximum up to three sub-levels and creates side navigation bar entries. This automation can be controlled by editing the *menu_info.json*.

The automation process will replace "_" with " " when creating menu entries from folder or MarkDown file names.

The entry ```side_bar_name``` corresponds to the side navigation bar file name. This is a **must** have entry in here and it needs to be **unique** so that it will not overwrite another side navigation bar file. **Do not create any names starting with ```remote``` since that is a keyword reserved to identify automated files in this documentation.**

The entry ```landing_page``` needs to be refering to a Markdown file. This is a **must** have entry. This content page will be shown when navigated to this entry.

The entry ```additional_menu_options``` is a **must** have entry. Under their, ```product``` entry should be added which reflects the side bar navigation title and ```title``` should be with value ```sidebar```.

**The entry ```custom_entries``` is optional. This is used in case if there is a need to have a specific order in the side bar menu entries. All the sub-folders, mark down files which are not in the ```custom_enties``` will be put in the side bar after in the alphabatical order after the entries which are found in the list ```custom_entries```.**

## Custom entries

Please refer to [Custom entries](Custom_entries.html)

## Example
```json
{
    "side_bar_name": "kratos_for_users",
    "landing_page": "General/Overview.md",
    "additional_menu_options": {
        "title": "sidebar",
        "product": "For Users"
    },
    "custom_entries": [
        "General",
        "How_to_get_Kratos",
        "Getting_started"
    ]
}
```