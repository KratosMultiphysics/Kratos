---
title: Documentation Pages
keywords: 
tags: [documentation, adding latex, adding images, page structure]
sidebar: documentation_guide
summary: 
---

## Page structure

The documentation pages are written in MarkDown language. The initial section of the *.md* file should contain following block.
```
---
title: Documentation Pages
keywords:
tags: [Pages.md]
sidebar: documentation_guide
summary:
---
```

The ```title``` entry refers to the title of the page. This will be used as the title of the page as well as the side navigation bar entry title.

The ```tags``` entry are used in the search to identify content in the corresponding page. So having tags corresponding to that page will help documentation user to find specific content easily with the help of the search provided in the documentation.

The ```sidebar``` entry will always be replaced with the correct side bar for the content page depending on where the conten MarkDown file is located.

If no initial section is found in a given *.md* file, then it will be automatically added to that specific file. In there, ```title``` will be the name of the file name where "_" are replaced with " ", and ```tags``` will be the name of the MarkDown file.

## Latex support

MarkDown language is used to create documentation guide. This is enabled with latex support to illustrate equations if necessary. Following is an example latex code.
```latex
<p align="center">$$ \underline{r} = \sum_{i=1}^N{\underline{x}_i}$$</p>
```

This is showed as following in the respective content page.
<p align="center">$$ \underline{r} = \sum_{i=1}^N{\underline{x}_i}$$</p>

## Adding images

It is allowed to create subfolders to store page data such as ```images```. If they do not have any MarkDown files within the given folder, then those folders will **not be navigated to find side bar navigation entries.** Images can be added by using the following MarkDown code block.
```markdown
<p align="center">
    <img src="images/vertex_morphing_filtering.png" alt="Vertex morphing filtering"/>
</p>
<p align="center">Figure 2: Vertex morphing filtering</p>
```


