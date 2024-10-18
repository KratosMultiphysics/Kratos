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

## Python snippets

Python snippets can be added to markdown as usual. In this case, if one requires to generate the output of the python snippet automatically and append after the python snippet, then following steps can be followed.
1. Add `## POST_PROCESS_PAGES_PYTHON_OUTPUT_GENERATION` anywhere in the python snippet. It should not have any leading or trailing spaces/tabs.
2. Run the local build in your computer [Make sure to initialize kratos environment or any other library environments in the terminal which is being used by the snippet]. This will run the snippet and capture the output. It will also remove the `## POST_PROCESS_PAGES_PYTHON_OUTPUT_GENERATION` tag line from the python snippet.
3. Now you will see some block after the python snippet which contains the python output.

This generation is only done if ```process_pages.py``` is passed with `-t local` flag.
Following is an example:
```python-example
print(1)
## POST_PROCESS_PAGES_PYTHON_OUTPUT_GENERATION
```

This will generate the following markdown
```
Expected output:
1
```


