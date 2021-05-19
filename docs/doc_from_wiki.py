import os
import re
import sys

def GenerateHeader(title, keywords, tags, sidebar, summary):
    buff = []
    
    buff.append("---\n")
    
    buff.append("title: {}\n".format(title))
    buff.append("keywords: {}\n".format(keywords))
    buff.append("tags: [{}]\n".format(tags))
    buff.append("sidebar: {}\n".format(sidebar))
    buff.append("summary: {}\n".format(summary))

    buff.append("---\n\n")

    return buff

def GenerateUrlsMap(wiki_file):
    urls_map = {}

    for l in wiki_file:
        map_url = re.match("(\[.+\]):(.+)", l)
        if map_url:
            map_key = map_url.group(1)[1:-1]
            map_url = map_url.group(2)[1:]

            urls_map[map_key] = map_url

    wiki_file.seek(0)

    return urls_map

def GenerateSideBar(wiki_file, urls_map):

    tk = "  "       # Tag Token
    nest_level = 0

    buff = []

    buff.append("entries:\n")
    buff.append("- title: sidebar\n")
    buff.append(tk + "product: Kratos Multiphysics\n")
    buff.append(tk + "version: 8.0\n")
    buff.append(tk + "folders:\n")

    for l in wiki_file:
        if l.startswith("##"):
            buff.append("\n")
            buff.append(1 * tk + "- title: {}\n".format(l.replace('## ','').replace("\n","")))
            buff.append(2 * tk + "output: web\n")
            buff.append(2 * tk + "folderitems:\n")

        if l.startswith("**"):
            buff.append("\n")
            buff.append(1 * tk + "- title: {}\n".format(l.replace('**','')))
            buff.append(2 * tk + "output: web\n")
            buff.append(2 * tk + "folderitems:\n")

        if l.startswith("* "):
            page_url = re.match(".*(\[.+\])(\(.+\)|\[.+\]).*", l)
            if page_url:
                page_name = page_url.group(1)[1:-1]
                page_name = page_name.replace(":", "")
                page_url = page_url.group(2)[1:-1]
                if page_url in urls_map:
                    page_url = urls_map[page_url]
                if page_url.startswith("http"):
                    page_url = page_url.split("/")[-1]
                
            if page_url:
                page_url = "/pages/" + page_url + ".html"

                buff.append(2 * tk + "- title: {}\n".format(page_name))
                buff.append(3 * tk + "url: {}\n".format(page_url))
                buff.append(3 * tk + "output: web\n")

        if l.startswith(" "):
            nest_level += 1
            # Count the number

    return buff

src_path = "wiki/"
dst_path = "pages/"

toc = ""

# Generate Doc Pages
for wiki_page in os.listdir(src_path):

    with open(os.path.join(src_path, wiki_page), "r") as wiki_file, open(os.path.join(dst_path, wiki_page), "w+") as docu_file:
        
        page_title = os.path.splitext(wiki_page)[0]
        
        page_title = page_title.replace("Tutorial:", "")         # Remove sections
        page_title = page_title.replace("How-to-use-", "")       # Remove sections
        page_title = page_title.replace("-"," ")                 # Remove odd symbols

        header = GenerateHeader(
            title=page_title,
            keywords="",
            tags=wiki_page,
            sidebar="kratos_sidebar",
            summary=""
        )

        docu_file.writelines(header)
        docu_file.writelines([l for l in wiki_file])

# Generate Side Menu
with open(os.path.join(src_path, "_Sidebar.md"), "r") as wiki_file, open(os.path.join("_data", "sidebars", "kratos_sidebar.yml"), "w+") as docu_file:
    url_map = GenerateUrlsMap(wiki_file)
    sidebar = GenerateSideBar(wiki_file, url_map)

    docu_file.writelines(sidebar)
