template = '''---
title: {title}
keywords: {keywords}
tags: {tags}
sidebar: {sidebar}
summary: {summary}
---'''

with open('process_list') as pl:

    summary = ''
    sidebar = 'kratos_core_processes'
    keywords = 'process core'

    for l in pl:
        l = l[:-3]
        with open(f'{l.strip()}.md', 'w') as f:
            header = template.format(
                title=f'{" ".join(l.strip().split("_")[:-1]).title()}', 
                keywords=keywords, 
                tags=f'[{" ".join(l.strip().split("_"))}]', 
                sidebar=sidebar, 
                summary=summary
            )
            f.write(header)
            f.write('\n\n')
            f.write(f'# {" ".join(l.strip().split("_")[:-1]).title()}')
            f.write('\n\n')
            f.write(f'## Description')
            f.write('\n\n')
            f.write(f'## Parameters & Defaults')