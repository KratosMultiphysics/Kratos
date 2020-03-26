# Guide to the documentation

This website is a collection of all MarkDown (MD) files found in the CoSimulationApplication folder and its subfolders.
These files are converted to a static website using [MkDocs][1], based on the outline defined in the `mkdocs.yml` file.



## Documentation in MarkDown 

MD files are easy to write and read in raw form, and can easily be converted to HTML for nicer rendering. It is easy to learn, see e.g. this [tutorial](https://www.markdowntutorial.com/). Several flavors of MD exist, which differ only in some details. GitHub has its own [GitHub Flavored Markdown (GFM)](https://github.github.com/gfm/), which is used to render MD documents found in repositories on the run. PyCharm comes with an automatic rendering tool built in, and I assume many other IDE's too. 

The MD documents for the CoCoNuT documentation must be written in the [Python-Markdown](https://python-markdown.github.io/) flavor, because they are processed by [MkDocs][1]. This flavor stays very close to the [reference implementation](reference implementation) of MD. Pay attention to lists: they should be indented with 4 spaces (contrary to GFM). 

Some important rules for writing MD files for this documentation website:

*   The name of the file must be unique: use names like `fluent.md` or `test_examples.md`, not `readme.md`.
*   Links to other MD files: I have yet to test this... probably only works on website, not on GitHub because paths change. Possible fix: automatically replace all names `example.md` with a proper link. 

A note on math in MD. There is no out-of-the-box LaTeX compatibiliy sadly, but basic equations can be typeset directly in HTML if necessary. Some usefull symbols:

*   lowercase Greek: αβγδεζηϑθικλμνξοπρστυφϕχψω
*   uppercase Greek: ΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩ
*   super and subscripts: `<sup></sup>` and `<sub></sub>`
*   operators: + - − · / × √ ∘ ∗
*   other symbols: ∂ Δ	∑ ≤ ≥ ∈
*   [more math symbols](http://www.unics.uni-hannover.de/nhtcapri/mathematics.html)

It is also possible to [render LaTeX equations on-the-run](https://alexanderrodin.com/github-latex-markdown/) and importing them as images, but that is quite cumbersome.




## Creating a static website with MkDocs

[MkDocs][1] can be installed using pip: 
```bash
pip install mkdocs
```
Furthermore, the [Material](https://squidfunk.github.io/mkdocs-material/) theme must be installed:
```bash
pip install mkdocs-material
```
Both can be installed locally by adding the argument `-t <install_dir>`. 

The structure/outline of the website is dictated by the `nav` variable in `mkdocs.yml`. This is the only variable that must be adjusted when new MD files are added to the code.

The complete process to create the documentation website is automated by `run_mkdocs.py`. This does the following things:

*   Copy all MD files in the CoSimulationApplication folder and its subfolders to folder `docs`. 
*   Check if there are duplicate filenames: these overwrite each other! If duplicates are found, a warning is given with the original paths.
*   Check if each MD file is mentioned in `mkdocs.yml`. If a file is not mentioned, a warning is given.
*   Build static HTML website using `mkdocs build`. 
*   Deploy website on GitHub Pages using `mkdocs gh-deploy`.

To activate the website, the user must have administrator privileges in the repository. 


[1]:    https://www.mkdocs.org/
