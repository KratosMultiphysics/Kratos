# mmg - Surface and volume remeshers
mmg is an open source software for bidimensional and tridimensional surface and volume remeshing.

It provides 3 applications and 4 libraries:
  * the **mmg2d** application and library: adaptation and optimization of a bidimensionnal triangulation;
  * the **mmgs** application and library: adaptation and optimization of a surface triangulation and isovalue discretization;
  * the **mmg3d** application and library: adaptation and optimization of a tetrahedral mesh and implicit domain meshing;
  * the **mmg** library gathering the **mmg2d**, **mmgs** and **mmg3d** libraries.

[//]: # ( comment )

## Get and compile the mmg project
  1. Get the repository:  
      ```Shell
      wget https://github.com/MmgTools/mmg/archive/master.zip
      ```
     or
      ```Shell
      git clone https://github.com/MmgTools/mmg.git
      ```

    The project sources are available under the **_src/_** directory, see:
      * **_src/mmg2d/_**   for files related to the mmg2d application;
      * **_src/mmgs/_**   for files related to the mmgs application;
      * **_src/mmg3d/_**  for files related to the mmg3d application;
      * **_src/common/_** for files related to the both.

  2. Fast compilation (build both **mmg2d**, **mmgs**, **mmg3d**, the mmg2d static library (**libmmg3d.a**), the mmgs static library (**libmmgs.a**), the mmg3d static library (**libmmg3d.a**) and the mmg static library (**libmmg.a**)):  
      ```Shell
      cd mmg  
      mkdir build  
      cd build  
      cmake ..  
      make  
      make install
      ```
    If the `make install` command fail, try to run the `sudo make install` command.
    If you don't have root access, please refers to the [Installation] section(https://github.com/MmgTools/Mmg/wiki/Setup-guide#iii-installation) of the [setup guide](https://github.com/MmgTools/Mmg/wiki/Setup-guide#setup-guide).

    The **mmg2d**, **mmgs** and **mmg3d** applications are available under the `mmg2d_O3`, `mmgs_O3` and `mmg3d_O3` commands. 

## Documentation
### Project's web page
Project's actualities and software tutorials can be found on the [mmgtools](http://www.mmgtools.org) web page.

### **Mmg**'s forum
Share your comments and issues with other members of the Mmg community on the [Mmg forum](https://forum.mmgtools.org/).

### GitHub's Wiki
More detailed informations about the compilation and configuration of the mmg's applications are available on the project [wiki](https://github.com/MmgTools/mmg/wiki).

### Man-pages
Man pages are available inside the **_doc/man_** directory:
  * To see the **mmg2d** man page, just tap `man ./doc/man/mmg2d.1.gz`
  * To see the **mmgs** man page, just tap `man ./doc/man/mmgs.1.gz`
  * To see the **mmg3d** man page, just tap `man ./doc/man/mmg3d.1.gz`

### Code documentation
Run the `make doc` command to build the Doxygen documentation.
  * To see the **mmg2d** documentation, open up the **_mmg/doc/mmg2d/html/index.html_** file;
  * To see the **mmgs** documentation, open up the **_mmg/doc/mmgs/html/index.html_** file;
  * To see the **mmg3d** documentation, open up the **_mmg/doc/mmg3d/html/index.html_** file.

## Platforms
The **mmg** applications are validated on OS X and on most of the Linux platforms. 

## Contributing
Your contributions to the **mmg** project are welcomed. You can help us to improve
our code by many means:
  * pull requests: please follow the [wiki's guideline](https://github.com/MmgTools/Mmg/wiki/Developers-wiki#pull-requests);
  * feature requests: please use the [Mmg forum](https://forum.mmgtools.org/);
  * bug reports: please use the [GitHub issue tracker](https://github.com/MmgTools/mmg/issues/new);

## About the team
mmg's current developers and maintainers are Charles Dapogny, Cécile Dobrzynski, Pascal Frey and Algiane Froehly.

Contact: contact@mmgtools.org

## License and copyright
Code is under the [terms of the GNU Lesser General Public License](https://raw.githubusercontent.com/MmgTools/mmg/master/LICENSE).

Copyright © Bx INP/Inria/UBordeaux/UPMC, 2004- .

## Reference
[Three-dimensional adaptive domain remeshing, implicit domain meshing, and applications to free and moving boundary problems - _C. Dapogny, C. Dobrzynski and P. Frey_ - April 1, 2014 - _JCP_](http://www.sciencedirect.com/science/article/pii/S0021999114000266)
