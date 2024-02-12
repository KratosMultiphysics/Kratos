---
title: Mesh node ordering
keywords: 
tags: [Mesh Node Ordering]
sidebar: kratos_for_developers
summary: 
---

# Overview

The node ordering followed by *Kratos* is the very same considered by the pre/post-processor [*GiD*](www.gidhome.com). You can see the notation considered by them [here](http://www-opale.inrialpes.fr/Aerochina/info/en/html-version/gid_11.html).

# Node ordering 

Following the documentation structure used in [*GMSH*](http://gmsh.info//doc/texinfo/gmsh.html#MSH-ASCII-file-format) we represent the node ordering considered in *Kratos*.

Additionally everything related with geometry implementation can be directly consulted on the geometries [code](https://github.com/KratosMultiphysics/Kratos/tree/master/kratos/geometries).

## Line 

```html
Line2D2/                   Line2D3/   
Line3D2:                   Line3D3:   
                                                
0----------1 --> u      0-----2----1

```

## Triangle

```html
Triangle2D3/               Triangle2D6/     
Triangle3D3:               Triangle3D6:     

v                                                              
^                                                               
|                                                              
2                       2                    
|`\                     |`\              
|  `\                   |  `\           
|    `\                 5    `4           
|      `\               |      `\          
|        `\             |        `\          
0----------1 --> u      0-----3----1   
             
```

## Quadrilateral

```html
Quadrilateral2D4/      Quadrilateral2D8/       Quadrilateral2D9/
Quadrilateral3D4:      Quadrilateral3D8:       Quadrilateral3D9:

      v
      ^
      |
3-----------2          3-----6-----2           3-----6-----2 
|     |     |          |           |           |           | 
|     |     |          |           |           |           | 
|     +---- | --> u    7           5           7     8     5 
|           |          |           |           |           | 
|           |          |           |           |           | 
0-----------1          0-----4-----1           0-----4-----1 
```

## Tetrahedra

```html
Tetrahedra3D4:                          Tetrahedra3D10:

                   v
                 .
               ,/
              /
           2                                     3                              
         ,/|`\                                 ,/|`\                          
       ,/  |  `\                             ,/  |  `\       
     ,/    '.   `\                         ,7    '.   `9     
   ,/       |     `\                     ,/       8     `\   
 ,/         |       `\                 ,/         |       `\ 
0-----------'.--------1 --> u         0--------6--'.--------2
 `\.         |      ,/                 `\.         |      ,/ 
    `\.      |    ,/                      `\.      |    ,5   
       `\.   '. ,/                           `4.   '. ,/     
          `\. |/                                `\. |/       
             `3                                    `1        
                `\.
                   ` w
```

## Hexahedron

```html
Hexahedron3D8:         Hexahedron3D20:       Hexahedron3D27:

       v
3----------2            3----10----2           3----10----2     
|\     ^   |\           |\         |\          |\         |\    
| \    |   | \          | 15       | 14        |15    23  | 14  
|  \   |   |  \         9  \       11 \        9  \ 20   11  \  
|   7------+---6        |   7----18+---6       |   7----18+---6 
|   |  +-- |-- | -> u   |   |      |   |       |22 |  26  | 24| 
0---+---\--1   |        0---+-8----1   |       0---+-8----1   | 
 \  |    \  \  |         \  17      \  19       \ 17    25 \  19
  \ |     \  \ |         12 |        13|        12 |  21    13| 
   \|      w  \|           \|         \|          \|         \| 
    4----------5            4----16----5           4----16----5 
```

## Prism

```html
Prism3D6:                      Prism3D15:              

           w
           ^
           |
           3                       3                          
         ,/|`\                   ,/|`\                    
       ,/  |  `\               12  |  14               
     ,/    |    `\           ,/    |    `\          
    4------+------5         4------13-----5       
    |      |      |         |      9      |        
    |    ,/|`\    |         |      |      |        
    |  ,/  |  `\  |         |      |      |       
    |,/    |    `\|         |      |      |       
   ,|      |      |\        10     |      11       
 ,/ |      0      | `\      |      0      |        
u   |    ,/ `\    |    v    |    ,/ `\    |       
    |  ,/     `\  |         |  ,6     `8  |        
    |,/         `\|         |,/         `\|       
    1-------------2         1------7------2        

```


