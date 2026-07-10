---
title: mdpa files
keywords: mpm mdpa
tags: [mpm mpdpa]
sidebar: mpm_application
summary: 
---

In Kratos, information about the geometry of the problem is contained in files having the `.mdpa` extension.

A detailed description of the `mdpa` file format is available [here](../../../Kratos/For_Users/Crash_Course/Input_Output_and_Visualization/Input_Data).

The `MPMApplication` requires two `mdpa` files, typically named `file_name_Grid.mdpa` and `file_name_Body.mdpa`, with the suffixes `_Grid` and `_Body` highlighting their distinct roles.

* `file_name_Grid.mdpa`: contains information about the nodes, elements and conditions of the **background grid**. This grid is used for computing the Finite Element solution and provides the domain where the material points move.
* `file_name_Body.mdpa`: defines the initial mesh discretisation of the continuum to be simulated by means of the Material Point Method. Indeed, in the `MPMApplication`, the body is first discretised with a mesh and then a certain number of material points (based on the `MATERIAL_POINTS_PER_ELEMENT` parameter defined in the [`ParticleMaterials.json` file](./json)) are placed within each element of the mesh. Finally, the mesh initially discretising the body and contained in the file `file_name_Body.mdpa` is disregarded, therefore serving only to define the initial position of the material points.

## `file_name_Grid.mdpa`

```
Begin ModelPartData
//  VARIABLE_NAME value
End ModelPartData

Begin Properties 0
End Properties

Begin Nodes
     1 -7.13647000000000000 -2.77405000000000000 0.0
     2 -7.136469999999999 0 -1.04399666666666650 0.0
     3 -5.148614285714286 0 -2.77405000000000000 0.0
     4 -5.404419192431123 0 -0.17896999999999963 0.0
     5 -7.13647000000000000 +0.68605666666666700 0.0
     6 -3.70035665111947360 -0.90015846792814620 0.0
     7 -3.16075857142857330 -2.77405000000000000 0.0
    11 -3.91619588299583340 +0.88743014490059560 0.0
    12 -7.13647000000000000 +2.41611000000000000 0.0
    13 -2.16683071428571640 -1.04199919243112230 0.0
    14 -5.148614285714286 0 +2.41611000000000000 0.0
    17 -1.17290285714285810 -2.77405000000000000 0.0
    18 -2.16683071428571640 +0.68405919243112300 0.0
    20 -3.16075857142857330 +2.41611000000000000 0.0
    21 -0.17897500000000166 -1.04199919243112230 0.0
    22 -0.17897500000000166 +0.68405919243112300 0.0
    23 -1.17290285714285810 +2.41611000000000000 0.0
    24 0.814952857142855300 -2.77405000000000000 0.0
    25 1.808880714285713000 -1.04199919243112320 0.0
    26 0.814952857142855300 +2.41611000000000000 0.0
    27 1.808880714285713000 +0.68405919243112300 0.0
    30 2.802808571428570000 -2.77405000000000000 0.0
    32 3.317802552771939000 -0.87099133333333360 0.0
    34 2.802808571428570000 +2.41611000000000000 0.0
    35 3.649111838486223000 +0.85906200000000020 0.0
    38 4.790664285714286000 -2.77405000000000000 0.0
    39 5.046469192431123000 -0.17896999999999963 0.0
    40 4.790664285714286000 +2.41611000000000000 0.0
    41 6.778520000000000000 -2.77405000000000000 0.0
    42 6.778520000000000000 -1.04399666666666650 0.0
    43 6.778520000000001000 +0.68605666666666700 0.0
    44 6.778520000000000000 +2.41611000000000000 0.0
End Nodes

Begin Elements Element2D3N
        1          0    12     5    14
        2          0     2     1     3
        3          0    42    43    39
        4          0     5     2     4
        5          0    41    42    38
        6          0    43    44    40
        7          0    39    43    40
        8          0    42    39    38
        9          0     4     2     3
       10          0     5     4    14
       11          0    23    20    18
       12          0     3     7     6
       13          0    30    38    32
       14          0    34    26    27
       15          0    17    24    21
       16          0    40    34    35
       17          0    20    14    11
       18          0     7    17    13
       19          0    26    23    22
       20          0    24    30    25
       21          0    18    20    11
       22          0    18    11     6
       23          0     6    11     4
       24          0     4    11    14
       25          0     6     4     3
       26          0    21    24    25
       27          0    21    25    27
       28          0    27    25    32
       29          0    32    25    30
       30          0    22    23    18
       31          0    22    18    13
       32          0    13    18     6
       33          0    13     6     7
       34          0    32    38    39
       35          0    32    39    35
       36          0    35    39    40
       37          0    32    35    27
       38          0    27    35    34
       39          0    27    26    22
       40          0    27    22    21
       41          0    21    22    13
       42          0    21    13    17
End Elements

Begin SubModelPart Grid
    Begin SubModelPartNodes
            1
            2
            3
            4
            5
            6
            7
           11
           12
           13
           14
           17
           18
           20
           21
           22
           23
           24
           25
           26
           27
           30
           32
           34
           35
           38
           39
           40
           41
           42
           43
           44
    End SubModelPartNodes
    Begin SubModelPartElements
            1
            2
            3
            4
            5
            6
            7
            8
            9
           10
           11
           12
           13
           14
           15
           16
           17
           18
           19
           20
           21
           22
           23
           24
           25
           26
           27
           28
           29
           30
           31
           32
           33
           34
           35
           36
           37
           38
           39
           40
           41
           42
    End SubModelPartElements
End SubModelPart

Begin SubModelPart DISPLACEMENT_BCs
    Begin SubModelPartNodes
            1
            3
            7
           17
           24
           30
           38
           41
    End SubModelPartNodes
End SubModelPart
```

## `file_name_Body.mdpa`

```
Begin ModelPartData
//  VARIABLE_NAME value
End ModelPartData

Begin Properties 0
End Properties

Begin Nodes
     8 -3.8769883038607390 -0.42505600000022326 0.0
     9 -2.9082800012962800 -0.98434000000044710 0.0
    10 -3.8769883032125985 +0.69351200112283560 0.0
    15 -1.9395716961392613 -0.42505600000022390 0.0
    16 -2.9082800000000000 +1.25279600000044700 0.0
    19 -1.9395716954911204 +0.69351199887761120 0.0
    28 2.34508276758495500 -0.94395950357364110 0.0
    29 2.34508276818629870 +0.09384750461520013 0.0
    31 3.24384999879731100 -1.46286300714728260 0.0
    33 3.24385000000000000 +0.61275100714728270 0.0
    36 4.14261723241504500 -0.94395950357364180 0.0
    37 4.14261723301638900 +0.09384750253208292 0.0
End Nodes

Begin Elements MPMUpdatedLagrangian2D3N // Body1
       43          0    16    10     9
       44          0     8     9    10
       45          0    15    19    16
       46          0     9    15    16
End Elements

Begin Elements MPMUpdatedLagrangian2D3N // Body2
       47          0    33    29    31
       48          0    28    31    29
       49          0    36    37    33
       50          0    31    36    33
End Elements

Begin SubModelPart Body1
    Begin SubModelPartNodes
            8
            9
           10
           15
           16
           19
    End SubModelPartNodes
    Begin SubModelPartElements
           43
           44
           45
           46
    End SubModelPartElements
End SubModelPart

Begin SubModelPart Body2
    Begin SubModelPartNodes
           28
           29
           31
           33
           36
           37
    End SubModelPartNodes
    Begin SubModelPartElements
           47
           48
           49
           50
    End SubModelPartElements
End SubModelPart
```
