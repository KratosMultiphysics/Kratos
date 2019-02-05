GiD Post Results File 1.0
GaussPoints "tri1_element_gp" ElemType Triangle
Number Of Gauss Points: 1
Natural Coordinates: Internal
End GaussPoints
GaussPoints "tet1_element_gp" ElemType Tetrahedra
Number Of Gauss Points: 1
Natural Coordinates: Internal
End GaussPoints
Result "VELOCITY" "Kratos" 0 Vector OnNodes
Values
1 0.707107 0 -0.707107
2 0.707107 0 0.707107
3 0.707107 0 0.707107
4 0.707107 0 -0.707107
5 -1 0 0
End Values
Result "VORTICITY" "Kratos" 0 Vector OnGaussPoints "tri1_element_gp"
Values
1 0 0 0
3 0 0 0
End Values
Result "VORTICITY" "Kratos" 0 Vector OnGaussPoints "tet1_element_gp"
Values
1 0 -3.12132 0
End Values
Result "NORMAL" "Kratos" 0 Vector OnGaussPoints "tri1_element_gp"
Values
1 -0.5 0 0.25
3 0.5 -0 0.25
End Values
Result "NORMAL" "Kratos" 0 Vector OnGaussPoints "tet1_element_gp"
Values
1 0 0 0
End Values
Result "ACTIVE" "Kratos" 0 Scalar OnGaussPoints "tri1_element_gp"
Values
1 0
2 0
3 0
4 0
End Values
Result "ACTIVE" "Kratos" 0 Scalar OnGaussPoints "tet1_element_gp"
Values
1 0
2 0
End Values
