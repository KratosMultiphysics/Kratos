lc = 0.4;
lc2=0.04;

Point(1) = {0, 0, 0, lc};
Point(2) = {1.62, 0, 0, lc};
Point(3) = {0, 4, 0, lc};
Point(4) = {1.62, 4, 0, lc};

Point(5) = {1.62, 1, 0, lc2};
Point(6) = {1.62, 3, 0, lc2};

Line(1)={1,2};
Line(2)={1,3};
Line(3)={3,4};
Line(5)={2,5};
Line(6)={5,6};
Line(7)={6,4};

//+
Curve Loop(1) = {2, 3, -7, -6, -5, -1};
//+
Plane Surface(1) = {1};
//+
Physical Curve("column_left_boundary", 8) = {2};
//+
Physical Curve("column_right_boundary", 9) = {7, 6, 5};
//+
Physical Curve("bottom", 10) = {1};
//+
Physical Surface("porous", 11) = {1};

Mesh.ElementOrder = 2;
Mesh 2;
