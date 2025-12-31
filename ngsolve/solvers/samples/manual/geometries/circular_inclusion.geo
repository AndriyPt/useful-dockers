// Gmsh project created on Wed Dec 31 14:13:33 2025
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, -0.2, 0, 1.0};
//+
Point(3) = {0, 0.2, 0, 1.0};
//+
Circle(1) = {2, 1, 3};
//+
Circle(2) = {3, 1, 2};
//+
Curve Loop(1) = {2, 1};
//+
Plane Surface(1) = {1};
//+
Physical Surface("inclusion_domain", 3) = {1};
//+
Physical Curve("inclusion", 4) = {1, 2};
