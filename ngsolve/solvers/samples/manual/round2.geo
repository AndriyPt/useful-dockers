// Gmsh project created on Tue Dec 30 16:25:35 2025
//+
Point(1) = {0.5, 0.5, 0, 1.0};
//+
Point(2) = {0.5, 0.3, 0, 1.0};
//+
Point(3) = {0.5, 0.7, 0, 1.0};
//+
Circle(1) = {2, 1, 3};
//+
Circle(2) = {3, 1, 2};
//+
Curve Loop(1) = {2, 1};
//+
Plane Surface(1) = {1};
//+
Physical Curve("inclusion", 3) = {1, 2};
//+
Physical Surface("inclusion_domain", 4) = {1};
