// Gmsh project created on Thu Jan  1 19:39:42 2026
//+
Point(1) = {-0.5, -0.5, 0, 1.0};
//+
Point(2) = {0.5, -0.5, 0, 1.0};
//+
Point(3) = {0.5, 0.5, 0, 1.0};
//+
Point(4) = {-0.5, 0.5, 0, 1.0};
//+
Line(1) = {4, 3};
//+
Line(2) = {3, 2};
//+
Line(3) = {2, 1};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {1, 2, 3, -4};
//+
Plane Surface(1) = {1};
//+
Physical Curve("left", 5) = {4};
//+
Physical Curve("right", 6) = {2};
//+
Physical Curve("top", 7) = {1};
//+
Physical Curve("bottom", 8) = {3};
//+
Physical Surface("body", 9) = {1};
