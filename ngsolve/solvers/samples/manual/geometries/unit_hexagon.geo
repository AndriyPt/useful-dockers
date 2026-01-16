// Gmsh project created on Wed Jan 14 21:58:47 2026
//+
// CAUTION: Order of points in border should be counterclockwise
//+
Point(1) = {-1, 0, 0.0};
//+
Point(2) = {-0.5, -0.5, 0.0};
//+
Point(3) = {0.5, -0.5, 0.0};
//+
Point(4) = {1, 0, 0.0};
//+
Point(5) = {0.5, 0.5, 0.0};
//+
Point(6) = {-0.5, 0.5, 0.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 1};
//+
Curve Loop(1) = {6, 1, 2, 3, 4, 5};
//+
Plane Surface(1) = {1};
//+
Physical Curve("left", 7) = {6, 1};
//+
Physical Curve("bottom", 8) = {2};
//+
Physical Curve("right", 9) = {3, 4};
//+
Physical Curve("top", 10) = {5};
//+
Physical Surface("body", 11) = {1};
