// -----------------------------
// 2D Human Bladder Geometry
// Anatomical approximation
// -----------------------------

SetFactory("OpenCASCADE");

// -----------------------------
// Bladder outline points
// (counter-clockwise)
// -----------------------------

Point(1)  = {  0,  40, 0};   // top
Point(2)  = { 25,  35, 0};
Point(3)  = { 45,  15, 0};
Point(4)  = { 50,   0, 0};
Point(5)  = { 45, -15, 0};
Point(6)  = { 25, -30, 0};
Point(7)  = { 10, -40, 0};
Point(8)  = {  5, -55, 0};  // bladder neck
Point(9)  = {  0, -65, 0};  // urethral outlet
Point(10) = { -5, -55, 0};
Point(11) = { -10, -40, 0};
Point(12) = { -25, -30, 0};
Point(13) = { -45, -15, 0};
Point(14) = { -50,   0, 0};
Point(15) = { -45,  15, 0};
Point(16) = { -25,  35, 0};

// -----------------------------
// Spline boundary
// -----------------------------

Spline(1) = {
 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1
};

// -----------------------------
// Surface definition
// -----------------------------

Curve Loop(1) = {1};
Plane Surface(1) = {1};

// -----------------------------
// Physical groups
// -----------------------------

Physical Curve("top") = {1};
Physical Surface("body") = {1};
