// Gmsh project created on Wed Mar 15 17:31:43 2023
SetFactory("OpenCASCADE");
//+
Sphere(1) = {0, 0, 0, 1, -Pi/2, Pi/2, 2*Pi};
//+
Transfinite Curve {2} = 30 Using Progression 0.95;
//+
Transfinite Curve {2} = 15 Using Progression 0.95;
//+
Transfinite Curve {2} = 15 Using Progression 0.7;
//+
Transfinite Curve {2, 2} = 10 Using Progression 0.8;
//+
Transfinite Curve {2} = 10 Using Progression 0.9;
