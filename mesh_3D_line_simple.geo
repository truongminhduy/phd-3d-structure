Include "parameters_gmsh_getdp.dat";
SetFactory("OpenCASCADE");

mid_lc    = lambda_m*2;
source_lc = lambda_m*0.5;
scat_lc   = lambda_m*1;

R_in  = D_pml_in/2; 
R_dif = pml_size;
R_out = D_pml_in/2+ pml_size;

Point(1) = {0,0,0};

Box(2) = {-R_in ,-R_in ,-R_in , R_in *2, R_in *2, R_in *2};
bnd_in()  = Boundary{Volume{2};};


Sphere(60) = {0,0,0,r_ellipse_2}; // scatterer

Point(90) = {xS         , yS+r_source, zS         };
Point(91) = {xS+r_source, yS         , zS         };
Point(92) = {xS+r_source, yS+r_source, zS         };
Point(93) = {xS         , yS         , zS         };
Point(94) = {xS         , yS         , zS+r_source};
Point(95) = {xS+r_source, yS         , zS+r_source};
Point(96) = {xS         , yS+r_source, zS+r_source};
Point(97) = {xS+r_source, yS+r_source, zS+r_source};
Line(637) = {93, 94};
Line(638) = {94, 95};
Line(639) = {95, 91};
Line(640) = {91, 93};
Line(641) = {93, 90};
Line(642) = {90, 92};
Line(643) = {92, 97};
Line(644) = {96, 97};
Line(645) = {96, 90};
Line(646) = {96, 94};
Line(647) = {97, 95};
Line(648) = {92, 91};
Line Loop(649) = {646, 638, -647, -644};
Plane Surface(650) = {649};
Line Loop(651) = {644, -643, -642, -645};
Plane Surface(652) = {651};
Line Loop(653) = {646, -637, 641, -645};
Plane Surface(654) = {653};
Line Loop(655) = {647, 639, -648, 643};
Plane Surface(656) = {655};
Line Loop(657) = {642, 648, 640, 641};
Plane Surface(658) = {657};
Line Loop(659) = {638, 639, 640, 637};
Plane Surface(660) = {659};
Surface Loop(661) = {650, 654, 660, 656, 658, 652}; // source
Volume(61) = {661};

// BooleanDifference(50) = { Volume{1}; Delete;}{ Volume{2};};
BooleanDifference(51) = { Volume{2}; Delete;}{ Volume{60,61};};

Physical Volume(1007) = {51};	// Vacuum no source
Physical Volume(1008) = {60};	// Object
Physical Volume(1009) = {61};	// source
Physical Point(1030) = {1};		// Printpoint
Physical Surface(111) = {bnd_in()}; // boundary in

Physical Line(1011)  = {637};	// source line z
Physical Line(1012)  = {641}; 	// source edge y
Physical Line(1013)  = {637}; 	// source edge z
Physical Point(1014) = {93};  	// source point

simple() = BooleanFragments{ Volume{:}; Delete; }{ };

Characteristic Length{ PointsOf{ Volume{:}; } } = mid_lc;
// Characteristic Length{ PointsOf{ Volume{50}; } } = mid_lc;
Characteristic Length{ PointsOf{ Volume{60}; } } = scat_lc;
Characteristic Length{ PointsOf{ Volume{61}; } } = source_lc;

