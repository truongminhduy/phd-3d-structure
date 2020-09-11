Include "parameters_gmsh_getdp.dat";
SetFactory("OpenCASCADE");

// mid_lc    = lambda_m*0.3;
// PML_lc    = lambda_m*0.7;
// source_lc = lambda_m*0.1;
// scat_lc   = lambda_m*0.15;

PML_lc    = lambda_m*3.5;
mid_lc    = lambda_m*2;
source_lc = lambda_m*0.5;
scat_lc   = lambda_m*1;


R_in  = D_pml_in/2; 
R_dif = pml_size;
R_out = D_pml_in/2+ pml_size;

Point(1) = {0,0,0};
Box(1) = {-R_out,-R_out,-R_out, R_out*2, R_out*2, R_out*2};
Box(2) = {-R_in ,-R_in ,-R_in , R_in *2, R_in *2, R_in *2};
bnd_in()  = Boundary{Volume{2};};
bnd_out() = Boundary{Volume{1};};

// Characteristic Length{ PointsOf{ Volume{1}; } }  = mid_lc*1.2;//PML_lc;

Box(3) = {-R_out,-R_in ,-R_in , R_dif, 2*R_in , 2*R_in };
Box(4) = { R_in ,-R_in ,-R_in , R_dif, 2*R_in , 2*R_in };
Box(5) = {-R_in ,-R_out,-R_in , 2*R_in , R_dif, 2*R_in };
Box(6) = {-R_in , R_in ,-R_in , 2*R_in , R_dif, 2*R_in };
Box(7) = {-R_in ,-R_in ,-R_out, 2*R_in , 2*R_in , R_dif};
Box(8) = {-R_in ,-R_in , R_in , 2*R_in , 2*R_in , R_dif};

// edge box
Box(11) = {-R_in ,-R_out,-R_out, 2*R_in , R_dif, R_dif};
Box(12) = {-R_in , R_in ,-R_out, 2*R_in , R_dif, R_dif};
Box(13) = {-R_in ,-R_out, R_in , 2*R_in , R_dif, R_dif};
Box(14) = {-R_in , R_in , R_in , 2*R_in , R_dif, R_dif};

Box(15) = {-R_out,-R_in ,-R_out, R_dif, 2*R_in, R_dif };
Box(16) = { R_in ,-R_in ,-R_out, R_dif, 2*R_in, R_dif };
Box(17) = {-R_out,-R_in , R_in , R_dif, 2*R_in, R_dif };
Box(18) = { R_in ,-R_in , R_in , R_dif, 2*R_in, R_dif };

Box(19) = {-R_out,-R_out,-R_in , R_dif, R_dif, 2*R_in };
Box(20) = { R_in ,-R_out,-R_in , R_dif, R_dif, 2*R_in };
Box(21) = {-R_out, R_in ,-R_in , R_dif, R_dif, 2*R_in };
Box(22) = { R_in , R_in ,-R_in , R_dif, R_dif, 2*R_in };

// corner box
v() = BooleanDifference{ Volume{1}; }{ Volume{2:8,11:22} ;};


// sphere
Sphere(60) = {0,0,0,r_ellipse_2}; // scatterer
// Box(60) = {-r_ellipse_2 ,-r_ellipse_2 ,-r_ellipse_2 , r_ellipse_2 *2, r_ellipse_2 *2, r_ellipse_2 *2};

Point(290) = {xS         , yS+r_source, zS         };
Point(291) = {xS+r_source, yS         , zS         };
Point(292) = {xS+r_source, yS+r_source, zS         };
Point(293) = {xS         , yS         , zS         };
Point(294) = {xS         , yS         , zS+r_source};
Point(295) = {xS+r_source, yS         , zS+r_source};
Point(296) = {xS         , yS+r_source, zS+r_source};
Point(297) = {xS+r_source, yS+r_source, zS+r_source};
Line(637) = {293, 294};
Line(638) = {294, 295};
Line(639) = {295, 291};
Line(640) = {291, 293};
Line(641) = {293, 290};
Line(642) = {290, 292};
Line(643) = {292, 297};
Line(644) = {296, 297};
Line(645) = {296, 290};
Line(646) = {296, 294};
Line(647) = {297, 295};
Line(648) = {292, 291};
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

BooleanDifference(50) = { Volume{1}; Delete;}{ Volume{2};};
BooleanDifference(51) = { Volume{2}; Delete;}{ Volume{60,61};};
// BooleanDifference(51) = { Volume{2}; Delete;}{ Volume{60};};


Physical Volume(1000) = {v()};   // PMLxyz
Physical Volume(1001) = {15:18}; // PMLxz
Physical Volume(1002) = {11:14}; // PMLyz
Physical Volume(1003) = {19:22}; // PMLxy
Physical Volume(1004) = {7,8};   // PMLz
Physical Volume(1005) = {5,6};   // PMLy
Physical Volume(1006) = {3,4};   // PMLx
Physical Volume(1007) = {51};	// Vacuum no source
Physical Volume(1008) = {60};	// Object
Physical Volume(1009) = {61};	// source
Physical Point(1030) = {1};		// Printpoint
Physical Surface(111) = {bnd_in()}; // boundary in
Physical Surface(112) = {bnd_out()}; // boundary out

Physical Line(1011)  = {637};	// source line z
Physical Line(1012)  = {641}; 	// source edge y
Physical Line(1013)  = {637}; 	// source edge z
Physical Point(1014) = {293};  	// source point

simple() = BooleanFragments{ Volume{:}; Delete; }{ };

Characteristic Length{ PointsOf{ Volume{:}; } }  = PML_lc;
Characteristic Length{ PointsOf{ Volume{51}; } } = mid_lc;
Characteristic Length{ PointsOf{ Volume{60}; } } = scat_lc;
Characteristic Length{ PointsOf{ Volume{61}; } } = source_lc;

