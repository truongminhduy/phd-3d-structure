//////////////////////////////////////////////////////////////////////////////
Include "parameters_gmsh_getdp.dat";
//////////////////////////////////////////////////////////////////////////////
Group {
  If (PML_or_not == 1)
    // with PML
    PMLxyz    = Region[1000];
    PMLxz     = Region[1001];
    PMLyz     = Region[1002];
    PMLxy     = Region[1003];
    PMLz      = Region[1004];
    PMLy      = Region[1005];
    PMLx      = Region[1006];
    vac_nosrc = Region[1007];
    scat      = Region[1008];
    source    = Region[1009];
    source_Lx    = Region[1011];
    source_Ly    = Region[1012];
    source_Lz    = Region[1013];
    source_point = Region[1014];
    source_line  = Region[{source_Lz}];
    SurfDirichlet = Region[112];
    Point_print = Region[1030];
    mid         = Region[{source,vac_nosrc}];
    PMLs        = Region[{PMLxyz,PMLxy,PMLxz,PMLyz,PMLx,PMLy,PMLz}]; 
    Omega_nopml = Region[{scat,mid}];
    Omega_out   = Region[{mid,PMLs}];
    Omega_in    = Region[{scat}];
    Omega       = Region[{Omega_in,Omega_out}];
    Omega_plot  = Region[{Omega,-source}];
  Else
    // with out PML
    source_Lx    = Region[1011];
    source_Ly    = Region[1012];
    source_Lz    = Region[1013];
    source_point = Region[1014];
    source_line  = Region[{source_Lz}];
    Point_print = Region[1030];
    vac_nosrc   = Region[1007];
    scat        = Region[1008];
    source      = Region[1009];
    SurfDirichlet = Region[111];
    mid         = Region[{source,vac_nosrc}];
    Omega_nopml = Region[{scat,mid}];
    Omega_out   = Region[{mid}];
    Omega_in    = Region[{scat}];
    Omega       = Region[{Omega_in,Omega_out}];
    Omega_plot  = Region[{Omega,-source}];

    // Point_print = Region[1030];
    // vac_nosrc   = Region[1007];
    // scat        = Region[1008];
    // SurfDirichlet = Region[111];
    // mid         = Region[{vac_nosrc}];
    // Omega_nopml = Region[{scat,mid}];
    // Omega_out   = Region[{mid}];
    // Omega_in    = Region[{scat}];
    // Omega       = Region[{Omega_in,Omega_out}];
    // Omega_plot  = Region[{Omega}];
  EndIf
}
///////////////////////////////////////////////////
Function{
  If (PML_or_not == 1)
    // a_pml = Cos[angle];
    // b_pml = Sin[angle];
    a_pml = 1;
    b_pml = PML;
    sx[Omega_nopml] = 1.;
    sy[Omega_nopml] = 1.;
    sz[Omega_nopml] = 1.;
    sx[PMLxyz] = Complex[a_pml,b_pml];
    sy[PMLxyz] = Complex[a_pml,b_pml];
    sz[PMLxyz] = Complex[a_pml,b_pml];
    sx[PMLxz]  = Complex[a_pml,b_pml];
    sy[PMLxz]  = 1.0;
    sz[PMLxz]  = Complex[a_pml,b_pml];
    sx[PMLyz]  = 1.0;
    sy[PMLyz]  = Complex[a_pml,b_pml];
    sz[PMLyz]  = Complex[a_pml,b_pml];
    sx[PMLxy]  = Complex[a_pml,b_pml];
    sy[PMLxy]  = Complex[a_pml,b_pml];
    sz[PMLxy]  = 1.0;
    sx[PMLx]   = Complex[a_pml,b_pml];
    sy[PMLx]   = 1.0;
    sz[PMLx]   = 1.0;
    sx[PMLy]   = 1.0;
    sy[PMLy]   = Complex[a_pml,b_pml];
    sz[PMLy]   = 1.0;
    sx[PMLz]   = 1.0;
    sy[PMLz]   = 1.0;
    sz[PMLz]   = Complex[a_pml,b_pml];
    pmltens[]  = TensorDiag[ sy[]*sz[]/sx[] , sz[]*sx[]/sy[] , sx[]*sy[]/sz[] ];
    epsilon[PMLs] = Complex[eps_mil_re,eps_mil_im]   * pmltens[] ;  
    epsilon1[PMLs] = Complex[eps_mil_re,eps_mil_im] * pmltens[] ;
    mur[PMLs] = pmltens[] ;  
  EndIf
 
  // epsilon[scat] = Complex[eps_diff_re,eps_diff_im] * TensorDiag[1,1,1];
  // epsilon[mid]  = Complex[eps_mil_re,eps_mil_im]   * TensorDiag[1,1,1];
  epsilon1[scat] = Complex[eps_mil_re,eps_mil_im] * TensorDiag[1,1,1];
  epsilon1[mid]  = Complex[eps_mil_re,eps_mil_im] * TensorDiag[1,1,1];
  mur[scat] = TensorDiag[1,1,1];
  mur[mid]  = TensorDiag[1,1,1];

  // epsilon[source_line] = TensorDiag[1,1,1];
  // mur[source_line]     = TensorDiag[1,1,1]; 

  EigFilter[] = (Norm[$EigenvalueReal] > 1e-3);
}

//// FUNCTION //////////////////////////////////////////////////////////////////
Constraint {
  {Name Dirichlet; Type Assign;
    Case { { Region SurfDirichlet ; Value 0.; }}
  }
}

Jacobian {
  { Name JVol ; Case { { Region All ; Jacobian Vol ; } } }
  { Name JSur ; Case { { Region All ; Jacobian Sur ; } } }
  { Name JLin ; Case { { Region All ; Jacobian Lin ; } } }
}

Integration {
  { Name Int_1 ;
    Case { 
      { Type Gauss ;
        Case { 
          { GeoElement Point       ; NumberOfPoints   1 ; }
          { GeoElement Line        ; NumberOfPoints   4 ; }
          { GeoElement Triangle    ; NumberOfPoints   6 ; }
          { GeoElement Tetrahedron ; NumberOfPoints  15 ; }
        }
      }
    }
  }
}
FunctionSpace {
  { Name Hcurl; Type Form1;
    BasisFunction {
      // { Name sn;  NameOfCoef un ; Function BF_Edge      ; Support Region[{Omega}]; Entity EdgesOf[All]; }
      // { Name sn2; NameOfCoef un2; Function BF_Edge_2E   ; Support Region[{Omega}]; Entity EdgesOf[All]; }
      { Name sn;  NameOfCoef un ; Function BF_Edge      ; Support Region[{Omega,source_line}]; Entity EdgesOf[All]; }
      { Name sn2; NameOfCoef un2; Function BF_Edge_2E   ; Support Region[{Omega,source_line}]; Entity EdgesOf[All]; }
    }
   Constraint {
     { NameOfCoef un;  EntityType EdgesOf  ; NameOfConstraint Dirichlet; }
     { NameOfCoef un2; EntityType EdgesOf  ; NameOfConstraint Dirichlet; }
     // { NameOfCoef un3; EntityType FacetsOf ; NameOfConstraint Dirichlet; }
     // { NameOfCoef un4; EntityType FacetsOf ; NameOfConstraint Dirichlet; }
     // { NameOfCoef un5; EntityType EdgesOf  ; NameOfConstraint Dirichlet; }
   }
  }
}
Include "spectral_3D_core_drude.pro";
////////////////////////////////////////////////////////////////////////////////
PostProcessing {
  { Name postpro_norm; NameOfFormulation Ez;
    Quantity {
      { Name norm1 ; Value { Integral { [ {u}*{u} ]; In Omega_in; Integration Int_1; Jacobian JVol; } } }
      { Name norm2 ; Value { Integral { [ epsilon1[]*{u}*{u} ]; In Omega_out; Integration Int_1; Jacobian JVol; } } }
      { Name norm3 ; Value { Integral { [ Inv[mur[]]*{d u}*{d u} ]; In Omega; Integration Int_1; Jacobian JVol; } } }
    }
  }
}

PostOperation {
  { Name postop_norm; NameOfPostProcessing postpro_norm ;
    Operation {
      Print [  norm1[Omega_in], OnElementsOf Point_print, Format TimeTable, File "norm1.txt"];
      Print [  norm2[Omega_out], OnElementsOf Point_print, Format TimeTable, File "norm2.txt"];
      Print [  norm3[Omega], OnElementsOf Point_print, Format TimeTable, File "norm3.txt"];
    }
  }
}
