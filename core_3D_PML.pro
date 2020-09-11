/////////////////////////////////////////////////////////////////////////////////////////////////
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
    // // with out PML
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
/////////////////////////////////////////////////////////////////////////////
Function{
  Freq   = cel/lambda0;
  omega0 = 2.*Pi*cel/lambda0;
  k0     = 2.*Pi/lambda0;   
  
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

  epsilon[scat] = Complex[eps_diff_re,eps_diff_im] * TensorDiag[1,1,1];
  epsilon[mid]  = Complex[eps_mil_re,eps_mil_im]   * TensorDiag[1,1,1];
  epsilon1[scat] = Complex[eps_mil_re,eps_mil_im] * TensorDiag[1,1,1];
  epsilon1[mid]  = Complex[eps_mil_re,eps_mil_im] * TensorDiag[1,1,1];
  mur[scat] = TensorDiag[1,1,1];
  mur[mid]  = TensorDiag[1,1,1];

  // epsilon[source_line] = TensorDiag[1,1,1];
  // mur[source_line] = TensorDiag[1,1,1]; 

  If (type_source == 0)
    // source[source_line] = 1/r_source;
    source[source_line] = Vector[0,0,1]*1/r_source;
    // source[Omega] = 0;
    // source[source_line] = 1/r_source;
    // source[source_line] = TensorDiag[0,0,1]*1/r_source;
    // source[Omega] = 0;
    // source[source] = Vector[0,0,1]*( 1/( 4/3*Pi*(r_source)^3 ) );
    // source[source] = 1/( 4/3*Pi*(r_source)^3 );
    // source[source] = TensorDiag[0,0,1]*( 1/( 4/3*Pi*(r_source)^3 ) );
  Else
    Ae     = 1.0;
    theta0 = 0;
    phi0   = 0;
    psi0   = 0;
    k_mid  = k0*Sqrt[eps_mil_re];
    alpha0 = -k_mid*Sin[theta0]*Cos[phi0];
    beta0  = -k_mid*Sin[theta0]*Sin[phi0];
    gamma0 = -k_mid*Cos[theta0];
    Prop[] = Complex[ Cos[alpha0*X[]+beta0*Y[]+gamma0*Z[]] , Sin[alpha0*X[]+beta0*Y[]+gamma0*Z[]] ];
    Ex0    =  Ae * (Cos[psi0]*Cos[theta0]*Cos[phi0] - Sin[psi0]*Sin[phi0]);
    Ey0    =  Ae * (Cos[psi0]*Cos[theta0]*Sin[phi0] + Sin[psi0]*Cos[phi0]);
    Ez0    =  Ae * (-Cos[psi0]*Sin[theta0]);
    Einc[] = Vector[Ex0*Prop[],Ey0*Prop[],Ez0*Prop[]];
    E0[]   = Vector[Ex0,Ey0,Ez0];
    source[] = (omega0/cel)^2*(epsilon1[]-epsilon[])*Einc[];
  EndIf

  If (E_or_H == 0)
    chi[] = epsilon[];
    xsi[] = Inv[mur[]];  
  Else
    chi[] = mur[];  
    xsi[] = Inv[epsilon[]];
  EndIf

  EigFilter[] = (Norm[$EigenvalueReal] > 1e-5);
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
          { GeoElement Point       ; NumberOfPoints  1 ; }
          { GeoElement Line        ; NumberOfPoints  3 ; }
          { GeoElement Triangle    ; NumberOfPoints  4 ; }
          { GeoElement Quadrangle  ; NumberOfPoints  4 ; }
          { GeoElement Tetrahedron ; NumberOfPoints  15; }
          // { GeoElement Hexahedron  ; NumberOfPoints  6 ; }
          // { GeoElement Prism       ; NumberOfPoints  6 ; } 
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
   }
  }
}