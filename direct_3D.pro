/////////////////////////////////////////////////////////////////////////
Include "parameters_gmsh_getdp.dat";
Include "core_3D.pro";
////////////////////////////////////////////////////////////////////////
Formulation {
  {Name F_main; Type FemEquation;
    Quantity {
      { Name u; Type Local; NameOfSpace Hcurl;}
    }
    Equation {
      Galerkin { [k0^2*epsilon[]*Dof{u}   , {u}   ]; In Omega; Jacobian JVol; Integration Int_1;  }
      Galerkin { [-Inv[mur[]]   *Dof{d u} , {d u} ]; In Omega; Jacobian JVol; Integration Int_1;  }
      If (type_source == 0)
        // Galerkin { [-source[]           , {u} ]; In source_line; Jacobian JLin; Integration Int_1; }
        Galerkin { [((CompZ[Tangent[]]<0) ? -1 : 1)*source[], {u} ]; In source_line; Jacobian JLin; Integration Int_1; }
      Else
        Galerkin { [-source[]           , {u}   ]; In Omega_in; Jacobian JVol; Integration Int_1;  }
      EndIf
    }
  } 
}

Resolution {
  { Name Scattering;
    System {
      { Name S; NameOfFormulation F_main; Type ComplexValue; Frequency Freq;}
    }
    Operation {
      Generate[S];
      Solve[S];
      SaveSolutions[S];
    }
  }
}

////// DATA PROCESS ////////////////////////////////////////////////////////////
PostProcessing {
  { Name get_scat; NameOfFormulation F_main;
    Quantity {
      { Name us  ; Value { Local { [ {u}  ] ; In Omega; Jacobian JVol; } } }   
      { Name E2 ; Value { Integral { [ Norm[{u}]] ; In Omega; Integration Int_1 ; Jacobian JVol ; } } } 
      // { Name ss  ; Value { Local { [ source[]   ] ; In Omega; Jacobian JVol; } } }
      // { Name source  ; Value {  Integral { [ source[] ]; In source_line; Integration Int_1 ;Jacobian JLin; } } }   
    If (type_source == 1)   
      { Name ut  ; Value { Local { [ {u} + Einc[] ] ; In Omega; Jacobian JVol; } } } 
      { Name E2t ; Value { Integral { [ Norm[{u}+ Einc[] ]] ; In Omega; Integration Int_1 ; Jacobian JVol ; } } } 
    EndIf      
    }
  }
}

PostOperation {
  { Name postop_scat; NameOfPostProcessing get_scat ;
    Operation {
      Print[ us, OnElementsOf Omega, File "us.pos" ];
      // Print[ source[source_line], OnGlobal, File "S2.txt" , Format Table  ];
      // Print[ ss, OnElementsOf Omega, File "ss.pos" ];
      Print[ E2[Omega_in], OnGlobal, File "E2.txt" , Format Table  ];
      Print[ E2[Omega_nopml], OnGlobal, File "E2_all.txt" , Format Table  ];
      Print[ us, OnPoint {0,0,0}, Format TimeTable, File "E0_0.txt"];
      Print[ us, OnPoint {xD1,yD1,zD1}, Format TimeTable, File "E0_1.txt"];
      Print[ us, OnPoint {xD2,yD2,zD2}, Format TimeTable, File "E0_2.txt"];

      If (type_source == 1) 
      Print[ E2t[Omega_in], OnGlobal, File "E2t.txt" , Format Table  ];
      Print[ E2t[Omega_nopml], OnGlobal, File "E2_allt.txt" , Format Table  ];
      Print[ ut, OnPoint {0,0,0}, Format TimeTable, File "E0_0t.txt"];
      Print[ ut, OnPoint {xD1,yD1,zD1}, Format TimeTable, File "E0_1t.txt"];
      Print[ ut, OnPoint {xD2,yD2,zD2}, Format TimeTable, File "E0_2t.txt"];
      EndIf 
    }
  }
}
///////////////////////////////////////////////////////////// END //////////////
////////////////////////////////////////////////////////////////////////////////
