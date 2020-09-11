//////////////////////////////////////////////////////////////////////////////
Include "parameters_gmsh_getdp.dat";
Include "Pns.dat";
Include "core_3D.pro";
//////////////////////////////////////////////////////////////////////////////
Formulation {
  {Name Ez; Type FemEquation;
    Quantity {{ Name u   ; Type Local; NameOfSpace Hcurl;}}
    Equation {
      Galerkin { [ Dof{u}       , {u} ]; In Omega; Jacobian JVol; Integration Int_1;}
      Galerkin { [-Field[XYZ[]] , {u} ]; In Omega; Jacobian JVol; Integration Int_1;}
    }
  }
}

Resolution {
  { Name Projection;
    System{
      { Name S; NameOfFormulation Ez ; Type ComplexValue; }
    }
    Operation{
      GmshRead["up.pos"];
      // Generate[S];
      // Solve[S];
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
PostProcessing {
  { Name postpro_modal; NameOfFormulation Ez;
    Quantity {        
      { Name uss; Value { Local { [ Field[XYZ[]]] ; In Omega; Jacobian JVol; } } }      
      { Name E2 ; Value { Integral { [ Norm[Field[XYZ[]]]] ; In Omega; Integration Int_1 ; Jacobian JVol ; } } }  
    If (type_source == 1)   
      { Name ut  ; Value { Local { [ Field[XYZ[]] + Einc[] ] ; In Omega; Jacobian JVol; } } } 
      { Name E2t ; Value { Integral { [ Norm[Field[XYZ[]]+ Einc[] ]] ; In Omega; Integration Int_1 ; Jacobian JVol ; } } } 
    EndIf   

    }
  }
}

PostOperation {
  { Name postop_modal; NameOfPostProcessing postpro_modal ;
    Operation {
      Print[ uss, OnElementsOf Omega, File "uss.pos" ];
      Print[ E2[Omega_in], OnGlobal, File "E2.txt" , Format Table  ];
      Print[ E2[Omega], OnGlobal, File "E2_all.txt" , Format Table  ];
      Print[ uss, OnPoint {0,0,0}, Format TimeTable, File "E0_0.txt"];
      Print[ uss, OnPoint {xD1,yD1,zD1}, Format TimeTable, File "E0_1.txt"];
      Print[ uss, OnPoint {xD2,yD2,zD2}, Format TimeTable, File "E0_2.txt"];

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
