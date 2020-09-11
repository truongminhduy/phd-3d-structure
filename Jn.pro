//////////////////////////////////////////////////////////////////////////////
Include "parameters_gmsh_getdp.dat";
Include "core_3D.pro";
Include "spectral_3D_core_drude.pro";
////////////////////////////////////////////////////////////////////////////////
PostProcessing {
  { Name postpro_Jn; NameOfFormulation Ez;
    Quantity {
      If (type_source == 0)
        { Name Jn ; Value { Integral { [ ((CompZ[Tangent[]]<0) ? 1 : -1)*source[]*{u} ]; In source_line; Integration Int_1; Jacobian JLin; } } }
        // { Name Jn ; Value { Integral { [ source[]*{u} ]; In source_line; Integration Int_1; Jacobian JLin; } } }
      Else
        { Name Jn ; Value { Integral { [ source[]*{u} ]; In Omega_in; Integration Int_1; Jacobian JVol; } } }  
      EndIf
    }
  }
}

PostOperation {
  { Name postop_Jn; NameOfPostProcessing postpro_Jn ;
    Operation {
      If (type_source == 0)
        Print [  Jn[source_line], OnElementsOf Point_print, Format TimeTable, File "Jns.txt"];
      Else
        // Print [  Jn[Omega], OnElementsOf Point_print, Format TimeTable, File "Jns.txt"];
        Print [  Jn[Omega_in], OnElementsOf Point_print, Format TimeTable, File "Jns.txt"];
      EndIf    
    }
  }
}
