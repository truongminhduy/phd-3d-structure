//////////////////////////////////////////////////////////////////////////////
Include "parameters_gmsh_getdp.dat";
Include "Pns.dat";
Function{
  For n In {0:neig-1}
    P~{n}[] = Complex[ Pns_re~{n}, Pns_im~{n} ];
  EndFor
}
Include "core_3D.pro";
Include "spectral_3D_core_drude.pro";
////////////////////////////////////////////////////////////////////////////////
PostProcessing {
  { Name postpro_modal; NameOfFormulation Ez;
    Quantity {
      { Name up   ;
        Value {
            For i In {0:neig-1}
              Local { [ P~{i}[] * {u}[neig-1-i] ]; In Omega; Jacobian JVol; }
            EndFor
        }
      }
    }
  }
}

PostOperation {
    { Name postop_modal; NameOfPostProcessing postpro_modal ;
          Operation {
              Print [ up, OnElementsOf Omega, File "up.pos", LastTimeStepOnly ];
          }
    }
}
