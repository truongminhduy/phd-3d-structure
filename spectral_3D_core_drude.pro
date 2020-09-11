Formulation {
  {Name Ez; Type FemEquation;
    Quantity {{ Name u   ; Type Local; NameOfSpace Hcurl;}}
    Equation {
      Galerkin { Eig[ 1/mur[]            * Dof{d u}, {d u} ]; Rational 1; In Omega    ; Jacobian JVol; Integration Int_1;}
      Galerkin { Eig[ epsilon1[]/cel^2   * Dof{u}  , {u}   ]; Rational 2; In Omega_in ; Jacobian JVol; Integration Int_1;}
      Galerkin { Eig[ epsilon1[]/cel^2   * Dof{u}  , {u}   ]; Rational 3; In Omega_out; Jacobian JVol; Integration Int_1;}
    }
  }
}
Resolution {
  { Name Projection;
    System{{ Name M1; NameOfFormulation Ez; Type ComplexValue; }}
    Operation{
      GenerateSeparate[M1];
      // EigenSolve[M1,neig,1,30,EigFilter[],
      EigenSolve[M1,neig,zone1,zone2,EigFilter[],
          { {1}, {-epsilon_oo,gamma*epsilon_oo,-omega_p^2,0}, {1,0,0} } ,
          { {1}, {-1,gamma}                                 , {1} } ];
      SaveSolutions[M1];
    }
  }
}
