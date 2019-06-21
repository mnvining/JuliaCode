using SymPy
using LinearAlgebra
using GenericSVD

function NCalc(d,w,FCMatC,FCMatS,C_C,C_S,IDOC,IDOS,m,mm)
  include("MP.jl")
  include("All.jl")
  include("Spaces.jl")
  (yc,ys,Gc,Gs)=MP(d,w,FCMatC,FCMatS,C_C,C_S,IDOC,IDOS)
  setprecision(230);
  D=BigFloat(350)/BigFloat(100);
  A=BigFloat(125)/BigFloat(100);
  n=10;
  F=10;
  x=Sym("x")

  f(x)=fmaker(d,n);
  xx=linspace(BigFloat(-1),BigFloat(1),Int(F*(n-1)+1));
  y=zeros(BigFloat,size(xx))
  y=evalasarray(f(x),xx)

  CAcc=maximum(abs.(C_C*yc.-y));
  SAcc=maximum(abs.(C_S*ys.-y));

  uc=C_C*IDOC*yc;
  us=C_S*IDOS*ys;
  CStab=norm(uc[1:10:end])/norm(y[1:10:end])
  SStab=norm(us[1:10:end])/norm(y[1:10:end])
  return CAcc,SAcc,CStab,SStab
end
