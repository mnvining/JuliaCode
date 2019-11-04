module AllModule

using SymPy
using LinearAlgebra
using GenericSVD

include("ConvertToBig.jl")
export ConvertToBig
include("hsol.jl")
export hsol,evalasarray
include("Spaces.jl")
export linspace,logspace
include("myhouse.jl")
export MyHouse,HRQR,BackSolve,SolveViaQR
include("coeffcalc.jl")
export coeffcalc,myvander,fliplr
include("All.jl")
export fmaker
include("CosCont.jl")
export CosCont
include("SinCont.jl")
export SinCont
include("DiffOp.jl")
export DiffOp
include("CosMtx.jl")
export CosMtx
include("SinMtx.jl")
export SinMtx
include("DiffOp.jl")
export DiffOp
include("MP.jl")
export MP
include("NCalc.jl")
export NCalc
include("IDOCalc.jl")
export IDOCalc
include("TotalCalc.jl")
export TotalCalc
include("Storage.jl")
export FCMat,AMat,FCMatS,AMatS,C_C,C_S,m,mm,IDOC,IDOS



end #module
