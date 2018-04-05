function dftmtx(n)
f=2*pi/n
w=(collect(0:f:2*pi-f/2)*1im)*BigFloat(1)
x=collect(0:n-1)'
    return D=exp.(-w*x)

end


function createFWDcont(size_original_grid,vec_of_modes,cl)
    ZZ=cl*size_original_grid/2+size_original_grid
    A=dftmtx(ZZ)
    # Pick specific modes, restrict
    return C=A[1:size_original_grid,vec_of_modes]
end

function createBackDoubleDomain(size_original_grid,vec_of_modes,cl)
    # Create dftmtx for extended period
    ZZ=cl*size_original_grid/2+size_original_grid
    A=dftmtx(ZZ)
    # Pick specific modes, restrict
    return B=A[:,vec_of_modes]
end

function createId2Derivative(NG,M,Alpha,cl)
    fd=cl+cl/(NG-1);
    ZZ=round(Int,cl*NG/2+NG)
    LL=round(Int,ZZ/2)
K=-pi*1im*vcat(collect(0:LL-1),0,collect(-(LL-1):-1))*BigFloat(1)
# derivative coefficients pertaining to the FC
    n=vcat(K[1:M+1],K[(ZZ)-M+1:(ZZ)])/fd;
    return eye(length(n))-Alpha*diagm(n.*n)
end



    



