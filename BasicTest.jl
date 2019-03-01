include("dftmtx.jl")
maxa=zeros(5,1)
maxpa=zeros(5,1)
maxb=zeros(5,1)
maxpb=zeros(5,1)
maxd=zeros(5,1)
maxdd=zeros(5,1)

for i=1:5
Alpha=1
N=2^(i+2)
M=round(Int,N/4)
j=collect(0:N-1)
x=j/N
y=x
xx=collect(0:2*N-1)/N
yy=xx;
f=2*pi/(2*N)
w=(collect(0:f:2*pi-f/2)*1im)
tt=collect(0:(2*N)-1)'
A=exp.(-w*tt)


    
# Analysis of Cont Matrix and pseudo inverses via svd
Ua,Sa,Va=svd(A)
maxa[i]=maximum(Sa)
PA=pinv(A)
Uaa,Saa,Vaa=svd(PA)
maxpa[i]=maximum(Saa)

B=A[1:N,vcat(1:M+1,2*N-M+1:2*N)]
C=A[:, vcat(1:M+1,2*N-M+1:2*N)]
dUb,Sb,Vb=svd(B)
maxb[i]=maximum(Sb)
B_dag=pinv(B)
Upb,Spb,Vpb=svd(B_dag)
maxpb[i]=maximum(Spb)

# construct derivative operator L=I-alphad^2/dx^2

K=-pi*1im*vcat(collect(0:N-1),0,collect(-(N-1):-1))
# derivative coefficients pertaining to the FC
    n=vcat(K[1:M+1],K[(2*N)-M+1:(2*N)])
    D=eye(length(n))-Alpha*diagm(n.*n)
    Ud,Sd,Vd=svd(D)
    maxd[i]=maximum(Sd)
    D_dag=pinv(D)
    Udd,Sdd,Vdd=svd(D_dag)
    maxdd[i]=maximum(Sdd)
    end
    
