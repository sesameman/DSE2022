

A1=Array{Float64}(undef,Pstep,kstep,zstep);
B1=Array{Float64}(undef,Pstep,kstep,zstep);
A2=Array{Float64}(undef,Pstep,kstep,zstep);
B2=Array{Float64}(undef,Pstep,kstep,zstep);
for i=1:Pstep
for j=1:kstep
for m=1:zstep
    A1[i,j,m]=AA((P2[i]/4+meshk[j]+sqrt(P2[i]*meshk[j])*meshz[m]))
    B1[i,j,m]=BB((P2[i]/4+meshk[j]+sqrt(P2[i]*meshk[j])*meshz[m]))
    A2[i,j,m]=AA((P2[i]/4+meshk[j]-sqrt(P2[i]*meshk[j])*meshz[m]))
    B2[i,j,m]=BB((P2[i]/4+meshk[j]-sqrt(P2[i]*meshk[j])*meshz[m]))
end
end
end

# q2=meshk

```
宏定义
```
macro genvar(var::String)
    return quote
    eval(Meta.parse(string($var,"Pqzq=zeros(Pstep,kstep,zstep)")))
    eval(Meta.parse(string($var,"Pq=zeros(Pstep,kstep)")))
    eval(Meta.parse(string($var,"P=zeros(Pstep)")))
    eval(Meta.parse(string("global ",$var,"Pqzq")))
    eval(Meta.parse(string("global ",$var,"Pq")))
    eval(Meta.parse(string("global ",$var,"P")))
    end
end



@genvar("f_gamma5mu")
@genvar("f_gamma5")
@genvar("vectordc")

Threads.@threads for i=1:Pstep
    for j=1:kstep
        # i for P2[i], j for q2[j]
        for m=1:zstep
            P2i=P2[i]
            q2i=meshk[j]
            zq=meshz[m]
            A1in=A1[i,j,m]
            A2in=A2[i,j,m]
            B1in=B1[i,j,m]
            B2in=B2[i,j,m]
            F1=F1k[i,j,m]
            F2=F2k[i,j,m]
            F3=F3k[i,j,m]
            F4=F4k[i,j,m]
            if (dataset["readsetting"]["readmode"] == 3)|(dataset["readsetting"]["readmode"] == 4)
            F5=F5k[i,j,m]
            F6=F6k[i,j,m]
            F7=F7k[i,j,m]
            F8=F8k[i,j,m]
            end
            # mF1=mF1k[i,j,m]
            # mF2=mF2k[i,j,m]
            # mF3=mF3k[i,j,m]
            # mF4=mF4k[i,j,m]
            # F1=0.
            # F2=0.
            # F3=0.
            # F4=0.
            forget_div=1/(A1in^2*(P2i/4+q2i+sqrt(P2i*q2i)*zq)+B1in^2)/(A2in^2*(P2i/4+q2i-sqrt(P2i*q2i)*zq)+B2in^2)

            f_gamma5Pqzq[i,j,m]=(
                2*A2in*B1in*F1 + 2*A1in*B2in*F1 - 4*B1in*B2in*F2 + A1in*A2in*F2*P2i + 4*(A1in*A2in*F2 + A2in*B1in*F4 + A1in*B2in*F4)*q2i + A1in*A2in*F3*P2i*q2i*zq^2 - (4*sqrt(P2i*q2i)*zq*((A2in*B1in - A1in*B2in)*F1 + sqrt(P2i*q2i)*(2*A1in*A2in*F2 + B1in*B2in*F3 + A2in*B1in*F4 + A1in*B2in*F4 + A1in*A2in*F3*q2i)*zq))/P2i
            )*forget_div
            f_gamma5muPqzq[i,j,m]=(
                4*F1*(B1in*B2in + A1in*A2in*q2i) + P2i*(-(A1in*A2in*F1) + 2*(A2in*B1in + A1in*B2in)*F2 - 4*A1in*A2in*F4*q2i) - 4*(A2in*B1in - A1in*B2in)*sqrt(P2i*q2i)*(F2 + F3*q2i)*zq + 2*(A2in*B1in*F3 + A1in*B2in*F3 + 2*A1in*A2in*F4)*P2i*q2i*zq^2
            )*forget_div
            if (dataset["readsetting"]["readmode"] == 3)|(dataset["readsetting"]["readmode"] == 4)
                vectordcPqzq[i,j,m]=P2i*(3*A1in*A2in*F1 - 2*(4*A1in*A2in*F4 + (A2in*B1in + A1in*B2in)*(4*F7 + F8))*q2i) + 16*A1in*A2in*F2*q2i^2*zq^4 - 4*q2i*zq^2*(2*A1in*A2in*F1 + A2in*B1in*F5 + A1in*B2in*F5 + 8*A1in*A2in*F2*q2i + (A2in*B1in - A1in*B2in)*(4*(F6 + F7) + F8)*sqrt(P2i*q2i)*zq + 2*A1in*A2in*F3*P2i*q2i*zq^2) + 2*(-6*B1in*B2in*F1 + 8*A1in*A2in*F2*q2i^2 + (4*A1in*A2in*F4 + (A2in*B1in + A1in*B2in)*(4*F7 + F8))*P2i*q2i*zq^2 + 2*q2i*(-(A1in*A2in*F1) + A2in*B1in*F5 + A1in*B2in*F5 + (A2in*B1in - A1in*B2in)*(4*(F6 + F7) + F8)*sqrt(P2i*q2i)*zq + 2*A1in*A2in*F3*P2i*q2i*zq^2)
                )*forget_div
            end
            vectordcPq[i,j] += weightz[m]*vectordcPqzq[i,j,m]
            f_gamma5Pq[i,j] += weightz[m]*f_gamma5Pqzq[i,j,m]
            f_gamma5muPq[i,j] += weightz[m]*f_gamma5muPqzq[i,j,m]

            # f_gammaPq[i,j] += weightz[m]*(P2[i]*A2[i,j,m]*((F1k[i,j,m] + 4*F4k[i,j,m]*q2[j])*A1[i,j,m] - 2*F2k[i,j,m]*B1[i,j,m]) - 2*F2k[i,j,m]*P2[i]*A1[i,j,m]*B2[i,j,m] + 4*sqrt(P2[i]*q2[j])*(F2k[i,j,m] + F3k[i,j,m]*q2[j])*meshz[m]*(A2[i,j,m]*B1[i,j,m] - A1[i,j,m]*B2[i,j,m]) - 2*P2[i]*q2[j]*meshz[m]^2*(2*F4k[i,j,m]*A1[i,j,m]*A2[i,j,m] + F3k[i,j,m]*A2[i,j,m]*B1[i,j,m] + F3k[i,j,m]*A1[i,j,m]*B2[i,j,m]) - 4*F1k[i,j,m]*(q2[j]*A1[i,j,m]*A2[i,j,m] + B1[i,j,m]*B2[i,j,m]))
        end
    end
end

for i=1:Pstep
    for j=1:kstep
        vectordcP[i]+=weightk[j]*vectordcPq[i,j]*meshk[j]
        f_gamma5muP[i]+=weightk[j]*f_gamma5muPq[i,j]*meshk[j]
        f_gamma5P[i]+=weightk[j]*f_gamma5Pq[i,j]*meshk[j]
    end
end

vectordcP=vectordcP *(2*pi)^-3*z2/2*3
f_gamma5P=f_gamma5P *(2*pi)^-3*z4/2*3
f_gamma5muP=f_gamma5muP *(2*pi)^-3*z2/2*3

