# 6-19
# kangjiayin
workdir = "/Users/kangjiayin/Desktop/program/julia/DSE2022"
cd(workdir)
using TOML
using LinearAlgebra
using Dierckx
using JLD2
using Gaussquad # View in github: https://github.com/kangjiayin/Gaussquad.jl
using FastGaussQuadrature
using ChebyshevFun 
using CSV
using DataFrames
include(joinpath(pwd(),"module/inter.jl"))
dataset = TOML.parsefile("config.toml")
hashname = hash(dataset)
# quark system
include(joinpath(pwd(),"src/equations/quark_gap.jl"))
quarkrepoint = dataset["readsetting"]["repoint"]
quarkintstep = dataset["readsetting"]["quarkintstep"]
quarkm = dataset["readsetting"]["quarkmass"]

# data
logofcutoff = dataset["readsetting"]["logofcutoff"]
cutup = 10. ^logofcutoff
cutdown = 10. ^(-4)

kstep = dataset["readsetting"]["kstep"]
zstep = dataset["readsetting"]["zstep"]
Pstep = dataset["readsetting"]["Pstep"]
epsilon = dataset["readsetting"]["epsilon"]
# 取点方式
meshk,weightk= gausslegendremesh(cutdown,cutup,kstep,2);
meshz,weightz= gausschebyshev(zstep,2);

# using SPMinterpolation
Ppoints=Pstep
# kstep=128
# k=gausslegendremesh(10. ^-4,10. ^4,32,2);
P2=Array{Float64}(undef,Ppoints,1)
F1k=Array{Float64}(undef,Ppoints,kstep,zstep)
F2k=Array{Float64}(undef,Ppoints,kstep,zstep)
F3k=Array{Float64}(undef,Ppoints,kstep,zstep)
F4k=Array{Float64}(undef,Ppoints,kstep,zstep)
if (dataset["readsetting"]["readmode"] == 3)|(dataset["readsetting"]["readmode"] == 4)
F5k=Array{Float64}(undef,Ppoints,kstep,zstep)
F6k=Array{Float64}(undef,Ppoints,kstep,zstep)
F7k=Array{Float64}(undef,Ppoints,kstep,zstep)
F8k=Array{Float64}(undef,Ppoints,kstep,zstep)
end
# mP2=Array{Float64}(undef,Ppoints,1)
# mF1k=Array{Float64}(undef,Ppoints,kstep,zstep)
# mF2k=Array{Float64}(undef,Ppoints,kstep,zstep)
# mF3k=Array{Float64}(undef,Ppoints,kstep,zstep)
# mF4k=Array{Float64}(undef,Ppoints,kstep,zstep)

# F1k=Array{Float64}(undef,Ppoints,kstep)
# F2k=Array{Float64}(undef,Ppoints,kstep)
# F3k=Array{Float64}(undef,Ppoints,kstep)
# F4k=Array{Float64}(undef,Ppoints,kstep)
kname=Array{String}(undef,Ppoints,1)
if dataset["readsetting"]["readmode"] == 1
    if ispath("data/pseudo_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint-$epsilon") == true
        print("导入--","data/pseudo_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint-$epsilon\n")
        for i=1:Ppoints
            # local a, b
            global P2, F1k, F2k, F3k, F4k
            P2[i], F1k[i,:,:], F2k[i,:,:], F3k[i,:,:], F4k[i,:,:] =load("data/pseudo_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint-$epsilon/P&F1-4_$i-$Pstep.jld2","P2", "F1", "F2", "F3", "F4")
        end #for i
    else #ispath
        print("不存在该文件","data/pseudo_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint-$epsilon","\n")
    end
elseif dataset["readsetting"]["readmode"] == 2
    if ispath("data/scalar_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint-$epsilon") == true
        print("导入","data/scalar_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint-$epsilon\n")
        for i=1:Ppoints
            # local a, b
            global P2, F1k, F2k, F3k, F4k
            P2[i], F1k[i,:,:], F2k[i,:,:], F3k[i,:,:], F4k[i,:,:] =load("data/scalar_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint-$epsilon/P&F1-4_$i-$Pstep.jld2","P2", "F1", "F2", "F3", "F4")
        end #for i
    else #ispath
        print("不存在该文件,","data/scalar_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint-$epsilon","\n")
    end
elseif dataset["readsetting"]["readmode"] == 3
    if ispath("data/vector_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint-$epsilon") == true
        print("导入","data/vector_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint-$epsilon\n")
        for i=1:Ppoints
            # local a, b
            global P2, F1k, F2k, F3k, F4k, F5k, F6k, F7k, F8k
            P2[i], F1k[i,:,:], F2k[i,:,:], F3k[i,:,:], F4k[i,:,:], F5k[i,:,:], F6k[i,:,:], F7k[i,:,:], F8k[i,:,:] =load("data/vector_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint-$epsilon/P&F1-4_$i-$Pstep.jld2","P2", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8")
        end #for i
    else #ispath
        print("不存在该文件,","data/vector_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint-$epsilon","\n")
    end
elseif dataset["readsetting"]["readmode"] == 4
    if ispath("data/av_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint-$epsilon") == true
        print("导入","data/av_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint-$epsilon\n")
        for i=1:Ppoints
            # local a, b
            global P2, F1k, F2k, F3k, F4k, F5k, F6k, F7k, F8k
            P2[i], F1k[i,:,:], F2k[i,:,:], F3k[i,:,:], F4k[i,:,:], F5k[i,:,:], F6k[i,:,:], F7k[i,:,:], F8k[i,:,:] =load("data/av_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint-$epsilon/P&F1-4_$i-$Pstep.jld2","P2", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8")
        end #for i
    else #ispath
        print("不存在该文件,","data/av_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint-$epsilon","\n")
    end
end # if

F1pq=zeros(Pstep,kstep);
F2pq=zeros(Pstep,kstep);
F3pq=zeros(Pstep,kstep);
F4pq=zeros(Pstep,kstep);
if (dataset["readsetting"]["readmode"] == 3)|(dataset["readsetting"]["readmode"] == 4)
    F5pq=zeros(Pstep,kstep);
    F6pq=zeros(Pstep,kstep);
    F7pq=zeros(Pstep,kstep);
    F8pq=zeros(Pstep,kstep);
end
for i=1:Pstep
    for j=1:kstep
        # i for P2[i], j for q2[j]
        for m=1:zstep
            F1pq[i,j] += weightz[m]*F1k[i,j,m]*2/pi
            F2pq[i,j] += weightz[m]*F2k[i,j,m]*2/pi
            F3pq[i,j] += weightz[m]*F3k[i,j,m]*2/pi
            F4pq[i,j] += weightz[m]*F4k[i,j,m]*2/pi
            if (dataset["readsetting"]["readmode"] == 3)|(dataset["readsetting"]["readmode"] == 4)
            F5pq[i,j] += weightz[m]*F5k[i,j,m]*2/pi
            F6pq[i,j] += weightz[m]*F6k[i,j,m]*2/pi
            F7pq[i,j] += weightz[m]*F7k[i,j,m]*2/pi
            F8pq[i,j] += weightz[m]*F8k[i,j,m]*2/pi
            end
            # f_gammaPq[i,j] += weightz[m]*(P2[i]*A2[i,j,m]*((F1k[i,j,m] + 4*F4k[i,j,m]*q2[j])*A1[i,j,m] - 2*F2k[i,j,m]*B1[i,j,m]) - 2*F2k[i,j,m]*P2[i]*A1[i,j,m]*B2[i,j,m] + 4*sqrt(P2[i]*q2[j])*(F2k[i,j,m] + F3k[i,j,m]*q2[j])*meshz[m]*(A2[i,j,m]*B1[i,j,m] - A1[i,j,m]*B2[i,j,m]) - 2*P2[i]*q2[j]*meshz[m]^2*(2*F4k[i,j,m]*A1[i,j,m]*A2[i,j,m] + F3k[i,j,m]*A2[i,j,m]*B1[i,j,m] + F3k[i,j,m]*A1[i,j,m]*B2[i,j,m]) - 4*F1k[i,j,m]*(q2[j]*A1[i,j,m]*A2[i,j,m] + B1[i,j,m]*B2[i,j,m]))
        end
    end
end

consta1=0.
consta2=0.
constb1=0.
constb2=0.
function Inport()
    global z2, z4, AA1, BB1, consta1, consta2, constb1, constb2
    local A, B, k, mk
    A, B, k, z2, z4=load("data/quark_gap/$hashname.jld2","A","B","k", "z2", "z4");
    consta1, constb1 = Main.interDSE.ExtraA(k,A)
    consta2, constb2 = Main.interDSE.ExtraB(k,B)
    AA1=Spline1D(k,A)
    BB1=Spline1D(k,B)
    mk=k[length(k)]
    return mk
end
maxofk=Inport()
function AA(x)
    if x<maxofk
        return AA1(x)
    else
        return consta1/log(x/0.234^2)^constb1
    end
end
function BB(x)
    if x<maxofk
        return BB1(x)
    else
        return consta2/log(x/0.234^2)^constb2
    end
end

# include("decayconst.jl")
# include("decayconst2.jl")
# print(F1pq[:,1])
function printdata()
    lengthoffunc = 4
    if (dataset["readsetting"]["readmode"] == 3)|(dataset["readsetting"]["readmode"] == 4)
        lengthoffunc = 8
    end
    for i=1:lengthoffunc
        print("dataF[",i,"]={")
        for  j=1:Pstep
            if i==1
                print("{",P2[j],",",1/F1pq[j,1],"}")
            elseif i==2
                print("{",P2[j],",",1/F2pq[j,1],"}")
            elseif i==3
                print("{",P2[j],",",1/F3pq[j,1],"}")
            elseif i==4
                print("{",P2[j],",",1/F4pq[j,1],"}")
            elseif i==5
                print("{",P2[j],",",1/F5pq[j,1],"}")
            elseif i==6
                print("{",P2[j],",",1/F6pq[j,1],"}")
            elseif i==7
                print("{",P2[j],",",1/F7pq[j,1],"}")
            elseif i==8
                print("{",P2[j],",",1/F8pq[j,1],"}")
            end
            if j==Pstep
                print("};")
            else
                print(",")
            end
        end
    print("\n")
    end
end