module mesonbse
using ProgressMeter
using ProgressBars
using TOML
using LinearAlgebra
using Dierckx
using JLD2
using Gaussquad # View in github: https://github.com/kangjiayin/Gaussquad.jl
using FastGaussQuadrature
# using ChebyshevFu
using DataFrames
using CSV

dataset = TOML.parsefile("config.toml")
hashname = hash(dataset)
# quark system
quarkrepoint = dataset["quarkDSE"]["repoint"]
quarkintstep = dataset["quarkDSE"]["quarkintstep"]
quarkm = dataset["quarkDSE"]["quarkmass"]

# data
logofcutoff = dataset["mesonBSE"]["logofcutoff"]
mt = dataset["data"]["mt"]
τ = dataset["data"]["tau"]
Λ = dataset["data"]["lambda"]
ω = dataset["data"]["omega"]
dd = dataset["data"]["dd"]
Nf = dataset["data"]["Nf"]
rm = dataset["data"]["rm"]
cutup = 10. ^logofcutoff
cutdown = 10. ^(-4)

epsilon = dataset["mesonBSE"]["epsilon"]
kstep = dataset["mesonBSE"]["kstep"]
zstep = dataset["mesonBSE"]["zstep"]
Pstep = dataset["mesonBSE"]["Pstep"]
dim = kstep * zstep
# 取点方式
plist=gausslegendremesh(10. ^-4,10. ^3,Pstep,2)[1]
if Pstep == 1
    plist = [0.0001]
end 

if Pstep == 10
    plist = [0.1*j for j=1:10]
end 

if Pstep == 20
    plist = [0.001+0.05*(j-1) for j=1:20]
end 

meshk,weightk= gausslegendremesh(cutdown,cutup,kstep,2);
meshz,weightz= gausschebyshev(zstep,2);

D(t::Float64)=8*pi^2*(dd*exp(-t/(ω^2))/ω^4+rm*((-expm1(-t/(4*mt)^2))/t)/log(τ+(1+t/Λ^2)^2))
D_infrared(t::Float64)=8*pi^2*(dd*exp(-t/(ω^2))/ω^4)
#切比雪夫展开
# include(joinpath(pwd(),"src/mesonfile/chebyshevD.jl"))

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
branchfunction(x::Float64)=(x*AA(x)^2+BB(x)^2)

print("参数导入完毕,开始计算mesonBSA\n")
if dataset["mesonBSE"]["mesonmode"] == 1
    println("将计算psmeson")
    include(joinpath(pwd(),"src/equations/psmeson.jl"))
elseif dataset["mesonBSE"]["mesonmode"] == 2
    println("将计算scalarmeson")
    include(joinpath(pwd(),"src/equations/scalarmeson.jl"))
elseif dataset["mesonBSE"]["mesonmode"] == 3
    println("将计算vectormeson")
    include(joinpath(pwd(),"src/equations/vectormeson.jl"))
elseif dataset["mesonBSE"]["mesonmode"] == 4
    println("将计算avmeson")
    include(joinpath(pwd(),"src/equations/avmeson.jl"))
end # if for mode

end # module
