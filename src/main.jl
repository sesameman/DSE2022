# 该文件旨在设定工作路径

wkdir="/Users/kangjiayin/Desktop/program/julia/DSE2022"
cd(wkdir)

# include(joinpath(pwd(),"module/testpkg.jl")) # 检测包完整性
include(joinpath(pwd(),"src/equations/quark_gap.jl"))
include(joinpath(pwd(),"module/inter.jl"))
include(joinpath(pwd(),"src/meson.jl"))