if ispath("data/scalar_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint-$epsilon") == false
    mkdir("data/scalar_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint-$epsilon")
    cp(joinpath(pwd(),"config.toml"),joinpath(pwd(),"data/scalar_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint-$epsilon/log.toml"))
    CSV.write(joinpath(pwd(),"data/scalar_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint-$epsilon/plist.CSV"),DataFrame(p2=plist))
    @showprogress for indexforp2=1:Pstep
        timetest1=time()
        global P2
        P2=plist[indexforp2]
        # 分配点与权重
        include(joinpath(pwd(),"src/mesonfile/setupkernel.jl"))
        # 计算kernel
        include(joinpath(pwd(),"src/mesonfile/scalarkernel.jl"))
        # 求解函数
        include(joinpath(pwd(),"src/mesonfile/solve_kernel.jl"))
        # 保存文件
        jldsave("data/scalar_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint-$epsilon/P&F1-4_$indexforp2-$Pstep.jld2";P2, F1, F2, F3, F4)
        # print("mesonBSA--$P2 for $indexforp2/$Pstep done, takes",round((time()-timetest1)*100)/100,"s \n")
    end # for
else # 
    print("已存在文件--",joinpath(pwd(),"data/scalar_BSE/meson-$kstep-$zstep-$Pstep-$quarkm-$logofcutoff-$quarkintstep-$quarkrepoint-$epsilon"),"\n")
end # ispath