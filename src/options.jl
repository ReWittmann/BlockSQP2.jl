function blockSQP_options()
    opts = Dict()
    opts["opttol"] = 1.0e-12
    opts["nlinfeastol"] = 1.0e-12
    opts["globalization"] = 0
    opts["hessUpdate"] = 2
    opts["fallbackUpdate"] = 2
    opts["hessScaling"] = 0
    opts["fallbackScaling"] = 0
    opts["hessLimMem"] = 1
    opts["hessMemsize"] = 20
    opts["maxConsecSkippedUpdates"] = 200
    opts["blockHess"] = 1
    opts["whichSecondDerv"] = 0
    opts["sparseQP"] = 0
    opts["printLevel"] = 2
    opts["debugLevel"] = 0
    opts["which_QPsolver"] = "qpOASES"
    return opts
end

function BSQP_options(param::Dict)
    cxx_opts = SQPoptions()
    opt_keys = keys(param)

    if "printLevel" in opt_keys
        set_printLevel(cxx_opts, Int32(param["printLevel"]))
    end
    if "printColor" in opt_keys
        set_printColor(cxx_opts, Int32(param["printColor"]))
    end
    if "debugLevel" in opt_keys
        set_debugLevel(cxx_opts, Int32(param["debugLevel"]))
    end
    if "eps" in opt_keys
        set_eps(cxx_opts, Float64(param["eps"]))
    end
    if "inf" in opt_keys
        set_inf(cxx_opts, Float64(param["inf"]))
    end
    if "opttol" in opt_keys
        set_opttol(cxx_opts, Float64(param["opttol"]))
    end
    if "nlinfeastol" in opt_keys
        set_nlinfeastol(cxx_opts, Float64(param["nlinfeastol"]))
    end
    if "sparseQP" in opt_keys
        set_sparseQP(cxx_opts, Int32(param["sparseQP"]))
    end
    if "globalization" in opt_keys
        set_globalization(cxx_opts, Int32(param["globalization"]))
    end
    if "restoreFeas" in opt_keys
        set_restoreFeas(cxx_opts, Int32(param["restoreFeas"]))
    end
    if "maxLineSearch" in opt_keys
        set_maxLineSearch(cxx_opts, Int32(param["maxLineSearch"]))
    end
    if "maxConsecReducedSteps" in opt_keys
        set_maxConsecReducedSteps(cxx_opts, Int32(param["maxConsecReducedSteps"]))
    end
    if "maxConsecSkippedUpdates" in opt_keys
        set_maxConsecSkippedUpdates(cxx_opts, Int32(param["maxConsecSkippedUpdates"]))
    end
    if "maxItQP" in opt_keys
        set_maxItQP(cxx_opts, Int32(param["maxItQP"]))
    end
    if "blockHess" in opt_keys
        set_blockHess(cxx_opts, Int32(param["blockHess"]))
    end
    if "hessScaling" in opt_keys
        set_hessScaling(cxx_opts, Int32(param["hessScaling"]))
    end
    if "fallbackScaling" in opt_keys
        set_fallbackScaling(cxx_opts, Int32(param["fallbackScaling"]))
    end
    if "maxTimeQP" in opt_keys
        set_maxTimeQP(cxx_opts, Int32(param["maxTimeQP"]))
    end
    if "iniHessDiag" in opt_keys
        set_iniHessDiag(cxx_opts, Float64(param["iniHessDiag"]))
    end
    if "colEps" in opt_keys
        set_colEps(cxx_opts, Float64(param["colEps"]))
    end
    if "colTau1" in opt_keys
        set_colTau1(cxx_opts, Float64(param["colTau1"]))
    end
    if "colTau2" in opt_keys
        set_colTau2(cxx_opts, Float64(param["colTau2"]))
    end
    if "hessDamp" in opt_keys
        set_hessDamp(cxx_opts, Int32(param["hessDamp"]))
    end
    if "hessDampFac" in opt_keys
        set_hessDampFac(cxx_opts, Float64(param["hessDampFac"]))
    end
    if "hessUpdate" in opt_keys
        set_hessUpdate(cxx_opts, Int32(param["hessUpdate"]))
    end
    if "fallbackUpdate" in opt_keys
        set_fallbackUpdate(cxx_opts, Int32(param["fallbackUpdate"]))
    end
    if "hessLimMem" in opt_keys
        set_hessLimMem(cxx_opts, Int32(param["hessLimMem"]))
    end
    if "hessMemsize" in opt_keys
        set_hessMemsize(cxx_opts, Int32(param["hessMemsize"]))
    end
    if "whichSecondDerv" in opt_keys
        set_whichSecondDerv(cxx_opts, Int32(param["whichSecondDerv"]))
    end
    if "skipFirstGlobalization" in opt_keys
        set_skipFirstGlobalization(cxx_opts, Int32(param["skipFirstGlobalization"]))
    end
    if "convStrategy" in opt_keys
        set_convStrategy(cxx_opts, Int32(param["convStrategy"]))
    end
    if "maxConvQP" in opt_keys
        set_maxConvQP(cxx_opts, Int32(param["maxConvQP"]))
    end
    if "maxSOCiter" in opt_keys
        set_maxSOCiter(cxx_opts, Int32(param["maxSOCiter"]))
    end

    return cxx_opts

end