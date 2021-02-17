using Distributed
using Base.Filesystem
using DataFrames
using CSV
using Query
using Statistics
using UnicodePlots
using ClusterManagers
using Dates
using DelimitedFiles

## load the packages by covid19abm

#using covid19abm

#addprocs(2, exeflags="--project=.")


#@everywhere using covid19abm

addprocs(SlurmManager(250), N=8, topology=:master_worker, exeflags="--project=.")
#addprocs(4)
@everywhere using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames
@everywhere include("codeabmBr.jl")



function run(myp::ModelParameters, nsims=1000, folderprefix="./")
    println("starting $nsims simulations...\nsave folder set to $(folderprefix)")
    dump(myp)
    myp.calibration && !myp.ignore_cal && error("can not run simulation, calibration is on.")
    # will return 6 dataframes. 1 total, 4 age-specific 
    cdr = pmap(1:nsims) do x                 
            runsim(x, myp)
    end      

    println("simulations finished")
    println("total size of simulation dataframes: $(Base.summarysize(cdr))")

    R0 = [cdr[i].R0 for i=1:nsims]
    R0_r = [cdr[i].R0_r for i=1:nsims]
    
    allag = vcat([cdr[i].a  for i = 1:nsims]...)
    ag1 = vcat([cdr[i].g1 for i = 1:nsims]...)
    ag2 = vcat([cdr[i].g2 for i = 1:nsims]...)
    ag3 = vcat([cdr[i].g3 for i = 1:nsims]...)
    ag4 = vcat([cdr[i].g4 for i = 1:nsims]...)
    ag5 = vcat([cdr[i].g5 for i = 1:nsims]...)
    ag6 = vcat([cdr[i].g6 for i = 1:nsims]...)
    mydfs = Dict("all" => allag, "ag1" => ag1, "ag2" => ag2, "ag3" => ag3, "ag4" => ag4, "ag5" => ag5, "ag6" => ag6)
    

    c1 = Symbol.((:LAT, :HOS, :ICU, :DED), :_INC)
    c2 = Symbol.((:LAT, :HOS, :ICU, :DED), :_PREV)
    for (k, df) in mydfs
        println("saving dataframe sim level: $k")
        # simulation level, save file per health status, per age group
        #for c in vcat(c1..., c2...)
        for c in vcat(c1...)
        #for c in vcat(c2...)
            udf = unstack(df, :time, :sim, c) 
            fn = string("$(folderprefix)/simlevel_", lowercase(string(c)), "_", k, ".dat")
            CSV.write(fn, udf)
        end
        println("saving dataframe time level: $k")
        # time level, save file per age group
        #yaf = compute_yearly_average(df)       
        #fn = string("$(folderprefix)/timelevel_", k, ".dat")   
        #CSV.write(fn, yaf)       
    end

    writedlm(string(folderprefix,"/R0.dat"),R0)
    writedlm(string(folderprefix,"/R0_r.dat"),R0_r)

end



function create_folder(ip::ModelParameters,heatmap=false)
    #strategy = ip.apply_vac_com == true ? "T" : "UT"
    #strategy = ip.vaccinating == true ? "$(ip.days_before)" : "NV"
    #n_strains = ip.ins_sec_strain ? 2 : 1
    #RF = string("heatmap/results_prob_","$(replace(string(ip.β), "." => "_"))","_vac_","$(replace(string(ip.vaccine_ef), "." => "_"))","_herd_immu_","$(ip.herd)","_$strategy","cov_$(replace(string(ip.cov_val)))") ## 
    main_folder = "/data/thomas-covid/Brazil"
    if ip.initialinf == 1
        RF = string(main_folder,"/results_prob_","$(replace(string(ip.β), "." => "_"))","_vac_","$(replace(string(ip.vac_efficacy), "." => "_"))","_herd_immu_","$(ip.herd)_redp_$(ip.reduction_protection)_fmild_$(ip.fmild)_rp_$(ip.red_risk_perc)_sev_$(ip.vac_ef_sev)_$(ip.prov)_$(ip.fd_1)") ##  
    else
        RF = string(main_folder,"/results_prob_","$(replace(string(ip.β), "." => "_"))","_vac_","$(replace(string(ip.vac_efficacy), "." => "_"))","_herd_immu_","$(ip.herd)_redp_$(ip.reduction_protection)_fmild_$(ip.fmild)_rp_$(ip.red_risk_perc)_sev_$(ip.vac_ef_sev)_$(ip.prov)_$(ip.fd_1)_$(ip.initialinf)") ##  
    end

    if !Base.Filesystem.isdir(RF)
        Base.Filesystem.mkpath(RF)
    end
    return RF
end


## now, running vaccine and herd immunity, focusing and not focusing in comorbidity, first  argument turns off vac
function run_param(b,prov="saopaulo",herd_im_v = [0],fs=0.0,fm=0.0,vaccinate = false,v_e = 0.0,red = 0.0,red_perc=0.0,ef_sev=0.0,nsims=1000)
    for h_i = herd_im_v
        #bd = Dict(5=>0.074, 10=>0.076, 20=>0.089)
        #b = bd[h_i]
        @everywhere ip = ModelParameters(β = $b,
        fsevere = $fs, fmild = $fm, vaccinating = $vaccinate,
        days_before = 0, vac_efficacy = $v_e, herd = $(h_i),
        reduction_protection = $red,red_risk_perc = $red_perc,prov=Symbol($prov),vac_ef_sev=$ef_sev)

        folder = create_folder(ip)

        #println("$v_e $(ip.vaccine_ef)")
        run(ip,nsims,folder)
    end
end


## now, running vaccine and herd immunity, focusing and not focusing in comorbidity, first  argument turns off vac
function run_param_calibration(b,prov="saopaulo",herd_im_v = [0],fs=0.0,fm=0.0,vaccinate = false,v_e = 0.0,red = 0.0,red_perc=0.0,nsims=1000)
    for h_i = herd_im_v
       
        #b = bd[h_i]
        @everywhere ip = ModelParameters(β = $b,
        fsevere = $fs, fmild = $fm, vaccinating = $vaccinate,
        days_before = 0, vac_efficacy = $v_e, herd = $(h_i),
        reduction_protection = $red,red_risk_perc = $red_perc, prov=Symbol($prov), modeltime = 50)

        folder = create_folder(ip)

        #println("$v_e $(ip.vaccine_ef)")
        run(ip,nsims,folder)

        R0 = readdlm(string(folder,"/R0.dat"),header=false)[:,1]
        m = mean(R0)
        sd = std(R0)
        println("mean R0: $(m) with std: $(sd)")
    end
end



## now, running vaccine and herd immunity, focusing and not focusing in comorbidity, first  argument turns off vac
function run_param_fixed(prov="saopaulo",herd_im_v = [0],fs=0.0,fm=0.0,vaccinate = false,v_e = 0.0,red = 0.0,red_perc=0.0,ef_sev=0.0,vac_rate=300,ii=1,vef=v_e/2.0,per=21,v_d = [14;7],nsims=1000)
    for h_i = herd_im_v
        bi = Dict(5=>1, 10=>2, 20=>3,30=>4)
        b = get_beta_city(prov,bi[h_i])
        #b = bd[h_i]
        @everywhere ip = ModelParameters(β = $b,
        fsevere = $fs, fmild = $fm, vaccinating = $vaccinate,
        days_before = 0, vac_efficacy = $v_e, herd = $(h_i),
        reduction_protection = $red,red_risk_perc = $red_perc,prov=Symbol($prov),vac_ef_sev=$ef_sev,fd_1=$vac_rate,modeltime=300,initialinf=$ii,vac_efficacy_fd = $vef,sec_dose_delay=$per,vac_period=$per,days_to_protection=$v_d)

        folder = create_folder(ip)

        #println("$v_e $(ip.vaccine_ef)")
        run(ip,nsims,folder)
    end
end

function get_beta_city(prov,i)
    pp = Symbol(prov)

    ret = @match pp begin
        :saopaulo => [0.0468;0.04695;0.053;0.0607]
         :para     => [0.0443;0.0451;0.054;0.0623]
        :bahia    => [0.0434;0.0442;0.0522;0.0635]
        :riograndedosul => [0.04255;0.0464;0.0562;0.0627]
        :goias =>  [0.0403;0.0442;0.0499;0.0607]
         _ => error("no state available")
    end
    return ret[i]
end