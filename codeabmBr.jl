using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames

@enum HEALTH SUS LAT PRE ASYMP MILD MISO INF IISO HOS ICU REC DED UNDEF

Base.@kwdef mutable struct Human
    idx::Int64 = 0 
    health::HEALTH = SUS
    swap::HEALTH = UNDEF
    sickfrom::HEALTH = UNDEF
    wentTo::HEALTH = UNDEF
    sickby::Int64 = -1
    nextday_meetcnt::Int16 = 0 ## how many contacts for a single day
    age::Int16   = 0    # in years. don't really need this but left it incase needed later
    ag::Int16   = 0
    tis::Int16   = 0   # time in state 
    exp::Int16   = 0   # max statetime
    dur::NTuple{4, Int8} = (0, 0, 0, 0)   # Order: (latents, asymps, pres, infs) TURN TO NAMED TUPS LATER
    doi::Int16   = 999   # day of infection.
    iso::Bool = false  ## isolated (limited contacts)
    isovia::Symbol = :null ## isolated via quarantine (:qu), preiso (:pi), intervention measure (:im), or contact tracing (:ct)    
    tracing::Bool = false ## are we tracing contacts for this individual?
    tracestart::Int16 = -1 ## when to start tracing, based on values sampled for x.dur
    traceend::Int16 = -1 ## when to end tracing
    tracedby::UInt32 = 0 ## is the individual traced? property represents the index of the infectious person 
    tracedxp::Int16 = 0 ## the trace is killed after tracedxp amount of days
    comorbidity::Int8 = 0 ##does the individual has any comorbidity?
    vac_status::Int8 = 0 ##
    vac_ef::Float16 = 0.0 
    got_inf::Bool = false
    herd_im::Bool = false
    hospicu::Int8 = -1
    ag_new::Int16 = -1
    ag_ct::Int16 = -1
    hcw::Bool = false
    days_vac::Int64 = -1
    vac_red::Float64 = 0.0
    first_one::Bool = false
    lethality::Float64 = 0.0
    vac_ef_sev::Float64 = 0.0
end



## default system parameters
@with_kw mutable struct ModelParameters @deftype Float64    ## use @with_kw from Parameters
    β = 0.0345       
    seasonal::Bool = false ## seasonal betas or not
    popsize::Int64 = 100000
    prov::Symbol = :saopaulo
    calibration::Bool = false
    calibration2::Bool = false 
    ignore_cal::Bool = false
    start_several_inf::Bool = false
    modeltime::Int64 = 300
    initialinf::Int64 = 1
    initialhi::Int64 = 0 ## initial herd immunity, inserts number of REC individuals
    τmild::Int64 = 0 ## days before they self-isolate for mild cases
    fmild::Float64 = 0.0  ## percent of people practice self-isolation
    fsevere::Float64 = 0.0 #
    eldq::Float64 = 0.0 ## complete isolation of elderly
    eldqag::Int8 = 5 ## default age group, if quarantined(isolated) is ag 5. 
    fpreiso::Float64 = 0.0 ## percent that is isolated at the presymptomatic stage
    tpreiso::Int64 = 0## preiso is only turned on at this time. 
    frelasymp::Float64 = 0.26 ## relative transmission of asymptomatic
    ctstrat::Int8 = 0 ## strategy 
    fctcapture::Float16 = 0.0 ## how many symptomatic people identified
    fcontactst::Float16 = 0.0 ## fraction of contacts being isolated/quarantined
    cidtime::Int8 = 0  ## time to identification (for CT) post symptom onset
    cdaysback::Int8 = 0 ## number of days to go back and collect contacts
    #vaccine_ef::Float16 = 0.0   ## change this to Float32 typemax(Float32) typemax(Float64)
    vac_com_dec_max::Float16 = 0.0 # how much the comorbidity decreases the vac eff
    vac_com_dec_min::Float16 = 0.0 # how much the comorbidity decreases the vac eff
    herd::Int8 = 0 #typemax(Int32) ~ millions
    
    hcw_vac_comp::Float64 = 0.9
    hcw_prop::Float64 = 0.03
    comor_comp::Float64 = 0.7
    eld_comp::Float64 = 0.80
    young_comp::Float64 = 0.22
    vac_period::Int64 = 21
    sec_dose_comp::Float64 = 0.7
    daily_cov::Float64 = 0.008 ####also run for 0.008 per day
    n_comor_comp::Float64 = 0.387
    vac_efficacy::Float64 = 0.52  #### 50:5:80
    vac_efficacy_fd::Float64 = vac_efficacy/2.0
    days_to_protection::Array{Int64,1} = [14;7]
    vaccinating::Bool = false
    days_before::Int64 = 0 ### six weeks of vaccination
    single_dose::Bool = false
    drop_rate::Float64 = 0.0
    fixed_cov::Float64 = 0.4
    vac_ef_sev::Float64 = 0.0
    vac_ef_sev_fd::Float64 = vac_ef_sev/2.0
    red_risk_perc::Float64 = 0.0
    reduction_protection::Float64 = 0.0
    fd_1::Int64 = 300
    fd_2::Int64 = 0
    sd1::Int64 = fd_1
    sec_dose_delay::Int64 = vac_period

    days_Rt::Array{Int64,1} = [50;100;150;200]
    priority::Bool = false
end


## constants 
const humans = Array{Human}(undef, 0) 
const p = ModelParameters()  ## setup default parameters
const agebraks = @SVector [0:4, 5:19, 20:49, 50:64, 65:100]
#const agebraks_ct = @SVector[0:4, 5:9, 10:14, 15:19, 20:24, 25:29, 30:34, 35:39, 40:44, 45:49, 50:54, 55:59, 60:64, 65:69, 70:74, 75:100]
const agebraks_dist = @SVector[0:4, 5:9, 10:14, 15:19, 20:24, 25:29, 30:34, 35:39, 40:44, 45:49, 50:54, 55:59, 60:64, 65:69, 70:74, 75:79, 80:100]

function runsim(simnum::Int64,ip::ModelParameters)
    hmatrix,hh1 = main(ip,simnum)  

    ags = [x.ag_new for x in humans] # store a vector of the age group distribution 
    all = _collectdf(hmatrix)
    spl = _splitstate(hmatrix, ags)
    ag1 = _collectdf(spl[1])
    ag2 = _collectdf(spl[2])
    ag3 = _collectdf(spl[3])
    ag4 = _collectdf(spl[4])
    ag5 = _collectdf(spl[5])
    ag6 = _collectdf(spl[6])
    insertcols!(all, 1, :sim => simnum); insertcols!(ag1, 1, :sim => simnum); insertcols!(ag2, 1, :sim => simnum); 
    insertcols!(ag3, 1, :sim => simnum); insertcols!(ag4, 1, :sim => simnum); insertcols!(ag5, 1, :sim => simnum);
    insertcols!(ag6, 1, :sim => simnum);

    R0 = zeros(Float64,size(hh1,1))

    for i = 1:size(hh1,1)
        if length(hh1[i]) > 0
            R0[i] = length(findall(k -> k.sickby in hh1[i] && k.wentTo == PRE,humans))/length(hh1[i])
        end
    end

    R0_r = zeros(Float64,size(hh1,1))

    for i = 1:size(hh1,1)
        if length(hh1[i]) > 0
            R0_r[i] = length(findall(k -> k.sickby in hh1[i],humans))/length(hh1[i])
        end
    end


    # calculating YLL
    
    pos = findall(y-> y == 11,hmatrix[:,end])

    vector_ded::Vector{Int64} = zeros(Int64,101)

    for i = pos
        x = humans[i]
        vector_ded[(x.age+1)] += 1
    end


    return (a=all, g1=ag1, g2=ag2, g3=ag3, g4=ag4, g5=ag5,g6=ag6, R0=R0, R0_r = R0_r, vector_dead = vector_ded)
end


function main(ip::ModelParameters,sim::Int64)

    Random.seed!(sim*726)

    lat_ct = Array{Int64,1}(undef,p.modeltime)
    mild_ct = Array{Int64,1}(undef,p.modeltime)
    miso_ct = Array{Int64,1}(undef,p.modeltime)
    inf_ct = Array{Int64,1}(undef,p.modeltime)
    infiso_ct = Array{Int64,1}(undef,p.modeltime)
    hos_ct = Array{Int64,1}(undef,p.modeltime)
    icu_ct = Array{Int64,1}(undef,p.modeltime)
    rec_ct = Array{Int64,1}(undef,p.modeltime)
    ded_ct = Array{Int64,1}(undef,p.modeltime)
    ## datacollection            
    # matrix to collect model state for every time step

    # reset the parameters for the simulation scenario
    reset_params(ip)  #logic: outside "ip" parameters are copied to internal "p" which is a global const and available everywhere. 

    p.popsize == 0 && error("no population size given")
    hmatrix = zeros(Int16, p.popsize, p.modeltime)
    initialize() # initialize population

    if p.start_several_inf
        herd_immu_dist()
        insert_infected(PRE, p.initialinf, 4)[1]
        
    else
        herd_immu_dist()
        insert_infected(PRE, 1, 4)[1]
        
    end
    h_init = findall(x->x.health  in (LAT,MILD,INF,PRE,ASYMP),humans)
    h_init = [h_init]
    grps = get_ag_dist()
    if p.vaccinating
        vac_ind2 = vac_selection()
        vac_ind = Array{Int64,1}(undef,length(vac_ind2))
        for i = 1:length(vac_ind2)
            vac_ind[i] = vac_ind2[i]
        end
        v1,v2 = vac_index_new(length(vac_ind))
        
        time_vac::Int64 = 1
        
        if p.days_before > 0
            tt = min(p.days_before,length(v1)-1)
            for time_vac_aux = 1:tt#p.days_before
               time_vac = time_vac_aux
                vac_ind2 = vac_time!(vac_ind,time_vac,v1,v2)
                resize!(vac_ind, length(vac_ind2))
                for i = 1:length(vac_ind2)
                    vac_ind[i] = vac_ind2[i]
                end
                for x in humans
                    vac_update(x)
                end
                if time_vac_aux%6 == 0
                    for x in humans
                        vac_update(x)
                    end
                end
            end
        end        
        for st = 1:p.modeltime
            # start of day
            #println("$st")
            if time_vac<=(length(v1)-1)
                #if st%7 > 0 #we are vaccinating everyday
                    vac_ind2 = vac_time!(vac_ind,time_vac,v1,v2)
                    #vac_ind = [vac_ind vac_ind2]
                    resize!(vac_ind, length(vac_ind2))
                    for i = 1:length(vac_ind2)
                        vac_ind[i] = vac_ind2[i]
                    end
                    time_vac += 1
                #end
            end
            #=if st == p.tpreiso ## time to introduce testing
            global  p.fpreiso = _fpreiso
            end=#
            _get_model_state(st, hmatrix) ## this datacollection needs to be at the start of the for loop
            dyntrans(st, grps)
            if st in p.days_Rt
                aux1 = findall(x->x.swap == LAT,humans)
                h_init = vcat(h_init,[aux1])
            end
            (lat_ct[st], mild_ct[st], miso_ct[st], inf_ct[st], infiso_ct[st], hos_ct[st], icu_ct[st], rec_ct[st], ded_ct[st]) = time_update()
            # end of day
        end
    else
        for st = 1:p.modeltime
            # start of day
            #println("$st")
            #=if st == p.tpreiso ## time to introduce testing
            global  p.fpreiso = _fpreiso
            end=#
            _get_model_state(st, hmatrix) ## this datacollection needs to be at the start of the for loop
            dyntrans(st, grps)
            if st in p.days_Rt
                aux1 = findall(x->x.swap == LAT,humans)
                h_init = vcat(h_init,[aux1])
            end
            (lat_ct[st], mild_ct[st], miso_ct[st], inf_ct[st], infiso_ct[st], hos_ct[st], icu_ct[st], rec_ct[st], ded_ct[st])  = time_update()
            # end of day
        end
    end

   

    return  hmatrix, h_init
end

function herd_immu_dist()
    
    age_group = [0:9,10:19,20:29,30:39,40:49,50:59,60:69,70:79,80:100]
    aux = herd_dist(p.prov)

    prop = p.herd/100

    aux = map(x-> Int(round(x*prop*p.popsize)),aux)

    for i = 1:length(age_group)
        pos = findall(x-> x.age in age_group[i],humans)
        for j in sample(pos,Int(aux[i]),replace=false)
            move_to_recovered(humans[j])
            humans[i].sickfrom = INF
            humans[i].herd_im = true
        end
    end
end
function herd_dist(prov)

    ret = @match prov begin
        :saopaulo => reverse([0.026;0.045;0.091;0.147;0.203;0.236;0.175;0.052;0.025])
        :bahia=> reverse([0.02450584913;0.04003630496;0.07513110125;0.1312020976;0.1981645825;0.2451593384;0.1795078661;0.07321500605;0.03307785397])
        :goias => reverse([0.01951902956;0.03780165433;0.07585777453;0.134037655;0.1924809824;0.2393476094;0.2071950975;0.06468821943;0.02907197792])
        :riograndedosul =>  reverse([0.02448613415;0.04601540822;0.09361775511;0.1471686999;0.1820251585;0.2251592752;0.1894463735;0.05909457302;0.03298662244])
        :para     => reverse([0.02151312882;0.04314254483;0.0822524524;0.1313377779;0.1944068899;0.2398681;0.1740141628;0.06871561312;0.04474933034])
         _=> error("no state available")
    end
    return ret   
end
## Data Collection/ Model State functions
function _get_model_state(st, hmatrix)
    # collects the model state (i.e. agent status at time st)
    for i=1:length(humans)
        hmatrix[i, st] = Int(humans[i].health)
    end    
end

function get_ag_dist() 
    # splits the initialized human pop into its age groups
    grps =  map(x -> findall(y -> y.ag_ct == x, humans), 1:length(agebraks)) 
    return grps
end
function _collectdf(hmatrix)
    ## takes the output of the humans x time matrix and processes it into a dataframe
    #_names_inci = Symbol.(["lat_inc", "mild_inc", "miso_inc", "inf_inc", "iiso_inc", "hos_inc", "icu_inc", "rec_inc", "ded_inc"])    
    #_names_prev = Symbol.(["sus", "lat", "mild", "miso", "inf", "iiso", "hos", "icu", "rec", "ded"])
    mdf_inc, mdf_prev = _get_incidence_and_prev(hmatrix)
    mdf = hcat(mdf_inc, mdf_prev)    
    _names_inc = Symbol.(string.((Symbol.(instances(HEALTH)[1:end - 1])), "_INC"))
    _names_prev = Symbol.(string.((Symbol.(instances(HEALTH)[1:end - 1])), "_PREV"))
    _names = vcat(_names_inc..., _names_prev...)
    datf = DataFrame(mdf, _names)
    insertcols!(datf, 1, :time => 1:p.modeltime) ## add a time column to the resulting dataframe
    return datf
end

function _get_incidence_and_prev(hmatrix)
    cols = instances(HEALTH)[1:end - 1] ## don't care about the UNDEF health status
    inc = zeros(Int64, p.modeltime, length(cols))
    pre = zeros(Int64, p.modeltime, length(cols))
    for i = 1:length(cols)
        inc[:, i] = _get_column_incidence(hmatrix, cols[i])
        pre[:, i] = _get_column_prevalence(hmatrix, cols[i])
    end
    return inc, pre
end

function _splitstate(hmatrix, ags)
    #split the full hmatrix into 4 age groups based on ags (the array of age group of each agent)
    #sizes = [length(findall(x -> x == i, ags)) for i = 1:4]
    matx = []#Array{Array{Int64, 2}, 1}(undef, 4)
    for i = 1:maximum(ags)#length(agebraks)
        idx = findall(x -> x == i, ags)
        push!(matx, view(hmatrix, idx, :))
    end
    return matx
end

function _get_column_prevalence(hmatrix, hcol)
    inth = Int(hcol)
    timevec = zeros(Int64, p.modeltime)
    for (i, c) in enumerate(eachcol(hmatrix))
        idx = findall(x -> x == inth, c)
        if idx !== nothing
            ps = length(c[idx])    
            timevec[i] = ps    
        end
    end
    return timevec
end

function _get_column_incidence(hmatrix, hcol)
    inth = Int(hcol)
    timevec = zeros(Int64, p.modeltime)
    for r in eachrow(hmatrix)
        idx = findfirst(x -> x == inth, r)
        if idx !== nothing 
            timevec[idx] += 1
        end
    end
    return timevec
end


function vac_selection()
    
    pos = findall(x-> humans[x].age>=20 && humans[x].age<65,1:length(humans))
    pos_hcw = sample(pos,Int(round(p.hcw_vac_comp*p.hcw_prop*p.popsize)),replace = false)
    
    for i in pos_hcw
        humans[i].hcw = true
    end

    pos_v_eld = findall(x-> humans[x].age>=80, 1:length(humans))
    pos_v_eld = sample(pos_v_eld,Int(round(p.eld_comp*length(pos_v_eld))),replace=false)

    pos_com = findall(x->humans[x].comorbidity == 1 && !(x in pos_hcw) && humans[x].age<60 && humans[x].age>=18, 1:length(humans))
    pos_com = sample(pos_com,Int(round(p.comor_comp*length(pos_com))),replace=false)


    pos_eld = findall(x-> humans[x].age>=60 && humans[x].age<=79, 1:length(humans))
    pos_eld = sample(pos_eld,Int(round(p.eld_comp*length(pos_eld))),replace=false)

    pos_n_com = findall(x->humans[x].comorbidity == 0 && !(x in pos_hcw) && humans[x].age<60 && humans[x].age>=18, 1:length(humans))
    pos_n_com = sample(pos_n_com,Int(round(length(pos_n_com))),replace=false)
    #pos_y = findall(x-> humans[x].age<18, 1:length(humans))
    #pos_y = sample(pos_y,Int(round(p.young_comp*length(pos_y))),replace=false)

    #= if p.priority
        aux1 = findall(y->humans[y].comorbidity==1,pos_eld)
        pos1 = shuffle([pos_com;pos_eld[aux1]])
        aux1 = findall(y->humans[y].comorbidity==0,pos_eld)
        pos1 = [pos1;pos_eld[aux1]]
    else =#
        #pos1 = shuffle([pos_com;pos_eld])
    
    #end
    pos1 = [pos_eld;pos_com]
    pos_h = shuffle([pos_hcw;pos_v_eld]) #in brasil the 80+ are vaccinated first
    #pos2 = shuffle([pos_n_com;pos_y])
    pos2 = shuffle(pos_n_com)

    v = [pos_h; pos1; pos2]
    
    if p.fixed_cov*p.popsize > length(v)
        error("general population compliance is not enough to reach the coverage.")
        exit(1)
    else
        aux = Int(round(p.fixed_cov*p.popsize))
        v = v[1:aux]
    end
   

    return v
end

function vac_index_new(l::Int64)


    v1 = Array{Int64,1}(undef,p.popsize);
    v2 = Array{Int64,1}(undef,p.popsize);
    n::Int64 = p.fd_2+p.sd1
    v1_aux::Bool = false
   
    
    v2_aux::Bool = false
    kk::Int64 = 2

    if p.single_dose
        for i = 1:p.popsize
            v1[i] = -1
            v2[i] = -1
           
        end
        v1[1] = 0
        while !v1_aux
            v1[kk] = v1[kk-1]+n
            if v1[kk] >= l
                v1[kk] = l
                v1_aux = true
            end
            kk += 1

        end
        a = findfirst(x-> x == l, v1)

        for i = (a+1):length(v1)
            v1[i] = -1
        end
        a = a+1
    else
        for i = 1:p.popsize
            v1[i] = -1
            v2[i] = -1
        end

        v1[1] = 0
        v2[1] = 0
        eligible::Int64 = 0
        for i = 2:(p.sec_dose_delay+1)
            v1[i] = (i-1)*p.fd_1
            v2[i] = 0
            if i > (p.vac_period+1)
                eligible = eligible+(v1[i-p.vac_period]-v1[i-p.vac_period-1])
            end
            if v1[i] >= l
                v1[i] = l
                
            end
        end

        kk = p.sec_dose_delay+2
       

        #eligible::Int64 = 0
        last_v2::Int64 = 0
        while !v1_aux || !v2_aux
        
            eligible = eligible+(v1[kk-p.vac_period]-v1[kk-p.vac_period-1])
            n1_a = v1_aux ? n : p.sd1
            v2_1 = min(n1_a,eligible-last_v2)

            v2[kk] = last_v2+v2_1
            last_v2 = v2[kk]
            n_aux = n-v2_1
           
            v1[kk] = v1[kk-1]+n_aux

            
            if v1[kk] >= l
                v1[kk] = l
                v1_aux = true
            end

            if v2[kk] >= l
                v2[kk] = l
                v2_aux = true
            end
            kk += 1

        end

        a = findfirst(x-> x == l, v1)

        for i = (a+1):length(v1)
            v1[i] = -1
        end


        a = findfirst(x-> x == -1, v2)
    end

    return v1[1:(a-1)],v2[1:(a-1)]
end 

function vac_time!(vac_ind::Array{Int64,1},t::Int64,n_1_dose::Array{Int64,1},n_2_dose::Array{Int64,1})
    
    ##first dose
    for i = (n_1_dose[t]+1):1:n_1_dose[t+1]
        x = humans[vac_ind[i]]
        if x.vac_status == 0
            if x.health in (MILD, MISO, INF, IISO, HOS, ICU, DED)
                pos = findall(k-> !(humans[vac_ind[k]].health in (MILD, MISO, INF, IISO, HOS, ICU, DED)) && k>n_1_dose[t+1],1:length(vac_ind))
                if length(pos) > 0
                    r = rand(pos)
                    aux = vac_ind[i]
                    vac_ind[i] = vac_ind[r]
                    vac_ind[r] = aux
                    x = humans[vac_ind[i]]
                    x.days_vac = 0
                    x.vac_status = 1
                end
            else
                x.days_vac = 0
                x.vac_status = 1
            end
            
        end
    end

    for i = (n_2_dose[t]+1):1:n_2_dose[t+1]
        x = humans[vac_ind[i]]

        if x.health in (MILD, MISO, INF, IISO, HOS, ICU, DED)
            
            if t != (length(n_2_dose)-1)
                vac_ind = [vac_ind; x.idx]
                n_2_dose[end] += 1
            end
            
        else
            if !x.hcw
                 drop_out_rate = [p.drop_rate;p.drop_rate;p.drop_rate] 
               #= drop_out_rate = [0;0;0] =#
                ages_drop = [17;64;999]
                age_ind = findfirst(k->k>=x.age,ages_drop)
                if rand() < (1-drop_out_rate[age_ind])#p.sec_dose_comp
                    x = humans[vac_ind[i]]
                    #red_com = p.vac_com_dec_min+rand()*(p.vac_com_dec_max-p.vac_com_dec_min)
                    #x.vac_ef = ((1-red_com)^x.comorbidity)*(p.vac_efficacy/2.0)+(p.vac_efficacy/2.0)
                    x.days_vac = 0
                    x.vac_status = 2
                end
            else
                #red_com = p.vac_com_dec_min+rand()*(p.vac_com_dec_max-p.vac_com_dec_min)
                #x.vac_ef = ((1-red_com)^x.comorbidity)*(p.vac_efficacy/2.0)+(p.vac_efficacy/2.0)
                x.days_vac = 0
                x.vac_status = 2
            end
        end
    end
    return vac_ind
end

function vac_update(x::Human)
    comm::Int64 = 0
    if x.age >= 65
        comm = 1
    else
        comm = x.comorbidity
    end

    if x.vac_status > 0 
        if x.days_vac == p.days_to_protection[x.vac_status]
            if x.vac_status == 1
                red_com = x.vac_red #p.vac_com_dec_min+rand()*(p.vac_com_dec_max-p.vac_com_dec_min)
                aux = p.single_dose ? ((1-red_com)^comm)*(p.vac_efficacy) : ((1-red_com)^comm)*p.vac_efficacy_fd
                x.vac_ef = aux
                x.vac_ef_sev = p.vac_ef_sev_fd
            else
                red_com = x.vac_red#p.vac_com_dec_min+rand()*(p.vac_com_dec_max-p.vac_com_dec_min)
                x.vac_ef = ((1-red_com)^comm)*(p.vac_efficacy-p.vac_efficacy_fd)+x.vac_ef
                x.vac_ef_sev = p.vac_ef_sev
            end
        end
        x.days_vac += 1
    end
end

function reset_params(ip::ModelParameters)
    # the p is a global const
    # the ip is an incoming different instance of parameters 
    # copy the values from ip to p. 
    for x in propertynames(p)
        setfield!(p, x, getfield(ip, x))
    end

    # resize the human array to change population size
    resize!(humans, p.popsize)
end


function comorbidity(ag::Int16)

    a = [4;17;24;29;34;39;44;49;54;59;999]
    g = findfirst(x->x>=ag,a)
    prob_f = [0.05; 0.1; 0.122;0.145;0.233;0.277;0.375;0.466;0.521;0.622;0.75]
    prob_m = [0.05; 0.1; 0.157;0.21;0.227;0.314;0.372;0.43;0.546;0.606;0.73]

    prob = rand() <= 0.5 ? prob_f : prob_m 
    com = rand() < prob[g] ? 1 : 0

    return com    
end


function initialize() 
    agedist = get_region_ag(p.prov)
    for i = 1:p.popsize
        #println(i)
        humans[i] = Human()              ## create an empty human       
        x = humans[i]
        x.idx = i 
        g = rand(agedist)
        x.age = rand(agebraks_dist[g])
        x.ag_ct = findfirst(y->x.age in y, agebraks) 
        a = [4;19;49;64;79;999]
        g = findfirst(y->y>=x.age,a)
        x.ag_new = g 
        x.ag = findfirst(y->x.age in y, agebraks) 
        x.exp = 999  ## susceptible people don't expire.
        x.dur = sample_epi_durations() # sample epi periods
        l = lethality(x.age)  
        x.lethality = (1-sev_prob[x.ag])*l/sev_prob[x.ag]+l
        if rand() < p.eldq && x.ag == p.eldqag   ## check if elderly need to be quarantined.
            x.iso = true   
            x.isovia = :qu         
        end
        x.comorbidity = comorbidity(x.age)
        x.vac_red = p.vac_com_dec_min+rand()*(p.vac_com_dec_max-p.vac_com_dec_min)
        # initialize the next day counts (this is important in initialization since dyntrans runs first)
        get_nextday_counts(x)
    end
end
## initialization functions 
function get_region_ag(prov) 
    ret = @match prov begin
       :saopaulo => Distributions.Categorical(@SVector[0.064838328175,0.069328176135,0.080580121287,0.080071059713,0.088160303817,0.091846486417,0.086191407298,0.077168015209,0.072237424864,0.066722619413,0.059039606687,0.048179472936,0.037249808232,0.026878984322,0.020750348279,0.014573338663,0.016184498553])
       :para     => Distributions.Categorical(@SVector[0.097170563818,0.103096259345,0.110297107881,0.103834811295,0.098799493632,0.093320174208,0.081423406860,0.067871328131,0.057356822952,0.046694976726,0.038808471279,0.030738086315,0.023009606452,0.017726565881,0.012576488405,0.008273918748,0.009001918072])
       :bahia    => Distributions.Categorical(@SVector[0.075614761203,0.084943710117,0.095567523960,0.094691439038,0.093056342106,0.093377311655,0.083529275291,0.070083155298,0.064466009831,0.055715219892,0.047352532720,0.038070170407,0.031157232559,0.024274401212,0.018762557158,0.012758307718,0.016580049834])
       :riograndedosul => Distributions.Categorical(@SVector[0.060217530900,0.067679428206,0.080586190539,0.081891043039,0.081439291396,0.083598740930,0.075584100100,0.069676823177,0.071102491890,0.072249965378,0.064827997268,0.054658021388,0.043481119054,0.032128416039,0.024504090124,0.017494786060,0.018879964511])
       :goias => Distributions.Categorical(@SVector[0.072931289379,0.078738123331,0.088437166669,0.088875556565,0.092298229051,0.092708470053,0.088656028494,0.078888861499,0.071788177730,0.061914911053,0.050862222317,0.040355855337,0.030953291489,0.023057109945,0.017297746023,0.011089498830,0.011147462236])
        _=> error("no state available")
    end       
    return ret  
end

function insert_infected(health, num, ag) 
    ## inserts a number of infected people in the population randomly
    ## this function should resemble move_to_inf()
    l = findall(x -> x.health == SUS && x.ag == ag, humans)
    if length(l) > 0 && num < length(l)
        h = sample(l, num; replace = false)
        @inbounds for i in h 
            x = humans[i]
            x.first_one = true
            if health == PRE 
                move_to_pre(x) ## the swap may be asymp, mild, or severe, but we can force severe in the time_update function
                x.wentTo = PRE
            elseif health == LAT 
                move_to_latent(x)
            elseif health == INF
                move_to_infsimple(x)
                x.wentTo = PRE
            elseif health == MILD
                move_to_mild(x)
                x.wentTo = PRE
            elseif health == ASYMP
                move_to_asymp(x)
                x.wentTo = ASYMP
            elseif health == REC 
                move_to_recovered(x)
            else 
                error("can not insert human of health $(health)")
            end       
            x.sickfrom = INF # this will add +1 to the INF count in _count_infectors()... keeps the logic simple in that function.    
        end
    end    
    return h
end


function time_update()
    # counters to calculate incidence
    lat=0; pre=0; asymp=0; mild=0; miso=0; inf=0; infiso=0; hos=0; icu=0; rec=0; ded=0;
    for x in humans 
        x.tis += 1 
        x.doi += 1 # increase day of infection. variable is garbage until person is latent
        if x.tis >= x.exp             
            @match Symbol(x.swap) begin
                :LAT  => begin move_to_latent(x); lat += 1; end
                :PRE  => begin move_to_pre(x); pre += 1; end
                :ASYMP => begin move_to_asymp(x); asymp += 1; end
                :MILD => begin move_to_mild(x); mild += 1; end
                :MISO => begin move_to_miso(x); miso += 1; end
                :INF  => begin move_to_inf(x); inf +=1; end    
                :IISO => begin move_to_iiso(x); infiso += 1; end
                :HOS  => begin move_to_hospicu(x); hos += 1; end 
                :ICU  => begin move_to_hospicu(x); icu += 1; end
                :REC  => begin move_to_recovered(x); rec += 1; end
                :DED  => begin move_to_dead(x); ded += 1; end
                _    => error("swap expired, but no swap set.")
            end
        end
        # run covid-19 functions for other integrated dynamics. 
        #ct_dynamics(x)

        # get the meet counts for the next day 
        get_nextday_counts(x)
        if p.vaccinating
            vac_update(x)
        end
    end
    return (lat, mild, miso, inf, infiso, hos, icu, rec, ded)
end


@inline _set_isolation(x::Human, iso) = _set_isolation(x, iso, x.isovia)
@inline function _set_isolation(x::Human, iso, via)
    # a helper setter function to not overwrite the isovia property. 
    # a person could be isolated in susceptible/latent phase through contact tracing
    # --> in which case it will follow through the natural history of disease 
    # --> if the person remains susceptible, then iso = off
    # a person could be isolated in presymptomatic phase through fpreiso
    # --> if x.iso == true from CT and x.isovia == :ct, do not overwrite
    # a person could be isolated in mild/severe phase through fmild, fsevere
    # --> if x.iso == true from CT and x.isovia == :ct, do not overwrite
    # --> if x.iso == true from PRE and x.isovia == :pi, do not overwrite
    x.iso = iso 
    x.isovia == :null && (x.isovia = via)
end

function sample_epi_durations()
    # when a person is sick, samples the 
    #lat_dist = Distributions.truncated(Gamma(3.122, 2.656),4,11.04) # truncated between 4 and 7
    lat_dist = Distributions.truncated(LogNormal(1.434, 0.661), 4, 7) # truncated between 4 and 7
    pre_dist = Distributions.truncated(Gamma(1.058, 5/2.3), 0.8, 3)#truncated between 0.8 and 3
    asy_dist = Gamma(5, 1)
    inf_dist = Gamma((3.2)^2/3.7, 3.7/3.2)

    latents = Int.(round.(rand(lat_dist)))
    pres = Int.(round.(rand(pre_dist)))
    latents = latents - pres # ofcourse substract from latents, the presymp periods
    asymps = Int.(ceil.(rand(asy_dist)))
    infs = Int.(ceil.(rand(inf_dist)))
    return (latents, asymps, pres, infs)
end


function move_to_latent(x::Human)
    ## transfers human h to the incubation period and samples the duration
    x.health = LAT
    x.doi = 0 ## day of infection is reset when person becomes latent
    x.tis = 0   # reset time in state 
    x.exp = x.dur[1] # get the latent period
    # the swap to asymptomatic is based on age group.
    # ask seyed for the references
    #asymp_pcts = (0.25, 0.25, 0.14, 0.07, 0.07)
    #symp_pcts = map(y->1-y,asymp_pcts) 
    #symp_pcts = (0.75, 0.75, 0.86, 0.93, 0.93) 
    
    #0-18 31 19 - 59 29 60+ 18 going to asymp
    symp_pcts = [0.7, 0.623, 0.672, 0.672, 0.812, 0.812] #[0.3 0.377 0.328 0.328 0.188 0.188]
        
    x.swap = rand() < (symp_pcts[x.ag_new])*(1-x.vac_ef) ? PRE : ASYMP
    x.wentTo = x.swap
    x.got_inf = true
    ## in calibration mode, latent people never become infectious.
   
end

function move_to_asymp(x::Human)
    ## transfers human h to the asymptomatic stage 
    x.health = ASYMP     
    x.tis = 0 
    x.exp = x.dur[2] # get the presymptomatic period
    x.swap = REC 
    # x.iso property remains from either the latent or presymptomatic class
    # if x.iso is true, the asymptomatic individual has limited contacts
end

function move_to_pre(x::Human)
    θ = (0.95, 0.9, 0.85, 0.6, 0.2)  # percentage of sick individuals going to mild infection stage
    x.health = PRE
    x.tis = 0   # reset time in state 
    x.exp = x.dur[3] # get the presymptomatic period
    
    if rand() < sev_prob[x.ag]*(1-x.vac_ef_sev)
        x.swap = INF
    else 
        x.swap = MILD
    end
    # calculate whether person is isolated
    rand() < p.fpreiso && _set_isolation(x, true, :pi)
end

function move_to_mild(x::Human)
    ## transfers human h to the mild infection stage for γ days
    x.health = MILD     
    x.tis = 0 
    x.exp = x.dur[4]
    x.swap = REC 
    # x.iso property remains from either the latent or presymptomatic class
    # if x.iso is true, staying in MILD is same as MISO since contacts will be limited. 
    # we still need the separation of MILD, MISO because if x.iso is false, then here we have to determine 
    # how many days as full contacts before self-isolation
    # NOTE: if need to count non-isolated mild people, this is overestimate as isolated people should really be in MISO all the time
    #   and not go through the mild compartment 
    aux = x.vac_status > 0 ? p.fmild*(1-p.red_risk_perc) : p.fmild
    if x.iso || rand() < aux#p.fmild
        x.swap = MISO  
        x.exp = p.τmild
    end
end


function move_to_miso(x::Human)
    ## transfers human h to the mild isolated infection stage for γ days
    x.health = MISO
    x.swap = REC
    x.tis = 0 
    x.exp = x.dur[4] - p.τmild  ## since tau amount of days was already spent as infectious
    _set_isolation(x, true, :mi) 
end

function move_to_iiso(x::Human)
    ## transfers human h to the sever isolated infection stage for γ days
    x.health = IISO   
    x.swap = REC
    x.tis = 0     ## reset time in state 
    x.exp = x.dur[4] - 1  ## since 1 day was spent as infectious
    _set_isolation(x, true, :mi)
end 

function move_to_infsimple(x::Human)
    ## transfers human h to the severe infection stage for γ days 
    ## simplified function for calibration/general purposes
    x.health = INF
    x.tis = 0 
    x.exp = x.dur[4]
    x.swap = REC 
    _set_isolation(x, false, :null) 
end


function move_to_inf(x::Human)
    ## transfers human h to the severe infection stage for γ days
    ## for swap, check if person will be hospitalized, selfiso, die, or recover
 
    # h = prob of hospital, c = prob of icu AFTER hospital    
    
    h = x.comorbidity == 1 ? 0.376 : 0.09
    c = x.comorbidity == 1 ? 0.33 : 0.25
    
   # death rate for severe cases.

    time_to_hospital = Int(round(rand(Uniform(3, 9)))) # duration symptom onset to hospitalization
    x.health = INF
    x.swap = UNDEF
    x.tis = 0 
    if rand() < h     # going to hospital or ICU but will spend delta time transmissing the disease with full contacts 
        x.exp = time_to_hospital    
        x.swap = rand() < c ? ICU : HOS        
    else ## no hospital for this lucky (but severe) individual 
        if rand() < x.lethality
            x.exp = x.dur[4]  
            x.swap = DED
        else 
            x.exp = x.dur[4]  
            x.swap = REC
            if x.iso || rand() < p.fsevere 
                x.exp = 1  ## 1 day isolation for severe cases     
                x.swap = IISO
            end  
        end
    end
    ## before returning, check if swap is set 
    x.swap == UNDEF && error("agent I -> ?")
end


function move_to_hospicu(x::Human)   
    #death prob taken from https://www.cdc.gov/nchs/nvss/vsrr/covid_weekly/index.htm#Comorbidities
    # on May 31th, 2020
    #= age_thres = [24;34;44;54;64;74;84;999]
    g = findfirst(y-> y >= x.age,age_thres) =#
    #= aux = [0:4, 5:19, 20:44, 45:54, 55:64, 65:74, 75:85, 85:99]
    mh = [0.001, 0.001, 0.0015, 0.0065, 0.02, 0.038, 0.0735, 0.1885]
    mc = [0.002,0.002,0.0022, 0.008, 0.022, 0.04, 0.08, 0.2] =#

    #gg = findfirst(y-> x.age in y,aux)

    psiH = Int(round(rand(Distributions.truncated(Gamma(4.5, 2.75), 8, 17))))
    psiC = Int(round(rand(Distributions.truncated(Gamma(4.5, 2.75), 8, 17)))) + 2
    muH = Int(round(rand(Distributions.truncated(Gamma(5.3, 2.1), 9, 15))))
    muC = Int(round(rand(Distributions.truncated(Gamma(5.3, 2.1), 9, 15)))) + 2

    swaphealth = x.swap 
    x.health = swaphealth ## swap either to HOS or ICU
    x.swap = UNDEF
    x.tis = 0
    _set_isolation(x, true) # do not set the isovia property here.  

    if swaphealth == HOS
        x.hospicu = 1 
        if rand() < x.lethality ## person will die in the hospital 
            x.exp = muH 
            x.swap = DED
        else 
            x.exp = psiH 
            x.swap = REC
        end        
    end
    if swaphealth == ICU
        x.hospicu = 2         
        if rand() < x.lethality ## person will die in the ICU 
            x.exp = muC
            x.swap = DED
        else 
            x.exp = psiC
            x.swap = REC
        end
    end 
    ## before returning, check if swap is set 
    x.swap == UNDEF && error("agent H -> ?")    
end


function move_to_dead(h::Human)
    # no level of alchemy will bring someone back to life. 
    h.health = DED
    h.swap = UNDEF
    h.tis = 0 
    h.exp = 999 ## stay recovered indefinitely
    h.iso = true # a dead person is isolated
    _set_isolation(h, true)  # do not set the isovia property here.  
    # isolation property has no effect in contact dynamics anyways (unless x == SUS)
end

function move_to_recovered(h::Human)
    h.health = REC
    h.swap = UNDEF
    h.tis = 0 
    h.exp = 999 ## stay recovered indefinitely
    h.iso = false ## a recovered person has ability to meet others
    _set_isolation(h, false)  # do not set the isovia property here.  
    # isolation property has no effect in contact dynamics anyways (unless x == SUS)
end


@inline function get_nextday_counts(x::Human)
    # get all people to meet and their daily contacts to recieve
    # we can sample this at the start of the simulation to avoid everyday    
    cnt = 0
    ag = x.ag_ct
    #if person is isolated, they can recieve only 3 maximum contacts
    if x.iso 
        cnt = rand(1:3) #rand() < 0.5 ? 0 : 
        #cnt = rand(nbs_shelter[ag]) ##using the contact average for shelter-in
    else 
        cnt = rand(Poisson(nbs[ag]))  # expensive operation, try to optimize
    end
    if x.health == DED 
        cnt = 0 
    end
    x.nextday_meetcnt = cnt
    return cnt
end


function dyntrans(sys_time, grps)
    totalmet = 0 # count the total number of contacts (total for day, for all INF contacts)
    totalinf = 0 # count number of new infected 
    ## find all the people infectious
    infs = findall(x -> x.health in (PRE, ASYMP, MILD, MISO, INF, IISO), humans)
    
    # go through every infectious person
    for xid in infs 
        x = humans[xid]
        xhealth = x.health
        cnts = x.nextday_meetcnt
        if cnts > 0  
            gpw = Int.(round.(cm[x.ag_ct]*cnts)) # split the counts over age groups
            for (i, g) in enumerate(gpw)
                # sample the people from each group
                meet = rand(grps[i], g)
                # go through each person
                for j in meet 
                    y = humans[j]
                    ycnt = y.nextday_meetcnt             
                    
                    if ycnt > 0 
                        y.nextday_meetcnt = y.nextday_meetcnt - 1 # remove a contact
                        totalmet += 1
                        # there is a contact to recieve
                         # tracing dynamics
                        
                        if x.tracing  
                            if y.tracedby == 0 && rand() < p.fcontactst
                                y.tracedby = x.idx
                                ct_data.totaltrace += 1 
                            end
                        end
                        
                    # tranmission dynamics
                        if  y.health == SUS && y.swap == UNDEF                  
                            beta = _get_betavalue(sys_time, xhealth)
                            if rand() < beta*(1-y.vac_ef*(1-p.reduction_protection))
                                totalinf += 1
                                y.swap = LAT
                                y.exp = y.tis   ## force the move to latent in the next time step.
                                y.sickfrom = xhealth ## stores the infector's status to the infectee's sickfrom
                                y.sickby = xid
                            end  
                        end
    
                    end
                    
                end
            end
        end        
    end
    return totalmet, totalinf
end


@inline function _get_betavalue(sys_time, xhealth) 
    bf = p.β ## baseline PRE
    # values coming from FRASER Figure 2... relative tranmissibilities of different stages.
    if xhealth == ASYMP
        bf = bf * p.frelasymp #0.26

    elseif xhealth == MILD || xhealth == MISO 
        bf = bf * 0.44

    elseif xhealth == INF || xhealth == IISO 
        bf = bf * 0.89
    elseif xhealth == PRE
        bf = bf
    else
        error("No value for beta: _get_betavalue")
    end
    return bf
end
function contact_distribution()
    
    #m_c = Array{Array{Float64,1},1}(undef,16)

   #=  m_c[1] = [0.285009075885926;0.0968481307401908;0.047287779698645;0.0319710045985053;0.0491021098576064;0.0813035055874084;0.102094611134498;0.0838073150878239;0.0508486779191276;0.0380450296078014;0.0396569205006603;0.0352981235828193;0.0234185808753338;0.0178129613645931;0.0115796954263125;0.00591647813274764]
    m_c[2] = [0.0909515812155396;0.468001246382323;0.0733506559351036;0.0216368802416714;0.0173701332995331;0.0439967678560112;0.0616140651076325;0.061971607833196;0.0536609840082811;0.029922682045757;0.0240937321756765;0.0183223482687961;0.0159003527895726;0.0103074761814226;0.00525751526290217;0.00364197139658195]
    m_c[3] = [0.0178871087370889;0.120105592489366;0.531326600647122;0.0489554515852244;0.0208374871164309;0.0247847734391823;0.0334721811079142;0.0490333276880978;0.0544452477046488;0.0356930574392555;0.0231019928307072;0.0142091048423455;0.00884644380652761;0.00734943992538906;0.00540206822040438;0.00455012242029545]
    m_c[4] = [0.0118542335451258;0.0236141759293866;0.155509163735534;0.470256898614389;0.0712956373485155;0.0388027326517064;0.0312410935385382;0.044216581249968;0.049169640464817;0.0477970854218185;0.0257218607016628;0.0136039532315842;0.00741742910243348;0.00489631385559543;0.00274767002701591;0.00185553058190873]
    m_c[5] = [0.026136367583263;0.0192199263898006;0.0232309760240073;0.160988446741435;0.26072690940203;0.123551790529315;0.0802140461412533;0.0695105813214408;0.0607370976649196;0.0702613868309756;0.0477632710807316;0.0284154125262377;0.012436627268191;0.00593929416177804;0.00599891865482513;0.00486894767979679]
    m_c[6] = [0.0524667386390997;0.0248867987542038;0.0132754188251924;0.0475147504653613;0.13710554317181;0.209120046393741;0.129899885827903;0.099986666977339;0.0823760740850051;0.0697543617106656;0.0664916545633888;0.0361952731198496;0.0178054103593703;0.00728332543707719;0.00363486676038617;0.00220318490960728]
    m_c[7] = [0.0554321526317877;0.0619423677068765;0.0456847433467982;0.0254151733562283;0.071090320277503;0.123470169831141;0.167713530064322;0.126614598594316;0.0954249361197941;0.0748732312188356;0.0609360225503205;0.0454621062746548;0.0243147832295427;0.0110525231561748;0.0056107101226935;0.0049626315190118]
    m_c[8] = [0.0528967883349672;0.075128099795076;0.064522821946784;0.0438930923311339;0.0472680663960901;0.0880864069456574;0.111809912329148;0.158390142721512;0.124264083843464;0.0811153049753075;0.0591479513247531;0.0362175083690677;0.0275615946143095;0.0162173537991792;0.00964872803844272;0.00383214423510793]
    m_c[9] = [0.0380061060861129;0.0612481872169174;0.0732482219528807;0.0619236776452987;0.0619515904729177;0.0788977782629785;0.104570136251816;0.118096072230954;0.149574298364556;0.102370605435434;0.0709723960743399;0.0290355435744051;0.0233785438436088;0.0134018789674983;0.00905802270491427;0.00426694091536707]
    m_c[10] = [0.028636659027944;0.0483050742054744;0.057737393775633;0.0871226372560299;0.0751180788975364;0.0798396801064653;0.0933903801845109;0.108836250454748;0.114494553912073;0.128377958191209;0.0860320204473161;0.041778474971005;0.0214601983165789;0.0110616373399833;0.00951782765922383;0.00829117525426996]
    m_c[11] = [0.039789538833428;0.0408127964126137;0.0543479688847952;0.0686487599023161;0.0781011200364344;0.0984341677529351;0.0886329366324754;0.0838056715620516;0.109414845555252;0.117331589464972;0.103153414191605;0.0635499885412668;0.0280001810041937;0.0111543834398679;0.00806490460057191;0.0067577331852213]
    m_c[12] = [0.0684582555869801;0.0677760435614284;0.0436287107992129;0.0583064399406976;0.0637328420154984;0.0988315554355947;0.106065193495515;0.080281177754568;0.0859540211000669;0.0747569815247282;0.09023114626485;0.0817543489060746;0.0450189525072465;0.0192294602936172;0.00947508009341298;0.00649979072050822]
    m_c[13] = [0.086730383165557;0.0801228353962362;0.0495602189517469;0.0507553546411704;0.0527656017562357;0.0801056846574907;0.0929533933187943;0.0968150175479138;0.0800858259073643;0.0668653143459611;0.0643468637617528;0.0724844379612646;0.0614646369820477;0.0351689438033478;0.0211947024076027;0.00858078539551409]
    m_c[14] = [0.075742243231936;0.106566687763789;0.0840857377496642;0.0460107964819059;0.0487844124451564;0.0650005900042552;0.0933410570254012;0.0887414076905565;0.0829450628563119;0.0471145216138226;0.0529835740431502;0.0560885459315333;0.056082586292592;0.0559705450804968;0.0275156529916792;0.0130265787977501]
    m_c[15] = [0.044952791668065;0.119269557337987;0.0921037912094059;0.080103142374512;0.0250598494255574;0.0501882180532716;0.0470223512434138;0.0828984461522972;0.0862083431217908;0.0671544115179381;0.050653868957725;0.0381428219843161;0.0672285241243526;0.0587726952377756;0.0619385620476334;0.0283026255439586]
    m_c[16] = [0.0642919880466826;0.0827469410006865;0.114924484109357;0.0925473488672617;0.0347696159592941;0.0329700561321326;0.0533911480838348;0.0768535718612446;0.0776915155675807;0.0797485159310261;0.078330069054638;0.0460364253119574;0.0352744013245568;0.0512129992327263;0.0462332916044098;0.0329776279126116]
 =#

    m_c = Array{Array{Float64,1},1}(undef,5) 
    m_c[1] = [0.285009075885926;0.176106915037341;0.405201249194266;0.098373624958813;0.035309134923653]
    m_c[2] = [0.035362349235548;0.643822691641948;0.25629375827145;0.04966993912862;0.014851261722434]
    m_c[3] = [0.042726265864658;0.163456704329509;0.652445996872462;0.119477660451184;0.021893372482187]
    m_c[4] = [0.060045136917024;0.169564663457276;0.529496374214793;0.202967142529434;0.037926682881474]
    m_c[5] = [0.063703866364509;0.266200691900669;0.388472172542373;0.161067428871984;0.120555840320466]
    

    return m_c
end
const cm = contact_distribution()

const nbs = [11.80601;21.3243189481628;15.8046374217745;12.7110727888636;6.33160612575525]
const sev_prob = map(y->1-y,[0.95, 0.9, 0.85, 0.6, 0.2])


function lethality(age::Int16)
    if p.prov == :brasil
        v=[]
        l=[]
    elseif p.prov == :para
        v=[999;
        79;
        69;
        59;
        49;
        39;
        29;
        19;
        9]
        l=[0.2420;
        0.2684;
        0.2274;
        0.1264;
        0.0733;
        0.0373;
        0.0142;
        0.0047;
        0.0064]
    elseif p.prov == :saopaulo
        v=[999;
        89;
        79;
        69;
        59;
        49;
        39;
        29;
        19;
        9]
        l=[0.379;
        0.308;
        0.177;
        0.081;
        0.027;
        0.01;
        0.004;
        0.001;
        0.001;
        0.001]
    elseif p.prov == :bahia
        v = [999;
        79;
        69;
        59;
        49;
        39;
        29;
        19;
        9;
        4;
        1]
        l = [0.2913;
        0.2455;
        0.2116;
        0.1242;
        0.0744;
        0.0344;
        0.0114;
        0.0034;
        0.0014;
        0.0003;
        0.0025]
    elseif p.prov == :riograndedosul
        v=[999;
        79;
        69;
        59;
        49;
        39;
        29;
        19;
        14]
        l=[0.2938;
        0.2874;
        0.2271;
        0.1062;
        0.0493;
        0.0225;
        0.0068;
        0.0015;
        0.0000]
    elseif p.prov == :goias
        v = [999;
        79;
        69;
        59;
        49;
        39;
        29;
        19;
        14;
        9]
        l = [0.2513;
        0.2606;
        0.2364;
        0.1267;
        0.0763;
        0.0334;
        0.0114;
        0.0017;
        0.0006;
        0.0015]
    else
        error("No region set")
    end

    g = findfirst(y-> y >= age,reverse(v))
    return reverse(l)[g]
end