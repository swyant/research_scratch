using Unitful

const Hf_a = 3.19u"Å"
const Hf_c = 5.05u"Å"

const kb = uconvert(u"eV*K^-1", Unitful.k)

function transition_rate(nu,Eij,T)
    return nu*exp(-Eij/(T*kb))
end

############ Trinkle

struct OT_data
    ν
    E_ot
    E_to
    E_oo
    E_tt
end

#Default to Trinkle data
function OT_data(; ν=17.47u"THz",
                   E_ot=1.97u"eV",
                   E_to=(1.97-0.91)u"eV",
                   E_oo=3.14u"eV",
                   E_tt=0.03u"eV")

    return OT_data(ν,E_ot, E_to, E_oo, E_tt)
end

TrinkleData = OT_data()

function λot(params::OT_data,T;units=u"ps^-1")
    nu = params.ν
    Eb = params.E_ot
    return uconvert(units,transition_rate(nu,Eb,T))
end

λot_tr(T;units=u"ps^-1") = λot(TrinkleData,T;units=units)

function λto(params::OT_data,T;units=u"ps^-1")
    nu = params.ν
    Eb = params.E_to
    return uconvert(units,transition_rate(nu,Eb,T))
end

λto_tr(T;units=u"ps^-1") = λto(TrinkleData,T;units=units)

function λtt(params::OT_data,T;units=u"ps^-1")
    nu = params.ν
    Eb = params.E_tt
    return uconvert(units,transition_rate(nu,Eb,T))
end

λtt_tr(T;units=u"ps^-1") = λtt(TrinkleData,T;units=units)

function λoo(params::OT_data,T;units=u"ps^-1")
    nu = params.ν
    Eb = params.E_oo
    return uconvert(units,transition_rate(nu,Eb,T))
end

λoo_tr(T;units=u"ps^-1") = λoo(TrinkleData,T;units=units)


function Dbasal(params::OT_data,T; units=u"Å^2*ps^-1")
    _λot = λot(params,T)
    _λto = λto(params,T)
    res = Hf_a^2*( (_λot*_λto)/(2*_λot + _λto) )
    u_res = uconvert(units,res)
    return u_res
end

Dbasal_tr(T;units=u"Å^2*ps^-1") = Dbasal(TrinkleData,T;units=units)


function Dc(params::OT_data,T; units=u"Å^2*ps^-1")
    _λot = λot(params,T)
    _λto = λto(params,T)
    _λtt = λtt(params,T)
    _λoo = λoo(params,T)
    num = _λot*(3*_λot*_λtt + 3*_λoo*_λto+2*_λoo*_λtt)
    den = (2*_λot+_λto)*(3*_λto+2*_λtt)
    res = 0.25*Hf_c^2*(num/den)
    u_res = uconvert(units,res)
    return u_res
end

Dc_tr(T;units=u"Å^2*ps^-1") = Dc(TrinkleData,T;units=units)

################ Demkov

struct OH_data
    ν_oh
    ν_ho
    ν_hh
    ν_oo
    E_oh
    E_ho
    E_hh
    E_oo
end

# Defaults are from Demkov paper
function OH_data(; ν_oh=2.5u"THz",
                   ν_ho=7.8u"THz",
                   ν_hh=5.1u"THz",
                   ν_oo=8.4u"THz",
                   E_oh=2.04u"eV",
                   E_ho=1.05u"eV",
                   E_hh=1.58u"eV", #this is very different from the Trinkle paper, likely the reason for isotropic behavior
                   E_oo=3.19u"eV",)    
    res = OH_data(ν_oh, ν_ho, ν_hh, ν_oo, E_oh, E_ho, E_hh, E_oo)
    return res
end

DemkovData = OH_data()                

function λoh(params::OH_data,T;units=u"ps^-1")
    nu = params.ν_oh
    Eb = params.E_oh
    return uconvert(units,transition_rate(nu,Eb,T))
end

λoh_dem(T;units=u"ps^-1") = λoh(DemkovData,T;units=units)

function λho(params::OH_data,T;units=u"ps^-1")
    nu = params.ν_ho
    Eb = params.E_ho
    return uconvert(units,transition_rate(nu,Eb,T))
end

λho_dem(T;units=u"ps^-1") = λho(DemkovData,T;units=units)

function λhh(params::OH_data,T;units=u"ps^-1")
    nu = params.ν_hh
    Eb = params.E_hh
    return uconvert(units,transition_rate(nu,Eb,T))
end

λhh_dem(T;units=u"ps^-1") = λhh(DemkovData,T;units=units)

function λoo(params::OH_data,T;units=u"ps^-1")
    nu = params.ν_oo
    Eb = params.E_oo
    return uconvert(units,transition_rate(nu,Eb,T))
end

λoo_dem(T;units=u"ps^-1") = λoo(DemkovData,T;units=units)

function Dbasal(params::OH_data,T; units=u"Å^2*ps^-1")
    _λoh = λoh(params,T) 
    _λho = λho(params,T)
    _λhh = λhh(params,T)
    _λoo = λoo(params,T)

    num = _λoh*(_λhh+2*_λho)
    den = 2*(_λho + _λoh)
    res = Hf_a^2*(num/den)
    u_res = uconvert(units,res)
    return u_res
end

Dbasal_dem(T;units=u"Å^2*ps^-1") = Dbasal(DemkovData,T;units=units)

function Dc(params::OH_data,T;units=u"Å^2*ps^-1")
    _λoh = λoh(params,T) 
    _λho = λho(params,T)
    _λhh = λhh(params,T)
    _λoo = λoo(params,T)

    num = 3*_λoh*(2*_λhh + _λho) + 2*_λho*_λoo
    den = 8*(_λho+_λoh)
    res = Hf_c^2*(num/den)    
    u_res = uconvert(units,res)
    return u_res
end

Dc_dem(T;units=u"Å^2*ps^-1") = Dc(DemkovData,T;units=units)



####### OHC data 

struct OHC_data
    ν_oo
    ν_oh
    ν_oc
    ν_ho
    ν_hc 
    ν_co
    ν_ch
    E_oo
    E_oh
    E_oc
    E_ho
    E_hc 
    E_co
    E_ch
end

function OHC_data(; ν_oo=11.4u"THz",
                    ν_oh=12.5u"THz",
                    ν_oc=29.7u"THz",
                    ν_ho=10.1u"THz",
                    ν_hc=16.2u"THz",
                    ν_co=19.4u"THz",
                    ν_ch=10.4u"THz",
                    E_oo=3.55u"eV",
                    E_oh=2.28u"eV",
                    E_oc=2.54u"eV",
                    E_ho=1.09u"eV",
                    E_hc=1.30u"eV",
                    E_co=0.18u"eV",
                    E_ch=0.13u"eV"
                  )

    res = OHC_data(ν_oo,ν_oh,ν_oc,ν_ho,
                   ν_hc,ν_co,ν_ch,E_oo,
                   E_oh,E_oc,E_ho,E_hc,
                   E_co,E_ch)
    return res 
end

DemkovOHCData = OHC_data()

λoo(params::OHC_data,T;units=u"ps^-1") = uconvert(units,
                                                  transition_rate(params.ν_oo,
                                                                  params.E_oo,
                                                                  T)
                                                  )
λoo_demOHC(T;units=u"ps^-1") = λoo(DemkovOHCData,T;units=units)

λoh(params::OHC_data,T;units=u"ps^-1") = uconvert(units,
                                                  transition_rate(params.ν_oh,
                                                                  params.E_oh,
                                                                  T)
                                                  )
λoh_demOHC(T;units=u"ps^-1") = λoh(DemkovOHCData,T;units=units)

λoc(params::OHC_data,T;units=u"ps^-1") = uconvert(units,
                                                  transition_rate(params.ν_oc,
                                                                  params.E_oc,
                                                                  T)
                                                  )
λoc_demOHC(T;units=u"ps^-1") = λoc(DemkovOHCData,T;units=units)

λho(params::OHC_data,T;units=u"ps^-1") = uconvert(units,
                                                  transition_rate(params.ν_ho,
                                                                  params.E_ho,
                                                                  T)
                                                  )
λho_demOHC(T;units=u"ps^-1") = λho(DemkovOHCData,T;units=units)

λhc(params::OHC_data,T;units=u"ps^-1") = uconvert(units,
                                                  transition_rate(params.ν_hc,
                                                                  params.E_hc,
                                                                  T)
                                                  )
λhc_demOHC(T;units=u"ps^-1") = λhc(DemkovOHCData,T;units=units)

λco(params::OHC_data,T;units=u"ps^-1") = uconvert(units,
                                                  transition_rate(params.ν_co,
                                                                  params.E_co,
                                                                  T)
                                                  )
λco_demOHC(T;units=u"ps^-1") = λco(DemkovOHCData,T;units=units)


λch(params::OHC_data,T;units=u"ps^-1") = uconvert(units,
                                                  transition_rate(params.ν_ch,
                                                                  params.E_ch,
                                                                  T)
                                                  )
λch_demOHC(T;units=u"ps^-1") = λch(DemkovOHCData,T;units=units)


function Dbasal(params::OHC_data,T; units=u"Å^2*ps^-1")
    _λoh = λoh(params,T) 
    _λoc = λoc(params,T)
    _λho = λho(params,T)
    _λhc = λhc(params,T)
    _λoo = λoo(params,T)
    
    res = Hf_a^2*(_λoh + 0.75*_λoc + 0.25*_λoh*_λhc/_λho + 0*_λoo)
    u_res = uconvert(units,res)
    return u_res
end

Dbasal_demOHC(T;units=u"Å^2*ps^-1") = Dbasal(DemkovOHCData,T;units=units)


function Dc(params::OHC_data,T;units=u"Å^2*ps^-1")
    _λoh = λoh(params,T) 
    _λoc = λoc(params,T)
    _λho = λho(params,T)
    _λhc = λhc(params,T)
    _λoo = λoo(params,T)

    res = Hf_c^2*(0.375*_λoh + 0*_λoc + 0.375*_λoh*_λhc/_λho + 0.25*_λoo)
    u_res = uconvert(units,res)
    return u_res
end

Dc_demOHC(T;units=u"Å^2*ps^-1") = Dc(DemkovOHCData,T;units=units)
