function md_activelearn(md_initialization::Function,init_param_dict,uq_metric,num_steps=1000)
    lmp, startstop = md_initialization(;init_param_dict...)

    uq_dict = initialize_uq_dict(uq_metric)
    thermo_dict = Dict("all_temps" => [],
                       "all_pes"   => [])
    uncertain_configs = []

    command(lmp, "thermo       1")
    command(lmp, "thermo_style custom step temp pe ke etotal press")
    command(lmp, "reset_timestep 0")
    command(lmp,"compute pe all pe")
    command(lmp,"compute temp all temp")

    command(lmp, "run 0")
    main_cfg_dict = extended_extract_ss_obs(lmp)
    update_thermo_record!(main_cfg_dict,thermo_dict)
    update_uq!(uq_metric,main_cfg_dict,uq_dict,uncertain_configs)

    for tstep in 1:num_steps       
        try
            run_command = "run 1"
            if !isnothing(startstop)
                run_command = run_command *""" start $(startstop[:start])"""*
                                           """ stop $(startstop[:stop])"""
            end
            command(lmp, run_command)
        catch e
            break 
        end

        main_cfg_dict = extended_extract_ss_obs(lmp)
        update_thermo_record!(main_cfg_dict,thermo_dict)
        update_uq!(uq_metric,main_cfg_dict,uq_dict,uncertain_configs) 
    end
    uq_dict,thermo_dict,uncertain_configs
end