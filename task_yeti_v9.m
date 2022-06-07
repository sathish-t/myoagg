function [] = task_yeti_v9(itask)

    %{
    
    master function that performs the ring simulation.
    Outputs are saved in mat files.
    
    arguments:
        itask: an index that can be used within the simulation to set
            parameter values. Here, all parameter values are the same,
            and itask just changes the seed of the random number generator.
            
    returns:
        Nothing.
    
    %}

    % seed the rng
    seed = itask;
    rng(seed);
    
    % set parameters
    run_pos_z_cyl
    parameters_v7
    
    %% decide which file to load
    a = dir;
    
    % name of the file to load
    prefix = 'prom_';
    suffix = 'min.mat';

    % save an initial file
    save(strcat(prefix,num2str(itask),'_0',suffix))

    %% run simulation
    totaltime = 160; % 160 minutes
    for t = 0:totaltime-1
    
        % load ring configuration
        load(strcat(prefix,num2str(itask),'_',num2str(t),suffix))
        
        % set timestep, number of steps per minute, and number
        % of saves per minute
        dt = .01;
        nstep = 100*60;
        nsave = 4;
        nstep = nstep / nsave;

        % set flags to turn off/not turn off turnover processes
        % and actin loss after 10 minutes of normal ring conditions.
        myoTurnoverOn = false;
        noActinLoss = false;
        actinTurnoverOn = false;
        
        if t > 10
        
            % after 10 minutes, implement turnover conditions as flags
            % dictate. Dissociation rates and actin severing rates
            % are greatly reduced in ghosts.
            
            if ~actinTurnoverOn
                rsev = 1.594 / 1.3 / 60 / 30;
                kofffor = 0; 
            end
            
            if noActinLoss
                rsev = 0;
            end
            
            if ~myoTurnoverOn
                koffmyo = 1/4500;
                koffmyp = 0* koffmyo;
            end
            
            koffx = 3.3;
            
            % set binding rates and actin synthesis rates to zero.
            
            if ~myoTurnoverOn
                binding_rate_myo2 = 0;
                binding_rate_myp2 = 0;
            end
            
            binding_rate_x = 0;
            
            if ~actinTurnoverOn
                binding_rate_for = 0; % 0.35 in dev cell
                vpol = 0;
            end
            
        end

        % run simulation and 
        % save data to files nsave times within one cycle
        
        tsave = t+1;
        for saveCount=1:nsave
            rk1;
            fracString = num2str(tsave);
            if saveCount < nsave
                fracString = strcat(num2str(tsave-1),'p',num2str(saveCount));
            end
            save(strcat(prefix,num2str(itask),'_',fracString,suffix), ...
            'rbead','rmyo','ifor','ipt','bancf','bancm','xmat','dbead_first','r','r_ring','fc','vmyo','vbead','seed')
        end
        clear t
        if r_ring < 0
            break
        end
    end
end
