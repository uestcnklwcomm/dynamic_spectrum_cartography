function  [Areset,resetinfo] = Factor_reset(A,sim,reset_mode,th)
    if ~iscell(A)
        A = {A};
    end



    [R2,R1] = size(sim);
    Rmax = max(R2,R1);

    sim_ext = zeros(Rmax,Rmax);
    sim_ext(1:R2,1:R1) = sim;

    for rr = 1:R2
        
        if sim_ext(rr,rr) >= th
            
            Areset{rr} = A{rr};

            resetinfo.flag = 0;
          
        else
            
           resetinfo.flag = 1;

            if reset_mode % 1 re-initialize
                [x,y] = size(A{1});
                Areset{rr} = randn(x,y);
            else % 0 reset to zero
                x = length(A{1});
                for xx = 1:x
                    Areset{rr}{xx} = 0;
                end
            end
        end
            

    end
        
end