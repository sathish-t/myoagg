%% use RK1 (Euler method) to evolve r and v

if ~exist('bmmat','var')
    bmmat = bm_nbr (rbead, rmyo, ipt, rcap);
end

for itimes = 1:nstep
    % add or remove components from the ring
    [rbead,rmyo,xmat,ipt,ifor,bancf,bancm,dbead_first] = update(rbead,rmyo,bmmat,vbead,vmyo,...
            xmat,bancf,bancm,ifor,ipt,kofffor,koffmyo,koffmyp,dt,d_for,d_myo,dbead,...
            dbead_first,rsev,koffx,rxbind,rcap,l_break_sq,r,r_ring,...
            wr,binding_rate_myo2,binding_rate_myp2,binding_rate_for,vpol,binding_rate_x);
    
    % calculate actin bead matrices - vector pointing from one bead to the next,
    % and identifying which myosins interact with which bead
    rdiff = get_rdiff(rbead,ifor,ipt);
    bmmat = bm_nbr (rbead, rmyo, ipt, rcap);
    
    % calculate forces, and one part of the drag matrix
    [force, gbeadfv, gppf, bmmatpt] = get_force(rbead,rmyo,bancm,rdiff,ifor,ipt,bmmat,nhead,...
            fhead,fheadmyp,fone,nsat,kcap,kcap_p,dbead,vmyo0,kod2,r_ring,kwall,kexc,rexc,xmat,rx0,kx);
    amat = get_a(rbead, rmyo, bmmat, gppf, rdiff, ga, ifor, ipt, bmmatpt);
    
    % constrict the ring, compute constraint related quantities.
    % NOTE: in the current simulation vs_abs = 0, so no constriction.
    r_ring = r_ring - vs_abs * dt;
    [pcpq, cvec, cdot] = get_pcpq_xy(rbead,rmyo,ifor,bancf,bancm,r_ring,dbead_sq,vs_abs,dbead_first,vpol,d_for,d_myo);
    
    % calculate the second part of the drag matrix and combine
    m = get_m(ga,ifor,ipt,rbead,rmyo,bancf,bancm,bmmat,rdiff,gbeadfv,gppf,gb,gf,gm,gmsol,bmmatpt);
    m = m - amat;
    
    % calculate component velocities
    tau = tau_o_dt * dt;
    [vbead, vmyo,fc] = velocity(pcpq,m,force,size(rbead,2),cvec,cdot,tau);
    
    % evolve positions
    rbead = rbead + dt * vbead;
    rmyo = rmyo + dt * vmyo;
end