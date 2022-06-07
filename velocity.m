function [vbead, vmyo,fc] = velocity(pcpq,m,force,nbead,cvec,cdot,tau)

    %{
    obtain constraint forces and component velocities
    
    arguments:
        pcpq - matrix of derivatives of constraint equations w.r.t. particle positions
        m - matrix used by constraint method
        force - forces on simulation components
        nbead - number of actin beads
        cvec - l.h.s. of constraint equations (constraint equations are cvec = 0)
        cdot - time derivative of constraints
        tau - timescale over which constraints are enforced
    
    returns: 
        vbead,vmyo - vector velocities of actin beads (incl. formins) and myosin
        fc - vector of constraint forces
    %}
   
    mm = [m, pcpq'];
    temp = sparse(size(pcpq,1),size(pcpq,1));
    temp2 = [pcpq, temp];
    mm = [mm; temp2];
    h = - cvec/tau - cdot;
    fvec = [force,h]';
    qvec = mm \ fvec;
    qdot = qvec(1:numel(force));
    lambda = qvec(numel(force)+1:end);
    fc = - pcpq' * lambda;
    
    % distribute the qdot vector into vbead and vmyo
    vcat = vec2mat(qdot,3)';
    vbead = vcat(:,1:nbead);
    vmyo = vcat(:,nbead+1:end);
end