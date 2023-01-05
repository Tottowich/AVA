function M = AngularMomentum(X,V,m)
%ANGULARMOMENTUM Calculate angular momentum given system of particles.
%   
%   Given the position and velocity for a system of particles over time
%   calculates the angular momentum of the system.
%
% INPUT
%
%   X - (mat) The positions of each particle at each time step.
%             Shape: (t_steps x NP x n_dims)
%
%   V - (mat) The velocities of the particles at each time step.
%             Shape: (t_steps x NP x n_dims)
%
%   m - (vec) The mass of each particle of the system. Must be aligned with
%             the order of X and V.
%             Shape: (NP x 1)
%
[t_steps,NP,n_dims] = size(X);
if n_dims==2
    x = X(:,:,1);
    y = X(:,:,2);
    vx = V(:,:,1);
    vy = V(:,:,2);
    Ms = x.*vy-y.*vx;
    M = Ms*m;
elseif n_dims==3
    % TODO
else
    error('Can only process 2 or 3 dimensions. Got: '+n_dims)
end

