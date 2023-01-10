function [X,X_marble,V,V_marble] = LeapFrogMarbleBounce(X_init,V_init,X_marble_init,V_marble_init,springs,M,M_marble,t_steps,dt);
%LeapFrogMarbleBounce Calculate the trajectory of the spring system using 'Leap frog'
% method
%
% Author: Theodor Jonsson
% Date: 10/01/2023
% 
%

% INPUT
%
%   X_init - (mat) Initial position of particles,
%           matrix of shape: (number of particles, number of dim)  or
%                             (NP x n_dims)
%   
%   X_marble_init - (mat) Initial position the marbles with radius of each
%             particle matrix of shape: (NM x n_dims+1)
%                            
%
%   V_init - (mat) Inital velocity of particles,
%            same size as X_init.
%
%   V_marble_init - (mat) Inital velocity of marbles,
%            same size as X_init.
%
%   t_steps - (int) Number of time steps to simulate
%
%   F - (function) Anonymous function using position matrix and velocity
%                  matrix. Return the force matrix by f=F(X,V)
%   
%   M - (mat) Diagonal matrix containing masses of corresponding particle.
%
%   M_marble - (mat) Diagonal matrix containing masses of each marble.
%
%   t_steps - (int) Number of time steps to simulate.
%
%   dt - (float) size of the time steps.
%
% OUTPUT
%
%   X - (mat) Positions of the particles for each timestep.
%             Shape: (t_steps x NP x n_dims)
%   
%   V - (mat) Velocities of the particles for each timestep.
%             Shape: (t_steps x NP x n_dims)
%
    NP = size(X_init,1); % Get the number of particles
    NM = size(X_marble_init,1); % Get the number of marbles.
    n_dims = size(X_init,2); % Get the number of dimensions.
    kd = springs.kd;
    ks = springs.ks;
    L = springs.L;
    % Initialize the position/velocity tensor.
    X = zeros(t_steps,NP,n_dims);
    X_marble = zeros(t_steps,NM,n_dims);
    V = zeros(t_steps,NP,n_dims);
    V_marble = zeros(t_steps,NP,n_dims);
    % Get the position and radii of the marbles
    X_marble(1,:,:) = X_marble_init(:,1:end-1);
    radii = X_marble_init(:,end); % Shape: (NM x 1)
    
    X(1,:,:) = X_init; % Set initial position of the net.
    F_mat = F(X_init,V_init); % Initial force.
    v = V_init-F_mat*dt/2; % Initialize with half Euler step.
    V(1,:,:) = v;
    M_inv = inv(M); % Compute inverse of diagonal mass matrix, i.e. 1./M.
    for n = 1:t_steps-1
        xs = squeeze(X(n,:,:)); % Remove singleton dimension.
        pos_marble = squeeze(X_marble(n,:,:));
        %F_mat = F(xs,v); % Calculate the force matrix.
        % ########
        R = xs - permute(xs, [3 2 1]); % Relative positions
        V_rels = v-permute(v, [3 2 1]); % Relative velocities
        rs = vecnorm(R,2,2); % Euclidian norm on the second channel to 
                             % get the length of each spring.
                             % This can then be used to construct r_bars.
        r_bars = R./rs; % Shape - (NP x n_dims x NP), will be anti symmetric.
        % Replace NaN with zeros. 
        r_bars(isnan(r_bars))=0;
        % We now want to compute the forces according to the formula (5) of the
        % lab instructions.
    
        % SPRING
        F_spring = ks.*(rs-L); % Shape (NP x 1 x NP), one spring from each 
                               % particle to another. ks couble be a matrix of
                               % shape (NP x 1 x NP) with different strengths
                               % for each spring.
        % DAMPING
        F_damping = kd.*dot(V_rels,R,2)./rs;
        F_damping(isnan(F_damping))=0;
        % Multiply with the unit vectors of each individual spring.
        F_tensor = -(F_spring+F_damping).*r_bars; % (NP x n_dims x NP)
        % The Entries F(i,:,i) should be zero since this corresponds to the
        % force asserted on particle i on particle i. Which is always zero.
        F_net = sum(F_tensor,3); % Sum along last channel.
                                 % Last channel corresponds to each
                                 % contribution from each spring
        % F_mat now have the correct shape of (NP x n_dims)
        
        % We must now compute the forces between the bouncing marble and the
        % net.
        marble_pos = X_marble(:,1:end-1);
        radii = X_marble(:,end);
        r = -(marble_pos - permute(X, [3 2 1])); % Shape (N_
        pos_diff = vecnorm(r,2,2);
        r_bars = r./pos_diff;
        inter = squeeze(pos_diff)<=radii; % Logical intersection matrix.
        [inter_circles,inter_particles] = find(inter==1);
        [inter_particles,id] = unique(inter_particles,'first');
        % This has only taken into account the spring system.
        % Now add gravity!
        if size(F_mat,2)==2
            F_g = ms*g*[0 -1]; % (NPx1)x(1x2) => (NPx2)
        else
            F_g = ms*g*[0 0 -1]; % (NPx1)x(1x3) => (NPx3)
        end
        F_mat = F_mat+F_g;
        % ########
        v = v+dt*M_inv*F_mat; % Calculate the next v(n+1/2).
        x_new=xs+dt*v;
        % Check if collision before proceeding.
        r = -(pos_marble - permute(x_new, [3 2 1]));
        pos_diff = vecnorm(r,2,2);
        r_bars = r./pos_diff;
        % pos_diff has Shape: (NP x N_circles)
        % Check which particles are intersecting
        inter = squeeze(pos_diff)<=radii;
        [inter_circles,inter_particles] = find(inter==1);
        % Only intersect with one circle. Using small dt should eliminate this
        % issue but for robustness select first cirle from left to right.
        [inter_particles,id] = unique(inter_particles,'first');
        inter_circles = inter_circles(id);
        % Unfortunatly I could only resort to a loop for this part...
        n_hat = zeros(length(id),n_dims);
        bounce = zeros(length(id),1);
        for i = 1:length(id) % Could not find pairwise indexing for multidimensional array
            n_hat(i,:) = r_bars(inter_circles(i),:,inter_particles(i));
            bounce(i) = 2*(radii(inter_circles(i))-pos_diff(inter_circles(i),:,inter_particles(i)));
        end
        v_par = dot(v(inter_particles,:),n_hat,2).*n_hat;
        v(inter_particles,:) = v(inter_particles,:)-2*v_par;
        x_new(inter_particles,:) = x_new(inter_particles,:)+bounce.*n_hat;
        X(n+1,:,:) = x_new;
        V(n+1,:,:) = v;
%         n_hats{n+1} = n_hat;
    end
end

function [F_net,F_marble] = ForceFunction

