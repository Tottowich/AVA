function [X,X_marble,V,V_marble] = LeapFrogMarbleBounce(X_init,V_init,X_marble_init,V_marble_init,fixed,springs,M,M_marble,g,t_steps,dt);
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
    X_fixed = X_init(fixed,:);
    % Initialize the position/velocity tensor.
    X = zeros(t_steps,NP,n_dims);
    X_marble = zeros(t_steps,NM,n_dims+1);
    V = zeros(t_steps,NP,n_dims);
    V(1,:,:) = V_init;
    V_marble = zeros(t_steps,NM,n_dims);
    V_marble(1,:,:) = V_marble_init;
    % Get the position and radii of the marbles
    X_marble(1,:,:) = X_marble_init;
    radii = X_marble_init(:,end); % Shape: (NM x 1)
    M_mat = diag(M);
    M_marble_mat = diag(M_marble);
    F = @(x_net,v_net,x_marble,v_marble) ForceFunction(x_net,v_net,x_marble,v_marble,M,M_marble,g,ks,kd,L);
    X(1,:,:) = X_init; % Set initial position of the net.
    [F_net,F_marble] = F(X_init,V_init,X_marble_init,V_marble_init); % Initial force.
    v_net = V_init-F_net*dt/2; % Initialize with half Euler step.
    v_marble = V_marble_init - F_marble*dt/2;
    v_net(fixed) = 0;
    M_inv = inv(M_mat); % Compute inverse of diagonal mass matrix, i.e. 1./M.
    M_marble_inv = inv(M_marble_mat);
    for n = 1:t_steps-1
        xs = squeeze(X(n,:,:)); % Remove singleton dimension.
        xs_marble = squeeze(X_marble(n,:,1:end-1));
        if NM==1
            xs_marble = xs_marble';
        end
        [F_net,F_marble] = F(xs,v_net,xs_marble,v_marble); % Calculate the force matricies.
        % These matricies are standalone from each other.
        % (Assuming no interaction)
        % ########
        % We now have the forces acting on each particle system.
        % Except for collision forces.
        % We use the fact that 
        
        % We must now compute the forces between the bouncing marble and the
        % net.
%         marble_pos = X_marble(:,1:end-1);
%         radii = X_marble(:,end);
%         r = -(marble_pos - permute(X, [3 2 1])); % Shape (N_
%         pos_diff = vecnorm(r,2,2);
%         r_bars = r./pos_diff;
%         inter = squeeze(pos_diff)<=radii; % Logical intersection matrix.
% %         [inter_circles,inter_particles] = find(inter==1);
% %         [inter_particles,id] = unique(inter_particles,'first');
%         [inter_marble,inter_nodes] = find(inter==1);
%         [inter_nodes,id] = unique(inter_nodes,'first');
%         % This has only taken into account the spring system.
%         inter_marble = inter_marble(id);
%         for i = 1:length(id) % Could not find pairwise indexing for multidimensional array
%             n_hat(i,:) = r_bars(inter_marble(i),:,inter_nodes(i));
%         end
%         v_col_marble = v_marble(inter_marble,:); % Velocity of the colliding marbles.
%         v_col_nodes = v_net(inter_nodes,:);
%         % Compute the force during the impulse.
%         % Since we know that momentum is conserved then.
%         % ########
%         v_net = v_net+dt*M_inv*F_mat; % Calculate the next v(n+1/2).
%         x_new=xs+dt*v_net;
%         % Check if collision before proceeding.
%         r = -(pos_marble - permute(x_new, [3 2 1]));
%         pos_diff = vecnorm(r,2,2);
%         r_bars = r./pos_diff;
%         % pos_diff has Shape: (NP x N_circles)
%         % Check which particles are intersecting
%         inter = squeeze(pos_diff)<=radii;
%         [inter_circles,inter_particles] = find(inter==1);
%         % Only intersect with one circle. Using small dt should eliminate this
%         % issue but for robustness select first cirle from left to right.
%         [inter_particles,id] = unique(inter_particles,'first');
%         inter_circles = inter_circles(id);
%         % Unfortunatly I could only resort to a loop for this part...
%         n_hat = zeros(length(id),n_dims);
%         bounce = zeros(length(id),1);
%         for i = 1:length(id) % Could not find pairwise indexing for multidimensional array
%             n_hat(i,:) = r_bars(inter_circles(i),:,inter_particles(i));
%             bounce(i) = 2*(radii(inter_circles(i))-pos_diff(inter_circles(i),:,inter_particles(i)));
%         end
%         v_par = dot(v_net(inter_particles,:),n_hat,2).*n_hat;
%         v_net(inter_particles,:) = v(inter_particles,:)-2*v_par;
%         x_new(inter_particles,:) = x_new(inter_particles,:)+bounce.*n_hat;
        % Check if in contact with the net.
        r = -(xs_marble - permute(xs, [3 2 1]));
        pos_diff = vecnorm(r,2,2);
        r_bars = r./pos_diff; % Shape (NM x n_dims x NP)
        inter = squeeze(pos_diff)<=radii;
        if NM==1 % Squeeze removes wrong dimension.
            inter = inter';
        end
        [inter_marbles,inter_nodes] = find(inter==1);
        [inter_nodes,id] = unique(inter_nodes,'first');
        inter_marbles = inter_marbles(id);
        n_hat = zeros(length(id),n_dims);
        for i = 1:length(id) % Could not find pairwise indexing for multidimensional array
            n_hat(i,:) = r_bars(inter_marbles(i),:,inter_nodes(i));
%             keyboard
            F_marble(inter_marbles(i),:) = F_marble(inter_marbles(i),:) + dot(F_net(inter_nodes(i),:),n_hat(i,:),2).*n_hat(i,:);
            F_net(inter_nodes(i),:) = F_marble(inter_marbles(i),:)-dot(F_net(inter_nodes(i),:),n_hat(i,:),2).*n_hat(i,:);
        end
  

        v_net = v_net+dt*M_inv*F_net; % Calculate the next v(n+1/2).
%         keyboard
        v_marble_new = v_marble + dt*M_marble_inv*F_marble;
        x_marble_new = xs_marble+dt*v_marble_new;
        xs(inter_nodes,:) = x_marble_new(inter_marbles,:)+(radii+eps(radii))*n_hat;
        x_new=xs+dt*v_net;
%         if ~isempty(inter_marbles) && length(inter_nodes)>2
%             % Each node can only collide with one marble at a time for
%             % simplicity.
%             disp("Marble: "+inter_marbles)
%             disp("Node: "+inter_nodes)
%             keyboard
%         end
%         x_marble_new(inter_marbles) = 
        
        X(n+1,:,:) = x_new;
        X(n+1,fixed,:) = X_fixed; % The edges should be still
        V(n+1,:,:) = v_net;
        V(n+1,fixed,:) = 0;
        X_marble(n+1,:,1:end-1) = x_marble_new;
        V_marble(n+1,:,:) = v_marble_new;
        v_marble = v_marble_new;
%         n_hats{n+1} = n_hat;
    end
end

function [F_net,F_marble] = ForceFunction(X,V,X_marble,V_marble,ms,ms_marble,g,ks,kd,L)
    % This is the force function of the current lab exercise.
    %
    % INPUT
    %   X - (mat) Positions of each particle at current time step.
    %             Shape: (NP x n_dims). Dimensions in order (x,y,z)
    %
    %   V - (mat) Velocities of each particle at current time step.
    %             Shape: (NP x n_dims)
    %   ms - (vec) The masses of each particle.
    %             Shape: (NP x 1)
    %   g - (float) gravitational acceleration constant.
    %
    %   ks - (mat/float) either matrix of shape (NP x 1 x NP) or float.
    %                    If matrix then ks(i,j) indicates coefficient of
    %                    the spring between particle i and particle j. The
    %                    matrix must be symmetric to make sense.
    %   kd - (mat/float) either matrix of shape (NP x 1 x NP) or float.
    %                    If matrix then ks(i,j) indicates damping coefficient 
    %                    of the spring between particle i and particle j. The
    %                    matrix must be symmetric to make sense.
    %
    %   L - (mat/float)  either matrix of shape (NP x NP) or float.
    %                    If matrix then ks(i,j) sym indicates length of
    %                    the spring between particle i and particle j at
    %                    rest.The matrix must be symmetric to make sense.
    %

    % Create a distance and relative velocity tensors of shape (NP x n_dims x NP)
    R = X - permute(X, [3 2 1]); % Relative positions
    V_rels = V-permute(V, [3 2 1]); % Relative velocities
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
%     centers = X_marble(:,1:end-1);
%     radii = X_marble(:,end);
%     r = -(centers - permute(X, [3 2 1])); % Shape (N_
%     pos_diff = vecnorm(r,2,2);
%     r_bars = r./pos_diff;
%     inter = squeeze(pos_diff)<=radii; % Logical intersection matrix.
%     [inter_circles,inter_particles] = find(inter==1);
%     [inter_particles,id] = unique(inter_particles,'first');
    % This has only taken into account the spring system.
    % Now add gravity!
    F_g = ms*g*[0 0 -1]; % (NPx1)x(1x3) => (NPx3)
    F_net = F_net+F_g;
    F_marble = ms_marble*g*[0 0 -1];
end

