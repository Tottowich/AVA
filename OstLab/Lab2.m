% Spring attached in both ends. 
% (1)---spring---(2)
% Make sure to generalize
% Use a struct or matrix.
% Symmetric Matrix, due to Newtown's 3rd law: F12 = -F21
% 
% No damping:
%   
%   F12 = ks*(r(1,2)-L), 
%   ks: Spring constant
%   L: Length of spring at rest.
%
% With damping:
%   
%   F12 =  ks*(r(1,2)-L) - kd*r'
%   kd: Damping constant
%   r': Relative velocity of (1) and (2)
%
% 2D No damping:
%
%   F12 = -ks(abs(r(1,2)-L)*r_bar, Must be a vector!
%   r_bar = r/abs(r): Unit vector of spring direction.
% 2D With damping:
%
%   F12 = -ks(abs(r(1,2)-L)*r_bar+kd*r'.r_bar*r_bar
<<<<<<< HEAD
%

% Exercies 2. Revloves around 2D spring system with no damping.

% -----GIVEN-----
% Initial position of the two particles.
x1 = [0 0]; % Refrensed as particle 1.
x2 = [1.8 0]; % - || -     particle 2.
% x3 = [0 1]; % REMOVE
masses = [1 1]; % The masses of particle 1 and two
L = 1; % Spring rest length.
ks = 10; % Spring constant.
kd = 0; % Damping coefficient, 0 since no damping.
g = 0; % NO GRAVITY.
% The particles are released from rest, i.e:
v1 = [0 0];
v2 = [0 0];
% v3 = [0 0]; % REMOVE
% ---------------

% Using the given values we can used the same methodology as in exercise 1.
% However the problem is that we must construct an accurate force function
% (recall F(X,V))

% As described in the theory above we can use the function
%      F(i,j) = -ks(abs(r(i,j))-L)*r_bar, Must be a vector!
% Note: r(i,j) is the distance between particle i and particle j
t_steps = 10;
dt = 0.01;
M = diag(masses);
X_init = [x1;x2];
V_init = [v1;v2];
F = @(X,V) ForceFunction(X,V,ks,L); % Anonymous function for LeapFrog
X = LeapFrog(X_init,V_init,F,M,t_steps,dt);
VisualizeSpringSystem(X)
function F_mat = ForceFunction(X,V,ks,L)
    % This is the force function of the current lab exercise.
    % No damping or gravity.
    % X has shape (NP x n_dims)
    % V has the same shape.
    % We want to return the force matrix of the same shape as X and V.
    
    % Create a distance tensor
    R = X - permute(X, [3 2 1]); % Shape of (NP x n_dims x NP)
    rs = vecnorm(R,2,2); % Euclidian norm on the second channel to 
                         % get the length of each spring.
                         % This can then be used to construct r_bars.
    r_bars = R./rs; % Shape - (NP x n_dims x NP), will be anti symmetric.
    % Replace NaN with zeros. 
    r_bars(isnan(r_bars))=0;
    % We now want to compute the forces according to the
    F_tensor = -ks*(rs-L).*r_bars;
    F_mat = sum(F_tensor,3); % Sum along last channel.
                                      % Last channel corresponds to each
                                      % contribution from each spring
end
=======
%
>>>>>>> 2330c1ce840d5fa5964b55a933699a3f1f30cacb
