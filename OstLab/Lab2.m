% This is the main file for exercise 2.
% 
% Author: Theodor Jonsson
% Date: 05/01/2023
%
% Theory of Lab 2.
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
%

% Exercies 2. Revloves around 2D spring system with no damping.

% -----GIVEN-----
% Initial position of the two particles.
x1 = [0 0]; % Refrensed as particle 1. (x,y)
x2 = [1.8 0]; % - || -     particle 2. (x,y)
x3 = [0 1]; % REMOVE
masses = [1; 1;]; % The masses of particle 1 and two
L = 1; % Spring rest length.
ks = 10; % Spring constant.
kd = 0.0; % Damping coefficient, 0 since no damping.
g = 0;%9.82; % NO GRAVITY.
% The particles are released from rest, i.e:
v1 = [0 0];
v2 = [0 0];
v3 = [0 0]; % REMOVE
% ---------------

% Using the given values we can used the same methodology as in exercise 1.
% However the problem is that we must construct an accurate force function
% (recall F(X,V))

% As described in the theory above we can use the function
%      F(i,j) = -ks(abs(r(i,j))-L)*r_bar, Must be a vector!
% Note: r(i,j) is the distance between particle i and particle j
t_steps = 100;
dt = 0.1;
M = diag(masses);
X_init = [x1;x2];
V_init = [v1;v2];
F = @(X,V) ForceFunction(X,V,ks,kd,L); % Anonymous function for LeapFrog
[X,V] = LeapFrog(X_init,V_init,F,M,t_steps,dt);
% a) Visualize the trajectory of the system.
%VisualizeSpringSystem(X)
% b) Calculate the energies.
[E,Ek,Es,Ep] = EnergyCalculation(X,V,masses,g,ks,L);
plot(linspace(0,dt*t_steps,t_steps),[E,Ek,Es,Ep])
legend("Total","Kinetic","Spring","Potential",Location="best")
xlabel("Time ( s )")
ylabel("Energy ( J )")
title("Energy over time in the coupled spring system.")
grid on

function F_mat = ForceFunction(X,V,ks,kd,L)
    % This is the force function of the current lab exercise.
    % No damping or gravity.
    % X has shape (NP x n_dims)
    % V has the same shape.
    % We want to return the force matrix of the same shape as X and V.
    
    % Create a distance tensor of shape (NP x n_dims x NP)
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
    F_mat = sum(F_tensor,3); % Sum along last channel.
                             % Last channel corresponds to each
                             % contribution from each spring
    % F_mat now have the correct shape of (NP x n_dims)
end