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

% Exercies 2. Revloves around 2D spring system
clear all
close all
% -----GIVEN-----
% Initial position of the two particles.
x1 = [0 0]; % Refrensed as particle 1. (x,y)
x2 = [1.8 0]; % - || -     particle 2. (x,y)
x3 = [0 1]; % REMOVE
masses = [1; 1]; % The masses of particle 1 and two
L = 1; % Spring rest length.
ks = 10; % Spring constant.
kd = 0.5; % Damping coefficient.w
g = 0; % NO GRAVITY.
% The particles are released from rest, i.e:
v = 0;
% The inital velocities for each particle.
v1 = [0 -v];
v2 = [0 v];
T = 4;
dt = 4e-2;
% ---------------
visualize = 0;
record = 0;
name = "path/to/save/name";

% As described in the theory above we can use the function
%      F(i,j) = -ks(abs(r(i,j))-L)*r_bar, Must be a vector!
% Note: r(i,j) is the distance between particle i and particle j

t_steps = ceil(T/dt); % Number of time steps in the simulation.
M = diag(masses); % The diagonal matrix of masses.
X_init = [x1;x2];
V_init = [v1;v2];
F = @(X,V) ForceFunction(X,V,ks,kd,L); % Anonymous function for LeapFrog
[X,V] = LeapFrog(X_init,V_init,F,M,t_steps,dt);
%timeit(@() LeapFrog(X_init,V_init,F,M,t_steps,dt))
% a) Visualize the trajectory of the system.
%
if visualize
    figure(1)
    VisualizeSpringSystem(X,[0 1;1 0],record,name)
end
% b) Calculate the energies.
[E,Ek,Es,Ep] = EnergyCalculation(X,V,masses,g,ks,L);

ts = linspace(0,T,t_steps);
% Plot energies over time.
figure(2);
PlotEnergies(E,Ek,Es,Ep,ts,kd)
% Plot displacement of the particles
%
figure(3)
subplot(2,2,[1 2])
plot(ts,X(:,:,1))
xlabel("Time ( s )")
ylabel("Displacement ( m )")
title("Displacement over time in x")
grid on;
legend("particle (1)","particle (2)",Location="best")
subplot(2,2,[3 4])
plot(ts,X(:,:,2))
legend("particle (1)","particle (2)",Location="best")
xlabel("Time ( s )")
ylabel("Displacement ( m )")
title("Displacement over time in y")
grid on
hold off

% The relative position of each particle is the length of the springs.
R = X - permute(X, [1 4 3 2]); % Relative positions, how long are the springs.
rs = vecnorm(R,2,3); % The magnitude of each spring.
spring_length = rs(:,1,2); % The spring. as seen from (1)->(2)
amps = abs(spring_length-L); % The amplitudes represent the strech of the spring.
%
% Calculate the difference in energy per occilation.
if kd==0
    [amp_peaks,t_peaks] = findpeaks(amps);
    E_diff = max(abs(E(t_peaks(2:end))-E(t_peaks(1:end-1)))./E(t_peaks(1:end-1))*100); % Difference in %
    fprintf("\nUsing %.4f results in maximum \nrelative differenceo of %.3f %%\n",dt,E_diff)
end
%

% Plot the length of the spring over time.
figure(4);
plot(ts,amps)
grid on;
title("Stretch over time")
xlabel("Time (s)")
ylabel("Length (m)")
hold on
% Analytical solution to the problem is as follows
r_analytical = @(t) exp(-0.5*t).*(0.8*cos(sqrt(19.75)*t)+0.090*sin(sqrt(19.75)*t));
plot(ts,abs(r_analytical(ts)))
legend(["Simulated","Analytical"],Location="best")
hold off
%
if kd>0
    % Find first peak which does not exceed 10% of the initial amplitude.
    % Initial energy:
    E_init = E(1);
    t_below_id = find((E<E_init*0.01)==1);
    if ~isempty(t_below_id)
        t_below = t_below_id(1)*dt;
        % Elapsed time can be found in amp_times
        fprintf("\nTime to reach 10%% of initial amplitude using the following values:\n")
        fprintf("Initial amplitude = %.3f\n",amps(1));
        fprintf("kd = %.3f\n",kd);
        fprintf("ks = %.3f\n",ks);
        fprintf("Time to reach = %.3f\n",t_below);
    else
        disp("Did not reach 10% of initial amplitude.")
    end
end
if abs(v)>0
    % d) Calculate the angular momentum
    ang_mom = AngularMomentum(X,V,masses);
    % We estimate the spring system as two point masses.
    % Find the moment of inertia at each time step.
    % Since we have equal masses we know that the axis of rotation is in the
    % middle of the spring. Resulting in two dual
    figure(4)
    plot(ts,ang_mom)
    grid on;
    title("Angular Momentum")
    xlabel("Time (s)")
    ylabel("kgm^2/s")
    hold off;
    I = sum((spring_length/2).^2.*masses',2);
    ang_freq = ang_mom./I;
    % Plotting the angular momentum vs spring length
    figure(5)
    plot(ts,ang_freq,DisplayName="Angular Frequency")
    hold on;
    plot(ts,spring_length,DisplayName="Spring Length")
    grid on;
    legend(Location="best")
    title("Angular frequency and spring length")
    xlabel("Time (s)")
    ylabel("rad/s & m")
    hold off;
end


function F_mat = ForceFunction(X,V,ks,kd,L)
    % This is the force function of the current lab exercise.
    % X has shape (NP x n_dims)
    % V has the same shape.
    % We want to return the force matrix of the same shape as X and V.
    
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
    F_mat = sum(F_tensor,3); % Sum along last channel.
                             % Last channel corresponds to each
                             % contribution from each spring
    % F_mat now have the correct shape of (NP x n_dims)
end