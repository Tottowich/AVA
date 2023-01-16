function PlotEnergies(E,Ek,Es,Ep,ts,kd)
%PLOTENERGIES Plot energies of a spring system over time.
%
%   E - (vec) Vector for the total energy of the system at each time step.
%
%   Ek - (vec) Vector for the kinetic energy of the system at each time step.
%
%   Es - (vec) Vector for the ts energy of the system at each time step.
%
%   Ep - (vec) Vector for the potential energy of the system at each time step.
%
%   ts - (vec) Each time step.
%
%   kd - (float) Damping coefficient of the system
%
hold on;
plot(ts,[E,Ek,Es,Ep])
legend("Total","Kinetic","Spring","Potential",Location="best")
xlabel("Time ( s )")
ylabel("Energy ( J )")
if kd>0
    title("Energy over time in the coupled damped spring system. kd = "+kd)
else
    title("Energy over time in the coupled non damped spring system.")
end
grid on;
hold off
end

