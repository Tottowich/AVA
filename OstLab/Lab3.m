% Modeling a sliding microscopic 2D object
% Author: Theodor Jonsson
% Date 5/1/2023
% Testing using vectorized implementation instead of nested loops first.
% The strategy is to construct an adjecency matrix and use that together
% with similar method as in exercise 2 to simulate the system.
% The force function will now be more complicated but it will probabily be
% significantly faster and more general.
%

% ------- GIVEN PROPERTIES -------
Nc = 3; % Number of columns of the 2D object.
Nr = 4; % Rows of columns of the 2D object.
masses = 1; % All particles have mass 1.
ks = 100;
kd = 0;
g = 1;
dt = 2e-3;
L = 1; % Evenly distributed particles => sqrt(2) on diagonal.
% --------------------------------
NP = Nc*Nr;

% ------- Set up the 2D object --------
n_dims = 2;
% The object is denoted by the matrix X
x = meshgrid(0:L:(Nc-1)/L,0:L:(Nr-1)/L);
y = flip(meshgrid(0:L:(Nr-1)/L,0:L:(Nc-1)/L)',1);
X = cat(3,x,y);
X = reshape(X,[NP n_dims]); % Flatten the matrix.
% X now has Shape (NP x n_dims) 

% Now we can construct the adjecency matrix for the list of particles.
% Which we can use as a mask to remove the forces from non connected
% particles.
A = zeros(NP,NP);
for r = 0:Nr-1
    for c = 0:Nc-1
        i = r*Nc+c+1

        if c > 0
            A(i-1,i) = 1;
            A(i,i-1) = 1;
        end
        if r > 0
            A(i-Nc,i) = 1;
            A(i,i-Nc) = 1;
        end
    end
end
image(A*256)
x = squeeze(X(:,1));
y = squeeze(X(:,2));
xy = [x y]; % xy should be (NP x n_dims)
% Update positions etc. and plot
gplot(A,xy,'b')
hold on
scatter(X(:,1),X(:,2))
hold off
axis padded





