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
for r = 0:1%Nr-1
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
figure(2)
gplot(A,xy,'b')
hold on
scatter(X(:,1),X(:,2))
hold off
axis padded
%figure(3)
%A = latticeAdjacencyMatrix(Nr,Nc)
%image(A*256)
figure(4)
mat = X;
[r, c,n_dims] = size(mat);                          % Get the matrix size
diagVec1 = repmat([ones(c-1, 1); 0], r, 1);  % Make the first diagonal vector
                                             %   (for horizontal connections)
diagVec1 = diagVec1(1:end-1);                % Remove the last value
diagVec2 = [0; diagVec1(1:(c*(r-1)))];       % Make the second diagonal vector
                                             %   (for anti-diagonal connections)
diagVec3 = ones(c*(r-1), 1);                 % Make the third diagonal vector
                                             %   (for vertical connections)
diagVec4 = diagVec2(2:end-1);                % Make the fourth diagonal vector
                                             %   (for diagonal connections)
adj = diag(diagVec1, 1)+...                  % Add the diagonals to a zero matrix
      diag(diagVec2, c-1)+...
      diag(diagVec3, c)+...
      diag(diagVec4, c+1);
adj = adj+adj.';
%image(adj*256)
gplot(adj,xy,'b')
function A = latticeAdjacencyMatrix(Nr,Nc)
  % Let Nr, Nc denote the size of the rectangular 2d grid
  % A - square adjacency matrix
  % Connect nodes (i,j) to (i+1,j)
  [i,j] = ndgrid(1:Nr-1,1:Nc);
  ind1 = sub2ind([Nr,Nc],i,j);
  ind2 = sub2ind([Nr,Nc],i+1,j);
  
  % Connect nodes (i,j) to (i,j+1)
  [i,j] = ndgrid(1:Nr,1:Nc-1);
  ind3 = sub2ind([Nr,Nc],i,j);
  ind4 = sub2ind([Nr,Nc],i,j+1);
  
  % build the global adjacency matrix
  totalnodes = Nr*(Nc-1) + (Nc-1)*Nc;
  A = sparse([ind1(:);ind3(:)],[ind2(:);ind4(:)],ones(totalnodes,1),Nr*Nc,Nr*Nc);
  
  % symmetrize, since the above computations only followed the edges in one direction.
  % that is to say, if a talks to b, then b also talks to a.
  A = A + A';
end



