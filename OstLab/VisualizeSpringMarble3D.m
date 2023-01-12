function  VisualizeSpringMarble3D(X,A,X_marble,record,name)
%VISUALIZECHEESE Visualize the simulated spring system'cheese'
%   Function used to animated the system of springs.
%
% Author: Theodor Jonsson
% Date: 05/01/2023
% INPUT
%
%   X - (mat) matrix of size (t_steps, NP, n_dims)
%
%   A - (mat) Adjacency matrix of the spring system.
%
%   circle_surface - (mat) Matrix of shape (N x n_dims+1). For each circle
%                          there is center (x,y,z) + radius.
%
[t_steps,NP,n_dims] = size(X);
NM = size(X_marble,2);
set(gcf,'Position',get(0,'Screensize')); % Maximize the window for quality
figure(1)
if record
    % Initialize animation file before main loop
    MOVE = VideoWriter(name+'.avi'); % Create an video struct, "MOVE.",
    % with the output file ExampleVideo.avi
    MOVE.Quality = 50; % Set the quality (0-100)
    MOVE.FrameRate = 25; % Set frames per seond to PAL standard (animation speed)
    open(MOVE);
end
set(gca,'nextplot','replacechildren');

% Set the limits to avoid window stuttering.
minminx = min([min(X_marble(1,:,1),[],'all'),min(X(1,:,1),[],'all')]);
minminy = min([min(X_marble(1,:,2),[],'all'),min(X(1,:,2),[],'all')]);
minminz = min([min(X_marble(1,:,3),[],'all'),min(X(1,:,3),[],'all')]);
maxmaxx = max([max(X_marble(1,:,1),[],'all'),max(X(1,:,1),[],'all')]);
maxmaxy = max([max(X_marble(1,:,2),[],'all'),max(X(1,:,2),[],'all')]);
maxmaxz = max([max(X_marble(1,:,3),[],'all'),max(X(1,:,3),[],'all')]);

xlim([minminx-2,maxmaxx+2]);
ylim([minminy-2,maxmaxy+2]);
zlim([minminz-2,maxmaxz+2]);
%xlim([-3,3]);
%ylim([-3,3]);
%zlim([-3,3]);
lx = xlim;
ly = ylim;
lz = zlim;
daspect([1,1,1]);
grid on
hold on
% Initialize the grid of springs.
x = squeeze(X(1,:,1));
y = squeeze(X(1,:,2));
z = squeeze(X(1,:,3));
springs_linespec = "k.-";
SPRINGS = plot(graph(A),springs_linespec,'XData',x,'YData',y,'ZData',z,'NodeLabel',{});
view(45,10)
[sx,sy,sz] = sphere(40); % Sphere for plotting
for n = 1:NM % Number of marbles
    % Create a cirlce as a rectangle.
    % Since the rectangle has position based bottom left, and top right.
    % circle_surface contains the center of the circle, find ["bottom
    % left","top right"] for the circle.
    pos = X_marble(n,n,1:n_dims);
    r = X_marble(n,n,end);
    spheres(n) = surf(sx*r+pos(1),sy*r+pos(2),sz*r+pos(3),'EdgeColor','k','EdgeAlpha',0.5,'FaceColor',[0.9,0.0,0.5],'FaceAlpha',0.2);
    % pos = [c(1)-r,c(2)-r,2*r,2*r];
    % spheres(n) = rectangle(Curvature=[1,1],Position=pos);
    % Create the rectangle which will be a circle.

end
% Begin main time loop
for t = 1:5:t_steps
    x = squeeze(X(t,:,1));
    y = squeeze(X(t,:,2));
    z = squeeze(X(t,:,3));
    % Update the grid.
    set(SPRINGS,'XData',x,'YData',y,'ZData',z)
    % Update the marbles.
    for n = 1:NM
        pos = X_marble(t,n,1:n_dims);
        r = X_marble(t,n,end);
        %xlim(2*(lx+pos(1)));
        %ylim(2*(ly+pos(2)));
        %zlim(2*(lz+pos(3)));

        set(spheres(n),'XData',sx*r+pos(1),'YData',sy*r+pos(2),'Zdata',sz*r+pos(3));
    end
    if record
        writeVideo(MOVE,getframe(gcf)); % Get a snapshot of the active figure frame
    else
        pause(0);
    end
    % End of main time loop
end
hold off
if record
    close(MOVE); % Close and save the avi-file
end

end % Function

