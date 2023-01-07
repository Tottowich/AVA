function PlotFrame(t,X,A,circle_surface)
[t_steps,NP,n_dims] = size(X);
minminx = min(X(:,:,1),[],'all');
minminy = min(X(:,:,2),[],'all');
maxmaxx = max(X(:,:,1),[],'all');
maxmaxy = max(X(:,:,2),[],'all');
% xlim([minminx-1,maxmaxx+1])
xlim([min(circle_surface(:,1))-1,max(circle_surface(:,1))+1])
ylim([minminy-1,maxmaxy+1])
daspect([1,1,1]);
grid on
hold on
% Initialize the grid of springs.
x = squeeze(X(t,:,1));
y = squeeze(X(t,:,2));
springs_linespec = "k.-";
set(gcf, 'Position',  [600, 100, 800, 800])
SPRINGS = plot(graph(A),springs_linespec,'XData',x,'YData',y,'NodeLabel',{});

% Display the circular floor.
N_cirlces = length(circle_surface);
for n = 1:N_cirlces
    % Create a cirlce as a rectangle.
    % Since the rectangle has position based bottom left, and top right.
    % circle_surface contains the center of the circle, find ["bottom
    % left","top right"] for the circle.
    c = circle_surface(n,1:n_dims);
    r = circle_surface(n,end);
    pos = [c(1)-r,c(2)-r,2*r,2*r];
    circles(n) = rectangle(Curvature=[1,1],Position=pos);
    % Create the rectangle which will be a circle.
end
% hold off % Doesnt call hold off such that it can be modified easily.