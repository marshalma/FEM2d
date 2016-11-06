function bounce2D(scene, nTiles)
% scene=1: square(default), scene=2: sphere
% nTiles: number of triangles, 8 by default


% FEM explicit triangles
if nargin < 1
    scene = 1;
end

if nargin < 2
    nTiles = 8; % set nTiles to 8 by default
end

global video;
video = [];

theta = pi/8;
dt = 1e-3; % time step
tEnd = 5.0; % end time
drawHz = 100; % refresh rate
grav = [0 -10]'; % gravity
rho = 5e0; % area density
damping = 0.0; % viscous damping
E = 10000; % Young's modulus
nu = -0.0; % Poisson's ratio


% Convert to lambda and mu
% ### TODO ###
lambda = E*nu / ((1+nu)*(1-2*nu));
mu = E / (2*(1+nu));

% Creates triangles from a regular grid of nodes
if scene == 1
    [nodes,tris] = createSquare(nTiles, theta);
else
    [nodes,tris] = createSphere(nTiles);
end
nNodes = length(nodes);
nTris = length(tris);


% Compute triangle mass and distribute to vertices
% ### TODO ###
for k = 1 : nTris
	tri = tris(k).nodes;
	Xa = nodes(tri(1)).X;
	Xb = nodes(tri(2)).X;
	Xc = nodes(tri(3)).X;
    L1 = norm(Xa-Xb, 2);
    L2 = norm(Xc-Xb, 2);
    L3 = norm(Xc-Xa, 2);
    LL = (L1 + L2 + L3) / 2;
    area = sqrt(LL*(LL-L1)*(LL-L2)*(LL-L3)); % Heron's formula for calculating the area of a triangle
	triMass = rho * area / 3;
	nodes(tri(1)).m = nodes(tri(1)).m + triMass;
	nodes(tri(2)).m = nodes(tri(2)).m + triMass;
	nodes(tri(3)).m = nodes(tri(3)).m + triMass;
    
    % compute the inverse for computing deformation matrix
    Xa = nodes(tris(k).nodes(1)).X;
    Xb = nodes(tris(k).nodes(2)).X;
    Xc = nodes(tris(k).nodes(3)).X;
    tris(k).inverseMatrix = inv([Xb(1)-Xa(1), Xc(1)-Xa(1); Xb(2)-Xa(2), Xc(2)-Xa(2)]); 
    
    % material parameters stored per triangle
    tris(k).mu = mu;
    tris(k).lambda = lambda;
%     if scene == 2 && (k == 1 || k == length(tris))
%         tris(k).mu = tris(k).mu * 4;
%         tris(k).lambda = tris(k).lambda * 4;
%     end
end


% Simulation loop
t0 = -inf;
for t = 0 : dt : tEnd
	% Draw scene
	if t - t0 > 1 / drawHz
		draw(t,nodes,tris);
		t0 = t;
	end

	% Gravity force
	for k = 1 : nNodes
		nodes(k).f = nodes(k).m*grav;
	end
	
	% FEM force
	% ### TODO ###
    for k = 1 : nTris
        triNodes = tris(k).nodes;
        nodeA = nodes(triNodes(1));
        nodeB = nodes(triNodes(2));
        nodeC = nodes(triNodes(3));
        xa = nodeA.x;
        xb = nodeB.x;
        xc = nodeC.x;
        F = [xb(1)-xa(1), xc(1)-xa(1); xb(2)-xa(2), xc(2)-xa(2)] ...
            * tris(k).inverseMatrix;
        epsilon = 0.5 * (F * F' - eye(2));
        P = F*(2*tris(k).mu*epsilon + tris(k).lambda*trace(epsilon)*eye(2));
        sigma = P*F/det(F);
        tris(k).stress = sigma;
        edge1 = xb-xa;
        edge2 = xc-xb;
        edge3 = xa-xc;
        f1 = sigma * [-edge1(2) edge1(1)]';
        f2 = sigma * [-edge2(2) edge2(1)]';
        f3 = sigma * [-edge3(2) edge3(1)]';
        nodes(triNodes(1)).f = nodes(triNodes(1)).f + f1 / 2 + f3 / 2;
        nodes(triNodes(2)).f = nodes(triNodes(2)).f + f1 / 2 + f2 / 2;
        nodes(triNodes(3)).f = nodes(triNodes(3)).f + f2 / 2 + f3 / 2;
        
    end
	
	% Integrate velocity and position
	% ### TODO ###
    for k = 1 : nNodes
        if nodes(k).fixed == 1
            continue
        end
        % using implicit damping
        nodes(k).v = (nodes(k).m * nodes(k).v + dt * nodes(k).f) / (nodes(k).m + dt * damping * nodes(k).m);
        nodes(k).x = nodes(k).x + nodes(k).v * dt;
        if nodes(k).x(2) < -2
            nodes(k).x(2) = -2;
            nodes(k).v(2) = 0; %-nodes(k).v(2);
        end
    end
	
end

if ~isempty(video)
	video.close();
end

end

%%
function draw(t,nodes,tris)

global video;

if t == 0
	clf;
	xlabel('X');
	ylabel('Y');
	axis equal;
 	axis([-1.5 1.5 -2 1]); % Change axis limits here
	grid on;
	view(2);
	colormap jet;
	caxis([0 25]); % Change color limits here (comment out for auto)
	cb = colorbar;
	ylabel(cb, 'stress')
	video = VideoWriter('output','MPEG-4');
	video.open();
end
cla;
hold on;

x = [nodes.x]'; % flattened positions
f = reshape([tris.nodes],3,length(tris))'; % flattened indices
stress = reshape([tris.stress],4,length(tris)); % flattened stress
col = round(max(abs(stress)))'; % max stress entry per triangle
patch('Faces',f,'Vertices',x,'FaceVertexCData',col,'FaceColor','flat');

str = sprintf('t = %.4f', t);
title(str);
drawnow;

frame = getframe(gcf);
video.writeVideo(frame);

end




%%
function [nodes,tris] = createSphere(n)
% n+1 nodes in total
x = zeros(n+1, 2);
x(1, :) = [0 0];
univec = [0 1]';
for i = 2:n+1
    theta = i * 2 * pi / n;
    x(i, :) = ([cos(theta) -sin(theta); sin(theta) cos(theta)] * univec)';
end

nodes = [];
for i = 1 : length(x)
	nodes(i).X = x(i,:)'; %#ok<*AGROW>
	nodes(i).x = nodes(i).X;
	nodes(i).v = [0 0]';
	nodes(i).m = 0;
	nodes(i).f = [0 0]';
	nodes(i).fixed = false;
end

tris = [];
for i = 2 : n
    tris(end+1).nodes = [1, i, i+1];
end
tris(end+1).nodes = [1, 2, n+1];

for k = 1 : length(tris)
	tris(k).stress = zeros(2);
    tris(k).mu = 0;
    tris(k).lambda = 0;
end

end

function [nodes,tris] = createSquare(n, theta)

% Regular grid with center points
x = linspace(-1,1,n+1);
y = linspace(-1,1,n+1);
[X,Y] = meshgrid(x,y);
x = reshape(X,(n+1)*(n+1),1);
y = reshape(Y,(n+1)*(n+1),1);
dx = 1/n;
dy = 1/n;
xc = linspace(-1+dx,1-dx,n);
yc = linspace(-1+dy,1-dy,n);
[Xc,Yc] = meshgrid(xc,yc);
xc = reshape(Xc,n*n,1);
yc = reshape(Yc,n*n,1);
x = [x;xc];
y = [y;yc];

nodes = [];
for i = 1 : length(x)
	nodes(i).X = [cos(theta) -sin(theta); sin(theta) cos(theta)] * [x(i),y(i)]'; %#ok<*AGROW>
	nodes(i).x = nodes(i).X;
	nodes(i).v = [0 0]';
	nodes(i).m = 0;
	nodes(i).f = [0 0]';
	nodes(i).fixed = false;
end

ng = (n+1)*(n+1);
tris = [];
for i = 1 : n
	% Index of the lower left node of the ith column
	ki = (i-1)*(n+1) + 1;
	% Index of the first center node of the ith column
	kci = ng + (i-1)*n + 1;
	for j = 1 : n
		% Index of the lower left node of the jth row of the ith column
		kij = ki + j - 1;
		% Index of the center node of the jth row of the ith column
		kcij = kci + j - 1;
		% Create the four triangles
		tris(end+1).nodes = [kij,kij+(n+1),kcij];
		tris(end+1).nodes = [kij+(n+1),kij+(n+1)+1,kcij];
		tris(end+1).nodes = [kij+(n+1)+1,kij+1,kcij];
		tris(end+1).nodes = [kij+1,kij,kcij];
	end
end
for k = 1 : length(tris)
	tris(k).stress = zeros(2);
end

end
