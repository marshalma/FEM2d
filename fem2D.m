function fem2D(scene)
% FEM explicit triangles

if nargin < 1
	scene = 0;
end

global video;
video = [];

switch(scene)
	case 0
		nTiles = 1; % Number of tiles in x and y
		dt = 1e-2; % time step
		tEnd = 1.0; % end time
		drawHz = 10; % refresh rate
		grav = [0 -9.81]'; % gravity
		rho = 1e0; % area density
		damping = 2.0; % viscous damping
		E = 1e2; % Young's modulus
		nu = 0.4; % Poisson's ratio
	case 1
end

% Convert to lambda and mu
% ### TODO ###
lambda = 0;
mu = 0;

% Creates triangles from a regular grid of nodes
[nodes,tris] = createSquare(nTiles);
nNodes = length(nodes);
nTris = length(tris);

% Find some nodes to fix
for k = 1 : nNodes
	if nodes(k).X(2) > 0.9
		nodes(k).fixed = true;
	else
		nodes(k).fixed = false;
	end
end

% Compute triangle mass and distribute to vertices
% ### TODO ###
for k = 1 : nTris
	tri = tris(k).nodes;
	Xa = nodes(tri(1)).X;
	Xb = nodes(tri(2)).X;
	Xc = nodes(tri(3)).X;
	triMass = 0;
	nodes(tri(1)).m = nodes(tri(1)).m + triMass;
	nodes(tri(2)).m = nodes(tri(2)).m + triMass;
	nodes(tri(3)).m = nodes(tri(3)).m + triMass;
end

% Simulation loop
t0 = -inf;
for t = 0 : dt : tEnd;
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
	
	% Integrate velocity and position
	% ### TODO ###
	
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
	axis(1.5*[-1 1 -1 1]); % Change axis limits here
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
function [nodes,tris] = createSquare(n)

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
	nodes(i).X = [x(i),y(i)]'; %#ok<*AGROW>
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
