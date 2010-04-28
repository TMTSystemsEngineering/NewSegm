function [Center_M1,Vertex_M1,nseg]=CreateBasePattern(a,Rin,Rout,Output);
% function [Center_M1,Vertex_M1,Connect,nseg]=CreateBasePattern(a,Rin,Rout,Output);
% Eric Ponslet, 2007, elp@ericandlucie.com
% Produces a base pattern (regular hexagonal grid) within sector A for a given set of geometric parameters:
%     a = side length of hexagonal cell
%     Rin = inner radius limit for center of cells 
%     Rout = outer radius limit for center of cells 
%     Output = 1 for graph and listing, 0 for no output

% Calculate vertex coordinates of hexagon of size a, in a local reference
angles=[0:5]*pi/3;
Vertices=[a*cos(angles);a*sin(angles);zeros(size(angles))];

s=a*cos(30*pi/180);  % distance from center to midpoint along edge of hexagon

ny=ceil(Rout/(cos(30*pi/180)*2*s));                 % number of possible hexagons along YM1
nt=ceil(Rout/(cos(30*pi/180)*sin(30*pi/180)*2*s));  % number of possible hexagons along TM1

% build list of center coordinates
nseg=0;
for i=0:ny  % count from center, up the YM1 axis
    for j=0:nt % count from YM1 axis, 
        center=[0;0;0]+i*[0;2*s;0]+j*[2*s*cos(30*pi/180);-2*s*sin(30*pi/180);0];   % possible center
        accept=1;
        if norm(center) < Rin, accept=0; end;                           % eliminate centers from inner hole
        if norm(center) >= Rout, accept=0; end;                         % trim M1 to outer diameter
        if center(2)*cos(30*pi/180)-center(1)*sin(30*pi/180)<a, accept=0; end;     % eliminate segments from sector F
        if accept   % if segment belongs in sector A
            nseg=nseg+1;
            Center_M1(:,nseg)=center;
            Vertex_M1(:,:,nseg)=center*ones(1,6)+Vertices;
        end;
    end;
end;

Rvtx_M1=sqrt(Vertex_M1(1,:,:).^2+Vertex_M1(2,:,:).^2);
Rmin=min(min(Rvtx_M1));
Rmax=max(max(Rvtx_M1));

if Output
    disp(['Base Pattern Creation']);
    disp(['  side length of hexagon = ' num2str(a)]);
    disp(['  min radius of array (before gaps) = ' num2str(Rmin)]);
    disp(['  max radius of array (before gaps) = ' num2str(Rmax)]);
    disp(['Created ' num2str(nseg) ' segments'])
    % 2D plot of array with segment numbers
    scrsz = get(0,'ScreenSize');
    h=figure('Position',[50 50 scrsz(3)*2/3 scrsz(4)*2/3],'Name','2D M1 Array','NumberTitle','off');
    plot(Center_M1(1,:),Center_M1(2,:),'.k');
    axis('equal');
    for i=1:nseg
        hold on;
        plot([Vertex_M1(1,:,i) Vertex_M1(1,1,i)],[Vertex_M1(2,:,i) Vertex_M1(2,1,i)],'-k');
        text(Center_M1(1,i)+0.1,Center_M1(2,i),num2str(i));
    end;
    xlabel('X_{M1} (m)');ylabel('Y_{M1} (m)');
    % plot 0/60/120 lines
    for theta=pi/6:2*pi/6:4.01*pi/6
        hold on;polar([theta theta],[0 Rmax],':k');
    end;
    % plot ID circle
    tt=[0:pi/200:pi/2];
    hold on; polar(tt,Rmin*ones(size(tt)),':k');
    % plot OD circle
    tt=[0:pi/200:pi/2];
    hold on; polar(tt,Rmax*ones(size(tt)),':k');
end;
