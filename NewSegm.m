function [] = NewSegm(GeneratePlots,debug);
% function [] = NewSegm(GeneratePlots,debug);
% Calculates M1 segmentation for TMT using the approved technique (scaling parameter = 0.165, clocking but not recentering the BFRH)
%
% INPUTS FLAGS:
%   * GeneratePlots: =1 to generate plots (on screen), =0 otherwise; note that this can take substantial CPU time for a large array such as TMT
%   * debug: =1 to operate in debug mode (highly curved mirror with 3 segments per sector, overwriting all user-defined geometry definition); =0 otherwise
%
% OUTPUT FILES:
%   * NewSegm.txt: ASCII file containing all results in printable form
%   * NewSegm.mat: Binary Matlab data file containing all results stored as specified by C. Baffes for comparison with his own code
%   * OneGComponents.txt: ASCII file containing components of gravity vector in PSACRS for all segments in all sectors as a function of telescope zenith angle
%
% EXAMPLES OF EXECUTIONS:
%   >> NewSegm(0,0): runs calculations for the master parameters set below, and produces files NewSegm.txt, NewSegm.mat, OneGComponents.txt, and a number of saved plots (.fig files)
%   >> NewSegm(1,0): same as above but also produces a number of diagnostic plots of the array and its statistics on screen (can take a while for a TMT array)
%   >> NewSegm(1,1): same as above but overwrites all master parameters with a standard set that produces a small array (16*3=8 segments), with high curvature, and wide gaps, for debugging purposes
%
% This program requires the following subroutines (M-files and functions); see ">> help fname" for more info
%   * CreateBasePattern.m
%   * OpticalSurface.m
%   * OpticalLocal.m
%   * FitHexNew.m
%   * HexNew.m
%   * circum.m
%   * circumradius.m
%   * Normal.m
%   * Intersect.m
%
% Execution time: without plotting, the program runs in about 1 second in debug mode, and about 35 seconds for a full TMT array, on a 2GHz Pentium laptop.
%
% Created: Eric Ponslet, HYTEC, Inc., January 2007
%
% Updated and corrected: Curt Baffes, TMT, January 2008, 626-395-1632, cbaffes@tmt.org
% Changes:
%   - corrected calculation of gapped vertex location so that gapped vertices exist in the optical surface
%   - adjusted dZcell to generate local AAP coordinates as given in segmentation database printout dated 26-NOV-2007 (EP-EW version)
%
% Modified and Extended: Eric Ponslet, February & March 2010, elp@ericandlucie.com
% Changes:
%     Modified code to make compatible with MATLAB Version 6.1.0.450 (R12.1)
%         removed all uses of cosd, sind, and tand functions (replaced with sin/cos/tan(arg*pi/180)
%         removed use of circshift function (replaced with explicit array construction)
%         rewrote best fit regular hexagon routines to avoid use of lsqcurvefit function (replaced with more explicit LSQ using fminsearch); this results in new m-files fithexnew.m and fitnew.m
%         modified the circumscribed circle routines circum.m and circumradius.m (different grammar of call to fminsearch in older version)
%     Made changes to make code consistent with C. Baffes REL07 of the code:
%         Changed dZcell to make Z of AAP/M1 Cell iinterface point = -388mm = CAD model; dZcell was dZcell=0.36222495, now is dZcell=0.39222699; As a result, Z_PSA was -358mm, now is -388mm
%         Conic constant was -1.000953; now is -1.00095348 per OAD
%         Replaced text output of sections 2, 3, and 4 with tab-delimited output
%     Changed calculation of gapping vectors ggI and gOptI, from 3D geometry to 2D geometry (the correct approach per definition of gapping process)
%     Added section to store key results into arrays and file "NewSegm.mat" for ease of comparison with CB codes
%     Added more comments throughout
%     Added calculation of reference frames and mating information for edge sensors 
%         Two new functions: Intersect.m: finds the intersection with the optical surface of a line along a given vector, originating from a given point (">>help Intersect" for more info) 
%                            Normal.m: calculates a vector in M1CRS that is normal to the optical surface at a point on the optical surface defined by X,Y (M1CRS) (">>help Normal" for more info)

%######################################## set master parameters ########################################
tstart=cputime;
dt=cputime-tstart;
disp(['T=0' '  Setting base parameters']);
ScalingParam=0.165;   % value of scaling rule parameter selected by project on 1/11/2007 (Larry, Ben, Jerry)
a=1.432/2;                              % segment size in base pattern (before scaling or extrusion into M1 surface)
NominalRadius=1.440/2;                  % nominal segment size for SSA design purposes
Gap=2.5/1000;                           % physical gap from glass to glass (nominal) 
Chamfer=0.5*cos(45*pi/180)/1000;        % width of chamfer at top surface of segment, in projection into the XY_PSA plane
t=0.045;                                % thickness of segment glass
Rin=1.5;              % inner radius used for cropping base pattern in routine CreateBasePattern (before scaling)
Rout=14.4;            % outer radius used for cropping base pattern in routine CreateBasePattern (before scaling)
Opt_k=60.0;           % paraxial radius of curvature
Opt_K=-1.00095348;    % conic constant

% locations of assembly tooling fiducials on optical surface of segments, in PSA coordinate system                      
Fiducial_Radius=0.589250;       % radial location of fiducials, in local PSA axes
Fiducial_Angles=[90 210 330];   % angular location of fiducials relative to XYZ_PSA, in degrees
% fiducials are located on the optical surface; Z coordinates will be calculated accordingly

% locations of Cell-PSA interface nodes, in PSA coordinate system
CellInterface_Radius=0.418375;      % radial location of Cell-PSA interface points, in local PSA axes
CellInterface_Angles=[90 210 330];  % angular location of Cell-PSA interface points relative to XYZ_PSA, in degrees
CellInterface_Z=-.388000;           % Z location of Cell-PSA interface points, in local PSA axes

% locations of origins of actuator coordinate systems, in PSA coordinate system
ActOrigin_Radius=0.531;           % radial location of Cell-PSA interface points, in local PSA axes
ActOrigin_Angles=[90 210 330];    % angular location of Cell-PSA interface points relative to XYZ_PSA, in degrees
ActOrigin_Z=0;                    % Z location of Cell-PSA interface points, in local PSA axes

% orientation of actuator output shaft, relative to PSA coordinate system
ActOS_PSA=[0 0 0;0 0 0;1 1 1];    % all three actuator output shafts are parallel to Z_PSA axis 

% location of actuation center of rotation, in PSA coordinate system
ActCOR_PSA=[0;0;-0.0557];   % Z location of segment rotation point, in local PSA axes

% Edge sensor numbering and relationships to vertices
ES_Vertex=[2 3;...
           3 2;...
           3 4;...
           4 3;...
           4 5;...
           5 4;...
           5 6;...
           6 5;...
           6 1;...
           1 6;...
           1 2;...
           2 1];  % line number is ES number; 1st column has vertex closest to ES; second column has second vertex on that edge

% Distance from near (ungapped) vertex to edge sensor
ES_Dist=0.100722;


% Exagerated M1 paramaters, for debugging purposes
if debug
    disp('WARNING: DEBUGGING MODE - your geometry settings are being overwritten!');
    Opt_k=11;
    a=5/2;
    Gap=a/10;
    Rout=10;
    ES_Dist=a/4;
    NominalRadius=a;
    Chamfer=Gap/4;
    Fiducial_Radius=0.8*a;
    CellInterface_Radius=0.5*a;
    CellInterface_Z=-0.5*a;
    ActOrigin_Radius=0.6*a;
end;

% ----------------------------------- Echo key settings ---------------------------------------
disp(['    M1 - Paraxial radius of curvature = ' num2str(Opt_k,8) ' m']);
disp(['    M1 - Conic constant = ' num2str(Opt_K,8)]);
disp(['    Diameter of base pattern hexagon = ' num2str(2*a,8) ' m']);
disp(['    Scaling rule parameter, alpha = ' num2str(ScalingParam,8)]);
disp(['    Segment Gaps = ' num2str(Gap*1000,8) ' mm']);
disp(['    Chamfer Width (projected into XY_PSA) = ' num2str(Chamfer*1000,8) ' mm']);
disp(['    Diameter of Nominal Segment = ' num2str(2*NominalRadius,8) ' m']);
disp(['    Assembly tooling fiducials, relative to PSA coordinate system (on optical surface): Radius = ' num2str(Fiducial_Radius,8) ' m;  ' ...
                                                                                              'Angles = ' num2str(Fiducial_Angles(1),4) ', ' num2str(Fiducial_Angles(2),4) ', ' num2str(Fiducial_Angles(3),4) ' degrees']); 
disp(['    Cell-PSA interface nodes, relative to PSA coordinate system: Radius = ' num2str(CellInterface_Radius,8) ' m;  ' ...
                                                                     'Angles = ' num2str(CellInterface_Angles(1),4) ', ' num2str(CellInterface_Angles(2),4) ', ' num2str(CellInterface_Angles(3),4) ' degrees;  ' ...
                                                                     'Z = ' num2str(CellInterface_Z,8) ' m']);

%######################################## Main Calculation Section ######################################## 
dt=cputime-tstart;
disp(['T=' num2str(dt) '  Calculating base geometry']);
s=a*cos(30*pi/180);  % radius of inscribed circle in base pattern (i.e. distance from center of hexagon to midpoint of edge)
SectorID='ABCDEF';  % sector identification characters (1=A,... 6=F)

% Calculate coordinates of Cell-PSA interface nodes in the PSA local system
CellInterface_PSA=[CellInterface_Radius*cos(CellInterface_Angles*pi/180); CellInterface_Radius*sin(CellInterface_Angles*pi/180); CellInterface_Z*ones(size(CellInterface_Angles))];   % each column contains coordinates of one SSA support point, expressed in the PSA system

% Calculate coordinates of Actuator origins in the PSA local system
ActOrigin_PSA=[ActOrigin_Radius*cos(ActOrigin_Angles*pi/180); ActOrigin_Radius*sin(ActOrigin_Angles*pi/180); ActOrigin_Z*ones(size(ActOrigin_Angles))];   % each column contains coordinates of one SSA support point, expressed in the PSA system

% Create base pattern
[BaseCenter_M1,BaseVertex_M1,n_segments]=CreateBasePattern(a,Rin,Rout,0);   % produces hexagonal pattern in sector A, with side length 'a', extending (centers of hexagons) from 'Rin' to 'Rout'
                                                                            % see ">>help CreateBasePattern" for more info

% Apply Scaling Rule: scaled Radial Position = Rj*(1 + Alpha*(Rmax/k)^2)  / (1 + Alpha*(Rj/k)^2)
Rvtx_M1=sqrt(BaseVertex_M1(1,:,:).^2+BaseVertex_M1(2,:,:).^2);  % radii to of vertices of base pattern, before scaling
Rctr_M1=sqrt(BaseCenter_M1(1,:).^2+BaseCenter_M1(2,:).^2);      % radii to centers of hexagons in base pattern, before scaling
RmaxBP=max(max(Rvtx_M1));  % maximum radius to any vertex, before scaling
scale=1+ScalingParam*(RmaxBP/Opt_k)^2;  % common multiplier for scaling rule
for i=1:n_segments  % loop on all segments
    for j=1:6  % loop on six vertices of each segment
        UngappedVertex_M1(1,j,i)=scale*BaseVertex_M1(1,j,i)./(1+ScalingParam*Rvtx_M1(1,j,i).^2/Opt_k^2);    % scale X coordinates of vertices
        UngappedVertex_M1(2,j,i)=scale*BaseVertex_M1(2,j,i)./(1+ScalingParam*Rvtx_M1(1,j,i).^2/Opt_k^2);    % scale Y coordinates of vertices
    end;
    % calculate new position of segment centers
    Center_M1(1,i)=scale*BaseCenter_M1(1,i)./(1+ScalingParam*Rctr_M1(1,i).^2/Opt_k^2);  % scale X coordinates of centers
    Center_M1(2,i)=scale*BaseCenter_M1(2,i)./(1+ScalingParam*Rctr_M1(1,i).^2/Opt_k^2);  % scale Y coordinates of centers
end;

% recalculate polar coordinates of centers and vertices (now after scaling)
Rctr_M1=sqrt(Center_M1(1,:).^2+Center_M1(2,:).^2);
Tctr_M1=atan2(Center_M1(2,:),Center_M1(1,:));
Rvtx_M1=sqrt(UngappedVertex_M1(1,:,:).^2+UngappedVertex_M1(2,:,:).^2);
Tvtx_M1=atan2(UngappedVertex_M1(2,:,:),UngappedVertex_M1(1,:,:));

% extrude centers and vertices into optical surface
[Center_M1(3,:),junk,dZdR]=OpticalSurface(Center_M1(1,:),Center_M1(2,:),Opt_k,Opt_K);   % dZdR now contains n_segments values of tan(radial slope) at the centers of the segments
                                                                                        % see "help OpticalSurface" for more info
[UngappedVertex_M1(3,:,:),junk,dZdR_V]=OpticalSurface(UngappedVertex_M1(1,:,:),UngappedVertex_M1(2,:,:),Opt_k,Opt_K);   % dZdR_V now contains n_segments by 6 values of tan(radial slope) at the vertices of the segments
                                                                                                                        % see "help OpticalSurface" for more info

% initialize empty arrays for BFRH-related segment statistics
Radius=[];
Clocking=[];
Irregularity=[];

% Preallocate array memory for speed of execution
Vertex_PSA=zeros(3,6,n_segments);

% -------------------------------------------- START MAIN LOOP ON SEGMENTS ------------------------------------------------
% this is the main calculation loop; it performs all key calculations for all segments in sector A
dt=cputime-tstart;
disp(['T=' num2str(dt) '  Entering main calculation loop - Sector A']);
for i=1:n_segments

    % calculate SEG local coordinate system, expressed in M1 axes
    XSeg=[cos(Tctr_M1(i)) sin(Tctr_M1(i)) dZdR(i)]';  % X_seg is in RZ_M1 plane, tangent to optical surface
    XSeg=XSeg/norm(XSeg);  % normalize
    ZSeg=[-dZdR(i)*cos(Tctr_M1(i)) -dZdR(i)*sin(Tctr_M1(i)) 1]';  % Z_seg is along the normal to the optical surface
    ZSeg=ZSeg/norm(ZSeg);  % normalize
    YSeg=cross(ZSeg,XSeg);  % Y=ZxX
    RSeg_M1=[XSeg YSeg ZSeg];  % 3x3 orthonormal matrix; projects vectors from SEG to M1
   
    % project ungapped segments into XY_SEG plane
    for j=1:6
        UngappedVertex_SEG(:,j,i)=RSeg_M1'*(UngappedVertex_M1(:,j,i)-Center_M1(:,i));
    end;
    
    % calculate TEMP coord system, expressed in M1 system (rotation about Z_SEG); 
    % ZPSA is the same as ZSeg; XTEMP lies in XM1-ZM1 plane, i.e. is orthogonal to YM1
    ZTemp=ZSeg;
    YM1=[0;1;0];  %cb note: this is only used to find a vector in the XZ plane (Xtemp)
    XTemp=cross(YM1,ZSeg);  % find a vector normal to YM1, and tangent to the optical surface (normal to ZSeg)
    XTemp=XTemp/norm(XTemp);  % normalize this vector into 1xTemp
    YTemp=cross(ZTemp,XTemp);  % 1y is the vector product of 1z and 1x 
    RTemp_M1=[XTemp YTemp ZTemp];   % 3x3 orthonormal matrix; expresses TEMP system in the M1 system
    
    % project ungapped segments into XYZ_Temp
    % build projection matrix st RTemp*XYZ_Temp=XYZ_M1, hence RTemp'*XYZ_M1=XYZ_Temp, since RTemp is orthonormal
    for j=1:6
        UngappedVertex_Temp(:,j,i)=RTemp_M1'*(UngappedVertex_M1(:,j,i)-Center_M1(:,i));
    end;
    
    % implement gaps in the local XY plane (2D geometry)
    V=UngappedVertex_Temp(1:2,:,i);   % copy X and Y vertex coordinates of current segment to 2*6 matrix (for programming convenience)
                                      % Note EP Jan 2010: this used to be incorrectly done in 3D space; fortunately the resulting error should be very small 
    UIm1=V-[V(:,6) V(:,1:5)];    % edge vectors from V_i-1 to V_i, not yet normalized (6 of those in the columns of a matrix)
    UI=[V(:,2:6) V(:,1)]-V;      % edge vectors from V_i to V_i+1, not yet normalized (6 of those in the columns of a matrix)
    for j=1:6  % loop on 6 vertices
        UIm1(:,j)=UIm1(:,j)/norm(UIm1(:,j));  % normalize edge vector i-1 to i
        UI(:,j)=UI(:,j)/norm(UI(:,j));        % normalize edge vector i to i+1
    end;
    gI=UI-UIm1;  % 2D vector in XY_TEMP plane in direction of vertex move, i.e. bisecting the angle between the two adjacent edges (not yet normalized)
    for j=1:6  % loop on 6 vertices
        ggI(:,j)=gI(:,j)*(Gap/2)/norm(gI(:,j))/cos(0.5*acos(dot(UI(:,j),UIm1(:,j))));  % 2D vector from ungapped vertex to vertex of nominal glass, i.e. NOT accounting for chamfers
        gOptI(:,j)=gI(:,j)*(Gap/2+Chamfer)/norm(gI(:,j))/cos(0.5*acos(dot(UI(:,j),UIm1(:,j))));  % same for vertex of actual optical surface, i.e. taking chamfers into account
    end;
    
    Vertex_Temp(1:2,:,i)=V+ggI;   % calculate in-plane vertex coordinates of gapped segments
    Vertex_Temp(3,:,i)=OpticalLocal(Vertex_Temp(1,:,i),Vertex_Temp(2,:,i),Opt_k,Opt_K,Center_M1(:,i),XTemp,YTemp,ZTemp);    % explicit calculation of Z_Temp coordinate of gapped vector (within optical surface)
                                                                                                                            % (replaces CB iterative correction code)
    
    VertexOpt_Temp(1:2,:,i)=V+gOptI;  % calculate vertex coordinates of optical surface (i.e. accounting for chamfers
    VertexOpt_Temp(3,:,i)=OpticalLocal(VertexOpt_Temp(1,:,i),VertexOpt_Temp(2,:,i),Opt_k,Opt_K,Center_M1(:,i),XTemp,YTemp,ZTemp);    % explicit calculation of Z_Temp coordinate of gapped vector (within optical surface)
    
    % Project gapped segments back into M1 frame and into XYZ_SEG
    for j=1:6
        % commented out CB 080104, replaced above
        Vertex_M1(:,j,i)=RTemp_M1*Vertex_Temp(:,j,i)+Center_M1(:,i);         % glass vertices - back from TEMP to M1
        Vertex_SEG(:,j,i)=RSeg_M1'*(Vertex_M1(:,j,i)-Center_M1(:,i));        % glass vertices - back from M1 to SEG
        VertexOpt_M1(:,j,i)=RTemp_M1*VertexOpt_Temp(:,j,i)+Center_M1(:,i);   % optical surface vertices - back from TEMP to M1
        VertexOpt_SEG(:,j,i)=RSeg_M1'*(VertexOpt_M1(:,j,i)-Center_M1(:,i));  % optical surface vertices - back from M1 to SEG
    end;
        
    % calculate Best Fit Regular Hexagon (BFRH); this hexagon is centered at (0,0)_TEMP, and has adjustable radius and clocking angle 
    [Radius(i),Clocking(i),Irregularity(i)]=FitHexNew(Vertex_Temp(1:2,:,i),a);  % New fitting technique, because old function not available in my older version of Matlab (Jan 2010)
                                                                                % see "help FitHexNew" for more info
        
    % define PSA frame in Temp frame and recalculate segment coordinates; PSA frame is obtained by rotating around Z_TEMP, following the BFRH
    ZPSA=[0; 0; 1];   % Z_PSA axis expressed in TEMP coordinate system
    XPSA=[  cos(Clocking(i)); sin(Clocking(i)); 0]; % X_PSA is rotated from X_TEMP by the BFRH clocking angle
    YPSA=[ -sin(Clocking(i)); cos(Clocking(i)); 0]; % Y_PSA is rotated from Y_TEMP by the BFRH clocking angle
    RPSA_Temp=[XPSA YPSA ZPSA];  % build orthonormal rotation matrix
    RPSA_M1(:,:,i)=RTemp_M1*RPSA_Temp;  % calculate and store the rotation matrix from PSA to M1 systems, for each segment type
    for j=1:6  % loop on 6 vertices
        Vertex_PSA(:,j,i)=RPSA_Temp'*Vertex_Temp(:,j,i);  % calculate vertex (gapped) coordinates in PSA system
        VertexOpt_PSA(:,j,i)=RPSA_Temp'*VertexOpt_Temp(:,j,i);  % calculate vertex (gapped & chamfered) coordinates in PSA system
        UngappedVertex_PSA(:,j,i)=RPSA_M1(:,:,i)'*(UngappedVertex_M1(:,j,i)-Center_M1(:,i));  % calculate ungapped vertex coordinates in PSA system
    end;
    
    % calculate vector from segment center to M1 center, in PSA frame, and the clocking angle of the projection of that vector into XY_PSA, relative to X_PSA
    ToCenter_PSA(:,i)=RPSA_M1(:,:,i)'*(-1*Center_M1(:,i));  % express the vector pointing from O_PSA to O_M1 in the PSA system
    ToCenterAngle(i)=atan2(ToCenter_PSA(2,i),ToCenter_PSA(1,i));  % calculate clocking angle of that vector relative to X_PSA, in projection in the XY_PSA plane
        
    % calculate local (PSA) coordinates of segment fiducials on optical surface
    Fiducial_PSA(1,:,i)=Fiducial_Radius*cos(Fiducial_Angles*pi/180);
    Fiducial_PSA(2,:,i)=Fiducial_Radius*sin(Fiducial_Angles*pi/180);
    Fiducial_PSA(3,:,i)=OpticalLocal(Fiducial_PSA(1,:,i),Fiducial_PSA(2,:,i),Opt_k,Opt_K,Center_M1(:,i),RPSA_M1(:,1,i),RPSA_M1(:,2,i),RPSA_M1(:,3,i));
    
    % calculate coordinates of Cell-PSA interface nodes in M1 frame
    CellInterface_M1(:,:,i)=RPSA_M1(:,:,i)*CellInterface_PSA+Center_M1(:,i)*[1 1 1];
    
    % calculate coordinates of actuator origins and segment rotation points, in M1 frame
    ActOrigin_M1(:,:,i)=RPSA_M1(:,:,i)*ActOrigin_PSA+Center_M1(:,i)*[1 1 1];
    ActCOR_M1(:,i)=RPSA_M1(:,:,i)*ActCOR_PSA+Center_M1(:,i);
    ActLOA_M1(:,:,i)=RPSA_M1(:,:,i)*ActOS_PSA;
    
    % Calculate location and orientation of Edge Sensor reference systems (per notes on drawing TMT.M1.M1CS-INT-001, Rev A). 
    % All defining calculations are done in the M1 system; results are then converted into SEG and PSA systems.
    for j=1:12 % loop on all 12 edge sensors for this segment
        % Step C - define datum E
        DatumE=UngappedVertex_M1(:,ES_Vertex(j,2),i)-UngappedVertex_M1(:,ES_Vertex(j,1),i);  % vector from vertex nearest edge sensor #j to other vertex on same edge (ungapped)
        DatumE=DatumE/norm(DatumE);  % normalize to unit length
        % Step D - define line normal to optical surface at near vertex (note that near vertex is in the optical surface)
        Normalnv=Normal(UngappedVertex_M1(1,ES_Vertex(j,1),i),UngappedVertex_M1(2,ES_Vertex(j,1),i),Opt_k,Opt_K);  % calculates a vector in M1CRS that is normal to the optical surface 
                                                                                                                   % at a point on the optical surface defined by its X,Y coordinates (M1CRS) (">>help Normal" for more info)
        % Step E - calculate unit vector normal to plane containing datum E and the normal to the optical surface at the near vertex
        PlaneNormal=cross(Normalnv,DatumE);  % normal to plane, points in same general direction as ESCRS_Y
        PlaneNormal=PlaneNormal/norm(PlaneNormal);  % normalize to unit length
        % Step F
        PointDatumE=UngappedVertex_M1(:,ES_Vertex(j,1),i)+DatumE*ES_Dist;  % coordinates of a point on datum E, at prescribed distance from near vertex
        ConstLine=cross(DatumE,PlaneNormal);  % vector in the construction plane, normal to datum E (points toward stars)
        ConstLine=ConstLine/norm(ConstLine);  % normalize to unit length
        [Pj,dist]=Intersect(PointDatumE,ConstLine,Opt_k,Opt_K);  % find point in optical surface, along ConstLine, from point PointDatumE (">> help Intersect" for more info)
        % Step G - define datum D, the normal to the optical surface at point Pj 
        DatumD=Normal(Pj(1),Pj(2),Opt_k,Opt_K);  % calculates a vector in M1CRS that is normal to the optical surface 
                                                 % at a point on the optical surface defined by its X,Y coordinates (M1CRS) (">>help Normal" for more info)
        % Step H - define datum F
        DatumF=cross(DatumD,DatumE);  % calculate normal to both DatumD and Datum E
        % define Edge Sensor Coordinate system in M1CRS
        ESOrigin_M1(:,i,j)=Pj;  % origin is at point Pj
        RES_M1(:,3,i,j)=DatumD;  % Z axis is along local normal to optical surface
        RES_M1(:,2,i,j)=DatumF;  % Y axis is along datum F
        RES_M1(:,1,i,j)=cross(RES_M1(:,2,i,j),RES_M1(:,3,i,j));  % X axis completes right handed system, i.e. equal to 1Yx1Z
        % calculate Edge Sensor Coordinate systems in PSACRS
        ESOrigin_PSA(:,i,j)=RPSA_M1(:,:,i)'*(ESOrigin_M1(:,i,j)-Center_M1(:,i));
        RES_PSA(:,:,i,j)=RPSA_M1(:,:,i)'*RES_M1(:,:,i,j);
    end;  % end loop on 12 edges sensors for this segment 
    
end;
% ---------------------------------------------- END OF MAIN LOOP ON SEGMENTS ---------------------------------------------


%######################################## Propagate some data to the other five sectors ########################################
dt=cputime-tstart;
disp(['T=' num2str(dt) '  Replicating data to sectors B through F']);
% Propagate some data to the other six sectors
% segment numbering across all sectors is: i + (j-1)*n_segments, where i is segment number in sector A, and j=1..6 for sectors A..F
for j=2:6 % loop on sector number 2 to 6
    Rangle=(j-1)*pi/3; % rotation angle from sector #1 (A) to sector #j 
    R=[cos(Rangle) -sin(Rangle) 0; sin(Rangle) cos(Rangle) 0; 0  0  1]; % rotation matrix from sector #1 (A) to sector #j, in M1_CRS (identity matrix if j=1) 
    for i=1:n_segments  % for each segment type, fill in data for sector j by applying rotation matrix to data from sector 1 (A)
        iseg=i+(j-1)*n_segments;  % extended segment number goes from 1 to 6*n_segments
        % Cell-PSA interface points
        CellInterface_M1(:,:,iseg)=R*CellInterface_M1(:,:,i);
        % Actuator locations and centers of rotation
        ActOrigin_M1(:,:,iseg)=R*ActOrigin_M1(:,:,i);
        ActCOR_M1(:,iseg)=R*ActCOR_M1(:,i);
        ActLOA_M1(:,:,iseg)=R*ActLOA_M1(:,:,i);
        % UNGapped Segment vertex coordinates (M1CRS)
        UngappedVertex_M1(:,:,iseg)=R*UngappedVertex_M1(:,:,i);
        % Gapped Segment vertex coordinates (M1CRS)
        Vertex_M1(:,:,iseg)=R*Vertex_M1(:,:,i);        
        % Segment center coordinates (M1CRS)
        Center_M1(:,iseg)=R*Center_M1(:,i);
        % Orientation of PSA coordinate systems
        RPSA_M1(:,:,iseg)=R*RPSA_M1(:,:,i);
        for k=1:12  % loop on 12 edges sensors for segment #i
            ESOrigin_M1(:,iseg,k)=R*ESOrigin_M1(:,i,k);
            RES_M1(:,:,iseg,k)=R*RES_M1(:,:,i,k);
        end;
    end;
end;

%########################################### Search through ES data for mating pairs ###########################################
% This search is done in two levels to reduce execution time: 
%     1) build a list of neigboring segments, then
%     2) search for edge sensor pairs only among neighboring segments; this avoids searching across all possible pairs of segments, which would results in unacceptable computing time

% Step 1 - build list of neighboring segments
dt=cputime-tstart;
disp(['T=' num2str(dt) '  Searching for neighboring segments']);
Neighbor=zeros(6*n_segments,6);  % intitialize list of neighbors to zero (any segment has at most 6 neighbors)
for i=1:6*n_segments
    k=0;  % intialize counter for number of neighbors of segment #i
    for j=1:6*n_segments
        dist=norm(Center_M1(:,i)-Center_M1(:,j));  % distance from center of segment #i to center of segment #j
        if dist<2.5*s & i~=j % neighboring segments are aproximately 2*s away from one another
            k=k+1;
            Neighbor(i,k)=j;
        end;
    end;
end;
% Step 2 - search among neighboring segments for pairs of mating edge sensors
dt=cputime-tstart;
disp(['T=' num2str(dt) '  Searching for edge segment pairs']);
ESmate=zeros(6*n_segments,12,2);  % initialize mating data to zero; 
                                  % for segment #i, edge sensor #j, the mating segment (#im) and edge sensor number (#jm) are found as ESmate(i,j,:)=[im jm]; 
for i=1:6*n_segments  % loop through all segments of all sectors
    nn=length(find(Neighbor(i,:)~=0));  % number of neighbors to segment #i
    for j=1:12  % loop through all edge sensors for each segment
        for k=1:nn  % loop through all neighboring segments of segment #i, looking for mating sensor
            im=Neighbor(i,k);  % index of segment neighboring segment i
            for jm=1:12  % loop through all edge sensors again, looking for mating sensor
                dist=norm(ESOrigin_M1(:,i,j)-ESOrigin_M1(:,im,jm));  % calculate distance between two edge sensors
                if dist<ES_Dist/10 & i~=im & j~=jm  % tolerance on coincidence of edge sensors (numerical roundoff); note that non-mating sensors are at least ~2*ES_Dist*sin(60deg) apart
                   ESmate(i,j,:)=[im jm];  % if sensor im,jm is closer to sensor i,j than the set tolerance, and sensors (i,j) is distinct from sensor (im,jm), then they are a mating pair;
                                           % otherwise, leave ESmate as initialized (i.e. zero), indicating no mate found.
                end;
            end;
        end;
    end;
end;

%############################## Store key results in MAT file for comparison with C.Baffes results #############################
%                                         (as specified by C. Baffes, Jan 2010)
dt=cputime-tstart;
disp(['T=' num2str(dt) '  Reorganizing data per C. Baffes and producing datafile NewSegm.mat ']);
% preallocate memory (for speed of execution)
Verts_ungapped=zeros(n_segments,5,7);
Verts_gapped=zeros(n_segments,9,7);
PSACRS=zeros(n_segments,10);
SEG_Clocking=zeros(n_segments,2);
Tooling_data=zeros(n_segments,5,3);
ES_M1CRS=zeros(n_segments,22,12,6);
ES_PSACRS=zeros(n_segments,19,12);
ES_Mates=zeros(12*6*n_segments,8);

% loop on all segments in sector A
for i=1:n_segments
    % Segment centers and UNGAPPED vertices, Sector A, in M1CRS polar and Cartesian coordinates
    % row = segment number - 1:82
    % columns = R (scaled, m), theta (rad), M1CRS X (m), M1CRS Y (m), M1CRS Z (m)
    % bin = vertex ID - 1:7, where 7 is the segment center
    Verts_ungapped(i,3:5,1:6)=UngappedVertex_M1(:,:,i);  % x,y,z coordinates of vertices in M1CRS
    Verts_ungapped(i,3:5,7)=Center_M1(:,i);              % x,y,z coordinates of center in M1CRS    
    Verts_ungapped(i,1,1:6)=sqrt(UngappedVertex_M1(1,:,i).^2+UngappedVertex_M1(2,:,i).^2);   % R coordinate of vertices in M1CRS
    Verts_ungapped(i,1,7)=sqrt(Center_M1(1,i)^2+Center_M1(2,i)^2);                           % R coordinate of center in M1CRS
    Verts_ungapped(i,2,1:6)=atan2(UngappedVertex_M1(2,:,i),UngappedVertex_M1(1,:,i));        % Theta coordinate of vertices in M1CRS
    Verts_ungapped(i,2,7)=atan2(Center_M1(2,i),Center_M1(1,i));                              % Theta coordinate of vertices in M1CRS
    
    % Segment centers and GAPPED vertices, sector A, in M1CRS and local Cartesian coordinates
    % row = segment number - 1:82
    % columns = 
        %1:3 gapped vertex location in M1CRS X,Y,Z (m)
        %4:6 gapped vertex locations in SCRS X,Y,Z (m)
        %7:9 gapped vertex locations in PSACRS X,Y,Z (m)
    % bin = vertex number - 1:7, where 7 is the segment center
    Verts_gapped(i,1:3,1:6)=Vertex_M1(:,:,i);  % Vertices, in M1CRS
    Verts_gapped(i,1:3,7)=Center_M1(:,i);  % Center, in M1CRS
    Verts_gapped(i,4:6,1:6)=Vertex_SEG(:,:,i);  % Vertices, in SCRS (aka SEG)
    Verts_gapped(i,4:6,7)=[0;0;0];  % Center, in SCRS
    Verts_gapped(i,7:9,1:6)=Vertex_PSA(:,:,i);  % Vertices, in PSACRS
    Verts_gapped(i,7:9,7)=[0;0;0];  % Center, in PSACRS
    
    % Unit vectors for PSA Coordinate System, Sector A, expressed in M1CRS Cartesian coordinates
    % row = segment number, 1:82
    % columns = 
        %1 segment number
        %2:4  PSACRS-X unit vector, expressed in M1CRS
        %5:7  PSACRS-Y unit vector, expressed in M1CRS  
        %8:10 PSACRS-Z unit vector, expressed in M1CRS
    PSACRS(i,1)=i;
    PSACRS(i,2:4)=RPSA_M1(:,1,i)';
    PSACRS(i,5:7)=RPSA_M1(:,2,i)';
    PSACRS(i,8:10)=RPSA_M1(:,3,i)';

    % Angle of arrow pointing to M1CRS center, Relative to PSACRS-X axis
	% rows = segment number -1:82
	% columns 
		%1 = segment number
		%2 = clocking angle
    SEG_Clocking(i,1)=i;
    SEG_Clocking(i,2)=ToCenterAngle(i);
    
    % Coordinates of Cell-PSA interface nodes
    % row = segment number, 1:82
    % column =
        %1 sector number 1:6; 1=A, 6=F
        %2 segment number 1:82
        %3 AAP point number 1:3
        %4:6 AAP coordinate in M1CRS X,Y,Z (m)
    % bin = AAP point number 1:3
    % bin2 = Sector 1:6; 1=A 6=F
    for j=1:6 % loop on sector number
        for k=1:3 % loop on AAP point number
            AAP_locations(i,1,k,j)=j;  % sector number
            AAP_locations(i,2,k,j)=i;  % segment number
            AAP_locations(i,3,k,j)=k;  % AAP number
            AAP_locations(i,4:6,k,j)=CellInterface_M1(:,k,i+(j-1)*n_segments);  % x,y,z coordinates in M1_CRS
        end;
    end;
    
    % Coordinates of and lines of action of M1CS segment position actuators 
    % row = segment number, 1:82
    % column =
        %1 sector number 1:6; 1=A, 6=F
        %2 segment number 1:82
        %3 Actuator number 1:3
        %4:6 Actuator coordinate in M1CRS X,Y,Z (m)
        %7:9 Actuator line of action, unit vectors expressed in M1CRS
        %10:12 Segment rotation point, in M1CRS
    % bin = Actuator number 1:3
    % bin2 = Sector 1:6; 1=A 6=F
    for j=1:6 % loop on sector number
        for k=1:3 % loop on actuator number
            Actuator_locations(i,1,k,j)=j;  % sector number
            Actuator_locations(i,2,k,j)=i;  % segment number
            Actuator_locations(i,3,k,j)=k;  % Actuator number
            Actuator_locations(i,4:6,k,j)=ActOrigin_M1(:,k,i+(j-1)*n_segments);  % x,y,z coordinates in M1_CRS
            Actuator_locations(i,7:9,k,j)=ActLOA_M1(:,k,i+(j-1)*n_segments);  % x,y,z components of line of action in M1_CRS
            Actuator_locations(i,10:12,k,j)=ActCOR_M1(:,i+(j-1)*n_segments);  % x,y,z coordinates in M1_CRS
        end;
    end;
        
    % Coordinates of points on the optical surface (in PSACRS) where assembly tooling will contact
    % row = segment number, 1:82
    % column =
        %1  segment number 1:82
        %2  post number 1:3
        %3:5  post contact point coordinate in PSACRS X,Y,Z (m)
    % bin = post number 1:3
    for k=1:3
        Tooling_data(i,1,k)=i;  % segment number
        Tooling_data(i,2,k)=k;  % contact point number
        Tooling_data(i,3:5,k)=Fiducial_PSA(:,k,i)';  % coordinates of tooling point in PSACRS
    end;
    
    % Edge sensor locations and orientations (in M1CRS Cartesian Coordinates), and mating
    % rows = segment number - 1:82
    % columns = 
        % 1:4 segment #, ES # , Vnear #, Vfar #
        % 5:7 Vnear coordinate in M1CRS X,Y,Z (m)
        % 8:10 P(ES#) coordinate in M1CRS X,Y,Z (m)
        % 11:13 ESCRS-X unit vector, expressed in M1CRS
        % 14:16 ESCRS-Y unit vector, expressed in M1CRS
        % 17:19 ESCRS-Z unit vector, expressed in M1CRS
        % 20:22 Sector, Segment, and ES# of mating edge sensor.  0=no mate
    % bin = edge sensor number 1:12
    % bin2 = Sector 1:6; 1=A 6=F
    for j=1:6 % loop on sector number
        for k=1:12 % loop on edge sensor #
            iseg=i+(j-1)*n_segments;  % generalized segment number (across all sectors)
            ES_M1CRS(i,1,k,j)=i;  % segment #
            ES_M1CRS(i,2,k,j)=k;  % ES #
            ES_M1CRS(i,3,k,j)=ES_Vertex(k,1);  % near vertex #
            ES_M1CRS(i,4,k,j)=ES_Vertex(k,2);  % far vertex #
            ES_M1CRS(i,5:7,k,j)=UngappedVertex_M1(:,ES_Vertex(k,1),iseg)';  % coordinates of near vertex in M1CRS X,Y,Z (m)
            ES_M1CRS(i,8:10,k,j)=ESOrigin_M1(:,iseg,k)';  % coordinates of origin of ESCRS in M1CRS X,Y,Z (m)
            ES_M1CRS(i,11:13,k,j)=RES_M1(:,1,iseg,k)';  % ESCRS-X unit vector, expressed in M1CRS
            ES_M1CRS(i,14:16,k,j)=RES_M1(:,2,iseg,k)';  % ESCRS-Y unit vector, expressed in M1CRS
            ES_M1CRS(i,17:19,k,j)=RES_M1(:,3,iseg,k)';  % ESCRS-Z unit vector, expressed in M1CRS
            if ESmate(iseg,k,1)~=0
                ES_M1CRS(i,20,k,j)=fix((ESmate(iseg,k,1)-1)/n_segments)+1;  %  sector # for mating sensor
                ES_M1CRS(i,21,k,j)=mod(ESmate(iseg,k,1)-1,n_segments)+1;  %  segment # for mating sensor
            else
                ES_M1CRS(i,20:21,k,j)=[0;0];
            end;
            ES_M1CRS(i,22,k,j)=ESmate(iseg,k,2);  %  edge sensor # for mating sensor
        end;
    end;

    % Edge sensor locations and orientations (in PSACRS Cartesian Coordinates)
    % rows = segment number 1:82
    % columns = 
        % 1:4 segment #, ES # , Vnear #, Vfar #
        % 5:7 Vnear coordinate in PSACRS X,Y,Z (m)
        % 8:10 P(ES#) coordinate in PSACRS X,Y,Z (m)
        % 11:13 ESCRS-X unit vector, expressed in PSACRS
        % 14:16 ESCRS-Y unit vector, expressed in PSACRS
        % 17:19 ESCRS-Z unit vector, expressed in PSACRS
    % bin = edge sensor number 1:12
    for k=1:12 % loop on edge sensor #
        ES_PSACRS(i,1,k)=i;  % segment #
        ES_PSACRS(i,2,k)=k;  % ES #
        ES_PSACRS(i,3,k)=ES_Vertex(k,1);  % near vertex #
        ES_PSACRS(i,4,k)=ES_Vertex(k,2);  % far vertex #
        ES_PSACRS(i,5:7,k)=UngappedVertex_PSA(:,ES_Vertex(k,1),i)';  % coordinates of near vertex in PSACRS X,Y,Z (m)
        ES_PSACRS(i,8:10,k)=ESOrigin_PSA(:,i,k)';  % coordinates of origin of ESCRS in PSACRS X,Y,Z (m)
        ES_PSACRS(i,11:13,k)=RES_PSA(:,1,i,k)';  % ESCRS-X unit vector, expressed in PSACRS
        ES_PSACRS(i,14:16,k)=RES_PSA(:,2,i,k)';  % ESCRS-Y unit vector, expressed in PSACRS
        ES_PSACRS(i,17:19,k)=RES_PSA(:,3,i,k)';  % ESCRS-Z unit vector, expressed in PSACRS
    end;

    % Edge sensor mating (connectivity) table
    % rows = Sequential edge sensor numbering, calculated as follows:
        % sequential number = k+12*(i-1)+82*12*(j-1)
        % where k=sensor #, i=segment #, and j=sector #
        % so sector A, segment 1, sensor 1 is #1, and we increment sensors first, then segments, then sectors
        % note: this assigns numbers to edge sensors at the edge of the array that will not really exist.  These sensors have no mates.
    % columns = 
        % 1 Sequential edge sensor number for edge sensor of interest
        % 2:4 sector, segment, sensor # for sensor of interest
        % 5 Sequential edge sensor number for mating sensor (0 = no mate)
        % 6:8 sector, segment, sensor # for mating sensor (0 = no mate)
    for j=1:6 % loop on sector number
        for k=1:12 % loop on edge sensor #
            ESnum=k+12*(i-1)+n_segments*12*(j-1);
            ES_Mates(ESnum,1)=ESnum;  % sequential edge sensor number for sensor of interest
            ES_Mates(ESnum,2)=j;  % sector # for sensor of interest
            ES_Mates(ESnum,3)=i;  % segment # for sensor of interest
            ES_Mates(ESnum,4)=k;  % sensor # for sensor of interest 
            ES_Mates(ESnum,8)=ESmate(i+(j-1)*n_segments,k,2);  % sensor # for mating sensor
            if ES_Mates(ESnum,8)~=0
                ES_Mates(ESnum,6)=fix((ESmate(i+(j-1)*n_segments,k,1)-1)/n_segments)+1;  % sector # for mating sensor
                ES_Mates(ESnum,7)=mod(ESmate(i+(j-1)*n_segments,k,1)-1,n_segments)+1;  % segment # for mating sensor
                ES_Mates(ESnum,5)=ES_Mates(ESnum,8)+12*(ES_Mates(ESnum,7)-1)+n_segments*12*(ES_Mates(ESnum,6)-1);  % sequential edge sensor number for mating sensor
            else
                ES_Mates(ESnum,6)=0;
                ES_Mates(ESnum,7)=0;
                ES_Mates(ESnum,5)=0;
            end;
        end;
    end;
end;

% save all this data into "NewSegm.mat" output file
save NewSegm Verts_ungapped Verts_gapped PSACRS SEG_Clocking Tooling_data AAP_locations Actuator_locations ES_M1CRS ES_PSACRS ES_Mates

%###################################################### SEGMENT STATISTICS ###########################################################
dt=cputime-tstart;
disp(['T=' num2str(dt) '  Calculating segment statistics']);
% calculate angular orientation of segment edges, in M1 system (i.e. in projection, as seen from the sky)
% this data is only used to calculate segment statistics related to diffraction from gaps - does not affect any other calculation
vtx_edge(1,1,:)=[1 2];  % connectivity, edge 1a
vtx_edge(1,2,:)=[5 4];  % connectivity, edge 1b (nearly parallel to 1a)
vtx_edge(2,1,:)=[2 3];  % connectivity, edge 2a
vtx_edge(2,2,:)=[6 5];  % connectivity, edge 2b (nearly parallel to 2a)
vtx_edge(3,1,:)=[3 4];  % connectivity, edge 3a
vtx_edge(3,2,:)=[1 6];  % connectivity, edge 3b (nearly parallel to 3a)
for i=1:n_segments
    for j=1:3 % loop on edge #
        dx=[UngappedVertex_M1(1,vtx_edge(j,:,2),i)-UngappedVertex_M1(1,vtx_edge(j,:,1),i)];  % delta X of edges a and b
        dy=[UngappedVertex_M1(2,vtx_edge(j,:,2),i)-UngappedVertex_M1(2,vtx_edge(j,:,1),i)];  % delta Y of edges a and b
        GapAngle(j,:,i)=atan2(dy,dx);  % angular orientation of gaps
    end;
end;

% Calculate center and radius of minimum bounding circle (in PSA frame, for segment statistics only)
for i=1:n_segments
    [CircumCenter(:,i),CircumRadii(i)]=circum(Vertex_PSA(1:2,:,i));  % see ">>help circum" for more info
end;

% calculate areas of hexagons (for segment statistics only)
% Note EP January 2010: this was incorrectly performed in 3D geometry; modified into 2D, in XY_Temp plane; very small difference since curvature of segments is small
for i=1:n_segments
    for j=1:6 % loop through 6 triangles formed by center of segment and two (gapped) vertices
        R1=Vertex_Temp(:,j,i); R1(3)=0;  % vector from center of segment to gapped vertex number "j" (set Z=0, i.e. work in projection in local XY plane)
        ROpt1=VertexOpt_Temp(:,j,i); ROpt1(3)=0;  % vector from center of segment to gapped & chamfered vertex number "j" (set Z=0, i.e. work in projection in local XY plane)
        if j==6
            k=1;
        else
            k=j+1;
        end
        R2=Vertex_Temp(:,k,i); R2(3)=0;  % vector from center of segment to gapped vertex number "j+1" (set Z=0, i.e. work in projection in local XY plane)
        ROpt2=VertexOpt_Temp(:,k,i); ROpt2(3)=0;  % vector from center of segment to gapped & chamfered vertex number "j+1" (set Z=0, i.e. work in projection in local XY plane)
        TriangleArea(j,i)=norm(cross(R1,R2))/2;  % area of triangle is half of the norm of the cross product of vectors from center to vertex j and from center to vertex j+1
        TriangleAreaOpt(j,i)=norm(cross(ROpt1,ROpt2))/2;  % same for area to edge of optical surface (i.e. accounting for chamfers)
    end;
    HexagonArea(i)=sum(TriangleArea(:,i));  % add 6 triangles to get area of hexagon
    HexagonAreaOpt(i)=sum(TriangleAreaOpt(:,i));  % same for area to edge of optical surface (i.e. accounting for chamfers)
end;    

% find extreme segments
[MinDiam,iMinDiam]=min(2*CircumRadii);
[MaxDiam,iMaxDiam]=max(2*CircumRadii);
[MinArea,iMinArea]=min(HexagonArea);
[MaxArea,iMaxArea]=max(HexagonArea);
[MinIrreg,iMinIrreg]=min(Irregularity);
[MaxIrreg,iMaxIrreg]=max(Irregularity);
[MinClocking,iMinClocking]=min(Clocking);
[MaxClocking,iMaxClocking]=max(Clocking);
[MinRadius,iMinRadius]=min(Radius);
[MaxRadius,iMaxRadius]=max(Radius);
MinGapAngle1=min(min(GapAngle(1,:,:)));
MaxGapAngle1=max(max(GapAngle(1,:,:)));
StdGapAngle1=std(reshape(GapAngle(1,:,:),2*n_segments,1));
MinGapAngle2=min(min(GapAngle(2,:,:)));
MaxGapAngle2=max(max(GapAngle(2,:,:)));
StdGapAngle2=std(reshape(GapAngle(2,:,:),2*n_segments,1));
MinGapAngle3=min(min(GapAngle(3,:,:)));
MaxGapAngle3=max(max(GapAngle(3,:,:)));
StdGapAngle3=std(reshape(GapAngle(3,:,:),2*n_segments,1));
StdGapAngle=sqrt( (1/3)*(StdGapAngle1^2+StdGapAngle2^2+StdGapAngle3^2) );
Rvtx_M1=sqrt(Vertex_M1(1,:,:).^2+Vertex_M1(2,:,:).^2);
RminG=min(min(Rvtx_M1));
RmaxG=max(max(Rvtx_M1));
RvtxOpt_M1=sqrt(VertexOpt_M1(1,:,:).^2+VertexOpt_M1(2,:,:).^2);
RminOpt=min(min(RvtxOpt_M1));
RmaxOpt=max(max(RvtxOpt_M1));

% AVERAGE AREA SEGMENT
% average segment area
Mean_Area=mean(HexagonArea);
% size of regular hexagon of that area (area of a hexagon = (3*sqrt(3)/2)*a^2
Mean_a=sqrt(2*Mean_Area/(3*sqrt(3)));

%###################################################### SEGMENTATION DATAFILE ###########################################################
dt=cputime-tstart;
disp(['T=' num2str(dt) '  Producing segmentation data file (ASCII)']);
fid=fopen('NewSegm.txt','wt');

fprintf(fid,'TMT M1 SEGMENTATION DATA - Eric Ponslet - HYTEC Inc. - %s\n',datestr(now));
fprintf(fid,'This file was generated by NewSegm.m Matlab program\n');
fprintf(fid,'Unless otherwise specified, all linear dimensions are in meters, and all angles are in degrees\n');
fprintf(fid,' \n');
fprintf(fid,'-----------------------------------------------  SECTION 1: CONTROL PARAMETERS AND STATISTICS ----------------------------------------------\n');
fprintf(fid,'1A: M1 GEOMETRY AND SEGMENTATION DATA \n');
fprintf(fid,'   M1 Radius of curvature: k = %7.3f m\n',Opt_k);
fprintf(fid,'   M1 Conic Constant: K = %12.9f \n',Opt_K);
fprintf(fid,'   Base pattern hex diameter: %7.4f m\n',2*a);
fprintf(fid,'   Inter-segment gap: %7.5f m (1/2 gap applied all around every segment (including outer edges of array)\n',Gap);
fprintf(fid,'   Segment chamfer width (projected into XY_PSA): %7.5f m\n',Chamfer);
fprintf(fid,'1B: SEGMENTATION PARAMETERS \n');
fprintf(fid,'   Scaling Parameter: alpha = %7.4f (radial scaling = (1+alpha*(Rmax/k)^2)/(1+alpha*(R/k)^2)\n',ScalingParam);
fprintf(fid,'1C: NOMINAL SEGMENT SIZE \n');
fprintf(fid,'   Nominal segment diameter = %9.6f m \n',2*NominalRadius);
fprintf(fid,'   Best Fit Regular Hexagon (BFRH) statistics:\n');
fprintf(fid,'     Min BFRH diameter = %9.6f m, or %7.2f mm smaller than nominal\n',2*min(Radius),2*1000*(NominalRadius-min(Radius)));
fprintf(fid,'     Mean BFRH diameter = %9.6f m \n',2*mean(Radius));
fprintf(fid,'     Max BFRH diameter = %9.6f m, or %7.2f mm larger than nominal\n',2*max(Radius),2*1000*(max(Radius)-NominalRadius));
fprintf(fid,'1E: SEGMENTATION STATISTICS \n');
fprintf(fid,'   M1 Inner Diameters\n');
fprintf(fid,'     Gapped segments (glass): min diameter = %9.5f m\n',RminG);
fprintf(fid,'             Optical surface: min diameter = %9.5f m\n',RminOpt);
fprintf(fid,'   M1 Outer Diameters\n');
fprintf(fid,'     Gapped segments (glass): max diameter = %9.5f m\n',RmaxG);
fprintf(fid,'             Optical surface: max diameter = %9.5f m\n',RmaxOpt);
fprintf(fid,'   Segment Irregularity (RMS): min = %9.3f mmRMS, at segment # %2i\n',1000*MinIrreg,iMinIrreg);
fprintf(fid,'                               max = %9.3f mmRMS, at segment # %2i\n',1000*MaxIrreg,iMaxIrreg);
fprintf(fid,'   Segment area: min = %9.4f m^2, at segment # %2i\n',MinArea,iMinArea);
fprintf(fid,'                 max = %9.4f m^2, at segment # %2i\n',MaxArea,iMaxArea);
fprintf(fid,'                 spread (max/min-1) = %6.2f percent\n',100*(MaxArea/MinArea-1));
fprintf(fid,'   Circumscribed diameter: min = %9.5f m, at segment # %2i\n',MinDiam,iMinDiam);
fprintf(fid,'                           max = %9.5f m, at segment # %2i\n',MaxDiam,iMaxDiam);
fprintf(fid,'                           spread (max/min-1) = %6.2f percent\n',100*(MaxDiam/MinDiam-1));
fprintf(fid,'   BFRH Clocking angle: min = %9.4f mrad, at segment # %2i\n',1000*MinClocking,iMinClocking);
fprintf(fid,'                        max = %9.4f mrad, at segment # %2i\n',1000*MaxClocking,iMaxClocking);
fprintf(fid,'   BFRH Radius: min = %9.5f m, at segment # %2i\n',MinRadius,iMinRadius);
fprintf(fid,'                max = %9.5f m, at segment # %2i\n',MaxRadius,iMaxRadius);
fprintf(fid,'                spread (max/min-1) = %6.2f percent\n',100*(MaxRadius/MinRadius-1));
fprintf(fid,'   Mean segment area = %9.2f m^2\n',Mean_Area);
fprintf(fid,'   Diameter of segment of mean area = %9.5f m\n',2*Mean_a);
fprintf(fid,'   Mean diameter of BFRH = %9.5f m\n',2*mean(Radius));
fprintf(fid,' \n');
fprintf(fid,'-----------------------------------------------  SECTION 2: DEFINITION OF PSA COORDINATE SYSTEMS -------------------------------------------\n');
fprintf(fid,'DEFINITION OF PSA COORDINATE SYSTEMS - SECTOR A\n');
fprintf(fid,'   Origin of PSA Coordinate System given as coordinates of segment center (ctr) expressed in the M1 Coordinate System, in meters\n');
fprintf(fid,'   Segment center lies in the M1 optical surface\n');
fprintf(fid,'   Orientation of PSA frame in M1 frame given as coordinates of 1xPSA, 1yPSA, and 1zPSA unit vectors, expressed in the M1 system\n');
fprintf(fid,'   For sectors B through F, PSA Coordinate Systems are rotated about Z_M1 by 60 to 300 degrees\n');
fprintf(fid,' \n');
fprintf(fid,'2A: SEGMENT CENTERS / ORIGINS OF PSA COORDINATE SYSTEMS (in meters)\n'); 
fprintf(fid,'seg#\t    X_M1(ctr)\t    Y_M1(ctr)\t    Z_M1(ctr)\n');
for j=1:6 % loop on sector number
    for i=1:n_segments
        fprintf(fid,' %c%-2i\t%13.9f\t%13.9f\t%13.9f\n',SectorID(j),i,Center_M1(:,i+(j-1)*n_segments));
    end;
end;
fprintf(fid,' \n');
fprintf(fid,'2B: ORIENTATIONS OF PSA COORDINATE SYSTEMS\n'); 
fprintf(fid,'seg#\t  X_M1(1xPSA)\t  Y_M1(1xPSA)\t  Z_M1(1xPSA)\t  X_M1(1yPSA)\t  Y_M1(1yPSA)\t  Z_M1(1yPSA)\t  X_M1(1zPSA)\t  Y_M1(1zPSA)\t  Z_M1(1zPSA)\n');
for j=1:6 % loop on sector number
    for i=1:n_segments
        fprintf(fid,' %c%-2i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',SectorID(j),i,RPSA_M1(:,:,i+(j-1)*n_segments));
    end;
end;
fprintf(fid,' \n');
fprintf(fid,'-------------------------------------------  SECTION 3: UNGAPPED SEGMENT VERTEX COORDINATES -------------------------------------------------\n');
fprintf(fid,'   Ungapped segment vertex coordinates in M1 coordinate system\n');
fprintf(fid,' \n');
fprintf(fid,'3A1: UNGAPPED SEGMENT VERTEX COORDINATES EXPRESSED IN M1 COORDINATE SYSTEM - VERTICES 1 THROUGH 3 (meters)\n');
fprintf(fid,'seg#\t   X_M1(vtx1)\t   Y_M1(vtx1)\t   Z_M1(vtx1)\t   X_M1(vtx2)\t   Y_M1(vtx2)\t   Z_M1(vtx2)\t   X_M1(vtx3)\t   Y_M1(vtx3)\t   Z_M1(vtx3)\n');
for j=1:6 % loop on sector number
    for i=1:n_segments
        fprintf(fid,' %c%-2i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',SectorID(j),i,UngappedVertex_M1(:,1:3,i+(j-1)*n_segments));
    end;
end;
fprintf(fid,' \n');
fprintf(fid,'3A2: UNGAPPED SEGMENT VERTEX COORDINATES EXPRESSED IN M1 COORDINATE SYSTEM - VERTICES 4 THROUGH 6 (meters)\n');
fprintf(fid,'seg#\t   X_M1(vtx4)\t   Y_M1(vtx4)\t   Z_M1(vtx4)\t   X_M1(vtx5)\t   Y_M1(vtx5)\t   Z_M1(vtx5)\t   X_M1(vtx6)\t   Y_M1(vtx6)\t   Z_M1(vtx6)\n');
for j=1:6 % loop on sector number
    for i=1:n_segments
        fprintf(fid,' %c%-2i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',SectorID(j),i,UngappedVertex_M1(:,4:6,i+(j-1)*n_segments));
    end;
end;
fprintf(fid,' \n');
fprintf(fid,'--------------------------------------------  SECTION 4: GAPPED SEGMENT VERTEX COORDINATES --------------------------------------------------\n');
fprintf(fid,'   Gapped segment vertex coordinates (before chamfering) expressed in PSA system (3A1 and 3A2) and in M1 coordinate system (3B1 and 3B2)\n');
fprintf(fid,'4A1: SEGMENT VERTEX COORDINATES EXPRESSED IN PSA COORDINATE SYSTEMS - VERTICES 1 THROUGH 3 (meters)\n');
fprintf(fid,'type#\t  X_PSA(vtx1)\t  Y_PSA(vtx1)\t  Z_PSA(vtx1)\t  X_PSA(vtx2)\t  Y_PSA(vtx2)\t  Z_PSA(vtx2)\t  X_PSA(vtx3)\t  Y_PSA(vtx3)\t  Z_PSA(vtx3)\n');
for i=1:n_segments    
    fprintf(fid,'%3i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',i,Vertex_PSA(:,1:3,i));  
end;
fprintf(fid,'\n');
fprintf(fid,'4A2: SEGMENT VERTEX COORDINATES EXPRESSED IN PSA COORDINATE SYSTEMS - VERTICES 4 THROUGH 6 (meters)\n');
fprintf(fid,'type#\t  X_PSA(vtx4)\t  Y_PSA(vtx4)\t  Z_PSA(vtx4)\t  X_PSA(vtx5)\t  Y_PSA(vtx5)\t  Z_PSA(vtx5)\t  X_PSA(vtx6)\t  Y_PSA(vtx6)\t  Z_PSA(vtx6)\n');
for i=1:n_segments    
    fprintf(fid,'%3i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',i,Vertex_PSA(:,4:6,i));
end;
fprintf(fid,' \n');
fprintf(fid,'4B1: SEGMENT VERTEX COORDINATES EXPRESSED IN M1 COORDINATE SYSTEM - VERTICES 1 THROUGH 3 (meters)\n');
fprintf(fid,'seg#\t   X_M1(vtx1)\t   Y_M1(vtx1)\t   Z_M1(vtx1)\t   X_M1(vtx2)\t   Y_M1(vtx2)\t   Z_M1(vtx2)\t   X_M1(vtx3)\t   Y_M1(vtx3)\t   Z_M1(vtx3)\n');
for j=1:6 % loop on sector number
    for i=1:n_segments
        fprintf(fid,' %c%-2i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',SectorID(j),i,Vertex_M1(:,1:3,i+(j-1)*n_segments));
    end;
end;
fprintf(fid,' \n');
fprintf(fid,'4B2: SEGMENT VERTEX COORDINATES EXPRESSED IN M1 COORDINATE SYSTEM - VERTICES 4 THROUGH 6 (meters)\n');
fprintf(fid,'seg#\t   X_M1(vtx4)\t   Y_M1(vtx4)\t   Z_M1(vtx4)\t   X_M1(vtx5)\t   Y_M1(vtx5)\t   Z_M1(vtx5)\t   X_M1(vtx6)\t   Y_M1(vtx6)\t   Z_M1(vtx6)\n');
for j=1:6 % loop on sector number
    for i=1:n_segments
        fprintf(fid,' %c%-2i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',SectorID(j),i,Vertex_M1(:,4:6,i+(j-1)*n_segments));
    end;
end;
fprintf(fid,' \n');
fprintf(fid,'-----------------------------------------------  SECTION 5: SEGMENT ORIENTATION MARKINGS ---------------------------------------------------\n');
fprintf(fid,'LOCATION OF SEGMENT CLOCKING MARK, IN PSA COORDINATE SYSTEM XY_PSA (in degrees, standard RH convention relative to PSA Coordinate System; angles from -180 to 180, measured from +X_PSA)\n');
fprintf(fid,'type#\t   A_PSA(M1center)\n');
for i=1:n_segments    
    fprintf(fid,'%3i\t%18.6f\n',i,180*ToCenterAngle(i)/pi);
end;
SensorLoc=zeros(24,2,3,n_segments);
fprintf(fid,' \n');
fprintf(fid,'------------------------------------------------------  SECTION 6: EDGE SENSORS ------------------------------------------------------------\n');
fprintf(fid,'6A: MATING AND LOCATIONS OF EDGE SENSOR RELATIVE TO M1CRS and PSACRS (in meters)\n'); 
fprintf(fid,'Current Edge Sensor      Mating Edge Sensor       Coordinates of Current Edge Sensor\n');
fprintf(fid,'   ES#   Seg  Sens        ES#   Seg   Sens\t         X_M1\t         Y_M1\t         Z_M1\t        X_PSA\t        Y_PSA\t        Z_PSA\n');
for j=1:6 % loop on sector
    for i=1:n_segments  % loop on segments
        for k=1:12  % loop on edge sensors
            iES=k+12*(i-1)+n_segments*12*(j-1);  
            iseg=i+n_segments*(j-1);
            iMateSensor=ESmate(iseg,k,2);
            if iMateSensor~=0
                iMateSector=fix((ESmate(iseg,k,1)-1)/n_segments)+1;  % number of mating sector
                iMateSegment=mod(ESmate(iseg,k,1)-1,n_segments)+1;  % number of mating segment
                iMateES=iMateSensor+12*(iMateSegment-1)+12*n_segments*(iMateSector-1);
                fprintf(fid,' %5i\t  %c%-2i\t%2i\t%5i\t %c%-2i\t%2i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',...
                          iES,SectorID(j),i,k,iMateES,SectorID(iMateSector),iMateSegment,iMateSensor,ESOrigin_M1(:,iseg,k),ESOrigin_PSA(:,i,k));
            else
                iMateSector=0;  % number of mating sector
                iMateSegment=0;  % number of mating segment
                iMateES=0;
                fprintf(fid,' %5i\t  %c%-2i\t%2i\t%5i\t %c%-2i\t%2i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',...
                          iES,SectorID(j),i,k,iMateES,' ',iMateSegment,iMateSensor,ESOrigin_M1(:,iseg,k),ESOrigin_PSA(:,i,k));
            end;          
        end;
    end;
end;
fprintf(fid,' \n');
fprintf(fid,'6B1: ORIENTATIONS OF ESCRS COORDINATE SYSTEMS, RELATIVE TO M1CRS\n'); 
fprintf(fid,'   ES#   Seg  Sens\t   X_M1(1xES)\t   Y_M1(1xES)\t   Z_M1(1xES)\t   X_M1(1yES)\t   Y_M1(1yES)\t   Z_M1(1yES)\t   X_M1(1zES)\t   Y_M1(1zES)\t   Z_M1(1zES)\n');
for j=1:6 % loop on sectors
    for i=1:n_segments  % loop on segments
        iseg=i+n_segments*(j-1);
        for k=1:12  % loop on edge sensors
            iES=k+12*(i-1)+n_segments*12*(j-1);  
            fprintf(fid,' %5i\t  %c%-2i\t%2i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',iES,SectorID(j),i,k,RES_M1(:,:,iseg,k));
        end;
    end;
end;
fprintf(fid,' \n');
fprintf(fid,'6B2: ORIENTATIONS OF ESCRS COORDINATE SYSTEMS, RELATIVE TO PSACRS\n'); 
fprintf(fid,'type# Sens\t  X_PSA(1xES)\t  Y_PSA(1xES)\t  Z_PSA(1xES)\t  X_PSA(1yES)\t  Y_PSA(1yES)\t  Z_PSA(1yES)\t  X_PSA(1zES)\t  Y_PSA(1zES)\t  Z_PSA(1zES)\n');
for i=1:n_segments  % loop on segments
    for k=1:12  % loop on edge sensors  
        fprintf(fid,'%2i\t%2i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',i,k,RES_PSA(:,:,i,k));
    end;
end;
fprintf(fid,' \n');
fprintf(fid,'-----------------------------------------------  SECTION 7: CELL-PSA INTERFACE NODES -------------------------------------------------------\n');
fprintf(fid,'Location of Cell to PSA interface nodes, expressed in the M1 Coordinate System\n'); 
fprintf(fid,'These points are the geometric centers of the bolt patterns designed to receive the PSA support posts\n');
fprintf(fid,'The faces in which the threaded hole patterns are machined are to be nominally orthogonal to the 1Z_PSA unit vector for the corresponding segment number, as listed in Section 2B\n');
fprintf(fid,'Locations and orientations in other sectors are obtained by rotating the data about the global Z_M1 axis, by +60 degrees from each sector to the next\n');
fprintf(fid,'Locations of AAP nodes in PSA system:\n');
fprintf(fid,'     Angles = %6.1f, %6.1f, %6.1f degrees\n',CellInterface_Angles);
fprintf(fid,'     Radius = %9.6f m\n',CellInterface_Radius);
fprintf(fid,'     Elevation = %9.6f m\n',CellInterface_Z);
fprintf(fid,' \n');
fprintf(fid,'LOCATIONS OF CELL TO PSA INTERFACE NODES, IN M1 COORDINATE SYSTEM XYZ_M1 (meters)\n');
fprintf(fid,'seg#\t  X_M1(intf1)\t  Y_M1(intf1)\t  Z_M1(intf1)\t  X_M1(intf2)\t  Y_M1(intf2)\t  Z_M1(intf2)\t  X_M1(intf3)\t  Y_M1(intf3)\t  Z_M1(intf3)\n');
for j=1:6 % loop on sector number
    for i=1:n_segments
        fprintf(fid,' %c%-2i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',SectorID(j),i,CellInterface_M1(:,:,i+(j-1)*n_segments));
    end;
end;
fprintf(fid,' \n');
fprintf(fid,'----------------------------------------------------  SECTION 8: ACTUATOR DATA -------------------------------------------------------------\n');
fprintf(fid,'Locations of actuator origins in PSA system:\n');
fprintf(fid,'     Angles = %6.1f, %6.1f, %6.1f degrees\n',ActOrigin_Angles);
fprintf(fid,'     Radius = %9.6f m\n',ActOrigin_Radius);
fprintf(fid,'     Elevation = %9.6f m\n',ActOrigin_Z);
fprintf(fid,'Location of center of rotation of segment for actuation, in PSA system:\n');
fprintf(fid,'     Radius = 0.0 m\n');
fprintf(fid,'     Elevation = %9.6f m\n',ActCOR_PSA(3));
fprintf(fid,'Components of actuator output shafts, in PSA system:\n');
for k=1:3 % loop on 3 actuators
    fprintf(fid,'     Actuator %1i: %9.6f, %9.6f, %9.6f\n',k,ActOS_PSA(:,k));
end;
fprintf(fid,' \n');
fprintf(fid,'8A: LOCATIONS OF CENTERS OF SEGMENT ROTATION, IN M1 COORDINATE SYSTEM XYZ_M1 (meters)\n');
fprintf(fid,'seg#\t    X_M1(COR)\t    Y_M1(COR)\t    Z_M1(COR)\n');
for j=1:6 % loop on sector number
    for i=1:n_segments
        fprintf(fid,' %c%-2i\t%13.9f\t%13.9f\t%13.9f\n',SectorID(j),i,ActCOR_M1(:,i+(j-1)*n_segments));
    end;
end;
fprintf(fid,' \n');
fprintf(fid,'8B: LOCATIONS OF ACTUATOR CRS ORIGINS, IN M1 COORDINATE SYSTEM XYZ_M1 (meters)\n');
fprintf(fid,'seg#\t   X_M1(act1)\t   Y_M1(act1)\t   Z_M1(act1)\t   X_M1(act2)\t   Y_M1(act2)\t   Z_M1(act2)\t   X_M1(act3)\t   Y_M1(act3)\t   Z_M1(act3)\n');
for j=1:6 % loop on sector number
    for i=1:n_segments
        fprintf(fid,' %c%-2i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',SectorID(j),i,ActOrigin_M1(:,:,i+(j-1)*n_segments));
    end;
end;
fprintf(fid,' \n');
fprintf(fid,'8C: COMPONENTS OF ACTUATOR LINES OF ACTION, IN M1 COORDINATE SYSTEM XYZ_M1\n');
fprintf(fid,'seg#\t  1X_M1(LOA1)\t  1Y_M1(LOA1)\t  1Z_M1(LOA1)\t  1X_M1(LOA2)\t  1Y_M1(LOA2)\t  1Z_M1(LOA2)\t  1X_M1(LOA3)\t  1Y_M1(LOA3)\t  1Z_M1(LOA3)\n');
for j=1:6 % loop on sector number
    for i=1:n_segments
        fprintf(fid,' %c%-2i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',SectorID(j),i,ActLOA_M1(:,:,i+(j-1)*n_segments));
    end;
end;
fprintf(fid,' \n');
fprintf(fid,'-----------------------------------------------  SECTION 9: SEGMENT TOOLING FIDUCIALS ------------------------------------------------------\n');
fprintf(fid,'Locations of tooling fiducials in PSA system:\n');
fprintf(fid,'     Angles = %6.1f, %6.1f, %6.1f degrees\n',Fiducial_Angles);
fprintf(fid,'     Radius = %9.6f m\n',Fiducial_Radius);
fprintf(fid,'LOCATIONS OF SEGMENT FIDUCIALS, IN LOCAL PSA COORDINATE SYSTEMS (meters)\n');
fprintf(fid,'type#\t  X_PSA(fdc1)\t  Y_PSA(fdc1)\t  Z_PSA(fdc1)\t  X_PSA(fdc2)\t  Y_PSA(fdc2)\t  Z_PSA(fdc2)\t  X_PSA(fdc3)\t  Y_PSA(fdc3)\t  Z_PSA(fdc3)\n');
for i=1:n_segments    
    fprintf(fid,'%3i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',i,Fiducial_PSA(1,1,i),Fiducial_PSA(2,1,i),Fiducial_PSA(3,1,i),Fiducial_PSA(1,2,i),Fiducial_PSA(2,2,i),Fiducial_PSA(3,2,i),Fiducial_PSA(1,3,i),Fiducial_PSA(2,3,i),Fiducial_PSA(3,3,i));
end;

fclose(fid);
    
%################################################### Components of 1g vector in PSA reference frame ####################################################
dt=cputime-tstart;
disp(['T=' num2str(dt) '  Calculating components of gravity vector in each PSACRS for various telescope zenith angles']);
% Prepare gravity vector data for each segment (for use in detailed analyses of segment deformations)
% calculate components of 1g (-1Z_AZS) unit vector, expressed in the PSA system, as a function
% of telescope zenith angle (to be used in evaluating RSS wavefront error of the entire M1), for all segments
OneG_AZS=[0;0;-1];  % unit gravity vector in AZS frame
XM1=[1;0;0];  % coordinates of XM1 in the AZS system (doe not change with azimuth angle since telescope rotates around XM1)
iTelZenith=0;  % initialize zenith angle step counter
for TelZenith=0:5:65.5  % loop on telescope Zenith angle
    iTelZenith=iTelZenith+1;
    YM1=[0;cos(TelZenith*pi/180);sin(TelZenith*pi/180)];   % coordinates of YM1 in the AZS system
    ZM1=[0;-sin(TelZenith*pi/180);cos(TelZenith*pi/180)];  % coordinates of ZM1 in the AZS system
    RM1_AZS=[XM1 YM1 ZM1];  % Rotation matrix: expresses M1 system in the AZS system
    for i=1:6*n_segments
        OneG_PSA(:,i,iTelZenith)=RPSA_M1(:,:,i)'*RM1_AZS'*OneG_AZS;
    end;
end
nZenith=iTelZenith;
        
fid=fopen('OneGComponents.txt','wt');
fprintf(fid,'TMT M1 GRAVITY ORIENTATION DATA - Eric Ponslet - HYTEC Inc. - %s\n',datestr(now));
for k=1:3  % loop on components in PSA
    fprintf(fid,'\nComponent # %1i in PSA\n',k);
    fprintf(fid,['seg\t         Z=0\t        Z=-5\t       Z=-10\t       Z=-15\t       Z=-20\t       Z=-25\t       Z=-30'...
                 '\t       Z=-35\t       Z=-40\t       Z=-45\t       Z=-50\t       Z=-55\t       Z=-60\t       Z=-65\n']);
    for i=1:6*n_segments
        fprintf(fid,'%3i',i);
        for iz=1:nZenith
            fprintf(fid,'\t%12.8f',OneG_PSA(k,i,iz));
        end;
        fprintf(fid,'\n');
    end;
end;
fclose(fid);

% end program if not plotting
if ~GeneratePlots
    dt=cputime-tstart;
    disp(['T=' num2str(dt) '  Execution completed']);
    return; 
end;

%############################################################### PLOT VARIOUS RESULTS ################################################################
dt=cputime-tstart;
disp(['T=' num2str(dt) '  Generating various plots of the results']);
scrsz = get(0,'ScreenSize');

% ----------------------------------------------------------------------------------------------------------------------------------------------------------------
% 2D plot of array with segment numbers and extreme segments identified
h=figure('Position',[50 50 scrsz(3)*2/3 scrsz(4)*2/3],'Name','2D M1 Array - Sector A - extreme segments','NumberTitle','off');
axis([-(a+a/3) cos(30*pi/180)*RmaxBP+a/3 -a/3 RmaxBP+a/3]);
axis('equal');
for i=1:n_segments
    hold on;
    plot(Center_M1(1,i),Center_M1(2,i),'.k');
    plot([Vertex_M1(1,:,i) Vertex_M1(1,1,i)],[Vertex_M1(2,:,i) Vertex_M1(2,1,i)],'-k');
    text(Center_M1(1,i)+0.1,Center_M1(2,i),num2str(i));
end;
hold on; i=iMinDiam;  d=0.9; h1=plot([d*Vertex_M1(1,:,i)+(1-d)*Center_M1(1,i) d*Vertex_M1(1,1,i)+(1-d)*Center_M1(1,i)],[[d*Vertex_M1(2,:,i)+(1-d)*Center_M1(2,i) d*Vertex_M1(2,1,i)+(1-d)*Center_M1(2,i)]],'--g','LineWidth',2);
hold on; i=iMaxDiam;  d=0.9; h2=plot([d*Vertex_M1(1,:,i)+(1-d)*Center_M1(1,i) d*Vertex_M1(1,1,i)+(1-d)*Center_M1(1,i)],[[d*Vertex_M1(2,:,i)+(1-d)*Center_M1(2,i) d*Vertex_M1(2,1,i)+(1-d)*Center_M1(2,i)]],'-g','LineWidth',2);
hold on; i=iMinArea;  d=0.8; h3=plot([d*Vertex_M1(1,:,i)+(1-d)*Center_M1(1,i) d*Vertex_M1(1,1,i)+(1-d)*Center_M1(1,i)],[[d*Vertex_M1(2,:,i)+(1-d)*Center_M1(2,i) d*Vertex_M1(2,1,i)+(1-d)*Center_M1(2,i)]],'--b','LineWidth',2);
hold on; i=iMaxArea;  d=0.8; h4=plot([d*Vertex_M1(1,:,i)+(1-d)*Center_M1(1,i) d*Vertex_M1(1,1,i)+(1-d)*Center_M1(1,i)],[[d*Vertex_M1(2,:,i)+(1-d)*Center_M1(2,i) d*Vertex_M1(2,1,i)+(1-d)*Center_M1(2,i)]],'-b','LineWidth',2);
hold on; i=iMinIrreg; d=0.7; h5=plot([d*Vertex_M1(1,:,i)+(1-d)*Center_M1(1,i) d*Vertex_M1(1,1,i)+(1-d)*Center_M1(1,i)],[[d*Vertex_M1(2,:,i)+(1-d)*Center_M1(2,i) d*Vertex_M1(2,1,i)+(1-d)*Center_M1(2,i)]],'--r','LineWidth',2);
hold on; i=iMaxIrreg; d=0.7; h6=plot([d*Vertex_M1(1,:,i)+(1-d)*Center_M1(1,i) d*Vertex_M1(1,1,i)+(1-d)*Center_M1(1,i)],[[d*Vertex_M1(2,:,i)+(1-d)*Center_M1(2,i) d*Vertex_M1(2,1,i)+(1-d)*Center_M1(2,i)]],'-r','LineWidth',2);
legend([h1 h2 h3 h4 h5 h6],'Min Circum. circle','Max Circum. circle','Min Area', 'Max Area','Min Irregularity','Max Irregularity');
xlabel('X_{M1} (m)');ylabel('Y_{M1} (m)');
title('Segment Outlines projected into XY_{M1} plane');
% plot 0/60/120 lines
for theta=pi/6:2*pi/6:4.01*pi/6
    hold on;polar([theta theta],[0 RmaxBP],':k');
end;
saveas(h,'2D_SectorA_Extremes.fig');

% ----------------------------------------------------------------------------------------------------------------------------------------------------------------
% 2D plot of complete array with segment numbers, vertex numbers, and edge sensor numbers
h=figure('Position',[60 60 scrsz(3)*2/3 scrsz(4)*2/3],'Name','2D M1 Array - All Sectors - with numbering','NumberTitle','off');
axis([-1*(cos(30*pi/180)*RmaxBP+a/3) (cos(30*pi/180)*RmaxBP+a/3) -1*(RmaxBP+a/3) (RmaxBP+a/3)]);
axis('equal');
for i=1:6  % loop on 6 sectors
    for j=1:n_segments  % loop on all segments
        iseg=j+(i-1)*n_segments;  % generalized segment number
        hold on;
        if mod(i,2)==0  % for every other segment:
            C=[1 1 1]*0.9;  % set color to gray
        else
            C=[1 0.8 0.8];  % set color to pale pink
        end;
        patch([Vertex_M1(1,:,iseg) Vertex_M1(1,1,iseg)],[Vertex_M1(2,:,iseg) Vertex_M1(2,1,iseg)],C);  % plot shaded segment with black edge
        plot(Center_M1(1,iseg),Center_M1(2,iseg),'.k');  % place dot at centers of segments
        ID=[SectorID(i) num2str(j)];  % segment number string
        text(Center_M1(1,iseg)+0.1,Center_M1(2,iseg),ID);  % label segments
    end;
    if debug % label vertices and edge sensors for segments #1 and 2 only in each sector (labeling all results in unbearably slow plots), unless in debugging mode (much smaller array)
        jm=n_segments; 
    else
        jm=2; 
    end;
    for j=1:jm 
        iseg=j+(i-1)*n_segments;  % generalized segment number
        d=0.15;
        plot(Vertex_M1(1,:,iseg),Vertex_M1(2,:,iseg),'.b');  % place dots at vertices
        text(Vertex_M1(1,:,iseg)-d*(Vertex_M1(1,:,iseg)-Center_M1(1,iseg)),Vertex_M1(2,:,iseg)-d*(Vertex_M1(2,:,iseg)-Center_M1(2,iseg)),num2str([1:6]'),...
                 'HorizontalAlignment','Center','VerticalAlignment','Middle','Color','b');  % label vertices
        for k=1:12  % loop on all edge sensors
            dR=ESOrigin_M1(1:2,iseg,k)-Center_M1(1:2,iseg);  % radial vector from segment center to edge sensor
            dR=dR/norm(dR);  % normalize to unit length
            plot(ESOrigin_M1(1,iseg,k)-Gap/2*dR(1),ESOrigin_M1(2,iseg,k)-Gap/2*dR(2),'.g');  % place dots for edges sensor markers (locations of dots are NOT representative)
            text(ESOrigin_M1(1,iseg,k)-(Gap/2+d*a)*dR(1),ESOrigin_M1(2,iseg,k)-(Gap/2+d*a)*dR(2),num2str(k),'HorizontalAlignment','Center','VerticalAlignment','Middle','Color','g');  % label edge sensors
        end;
    end;
end;
xlabel('X_{M1} (m)');ylabel('Y_{M1} (m)');
title('Segment Outlines projected into XY_{M1} plane - with segments (black), vertices (blue), and edge sensors (green) identified');
% plot 0/60/120 lines
for theta=pi/6:2*pi/6:4.01*pi/6
    hold on;polar([theta theta],[0 RmaxBP],':k');
end;
saveas(h,'2D_AllSectors_Numbered.fig');

% ----------------------------------------------------------------------------------------------------------------------------------------------------------------
% 2D plot of array with BFRH statistics: clocking angle, irregularity, and radius color code and corresponding histograms
h=figure('Position',[70 70 scrsz(3)*2/3 scrsz(4)*2/3],'Name','Statistics: BFRH Clocking, Irregularity & Radius','NumberTitle','off');
subplot(231);
  %axis([-(a+a/3) cos(30*pi/180)*RmaxBP+a/3 -a/3 RmaxBP+a/3]);
  plot(0,0,'.k');
  for i=1:n_segments
      hold on;
      fill([Vertex_M1(1,:,i) Vertex_M1(1,1,i)],[Vertex_M1(2,:,i) Vertex_M1(2,1,i)],1000*Clocking(i));
  end;
  axis('equal');
  xlabel('X_{M1} (m)');ylabel('Y_{M1} (m)');
  title(': Dist. of Clocking Angles (mrad)');
  colorbar('vert');
subplot(234);
  hist(Clocking*1000,15);
  xlabel('Clocking Angle (mrad)');ylabel('count');
  title('Dist. of Clocking Angles (mrad)');
subplot(232);
  %axis([-(a+a/3) cos(30*pi/180)*RmaxBP+a/3 -a/3 RmaxBP+a/3]);
  plot(0,0,'.k');
  for i=1:n_segments
      hold on;
      fill([Vertex_M1(1,:,i) Vertex_M1(1,1,i)],[Vertex_M1(2,:,i) Vertex_M1(2,1,i)],1000*Irregularity(i));
  end;
  axis('equal');
  xlabel('X_{M1} (m)');ylabel('Y_{M1} (m)');
  title('Dist. of RMS Irregularity (mm)');
  colorbar('vert');
subplot(235);
  hist(Irregularity*1000,15);
  xlabel('RMS Irregularity (mm)');ylabel('count');
  title('Dist. of RMS Irregularity (mm)');
subplot(233);
  %axis([-(a+a/3) cos(30*pi/180)*RmaxBP+a/3 -a/3 RmaxBP+a/3]);
  plot(0,0,'.k');
  for i=1:n_segments
      hold on;
      fill([Vertex_M1(1,:,i) Vertex_M1(1,1,i)],[Vertex_M1(2,:,i) Vertex_M1(2,1,i)],Radius(i));
  end;
  axis('equal');
  xlabel('X_{M1} (m)');ylabel('Y_{M1} (m)');
  title('Dist. of Best Fit Radius (m)');
  colorbar('vert');
subplot(236);
  hist(Radius,15);
  xlabel('Best Fit Radius (m)');ylabel('count');
  title('Dist. of Best Fit Radius (m)');
saveas(h,'Stats_Clocking_Irregularity_Radius.fig');

% ----------------------------------------------------------------------------------------------------------------------------------------------------------------
% 2D plot of array with area and circumscribed radius color code and corresponding histograms
h=figure('Position',[80 80 scrsz(3)*2/3 scrsz(4)*2/3],'Name','Statistics: Segment Area & Radius of Circumscribed Circle','NumberTitle','off');
subplot(221);
  %axis([-(a+a/3) cos(30*pi/180)*RmaxBP+a/3 -a/3 RmaxBP+a/3]);
  plot(0,0,'.k');
  for i=1:n_segments
      hold on;
      fill([Vertex_M1(1,:,i) Vertex_M1(1,1,i)],[Vertex_M1(2,:,i) Vertex_M1(2,1,i)],HexagonArea(i));
  end;
  axis('equal');
  xlabel('X_{M1} (m)');ylabel('Y_{M1} (m)');
  title('Dist. of Segment Area ({m^2})');
  colorbar('vert');
subplot(222);
  %axis([-(a+a/3) cos(30*pi/180)*RmaxBP+a/3 -a/3 RmaxBP+a/3]);
  plot(0,0,'.k');
  for i=1:n_segments
      hold on;
      fill([Vertex_M1(1,:,i) Vertex_M1(1,1,i)],[Vertex_M1(2,:,i) Vertex_M1(2,1,i)],CircumRadii(i));
  end;
  axis('equal');
  xlabel('X_{M1} (m)');ylabel('Y_{M1} (m)');
  title('Dist. of Radius of Circumscribed Circle (m)');
  colorbar('vert');
subplot(223);
  hist(HexagonArea,15);
  xlabel('Segment Area ({m^2})');ylabel('count');
  title('Dist. of Segment Area ({m^2})');
subplot(224);
  hist(CircumRadii,15);
  xlabel('Radius of Circumscribed Circle (m)');ylabel('count');
  title('Dist. of Radius of Circumscribed Circle (m)');
saveas(h,'Stats_Area_Circum.fig');

% ----------------------------------------------------------------------------------------------------------------------------------------------------------------
% 3D plot of segment outlines and PSA and Edge sensor coordinate systems - sector A
h=figure('Position',[90 90 scrsz(3)*2/3 scrsz(4)*2/3],'Name','3D Array - Sector A - Gapped and Ungapped, with PSA and ES coordinate systems','NumberTitle','off');
hold off;
plot3(Center_M1(1,1:n_segments),Center_M1(2,1:n_segments),Center_M1(3,1:n_segments),'.k');
axis([-(a/2+a/3) cos(30*pi/180)*RmaxBP+a/3 -a/3 RmaxBP+a/3]);
axis('equal');
% plot segment outlines
for i=1:n_segments
    hold on;
    plot3([UngappedVertex_M1(1,:,i) UngappedVertex_M1(1,1,i)],[UngappedVertex_M1(2,:,i) UngappedVertex_M1(2,1,i)],[UngappedVertex_M1(3,:,i) UngappedVertex_M1(3,1,i)],':k');  % plot sector A, ungapped contours
    plot3([Vertex_M1(1,:,i) Vertex_M1(1,1,i)],[Vertex_M1(2,:,i) Vertex_M1(2,1,i)],[Vertex_M1(3,:,i) Vertex_M1(3,1,i)],'k');  % plot sector A, gapped contours
    plot3([VertexOpt_M1(1,:,i) VertexOpt_M1(1,1,i)],[VertexOpt_M1(2,:,i) VertexOpt_M1(2,1,i)],[VertexOpt_M1(3,:,i) VertexOpt_M1(3,1,i)],'--k');  % plot sector A, chamfered contours
end;
% plot PSA coordinate systems
PSAaxis=0.3*a;
for i=1:n_segments
    for k=1:3  % loop on three axes of PSA coordinate system
        plot3([Center_M1(1,i) Center_M1(1,i)+PSAaxis*RPSA_M1(1,k,i)],[Center_M1(2,i) Center_M1(2,i)+PSAaxis*RPSA_M1(2,k,i)],[Center_M1(3,i) Center_M1(3,i)+PSAaxis*RPSA_M1(3,k,i)],'b');
    end;

end;
% plot edge sensor coordinate systems
ESaxis=0.3*a;
for i=1:n_segments
%i=1;
    for j=1:12  % loop on all edge sensors of segment 'i'
        plot3(ESOrigin_M1(1,i,j),ESOrigin_M1(2,i,j),ESOrigin_M1(3,i,j),'.g');  % center of ES coordinate system
        for k=1:3  % loop on three axes of ES coordinate system
            plot3([ESOrigin_M1(1,i,j) ESOrigin_M1(1,i,j)+ESaxis*RES_M1(1,k,i,j)],[ESOrigin_M1(2,i,j) ESOrigin_M1(2,i,j)+ESaxis*RES_M1(2,k,i,j)],[ESOrigin_M1(3,i,j) ESOrigin_M1(3,i,j)+ESaxis*RES_M1(3,k,i,j)],'g');
        end;
    end;
end;
% plot M1CRS Coordinate system
LM1CRS=0.8*a;
LM1CRStxt=0.9*a;
plot3(0,0,0,'.k');  % center of M1CRS
plot3([0 LM1CRS],[0 0],[0 0],'-k');  % X axis of M1CRS
text(LM1CRStxt,0,0,'X','HorizontalAlignment','Center','VerticalAlignment','Middle','Color','k');
plot3([0 0],[0 LM1CRS],[0 0],'-k');  % Y axis of M1CRS
text(0,LM1CRStxt,0,'Y','HorizontalAlignment','Center','VerticalAlignment','Middle','Color','k');
plot3([0 0],[0 0],[0 LM1CRS],'-k');  % Z axis of M1CRS
text(0,0,LM1CRStxt,'Z','HorizontalAlignment','Center','VerticalAlignment','Middle','Color','k');

xlabel('X_{M1} (m)');ylabel('Y_{M1} (m)');zlabel('Z_{M1} (m)');
title('Segment Outlines in XYZ_{M1} (black; dotted=ungapped, solid=gapped, dashed=chamfered), PSACRS (blue), and ESCRS (green)');
saveas(h,'3D_SectorA_PSA_ES.fig');

% ----------------------------------------------------------------------------------------------------------------------------------------------------------------
% 3D plot of entire array with supports and actuators shown
h=figure('Position',[100 100 scrsz(3)*2/3 scrsz(4)*2/3],'Name','3D Array - Complete - Gapped, with actuators & cell-PSA interface points','NumberTitle','off');
hold off
% plot segment centers
plot3(Center_M1(1,:),Center_M1(2,:),Center_M1(3,:),'.k');
axis([-1*ceil(cos(30*pi/180)*RmaxBP*1.1) ceil(cos(30*pi/180)*RmaxBP*1.1) -1*ceil(RmaxBP*1.1) ceil(RmaxBP*1.1)]);
axis equal;
ActL=0.4*a;  % length of line segments to represent actuators lines of action
for j=1:6  % loop on 6 sectors
    for i=1:n_segments
        hold on;
        iseg=i+(j-1)*n_segments;
        % plot segment outlines
        plot3([Vertex_M1(1,:,iseg) Vertex_M1(1,1,iseg)],[Vertex_M1(2,:,iseg) Vertex_M1(2,1,iseg)],[Vertex_M1(3,:,iseg) Vertex_M1(3,1,iseg)],'k');  % plot segment
        hold on;
        % plot SSA supports
        CIcenter_M1=1/3*(CellInterface_M1(:,1,iseg)+CellInterface_M1(:,2,iseg)+CellInterface_M1(:,3,iseg)); % center (mean) of cell interface points (for display only)
        for k=1:3
            plot3(CellInterface_M1(1,k,iseg),CellInterface_M1(2,k,iseg),CellInterface_M1(3,k,iseg),'.b');
            plot3([CIcenter_M1(1) CellInterface_M1(1,k,iseg)],[CIcenter_M1(2) CellInterface_M1(2,k,iseg)],[CIcenter_M1(3) CellInterface_M1(3,k,iseg)],':b'); % plot lines connecting each of the 3 cell interface points to the mean of those three points (for ease of visualization)
        end;
        % plot actuators
        plot3(ActCOR_M1(1,iseg),ActCOR_M1(2,iseg),ActCOR_M1(3,iseg),'.r');  % center of rotation
        for k=1:3
            plot3(ActOrigin_M1(1,k,iseg),ActOrigin_M1(2,k,iseg),ActOrigin_M1(3,k,iseg),'.r');  % ACTCRS origin
            plot3([ActOrigin_M1(1,k,iseg) ActOrigin_M1(1,k,iseg)+ActL*ActLOA_M1(1,k,iseg)],...
                  [ActOrigin_M1(2,k,iseg) ActOrigin_M1(2,k,iseg)+ActL*ActLOA_M1(2,k,iseg)],...
                  [ActOrigin_M1(3,k,iseg) ActOrigin_M1(3,k,iseg)+ActL*ActLOA_M1(3,k,iseg)],'-r');  % lines of action
            plot3([ActOrigin_M1(1,k,iseg) ActCOR_M1(1,iseg)],...
                  [ActOrigin_M1(2,k,iseg) ActCOR_M1(2,iseg)],...
                  [ActOrigin_M1(3,k,iseg) ActCOR_M1(3,iseg)],':r');  % plot lines connecting actuator center of rotation to each of the 3 actuator centers (for ease of visualization)
        end;
    end;
end;
% plot M1CRS Coordinate system
LM1CRS=0.8*a;
LM1CRStxt=0.9*a;
plot3(0,0,0,'.k');  % center of M1CRS
plot3([0 LM1CRS],[0 0],[0 0],'-k');  % X axis of M1CRS
text(LM1CRStxt,0,0,'X','HorizontalAlignment','Center','VerticalAlignment','Middle','Color','k');
plot3([0 0],[0 LM1CRS],[0 0],'-k');  % Y axis of M1CRS
text(0,LM1CRStxt,0,'Y','HorizontalAlignment','Center','VerticalAlignment','Middle','Color','k');
plot3([0 0],[0 0],[0 LM1CRS],'-k');  % Z axis of M1CRS
text(0,0,LM1CRStxt,'Z','HorizontalAlignment','Center','VerticalAlignment','Middle','Color','k');

xlabel('X_{M1} (m)');ylabel('Y_{M1} (m)');zlabel('Z_{M1} (m)');
title('Segment Outlines in XYZ_{M1} (black), Actuator Centers & Lines of Action (red), Actuation Center of Rotation (red), and Cell Interface Points (blue)');
saveas(h,'3D_AllSectors_Act_Interf.fig');

% ----------------------------------------------------------------------------------------------------------------------------------------------------------------
% 2D plot of segment outlines in TEMP frame - REDUCED by substracting a regular hexagon of size a_reduced;
h=figure('Position',[110 110 scrsz(3)*2/3 scrsz(4)*2/3],'Name','Segment Outlines in TEMP frame','NumberTitle','off');
angles=[0:5]*pi/3;
a_reduced=(100-50*(Rout/Opt_k)^2)/100*a;
DVx=a_reduced*cos(angles);
DVy=a_reduced*sin(angles);
axis([-1 1 -1 1]*(MaxDiam/2-a_reduced)*1.1);
axis('equal');
for i=1:n_segments
    hold on;
    plot([Vertex_Temp(1,:,i)-DVx Vertex_Temp(1,1,i)-a_reduced],[Vertex_Temp(2,:,i)-DVy Vertex_Temp(2,1,i)],'-y');
end;
hold on; i=iMinDiam;  h1=plot([Vertex_Temp(1,:,i)-DVx Vertex_Temp(1,1,i)-a_reduced],[Vertex_Temp(2,:,i)-DVy Vertex_Temp(2,1,i)],'--g','LineWidth',2);
hold on; i=iMaxDiam;  h2=plot([Vertex_Temp(1,:,i)-DVx Vertex_Temp(1,1,i)-a_reduced],[Vertex_Temp(2,:,i)-DVy Vertex_Temp(2,1,i)],'-g','LineWidth',2);
hold on; i=iMinArea;  h3=plot([Vertex_Temp(1,:,i)-DVx Vertex_Temp(1,1,i)-a_reduced],[Vertex_Temp(2,:,i)-DVy Vertex_Temp(2,1,i)],'--b','LineWidth',2);
hold on; i=iMaxArea;  h4=plot([Vertex_Temp(1,:,i)-DVx Vertex_Temp(1,1,i)-a_reduced],[Vertex_Temp(2,:,i)-DVy Vertex_Temp(2,1,i)],'-b','LineWidth',2);
hold on; i=iMinIrreg; h5=plot([Vertex_Temp(1,:,i)-DVx Vertex_Temp(1,1,i)-a_reduced],[Vertex_Temp(2,:,i)-DVy Vertex_Temp(2,1,i)],'--r','LineWidth',2);
hold on; i=iMaxIrreg; h6=plot([Vertex_Temp(1,:,i)-DVx Vertex_Temp(1,1,i)-a_reduced],[Vertex_Temp(2,:,i)-DVy Vertex_Temp(2,1,i)],'-r','LineWidth',2); 
legend([h1 h2 h3 h4 h5 h6],'Min Circum. circle','Max Circum. circle','Min Area', 'Max Area','Min Irregularity','Max Irregularity');
xlabel('X_{Temp} (m)');ylabel('Y_{Temp} (m)');
title('Segment Outlines projected into XY_{Temp} plane (reduced)');
saveas(h,'Outlines_Temp_top.fig');

% ----------------------------------------------------------------------------------------------------------------------------------------------------------------
% 2D plot of segment outlines in PSA frame - REDUCED by substracting a regular hexagon of size a_reduced;
h=figure('Position',[120 120 scrsz(3)*2/3 scrsz(4)*2/3],'Name','Segment Outlines in PSA frame','NumberTitle','off');
angles=[0:5]*pi/3;
DVx=a_reduced*cos(angles);
DVy=a_reduced*sin(angles);
axis([-1 1 -1 1]*(MaxDiam/2-a_reduced)*1.1);
axis('equal');
for i=1:n_segments
    hold on;
    plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-a_reduced],[Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)],'-y');
end;
hold on; i=iMinDiam;  h1=plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-a_reduced],[Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)],'--g','LineWidth',2);
hold on; i=iMaxDiam;  h2=plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-a_reduced],[Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)],'-g','LineWidth',2);
hold on; i=iMinArea;  h3=plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-a_reduced],[Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)],'--b','LineWidth',2);
hold on; i=iMaxArea;  h4=plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-a_reduced],[Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)],'-b','LineWidth',2);
hold on; i=iMinIrreg; h5=plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-a_reduced],[Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)],'--r','LineWidth',2);
hold on; i=iMaxIrreg; h6=plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-a_reduced],[Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)],'-r','LineWidth',2); 
legend([h1 h2 h3 h4 h5 h6],'Min Circum. circle','Max Circum. circle','Min Area', 'Max Area','Min Irregularity','Max Irregularity');
xlabel('X_{PSA} (m)');ylabel('Y_{PSA} (m)');
title('Segment Outlines projected into XY_{PSA} plane (reduced)');
saveas(h,'Outlines_PSA_top.fig');

% ----------------------------------------------------------------------------------------------------------------------------------------------------------------
% XY and XZ plots of segment outlines in PSA frame - REDUCED by substracting a
% regular hexagon of size a_reduced;
h=figure('Position',[130 130 scrsz(3)*2/3 scrsz(4)*2/3],'Name','Segment Profiles in PSA frame','NumberTitle','off');
angles=[0:5]*pi/3;
DVx=a_reduced*cos(angles);
DVy=a_reduced*sin(angles);
DVz=zeros(size(DVx));
subplot(211);
  for i=1:n_segments
      hold on;
      plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-DVx(1)],[Vertex_PSA(3,:,i)-DVz Vertex_PSA(3,1,i)],'-y');
  end;
  hold on; i=iMinDiam;  h1=plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-DVx(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'--g','LineWidth',2);
  hold on; i=iMaxDiam;  h2=plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-DVx(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'-g','LineWidth',2);
  hold on; i=iMinArea;  h3=plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-DVx(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'--b','LineWidth',2);
  hold on; i=iMaxArea;  h4=plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-DVx(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'-b','LineWidth',2);
  hold on; i=iMinIrreg; h5=plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-DVx(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'--r','LineWidth',2);
  hold on; i=iMaxIrreg; h6=plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-DVx(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'-r','LineWidth',2); 
  legend([h1 h2 h3 h4 h5 h6],'Min Circum. circle','Max Circum. circle','Min Area', 'Max Area','Min Irregularity','Max Irregularity');
  xlabel('X_{PSA} (m)');ylabel('Z_{PSA} (m)');
  title('Segment Outlines projected into XZ_{PSA} plane (reduced)');
subplot(212);
  for i=1:n_segments
      hold on;
      plot([Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)-DVy(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'-y');
  end;
  hold on; i=iMinDiam;  h1=plot([Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)-DVy(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'--g','LineWidth',2);
  hold on; i=iMaxDiam;  h2=plot([Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)-DVy(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'-g','LineWidth',2);
  hold on; i=iMinArea;  h3=plot([Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)-DVy(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'--b','LineWidth',2);
  hold on; i=iMaxArea;  h4=plot([Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)-DVy(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'-b','LineWidth',2);
  hold on; i=iMinIrreg; h5=plot([Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)-DVy(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'--r','LineWidth',2);
  hold on; i=iMaxIrreg; h6=plot([Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)-DVy(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'-r','LineWidth',2); 
  legend([h1 h2 h3 h4 h5 h6],'Min Circum. circle','Max Circum. circle','Min Area', 'Max Area','Min Irregularity','Max Irregularity');
  xlabel('Y_{PSA} (m)');ylabel('Z_{PSA} (m)');
  title('Segment Outlines projected into YZ_{PSA} plane (reduced)');
saveas(h,'Outlines_PSA_side.fig');

dt=cputime-tstart;
disp(['T=' num2str(dt) '  Execution completed']);
