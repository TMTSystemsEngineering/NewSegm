function [] = NewSegm(GeneratePlots,debugParameters);
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
%   - adjusted dNormalM1OSToCellNode to generate local AAP coordinates as given in segmentation database printout dated 26-NOV-2007 (EP-EW version)
%
% Modified and Extended: Eric Ponslet, February & March 2010, elp@ericandlucie.com
% Changes:
%     Modified code to make compatible with MATLAB Version 6.1.0.450 (R12.1)
%         removed all uses of cosd, sind, and tand functions (replaced with sin/cos/tan(arg*pi/180)
%         removed use of circshift function (replaced with explicit array construction)
%         rewrote best fit regular hexagon routines to avoid use of lsqcurvefit function (replaced with more explicit LSQ using fminsearch); this results in new m-files fithexnew.m and fitnew.m
%         modified the circumscribed circle routines circum.m and circumradius.m (different grammar of call to fminsearch in older version)
%     Made changes to make code consistent with C. Baffes REL07 of the code:
%         Changed dNormalM1OSToCellNode to make Z of AAP/M1 Cell iinterface point = -388mm = CAD model; dNormalM1OSToCellNode was dNormalM1OSToCellNode=0.36222495, now is dNormalM1OSToCellNode=0.39222699; As a result, Z_PSA was -358mm, now is -388mm
%         Conic constant was -1.000953; now is -1.00095348 per OAD
%         Replaced text output of sections 2, 3, and 4 with tab-delimited output
%     Changed calculation of gapping vectors ggI and gOptI, from 3D geometry to 2D geometry (the correct approach per definition of gapping process)
%     Added section to store key results into arrays and file "NewSegm.mat" for ease of comparison with CB codes
%     Added more comments throughout
%     Added calculation of reference frames and mating information for edge sensors 
%         Two new functions: Intersect.m: finds the intersection with the optical surface of a line along a given vector, originating from a given point (">>help Intersect" for more info) 
%                            Normal.m: calculates a vector in M1CRS that is normal to the optical surface at a point on the optical surface defined by X,Y (M1CRS) (">>help Normal" for more info)

% extended: Eric Ponslet, 2010
% Additions: 
%     Now computing and listing off axis segment prescription parameters

% Updated and extended: Eric Ponslet, June-September, 2012, based on e-mails from Eric Williams, dated 3/14/2012 and later
% Changes:
%     Several changes/additions/modifications in labeling and titling of the database file (for clarity)
%     Optical surface fiducials moved from radius 0.589250 to 1.115000m
%     Added section 10B, which lists the coordinates of optical surface probe points
%     Added section 10C, which lists the coordinates of subcell alignment targets
%     Added section 6C, which lists the coordinates of the origins of the edge sensor pocket coordinate systems, expressed in PSACRS
%     Added section 8, with coordinates of cell top chord centerline nodes
%
% Updated: Eric Ponslet, October-November 2012, based on e-mails with Eric Williams
%     New definition of cell-side AAP points based on top cell truss plates being parallel to Z_PSA
%     Changed radius to AAP points in PSA (re-optimized to minimize AAP adjustment range)
%     Optical fiducial angles and radius changed (M1S-001-01000-RevF)
%     Optical surface probes angles and radius changed
%
% Corrected and Updated: Eric Ponslet, Dec 2012, based on e-mails with Eric
% Williams
%     Fixed bug in clocking of cell top chord members from sector to sector
%     Cahnged definition of cell top chord: was using verties 1,3,5 in
%     sector A (triangles pointing to +X_M1); now using vertices 2,4,6 in
%     sector A (triangles pointing to +X_M1).
%
% Extended: Eric Ponslet, Feb 2013, based on e-mail with E. Williams 02/08/2013
%     Added computationn and listing of angle, measured about PSACRS_Z (or, in the PSACRS_XY plane), from the intersection between planes PSACRS_XY and M1CRS_RZ to PSACRS_X 
%
% Corrected: Eric Ponslet, August 2013, based on e-mail from Alan Tubb 08/01/13
%     Y_ESCRS offest from ESCRS origin to ESPCRS origin had wrong sign for even numbered edge sensors
%
% Extended: Eric Ponslet, Sept 2013, based on e-mails from Eric Williams 09/25-26/2013
%     Added section 10D "Fixed Frame Alignment Targets" section
% Modified: Eric Ponslet, Oct 15, 2013, based on e-mail from Eric Williams 10/10/2013
%     Removed section 10D and instead implemented same functionality in
%     section 10C with new coordinates
%
% Extended: Eric Ponslet, Nov 2014, based on e-mails from E.Williams
%     Added a new set of OS probe points (old set renamed SET1, new set named SET2)
%     Added section 10C to list coordinates of new OS Probe point set (SET2) into DB file
%     Renumbered following sections: 10C to 10D, 10C1 to 10D1, 10C2 to 10D2
%
% Extended: Eric Ponslet, Feb 26, 2015, based on e-mail from E. Williams dated 02/02/2015
%     Added section 11D containing B&W Zernike coefficients in PSA coordinates
%     Added section 12 (based on existence in REL13)
%
% Modified: Brian Sutin, 2016-01-27
%     Replaced Eric Ponslet's name with the TMT document number
%     Fixed 11C last column to correctly print Noll, not B&W Zernike's
%     Fixed 1E showing radii for range of M1 vertices, rather than diameter
%     Added note clarifying section 11 coordinate system
%     Added note to section 4 calling out the plane defining the hex sides
%
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

% locations of optical surface fiducials on optical surface of segments, in PSA coordinate system                      
OSFiducials_Radius=1.3400/2.0; % was 0.589250;       % radial location of fiducials, in local PSA axes (per M1S-001-01000-RevF)
OSFiducials_Angles=[0 120 240];   % angular location of fiducials relative to XYZ_PSA, in degrees (per M1S-001-01000-RevF)
% fiducials are located on the optical surface; Z coordinates will be calculated accordingly

% locations of optical surface probe ***set 1*** locations
OSProbes1_Radius=0.589250;       % radial location of probe locations, in local PSA axes
OSProbes1_Angles=[90 210 330];   % angular location of probe locations to XYZ_PSA, in degrees
% probe locations are located on the optical surface; Z coordinates will be calculated accordingly

% locations of optical surface probe ***set 2*** locations
OSProbes2_Radius=0.6700;       % radial location of probe locations, in local PSA axes
OSProbes2_Angles=[60 180 300];   % angular location of probe locations to XYZ_PSA, in degrees
% probe locations are located on the optical surface; Z coordinates will be calculated accordingly


% locations of subcell alignment targets
SubcellAlignmentTargets_Radius=0.6700;       % radial location of alignment targets, in local PSA axes
SubcellAlignmentTargets_Angles=[0 120 240];   % angular location of alignment targets relative to XYZ_PSA, in degrees
SubcellAlignmentTargets_Z = 0.00363000*[1 1 1];  % Z location of target centers
% data above from e-mail from Eric Williams, dated 10/10/2013

% locations of Cell-PSA interface nodes, in PSA coordinate system
PSASideAAPFace_Radius=0.418276;  %0.418171;      % radial location of AAP mounting faces, in local PSA axes (reoptimized on Nov 13, 2012)
PSASideAAPFace_Angles=[90 210 330];  % angular location of AAP mounting faces relative to XYZ_PSA, in degrees
PSASideAAPFace_Z=-.388000;           % Z location of AAP mounting faces, in local PSA axes
dZTopChordToAAPFace = -0.030;  % distance along normal to M1 surface from cell-side SSA
dNormalM1OSToCellNode=0.39222+dZTopChordToAAPFace;    % distance between optical surface and cell nodes (ADJUST BY TRIAL AND ERROR S.T. Z_PSA COORDINATE OF AAP CENTERS MATCHES CAD MODEL OF PSA)
dZAAPFaceToAAPCenter = -0.067375;       % distance along Z_PSA from AAP mounting face to nominal center of AAP ball joint

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

% locations of edge sensor pocket coordinate systems ESPCRS, relative to
% ESCRS (per E. Williams e-mail 5/7/2012 
% (!!!Y offset must be flipped for even-numbered edge sensors)
ESPOrigin_ESCRS = [0.0; 0.021325; -0.04000];

% Exagerated M1 paramaters, for debugging purposes
if debugParameters
    disp('WARNING: DEBUGGING MODE - your geometry settings are being overwritten!');
    Opt_k=11;
    a=5/2;
    %Gap=a/10;
    Gap=0.5*a;
    Rout=10;
    ES_Dist=a/4;
    NominalRadius=a;
    Chamfer=Gap/4;
    OSFiducial_Radius=0.8*a;
    PSASideAAPFace_Radius=0.5*a;
    PSASideAAPFace_Z=-0.5*a;
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
disp(['    Optical Surface Fiducials, relative to PSA coordinate system (on optical surface): Radius = ' num2str(OSFiducials_Radius,8) ' m;  ' ...
                                                                                              'Angles = ' num2str(OSFiducials_Angles(1),4) ', ' num2str(OSFiducials_Angles(2),4) ', ' num2str(OSFiducials_Angles(3),4) ' degrees']); 
disp(['    Optical Surface Probe SET1, relative to PSA coordinate system (on optical surface): Radius = ' num2str(OSProbes1_Radius,8) ' m;  ' ...
                                                                                              'Angles = ' num2str(OSProbes1_Angles(1),4) ', ' num2str(OSProbes1_Angles(2),4) ', ' num2str(OSProbes1_Angles(3),4) ' degrees']); 
disp(['    Optical Surface Probe SET2, relative to PSA coordinate system (on optical surface): Radius = ' num2str(OSProbes2_Radius,8) ' m;  ' ...
                                                                                              'Angles = ' num2str(OSProbes2_Angles(1),4) ', ' num2str(OSProbes2_Angles(2),4) ', ' num2str(OSProbes2_Angles(3),4) ' degrees']); 
disp(['    Subcell alignment targets, relative to PSA coordinate system (on optical surface): Radius = ' num2str(SubcellAlignmentTargets_Radius,8) ' m;  ' ...
                                                                                              'Angles = ' num2str(SubcellAlignmentTargets_Angles(1),4) ', ' num2str(SubcellAlignmentTargets_Angles(2),4) ', ' num2str(SubcellAlignmentTargets_Angles(3),4) ' degrees;  '...
                                                                                              'Z =' num2str(SubcellAlignmentTargets_Z(1),4) ', ' num2str(SubcellAlignmentTargets_Z(2),4) ', ' num2str(SubcellAlignmentTargets_Z(3),4) ' m']); 
                                                                                      
disp(['    PSA-side AAP mounting face locations, relative to PSA coordinate system: Radius = ' num2str(PSASideAAPFace_Radius,8) ' m;  ' ...
                                                                            'Angles = ' num2str(PSASideAAPFace_Angles(1),4) ', ' num2str(PSASideAAPFace_Angles(2),4) ', ' num2str(PSASideAAPFace_Angles(3),4) ' degrees;  ' ...
                                                                            'Z = ' num2str(PSASideAAPFace_Z,8) ' m']);
disp(['    Nominal distance from optical surface to cell nodes (MUST BE ADJUSTED IF DESIGN CHANGES) : ', num2str(dNormalM1OSToCellNode)]);
disp(['    Nominal distance from M1 Cell top chord centerline to AAP Post mounting face : ', num2str(dZTopChordToAAPFace)]);
disp(['    Nominal distance from AAP Post mounting face to AAP center: ', num2str(dZAAPFaceToAAPCenter)]);
disp(['    PSA-side AAP center locations, relative to PSA coordinate system: Radius = ' num2str(PSASideAAPFace_Radius,8) ' m;  ' ...
                                                                            'Angles = ' num2str(PSASideAAPFace_Angles(1),4) ', ' num2str(PSASideAAPFace_Angles(2),4) ', ' num2str(PSASideAAPFace_Angles(3),4) ' degrees;  ' ...
                                                                            'Z = ' num2str(PSASideAAPFace_Z + dZAAPFaceToAAPCenter,8) ' m']);


%######################################## Main Calculation Section ######################################## 
dt=cputime-tstart;
disp(['T=' num2str(dt) '  Calculating base geometry']);
s=a*cos(30*pi/180);  % radius of inscribed circle in base pattern (i.e. distance from center of hexagon to midpoint of edge)
SectorID='ABCDEF';  % sector identification characters (1=A,... 6=F)

% Calculate coordinates of Cell-PSA interface nodes in the PSA local system
PSASideAAPFace_PSA=[PSASideAAPFace_Radius*cos(PSASideAAPFace_Angles*pi/180); PSASideAAPFace_Radius*sin(PSASideAAPFace_Angles*pi/180); PSASideAAPFace_Z*ones(size(PSASideAAPFace_Angles))];   % each column contains coordinates of one SSA support point, expressed in the PSA system
PSASideAAPCenter_PSA=[PSASideAAPFace_Radius*cos(PSASideAAPFace_Angles*pi/180); PSASideAAPFace_Radius*sin(PSASideAAPFace_Angles*pi/180); (PSASideAAPFace_Z+dZAAPFaceToAAPCenter)*ones(size(PSASideAAPFace_Angles))];   % each column contains coordinates of one SSA support point, expressed in the PSA system

% Calculate coordinates of Actuator origins in the PSA local system
ActOrigin_PSA=[ActOrigin_Radius*cos(ActOrigin_Angles*pi/180); ActOrigin_Radius*sin(ActOrigin_Angles*pi/180); ActOrigin_Z*ones(size(ActOrigin_Angles))];   % each column contains coordinates of one SSA support point, expressed in the PSA system

% Create base pattern
[BaseCenter_M1,BaseVertex_M1,n_segments]=CreateBasePattern(a,Rin,Rout,0);   % produces hexagonal pattern in sector A, with side length 'a', extending (centers of hexagons) from 'Rin' to 'Rout'
                                                                            % see ">>help CreateBasePattern" for more info

% Apply Scaling Rule: scaled Radial Position = Rj*(1 + Alpha*(Rmax/k)^2)  / (1 + Alpha*(Rj/k)^2)
Rvtx_M1=sqrt(BaseVertex_M1(1,:,:).^2+BaseVertex_M1(2,:,:).^2);  % radii to vertices of base pattern, before scaling
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

superHexPoints = [];

% Preallocate array memory for speed of execution
Vertex_PSA=zeros(3,6,n_segments);






% Calculate coordinates of cell nodes: cell nodes are located at a distance dNormalM1OSToCellNode below vertices 1,3, and 5 of each segment, along the local normal to the optical surface
VnCnACE=[2 4 6];  % list of vertices to use for cell nodes in sectors A, C, and E
VnCnBDF=[1 3 5];  % list of vertices to use for cell nodes in sectors B, D, and F
for i=1:n_segments
    for j=1:3 % loop on 3 selected vertices
        VnACE=VnCnACE(j);
        VnBDF=VnCnBDF(j);
        % calculate local normal vector
        VZACE=[-dZdR_V(1,VnACE,i)*cos(Tvtx_M1(1,VnACE,i)) -dZdR_V(1,VnACE,i)*sin(Tvtx_M1(1,VnACE,i)) 1]';
        VZACE=VZACE/norm(VZACE);     % normalize unit normal to opt surface
        VZBDF=[-dZdR_V(1,VnBDF,i)*cos(Tvtx_M1(1,VnBDF,i)) -dZdR_V(1,VnBDF,i)*sin(Tvtx_M1(1,VnBDF,i)) 1]';
        VZBDF=VZBDF/norm(VZBDF);     % normalize unit normal to opt surface
        % Calculate position of cell "node"
        CellNode_M1_ACE(:,j,i)=UngappedVertex_M1(:,VnACE,i)-dNormalM1OSToCellNode*VZACE;
        CellNode_M1_BDF(:,j,i)=UngappedVertex_M1(:,VnBDF,i)-dNormalM1OSToCellNode*VZBDF;
        % and positions of cell top chord nodes
        %CellTopChordNode_M1(:,j,i)=UngappedVertex_M1(:,Vn,i)-(dNormalM1OSToCellNode-dNormalAAPFaceNodeToCellNode)*VZ;
    end;
end;
% calculate coordinates of cell-side interface nodes (2/3 along each cell top chord centerline) and direction of top chord centerlines
for i=1:n_segments
    for j=1:3  % loop on 3 SSA interface nodes
        k=j+1; 
        if k>3, k=1; end;
        CellTopChordThirdPoints_M1_ACE(:,j,i)=(2/3)*CellNode_M1_ACE(:,j,i) + (1/3)*CellNode_M1_ACE(:,k,i);
        CellTopChordVectors_M1_ACE(:,j,i)=CellNode_M1_ACE(:,k,i)-CellNode_M1_ACE(:,j,i);
        CellTopChordThirdPoints_M1_BDF(:,j,i)=(1/3)*CellNode_M1_BDF(:,j,i) + (2/3)*CellNode_M1_BDF(:,k,i);
        CellTopChordVectors_M1_BDF(:,j,i)=CellNode_M1_BDF(:,k,i)-CellNode_M1_BDF(:,j,i);
    end;
end;









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
    
    
    
    % Calculate coordinates of cell side AAP Centers
    % new assumption: cell top chord trusses are nominally aligned such that 
    %  * their top chord centerline goes from cell node to cell node
    %  * their "roll" angulation about the top chord centerline is such that the metal plate they are cut from is nominally parallel to Z_PSA
    dZTopChordToAAPFace_PSA = -1.0 * (dZTopChordToAAPFace);  % see below: vector used to offset points along -Z, so these numbers must be positive
    dZTopChordToAAPCenter_PSA = -1.0 * (dZTopChordToAAPFace + dZAAPFaceToAAPCenter);  % see below: vector used to offset points along -Z, so these numbers must be positive
        
    CellTopChordThirdPoints_PSA_ACE = RPSA_M1(:,:,i)'*(CellTopChordThirdPoints_M1_ACE(:,:,i)-Center_M1(:,i)*[1 1 1]);
    CellTopChordThirdPoints_PSA_BDF = RPSA_M1(:,:,i)'*(CellTopChordThirdPoints_M1_BDF(:,:,i)-Center_M1(:,i)*[1 1 1]);

    % build a unit vector along the nominal AAP post (orthogonal to top chord centerline and // to Z_PSA)
    % Sectors ACE
    for j = 1:3
        topChordVector_M1_ACE = CellTopChordVectors_M1_ACE(:,j,i);
        ZPSA_M1 = RPSA_M1(:,3,i);
        normalToTopChordFace_M1_ACE = cross(topChordVector_M1_ACE, ZPSA_M1);
        AAPPostVector_M1_ACE = cross(normalToTopChordFace_M1_ACE, topChordVector_M1_ACE);
        if AAPPostVector_M1_ACE(3) > 0 % make sure AAPPostVector points away from M1 surface
            AAPPostVector_M1_ACE = -1.0 * AAPPostVector_M1_ACE;
        end;
        % and normalize
        AAPPostVector_M1_ACE = AAPPostVector_M1_ACE / norm(AAPPostVector_M1_ACE);
        CellSideAAPFace_M1_ACE(:,j,i) = CellTopChordThirdPoints_M1_ACE(:,j,i) + AAPPostVector_M1_ACE * dZTopChordToAAPFace_PSA;
        CellSideAAPCenter_M1_ACE(:,j,i) = CellTopChordThirdPoints_M1_ACE(:,j,i) + AAPPostVector_M1_ACE * dZTopChordToAAPCenter_PSA;
    end;
    CellSideAAPFace_PSA_ACE(:,:,i)=RPSA_M1(:,:,i)'*(CellSideAAPFace_M1_ACE(:,:,i)-Center_M1(:,i)*[1 1 1]);
    CellSideAAPCenter_PSA_ACE(:,:,i)=RPSA_M1(:,:,i)'*(CellSideAAPCenter_M1_ACE(:,:,i)-Center_M1(:,i)*[1 1 1]);
    
    % Sectors BDF
    for j = 1:3
        topChordVector_M1_BDF = CellTopChordVectors_M1_BDF(:,j,i);
        ZPSA_M1 = RPSA_M1(:,3,i);
        normalToTopChordFace_M1_BDF = cross(topChordVector_M1_BDF, ZPSA_M1);
        AAPPostVector_M1_BDF = cross(normalToTopChordFace_M1_BDF, topChordVector_M1_BDF);
        if AAPPostVector_M1_BDF(3) > 0 % make sure AAPPostVector points away from M1 surface
            AAPPostVector_M1_BDF = -1.0 * AAPPostVector_M1_BDF;
        end;
        % and normalize
        AAPPostVector_M1_BDF = AAPPostVector_M1_BDF / norm(AAPPostVector_M1_BDF);
        CellSideAAPFace_M1_BDF(:,j,i) = CellTopChordThirdPoints_M1_BDF(:,j,i) + AAPPostVector_M1_BDF * dZTopChordToAAPFace_PSA;
        CellSideAAPCenter_M1_BDF(:,j,i) = CellTopChordThirdPoints_M1_BDF(:,j,i) + AAPPostVector_M1_BDF * dZTopChordToAAPCenter_PSA;
    end;
    CellSideAAPFace_PSA_BDF(:,:,i)=RPSA_M1(:,:,i)'*(CellSideAAPFace_M1_BDF(:,:,i)-Center_M1(:,i)*[1 1 1]);
    CellSideAAPCenter_PSA_BDF(:,:,i)=RPSA_M1(:,:,i)'*(CellSideAAPCenter_M1_BDF(:,:,i)-Center_M1(:,i)*[1 1 1]);
    
    % calculate coordinates of PSA-side Cell-PSA interface nodes in M1 frame
    PSASideAAPFace_M1(:,:,i)=RPSA_M1(:,:,i)*PSASideAAPFace_PSA+Center_M1(:,i)*[1 1 1];
    PSASideAAPCenter_M1(:,:,i)=RPSA_M1(:,:,i)*PSASideAAPCenter_PSA+Center_M1(:,i)*[1 1 1];
    
    % compute AAP adjustments for this segment
    for j=1:3 % loop on 3 AAPs per segment
        dX=CellSideAAPCenter_PSA_ACE(1,j,i) - PSASideAAPCenter_PSA(1,j);
        dY=CellSideAAPCenter_PSA_ACE(2,j,i) - PSASideAAPCenter_PSA(2,j);
        AAPRadialOffset(j,i) = sqrt(dX^2+dY^2);
        AAPVerticalOffset(j,i) = CellSideAAPCenter_PSA_ACE(3,j,i) - PSASideAAPCenter_PSA(3,j);
    end;
    
    
    
    % accumulate 2D (XY_PSA) projections of gapped vertices for calculation of superhex
    superHexPoints = [superHexPoints Vertex_PSA(1:2,:,i)];
    
    % calculate vector from segment center to M1 center, in PSA frame, and the clocking angle of the projection of that vector into XY_PSA, relative to X_PSA
    ToCenter_PSA(:,i)=RPSA_M1(:,:,i)'*(-1*Center_M1(:,i));  % express the vector pointing from O_PSA to O_M1 in the PSA system
    ToCenterAngle(i)=atan2(ToCenter_PSA(2,i),ToCenter_PSA(1,i));  % calculate clocking angle of that vector relative to X_PSA, in projection in the XY_PSA plane
        
    
    % calculate the angle, measured about PSACRS_Z (or, in the PSACRS_XY plane), from the intersection between planes PSACRS_XY and M1CRS_RZ to PSACRS_X  
    % the intersection in question can be computed as a vector which is
    % perpendicular to PSACRS_Z and to M1CRS_Theta
    M1R_M1 = Center_M1(:,i);
    M1Z_M1 = [0; 0; 1];
    M1Theta_M1 = cross(M1Z_M1,M1R_M1);
    M1Theta_PSA = RPSA_M1(:,:,i)'*M1Theta_M1;
    PSAZ_PSA = [0; 0; 1];
    M1RadialVector_PSA = cross(M1Theta_PSA, PSAZ_PSA);
    AngleFromXtoM1RPlane_PSA(i) = -atan2(M1RadialVector_PSA(2),M1RadialVector_PSA(1));

    
    
    % calculate local (PSA) coordinates of segment fiducials on optical surface
    OSFiducials_PSA(1,:,i)=OSFiducials_Radius*cos(OSFiducials_Angles*pi/180);
    OSFiducials_PSA(2,:,i)=OSFiducials_Radius*sin(OSFiducials_Angles*pi/180);
    OSFiducials_PSA(3,:,i)=OpticalLocal(OSFiducials_PSA(1,:,i),OSFiducials_PSA(2,:,i),Opt_k,Opt_K,Center_M1(:,i),RPSA_M1(:,1,i),RPSA_M1(:,2,i),RPSA_M1(:,3,i));
    
    OSProbes1_PSA(1,:,i)=OSProbes1_Radius*cos(OSProbes1_Angles*pi/180);
    OSProbes1_PSA(2,:,i)=OSProbes1_Radius*sin(OSProbes1_Angles*pi/180);
    OSProbes1_PSA(3,:,i)=OpticalLocal(OSProbes1_PSA(1,:,i),OSProbes1_PSA(2,:,i),Opt_k,Opt_K,Center_M1(:,i),RPSA_M1(:,1,i),RPSA_M1(:,2,i),RPSA_M1(:,3,i));
    
    OSProbes2_PSA(1,:,i)=OSProbes2_Radius*cos(OSProbes2_Angles*pi/180);
    OSProbes2_PSA(2,:,i)=OSProbes2_Radius*sin(OSProbes2_Angles*pi/180);
    OSProbes2_PSA(3,:,i)=OpticalLocal(OSProbes2_PSA(1,:,i),OSProbes2_PSA(2,:,i),Opt_k,Opt_K,Center_M1(:,i),RPSA_M1(:,1,i),RPSA_M1(:,2,i),RPSA_M1(:,3,i));
    
    SubcellAlignmentTargets_PSA(1,:,i)=SubcellAlignmentTargets_Radius*cos(SubcellAlignmentTargets_Angles*pi/180);
    SubcellAlignmentTargets_PSA(2,:,i)=SubcellAlignmentTargets_Radius*sin(SubcellAlignmentTargets_Angles*pi/180);
    SubcellAlignmentTargets_PSA(3,:,i)=SubcellAlignmentTargets_Z;
    
    % and transform into global M1CRS 
    OSFiducials_M1(:,:,i) = RPSA_M1(:,:,i)*OSFiducials_PSA(:,:,i)+Center_M1(:,i)*[1 1 1];
    OSProbes1_M1(:,:,i) = RPSA_M1(:,:,i)*OSProbes1_PSA(:,:,i)+Center_M1(:,i)*[1 1 1];
    OSProbes2_M1(:,:,i) = RPSA_M1(:,:,i)*OSProbes2_PSA(:,:,i)+Center_M1(:,i)*[1 1 1];
    SubcellAlignmentTargets_M1(:,:,i) = RPSA_M1(:,:,i)*SubcellAlignmentTargets_PSA(:,:,i)+Center_M1(:,i)*[1 1 1];
    
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
        
        if mod(j,2) % odd
            ESPOrigin_PSA(:,i,j) = ESOrigin_PSA(:,i,j) + RES_PSA(:,:,i,j)*(ESPOrigin_ESCRS);  % per e-mail from E. Williams, 6/24/2012 and 7/??/2012
        else % even
            ESPOrigin_PSA(:,i,j) = ESOrigin_PSA(:,i,j) + RES_PSA(:,:,i,j)*(ESPOrigin_ESCRS.*[1; -1; 1]);  % per e-mail from E. Williams, 6/24/2012 and 7/??/2012
        end;
        
        ESPOrigin_M1(:,i,j) = RPSA_M1(:,:,i)*ESPOrigin_PSA(:,i,j) + Center_M1(:,i);
                
    end;  % end loop on 12 edges sensors for this segment
    
    
    % Compute off-axis Zernike expansion of segment surfaces
    
    % first, compute cylindrical monomial expansion per Nelson & Temple-Raston 1982
    % copy key variables for ease of notation
    e = Rctr_M1(i) / Opt_k;
    k = Opt_k;
    K = Opt_K;
    r = NominalRadius;  % !!!!! NOTE that I am using the radius of the BFRH for the prescription calculation...
    prescriptionRadius(i) = r;  % store that radius
    monomials(i).a20 = (r^2/k) * (2.0 - K * e^2) / 4.0 / (1.0 - K * e^2)^(3.0/2.0);
    monomials(i).a22 = (r^2/k) * (K * e^2) / 4.0 / (1.0 - K * e^2)^(3.0/2.0);
    monomials(i).a31 = (r^3/k^2) * K * e * (1.0 - (K+1.0) * e^2)^0.5 * (4.0 - K * e^2) / 8.0 / (1.0 - K * e^2)^3.0;
    monomials(i).a33 = (r^3/k^2) * K^2 * e^3 * (1.0 - (K+1.0) * e^2)^0.5 / 8.0 / (1.0 - K * e^2)^3.0;
    monomials(i).a40 = (r^4/k^3) * ( 8.0*(1.0+K) - 24.0*K*e^2 + 3.0*K^2*e^4*(1.0-3.0*K) - K^3*e^6*(2.0-K)) / 64.0 / (1.0 - K * e^2)^(9.0/2.0);
    monomials(i).a42 = (r^4/k^3) * ( K * e^2 * ( 2.0*(1.0+3.0*K) - (9.0+7.0*K)*K*e^2 + (2.0+K)*K^2*e^4 ) ) / 16.0 / (1.0 - K * e^2)^(9.0/2.0);
    monomials(i).a44 = (r^4/k^3) * ( K^2 * e^4 * ( 1.0 + 5.0*K - K*e^2*(6.0+5.0*K) ) ) / 64.0 / (1.0 - K * e^2)^(9.0/2.0);
    % the following were not expanded in Nelson's paper so I am setting them to zero as a place holder
    monomials(i).a60 = 0.0;
    monomials(i).a80 = 0.0;
    monomials(i).a51 = 0.0;
    
    % now convert the cylindrical monomial coefficients into Zernike coefficients, using the BW normalization, per Nelson June 22, 2011 and Eric Williams (e-mail to E.Ponslet 06/27/2011) 
    BWZernike(i).c20 = ( monomials(i).a20/2.0 + monomials(i).a40/2.0 + 9.0*monomials(i).a60/20.0 + 2.0*monomials(i).a80/5.0 );
    BWZernike(i).c40 = ( monomials(i).a40/6.0 + monomials(i).a60/4.0 + 2.0*monomials(i).a80/7.0 );
    BWZernike(i).c22 = ( monomials(i).a22 + 3.0*monomials(i).a42/4.0 );
    BWZernike(i).c42 = ( monomials(i).a42/4.0 );
    BWZernike(i).c31 = ( monomials(i).a31/3.0 + 2.0*monomials(i).a51/5.0 );
    BWZernike(i).c33 = monomials(i).a33;
    
    % and Noll normalizations: ( per Nelson, June 22, 2011 )
    NollZernike(i).c20 = (1.0/sqrt(3.0)) * BWZernike(i).c20;
    NollZernike(i).c40 = (1.0/sqrt(5.0)) * BWZernike(i).c40;
    NollZernike(i).c22 = (1.0/sqrt(6.0)) * BWZernike(i).c22;
    NollZernike(i).c42 = (1.0/sqrt(10.0)) * BWZernike(i).c42;
    NollZernike(i).c31 = (1.0/sqrt(8.0)) * BWZernike(i).c31;
    NollZernike(i).c33 = (1.0/sqrt(8.0)) * BWZernike(i).c33;
    
    % and BW Zernike coefficients rotated into PSA coordinates
    angle = ToCenterAngle(i)+pi;
    
    m = 0;
    PSABWZernike(i).c20 = BWZernike(i).c20;
    PSABWZernike(i).c40 = BWZernike(i).c40;
    
    m = 1;
    R11 = cos(m*angle);
    R21 = -sin(m*angle);
    PSABWZernike(i).c31 = R11 * BWZernike(i).c31;
    PSABWZernike(i).c3m1 = R21 * BWZernike(i).c31;  
    
    m = 2;
    R11 = cos(m*angle);
    R21 = -sin(m*angle);
    PSABWZernike(i).c22 = R11 * BWZernike(i).c22;
    PSABWZernike(i).c2m2 = R21 * BWZernike(i).c22;  
    PSABWZernike(i).c42 = R11 * BWZernike(i).c42;
    PSABWZernike(i).c4m2 = R21 * BWZernike(i).c42; 
    
    m = 3;
    R11 = cos(m*angle);
    R21 = -sin(m*angle);
    PSABWZernike(i).c33 = R11 * BWZernike(i).c33;
    PSABWZernike(i).c3m3 = R21 * BWZernike(i).c33;  
    
    
end;
% ---------------------------------------------- END OF MAIN LOOP ON SEGMENTS ---------------------------------------------






disp(['AAP adjustments: Radial: Min = ',num2str(min(min(AAPRadialOffset)))]);
disp(['                         Max = ',num2str(max(max(AAPRadialOffset)))]);
disp(['                         Mean = ',num2str(mean(mean(AAPRadialOffset)))]);
disp(['                 Vertical: Min = ',num2str(min(min(AAPVerticalOffset)))]);
disp(['                           Max = ',num2str(max(max(AAPVerticalOffset)))]);
disp(['                           Mean = ',num2str(mean(mean(AAPVerticalOffset)))]);





% FOLLOWING CODE FOR VERIFICATION OF OPTIMAL PSA SUPPORT DIMENSIONS
% check the min/max values of the Z coordinate of the cell-side interface nodes in PSA
minZ = min(min(CellSideAAPCenter_PSA_ACE(3,:,:)));
maxZ = max(max(CellSideAAPCenter_PSA_ACE(3,:,:)));
meanZ = mean(mean(CellSideAAPCenter_PSA_ACE(3,:,:)));
disp(['  Min Z of all cell-side interface nodes, as expressed in PSA = ', num2str(minZ)]);
disp(['  Max Z of all cell-side interface nodes, as expressed in PSA = ', num2str(maxZ)]);
disp(['  Mid-range Z of all cell-side interface nodes, as expressed in PSA = ', num2str((maxZ+minZ)/2.0)]);
disp(['  Mean Z of all cell-side interface nodes, as expressed in PSA = ', num2str(meanZ)]);


% ############################################################# old code reintroduced to compute optimal radial location of AAP points #####################################################
%################################## reintroducing optimization of radial location of SSA-Cell interface nodes ##################
% Calculate coordinates of Cell-SSA interface nodes (SSA-side) in the PSA
% local system - in Sector A only
RSSASup = a/2/cos(30*pi/180);   % initial guess 

% following statement not accepted under this version of Matlab... replacing with brute force search below
%gamma = fminbnd(@(x) adjust(x,IntfNode_PSA,RSSASup,n_segments),0.5,1.5);  % iterative refinement (adjust multiplier gamma)
% brute force optimization of AAP radius
minAAPRange = 1e20;
minMult = 1.00;
maxMult = 1.02;
dMult = 1e-6;
for multiplier= minMult:dMult:maxMult
    f = adjust(multiplier,CellSideAAPCenter_PSA_ACE,RSSASup,n_segments);
    if f<minAAPRange
        minAAPRange = f;
        multiplierForMinAAPRange = multiplier;
    end;
end;

if (multiplierForMinAAPRange==minMult) | (multiplierForMinAAPRange==maxMult)
    disp('multiplier optimization incomplete - check limits');
    return;
end;

% calculate best value of RSSASup and ZSSASup (use mid-range value of Z_PSA of the cell-side interface nodes)
disp(['  Best value multiplier = ', num2str(multiplierForMinAAPRange)]);
RSSASup=RSSASup*multiplierForMinAAPRange;    % implement multiplier
ZSSASup=(min(min(CellSideAAPCenter_PSA_ACE(3,:,:)))+max(max(CellSideAAPCenter_PSA_ACE(3,:,:))))/2;
disp(['  Optimal Radius to PSA AAP = ' num2str(RSSASup,8) ' m']);
disp(['  Optimal Elevation of PSA AAP = ' num2str(ZSSASup,8) ' m']);













%######################################## Compute size of superhex ########################################
superHexRadius = MinContainingHexRadius(superHexPoints, 1e-6);  % to one micron


%######################################## Propagate some data to the other five sectors ########################################
dt=cputime-tstart;
disp(['T=' num2str(dt) '  Replicating data to sectors B through F']);
% Propagate some data to the other six sectors

CellNode_M1=CellNode_M1_ACE;
CellTopChordThirdPoints_M1=CellTopChordThirdPoints_M1_ACE;
% AAP points
CellSideAAPFace_M1 = CellSideAAPFace_M1_ACE;
CellSideAAPCenter_M1 = CellSideAAPCenter_M1_ACE;

% segment numbering across all sectors is: i + (j-1)*n_segments, where i is segment number in sector A, and j=1..6 for sectors A..F
for j=2:6 % loop on sector number 2 to 6
    Rangle=(j-1)*pi/3; % rotation angle from sector #1 (A) to sector #j 
    R=[cos(Rangle) -sin(Rangle) 0; sin(Rangle) cos(Rangle) 0; 0  0  1]; % rotation matrix from sector #1 (A) to sector #j, in M1_CRS (identity matrix if j=1) 
    for i=1:n_segments  % for each segment type, fill in data for sector j by applying rotation matrix to data from sector 1 (A)
        iseg=i+(j-1)*n_segments;  % extended segment number goes from 1 to 6*n_segments
        % Optical surface fiducials
        OSFiducials_M1(:,:,iseg)=R*OSFiducials_M1(:,:,i);
        % Optical surface probe points
        OSProbes1_M1(:,:,iseg)=R*OSProbes1_M1(:,:,i);
        OSProbes2_M1(:,:,iseg)=R*OSProbes2_M1(:,:,i);
        % subcell alignment targets
        SubcellAlignmentTargets_M1(:,:,iseg)=R*SubcellAlignmentTargets_M1(:,:,i);
        % Cell-PSA interface points
        PSASideAAPCenter_M1(:,:,iseg)=R*PSASideAAPCenter_M1(:,:,i);
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
        % cell top chord nodes
        if (j == 3 || j==5)
            CellNode_M1(:,:,iseg)=R*CellNode_M1_ACE(:,:,i);
            CellTopChordThirdPoints_M1(:,:,iseg)=R*CellTopChordThirdPoints_M1_ACE(:,:,i);
            % AAP points
            CellSideAAPFace_M1(:,:,iseg) = R*CellSideAAPFace_M1_ACE(:,:,i);
            CellSideAAPCenter_M1(:,:,iseg) = R*CellSideAAPCenter_M1_ACE(:,:,i);
        else
            CellNode_M1(:,:,iseg)=R*CellNode_M1_BDF(:,:,i);
            CellTopChordThirdPoints_M1(:,:,iseg)=R*CellTopChordThirdPoints_M1_BDF(:,:,i);
            % AAP points
            CellSideAAPFace_M1(:,:,iseg) = R*CellSideAAPFace_M1_BDF(:,:,i);
            CellSideAAPCenter_M1(:,:,iseg) = R*CellSideAAPCenter_M1_BDF(:,:,i);
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
        if dist<2.5*s && i~=j % neighboring segments are aproximately 2*s away from one another
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
                if dist<ES_Dist/10 && i~=im && j~=jm  % tolerance on coincidence of edge sensors (numerical roundoff); note that non-mating sensors are at least ~2*ES_Dist*sin(60deg) apart
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
            AAP_locations(i,4:6,k,j)=PSASideAAPCenter_M1(:,k,i+(j-1)*n_segments);  % x,y,z coordinates in M1_CRS
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
        
    % Coordinates of optical surface fiducials (in PSACRS)
    % row = segment number, 1:82
    % column =
        %1  segment number 1:82
        %2  post number 1:3
        %3:5  post contact point coordinate in PSACRS X,Y,Z (m)
    % bin = post number 1:3
    for k=1:3
        Tooling_data(i,1,k)=i;  % segment number
        Tooling_data(i,2,k)=k;  % contact point number
        Tooling_data(i,3:5,k)=OSFiducials_PSA(:,k,i)';  % coordinates of tooling point in PSACRS
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
filename=['NewSegm_' strrep(date,'-','') '.txt']
fid=fopen(filename,'wt');

fprintf(fid,'TMT M1 SEGMENTATION DATA - TMT.OPT.TEC.07.044.CCR15 - %s\n',datestr(now));
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
fprintf(fid,'   Edge sensor geometry defined in TMT.M1.M1CS-INT-001.  Edge sensor located %10.6f m from ungapped vertices\n',ES_Dist);
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
fprintf(fid,'     Gapped segments (glass): min diameter = %9.5f m\n',2*RminG);
fprintf(fid,'             Optical surface: min diameter = %9.5f m\n',2*RminOpt);
fprintf(fid,'   M1 Outer Diameters\n');
fprintf(fid,'     Gapped segments (glass): max diameter = %9.5f m\n',2*RmaxG);
fprintf(fid,'             Optical surface: max diameter = %9.5f m\n',2*RmaxOpt);
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
fprintf(fid,'1F: SUPERHEX \n');
fprintf(fid,'   The SuperHex is a minimum radius regular hexagon in the X_PSA plane, centered at O_PSA, with two edges parallel to X_PSA, that contains the gapped vertices of all segment types projected in the XY_PSA plane\n');
fprintf(fid,'   SuperHex Radius  = %13.6f m\n',superHexRadius);
fprintf(fid,'   SuperHex Diameter = %13.6f m\n',superHexRadius*2.0);
fprintf(fid,' \n');
fprintf(fid,'1G: AAP Travel Range Optimization\n');
fprintf(fid,'   The optimum radius for the fixed frame AAP Interface is determined in order to minimize the AAP travel range due to segmentation;\n');
fprintf(fid,'   Fixed frame AAP interface radius: current value = %13.6f m\n', PSASideAAPFace_Radius);
fprintf(fid,'                                     optimal value = %13.6f m\n', RSSASup);
fprintf(fid,'   Locations of PSA-Side AAP flanges, in PSA system:\n');
fprintf(fid,'        Angles = %6.1f, %6.1f, %6.1f degrees\n',PSASideAAPFace_Angles);
fprintf(fid,'        Radius = %9.6f m (this value calculated above in section 1G)\n',PSASideAAPFace_Radius);
fprintf(fid,'        Elevation = %9.6f m\n', PSASideAAPFace_Z);
fprintf(fid,'   Distance from AAP mounting flanges to AAP centers: %9.6f m\n',abs(dZAAPFaceToAAPCenter));
fprintf(fid,'   Locations of PSA-Side AAP Centers, in PSA system:\n');
fprintf(fid,'        Angles = %6.1f, %6.1f, %6.1f degrees\n',PSASideAAPFace_Angles);
fprintf(fid,'        Radius = %9.6f m (this value calculated above in section 1G)\n',PSASideAAPFace_Radius);
fprintf(fid,'        Elevation = %9.6f m\n', PSASideAAPCenter_PSA(3,1,1));
fprintf(fid,'   Maximum AAP Adjustment Ranges with current AAP interface radius:  radial = %13.6f\n',max(max(abs(AAPRadialOffset))));
fprintf(fid,'                                                                     vertical = %13.6f\n',max(max(abs(AAPVerticalOffset))));
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
fprintf(fid,'2B: ORIENTATIONS (UNIT VECTORS) OF PSA COORDINATE SYSTEMS\n'); 
fprintf(fid,'seg#\t  X_M1(1xPSA)\t  Y_M1(1xPSA)\t  Z_M1(1xPSA)\t  X_M1(1yPSA)\t  Y_M1(1yPSA)\t  Z_M1(1yPSA)\t  X_M1(1zPSA)\t  Y_M1(1zPSA)\t  Z_M1(1zPSA)\n');
for j=1:6 % loop on sector number
    for i=1:n_segments
        fprintf(fid,' %c%-2i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',SectorID(j),i,RPSA_M1(:,:,i+(j-1)*n_segments));
    end;
end;
fprintf(fid,' \n');
fprintf(fid,'-------------------------------------------  SECTION 3: UNGAPPED SEGMENT VERTEX COORDINATES -------------------------------------------------\n');
fprintf(fid,'   Ungapped segment vertex coordinates in M1 coordinate system (M1CRS)\n');
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
fprintf(fid,'   Sides of the hexed segment are defined by the plane passing through two vertex points and parallel to the PSACRS z-axis.\n');
fprintf(fid,'   This same plane (modulo the chamfer) therefore defines the edge of the segment aperture.\n');
fprintf(fid,'\n');
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
fprintf(fid,'6A: MATING AND LOCATIONS OF EDGE SENSOR RELATIVE TO M1CRS and PSACRS (in meters) (LOCATION OF TYPICAL POINT P8 ON TMT.M1.M1CS-INT-001)\n');
fprintf(fid,'Note: ES1, ES3, ES5, ES7, ES9, ES11 are "sense" side, and ES2, ES4, ES6, ES8, ES10, ES12 are "drive" side.\n');
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
fprintf(fid,'6C: ORIGIN OF EDGE SENSOR POCKET COORDINATE SYSTEM (ESPCRS, AS SHOWN ON POLISHED SEGMENT DRAWING M1S-001-01000) IN PSACRS\n'); 
fprintf(fid,'Edge sensor pocket center (origin of ESPCRS) is located at coordinates [%13.9f, +/-%13.9f,%13.9f] relative to the corresponding ESCRS coordinate system\n',ESPOrigin_ESCRS); 
fprintf(fid,'type# Sens\t  X_PSA\t  Y_PSA\t  Z_PSA\n');
for i=1:n_segments  % loop on segments
    for k=1:12  % loop on edge sensors  
        fprintf(fid,'%2i\t%2i\t%13.9f\t%13.9f\t%13.9f\n',i,k,ESPOrigin_PSA(:,i,k));
    end;
end;
fprintf(fid,' \n');

fprintf(fid,'----------------------------------------  SECTION 7: CELL-SIDE AAP FLANGE CENTER INTERFACE NODES -------------------------------------------\n');
fprintf(fid,'Location of Cell-Side Cell to PSA interface nodes, expressed in the M1 Coordinate System\n'); 
fprintf(fid,'These points are the theoretical geometric centers of the AAP post mounting flanges which are part of the top chord trusses.\n');
fprintf(fid,'Notes: - top chord truss plates are now assumed to be nominally parallel to Z_PSA for each segment\n');
fprintf(fid,'       - centerline of top chord members assumed to nominally span exactly between cell top chord nodes (whose coordinates are listed in the next section)\n');
fprintf(fid,'Distance from centerline of top chord member to AAP face = %9.6f m\n',abs(dZTopChordToAAPFace));
fprintf(fid,' \n');
fprintf(fid,'LOCATIONS OF CELL-SIDE AAP FLANGE CENTERS, IN M1 COORDINATE SYSTEM XYZ_M1 (meters)\n');
fprintf(fid,'seg#\t  X_M1(intf1)\t  Y_M1(intf1)\t  Z_M1(intf1)\t  X_M1(intf2)\t  Y_M1(intf2)\t  Z_M1(intf2)\t  X_M1(intf3)\t  Y_M1(intf3)\t  Z_M1(intf3)\n');
for j=1:6 % loop on sector number
    for i=1:n_segments
        fprintf(fid,' %c%-2i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',SectorID(j),i,CellSideAAPFace_M1(:,:,i+(j-1)*n_segments));
    end;
end;
fprintf(fid,' \n');

fprintf(fid,'----------------------------------------------------  SECTION 8: CELL TOP CHORD NODES ------------------------------------------------------\n');
fprintf(fid,'Location of the nodes of the CL of the top chord of the M1 cell, expressed in the M1 Coordinate System\n'); 
fprintf(fid,'These points are located along the normals to the M1 optical surface at the ungapped vertices, at a distance %9.6f m from the optical surface.\n',dNormalM1OSToCellNode);
fprintf(fid,' \n');
fprintf(fid,'LOCATIONS OF CELL TOP CHORD NODES, EXPRESSED IN M1 COORDINATE SYSTEM (meters)\n');
fprintf(fid,'(note that 3 cell nodes are listed for each segment; since most cell nodes are shared by three segments, most cell nodes are listed three times)\n');
fprintf(fid,'(listed as nodes 1, 2, 3, corresponding to vertex numbers 1,3,5 for each segment)\n');
fprintf(fid,'seg#\t   X_M1(nde1)\t   Y_M1(nde1)\t   Z_M1(nde1)\t   X_M1(nde2)\t   Y_M1(nde2)\t   Z_M1(nde2)\t   X_M1(nde3)\t   Y_M1(nde3)\t   Z_M1(nde3)\n');
for j=1:6 % loop on sector number
    for i=1:n_segments
        fprintf(fid,' %c%-2i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',SectorID(j),i,CellNode_M1(:,1:3,i+(j-1)*n_segments));
    end;
end;
fprintf(fid,' \n');
fprintf(fid,' \n');

fprintf(fid,'----------------------------------------------------  SECTION 9: ACTUATOR DATA -------------------------------------------------------------\n');
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
fprintf(fid,'9A: LOCATIONS OF CENTERS OF SEGMENT ROTATION, IN M1 COORDINATE SYSTEM XYZ_M1 (meters) (assumed to be the mid-plane of the SSA Lateral Guide Flexure)\n');
fprintf(fid,'seg#\t    X_M1(COR)\t    Y_M1(COR)\t    Z_M1(COR)\n');
for j=1:6 % loop on sector number
    for i=1:n_segments
        fprintf(fid,' %c%-2i\t%13.9f\t%13.9f\t%13.9f\n',SectorID(j),i,ActCOR_M1(:,i+(j-1)*n_segments));
    end;
end;
fprintf(fid,' \n');
fprintf(fid,'9B: LOCATIONS OF ACTUATOR CRS ORIGINS, IN M1 COORDINATE SYSTEM XYZ_M1 (meters)\n');
fprintf(fid,'seg#\t   X_M1(act1)\t   Y_M1(act1)\t   Z_M1(act1)\t   X_M1(act2)\t   Y_M1(act2)\t   Z_M1(act2)\t   X_M1(act3)\t   Y_M1(act3)\t   Z_M1(act3)\n');
for j=1:6 % loop on sector number
    for i=1:n_segments
        fprintf(fid,' %c%-2i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',SectorID(j),i,ActOrigin_M1(:,:,i+(j-1)*n_segments));
    end;
end;
fprintf(fid,' \n');
fprintf(fid,'9C: COMPONENTS OF ACTUATOR LINES OF ACTION, IN M1 COORDINATE SYSTEM XYZ_M1\n');
fprintf(fid,'seg#\t  1X_M1(LOA1)\t  1Y_M1(LOA1)\t  1Z_M1(LOA1)\t  1X_M1(LOA2)\t  1Y_M1(LOA2)\t  1Z_M1(LOA2)\t  1X_M1(LOA3)\t  1Y_M1(LOA3)\t  1Z_M1(LOA3)\n');
for j=1:6 % loop on sector number
    for i=1:n_segments
        fprintf(fid,' %c%-2i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',SectorID(j),i,ActLOA_M1(:,:,i+(j-1)*n_segments));
    end;
end;
fprintf(fid,' \n');
fprintf(fid,'-----------------------------------------------  SECTION 10: FIDUCIALS AND TARGETS ------------------------------------------------------\n');
fprintf(fid,' \n');
fprintf(fid,'10A: OPTICAL SURFACE FIDUCIALS\n');
fprintf(fid,'Locations of optical surface fiducials in PSA system:\n');
fprintf(fid,'     Angles = %6.1f, %6.1f, %6.1f degrees\n',OSFiducials_Angles);
fprintf(fid,'     Radius = %9.6f m\n',OSFiducials_Radius);
fprintf(fid,'\n');
fprintf(fid,'10A1: LOCATIONS OF OPTICAL SURFACE FIDUCIALS, IN LOCAL PSA COORDINATE SYSTEMS (meters)\n');
fprintf(fid,'type#\t  X_PSA(fdc1)\t  Y_PSA(fdc1)\t  Z_PSA(fdc1)\t  X_PSA(fdc2)\t  Y_PSA(fdc2)\t  Z_PSA(fdc2)\t  X_PSA(fdc3)\t  Y_PSA(fdc3)\t  Z_PSA(fdc3)\n');
for i=1:n_segments    
    fprintf(fid,'%3i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',i,OSFiducials_PSA(1,1,i),OSFiducials_PSA(2,1,i),OSFiducials_PSA(3,1,i),OSFiducials_PSA(1,2,i),OSFiducials_PSA(2,2,i),OSFiducials_PSA(3,2,i),OSFiducials_PSA(1,3,i),OSFiducials_PSA(2,3,i),OSFiducials_PSA(3,3,i));
end;
fprintf(fid,'\n');
fprintf(fid,'10A2: LOCATIONS OF OPTICAL SURFACE FIDUCIALS, IN M1CRS COORDINATE SYSTEMS (meters)\n');
fprintf(fid,'seg#\t   X_M1(fdc1)\t   Y_M1(fdc1)\t   Z_M1(fdc1)\t   X_M1(fdc2)\t   Y_M1(fdc2)\t   Z_M1(fdc2)\t   X_M1(fdc3)\t   Y_M1(fdc3)\t   Z_M1(fdc3)\n');
for j=1:6 % loop on sector number
    for i=1:n_segments
        k=i+(j-1)*n_segments;
        fprintf(fid,' %c%-2i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',SectorID(j),i,OSFiducials_M1(1,1,k),OSFiducials_M1(2,1,k),OSFiducials_M1(3,1,k),OSFiducials_M1(1,2,k),OSFiducials_M1(2,2,k),OSFiducials_M1(3,2,k),OSFiducials_M1(1,3,k),OSFiducials_M1(2,3,k),OSFiducials_M1(3,3,k));
    end;
end;
fprintf(fid,' \n');
fprintf(fid,'10B: OPTICAL SURFACE PROBE SET1 LOCATIONS\n');
fprintf(fid,'Locations of optical surface probe SET1 locations in PSA system:\n');
fprintf(fid,'     Angles = %6.1f, %6.1f, %6.1f degrees\n',OSProbes1_Angles);
fprintf(fid,'     Radius = %9.6f m\n',OSProbes1_Radius);
fprintf(fid,'\n');
fprintf(fid,'10B1: LOCATIONS OF OPTICAL SURFACE PROBE SET1 POINTS (S1p), IN LOCAL PSA COORDINATE SYSTEMS (meters)\n');
fprintf(fid,'type#\t  X_PSA(S1p1)\t  Y_PSA(S1p1)\t  Z_PSA(S1p1)\t  X_PSA(S1p2)\t  Y_PSA(S1p2)\t  Z_PSA(S1p2)\t  X_PSA(S1p3)\t  Y_PSA(S1p3)\t  Z_PSA(S1p3)\n');
for i=1:n_segments    
    fprintf(fid,'%3i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',i,OSProbes1_PSA(1,1,i),OSProbes1_PSA(2,1,i),OSProbes1_PSA(3,1,i),OSProbes1_PSA(1,2,i),OSProbes1_PSA(2,2,i),OSProbes1_PSA(3,2,i),OSProbes1_PSA(1,3,i),OSProbes1_PSA(2,3,i),OSProbes1_PSA(3,3,i));
end;
fprintf(fid,'\n');
fprintf(fid,'10B2: LOCATIONS OF OPTICAL SURFACE PROBE SET1 POINTS (S1p), IN M1CRS COORDINATE SYSTEMS (meters)\n');
fprintf(fid,'seg#\t   X_M1(S1p1)\t   Y_M1(S1p1)\t   Z_M1(S1p1)\t   X_M1(S1p2)\t   Y_M1(S1p2)\t   Z_M1(S1p2)\t   X_M1(S1p3)\t   Y_M1(S1p3)\t   Z_M1(S1p3)\n');
for j=1:6 % loop on sector number
    for i=1:n_segments    
        k=i+(j-1)*n_segments;
        fprintf(fid,'% c%-2i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',SectorID(j),i,OSProbes1_M1(1,1,k),OSProbes1_M1(2,1,k),OSProbes1_M1(3,1,k),OSProbes1_M1(1,2,k),OSProbes1_M1(2,2,k),OSProbes1_M1(3,2,k),OSProbes1_M1(1,3,k),OSProbes1_M1(2,3,k),OSProbes1_M1(3,3,k));
    end;
end;


fprintf(fid,' \n');
fprintf(fid,'10C: OPTICAL SURFACE PROBE SET2 LOCATIONS\n');
fprintf(fid,'Locations of optical surface probe SET2 locations in PSA system:\n');
fprintf(fid,'     Angles = %6.1f, %6.1f, %6.1f degrees\n',OSProbes2_Angles);
fprintf(fid,'     Radius = %9.6f m\n',OSProbes2_Radius);
fprintf(fid,'\n');
fprintf(fid,'10C1: LOCATIONS OF OPTICAL SURFACE PROBE SET2 POINTS (S2p), IN LOCAL PSA COORDINATE SYSTEMS (meters)\n');
fprintf(fid,'type#\t  X_PSA(S2p1)\t  Y_PSA(S2p1)\t  Z_PSA(S2p1)\t  X_PSA(S2p2)\t  Y_PSA(S2p2)\t  Z_PSA(S2p2)\t  X_PSA(S2p3)\t  Y_PSA(S2p3)\t  Z_PSA(S2p3)\n');
for i=1:n_segments    
    fprintf(fid,'%3i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',i,OSProbes2_PSA(1,1,i),OSProbes2_PSA(2,1,i),OSProbes2_PSA(3,1,i),OSProbes2_PSA(1,2,i),OSProbes2_PSA(2,2,i),OSProbes2_PSA(3,2,i),OSProbes2_PSA(1,3,i),OSProbes2_PSA(2,3,i),OSProbes2_PSA(3,3,i));
end;
fprintf(fid,'\n');
fprintf(fid,'10C2: LOCATIONS OF OPTICAL SURFACE PROBE SET2 POINTS, IN M1CRS COORDINATE SYSTEMS (meters)\n');
fprintf(fid,'seg#\t   X_M1(S2p1)\t   Y_M1(S2p1)\t   Z_M1(S2p1)\t   X_M1(S2p2)\t   Y_M1(S2p2)\t   Z_M1(S2p2)\t   X_M1(S2p3)\t   Y_M1(S2p3)\t   Z_M1(S2p3)\n');
for j=1:6 % loop on sector number
    for i=1:n_segments    
        k=i+(j-1)*n_segments;
        fprintf(fid,'% c%-2i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',SectorID(j),i,OSProbes2_M1(1,1,k),OSProbes2_M1(2,1,k),OSProbes2_M1(3,1,k),OSProbes2_M1(1,2,k),OSProbes2_M1(2,2,k),OSProbes2_M1(3,2,k),OSProbes2_M1(1,3,k),OSProbes2_M1(2,3,k),OSProbes2_M1(3,3,k));
    end;
end;


fprintf(fid,' \n');
fprintf(fid,'10D: SUBCELL ALIGNMENT TARGETS\n');
fprintf(fid,'Locations of subcell alignment targets in PSA system:\n');
fprintf(fid,'     Angles = %6.1f, %6.1f, %6.1f degrees\n',SubcellAlignmentTargets_Angles);
fprintf(fid,'     Elevations = %12.8f, %12.8f, %12.8f m\n',SubcellAlignmentTargets_Z);
fprintf(fid,'     Radius = %9.6f m\n',SubcellAlignmentTargets_Radius);
fprintf(fid,'\n');
fprintf(fid,'10D1: LOCATIONS OF SUBCELL ALIGNMENT TARGETS, IN LOCAL PSA COORDINATE SYSTEMS (meters)\n');
fprintf(fid,'type#\t  X_PSA(sat1)\t  Y_PSA(sat1)\t  Z_PSA(sat1)\t  X_PSA(sat2)\t  Y_PSA(sat2)\t  Z_PSA(sat2)\t  X_PSA(sat3)\t  Y_PSA(sat3)\t  Z_PSA(sat3)\n');
for i=1:n_segments    
    fprintf(fid,'%3i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',i,SubcellAlignmentTargets_PSA(1,1,i),SubcellAlignmentTargets_PSA(2,1,i),SubcellAlignmentTargets_PSA(3,1,i),SubcellAlignmentTargets_PSA(1,2,i),SubcellAlignmentTargets_PSA(2,2,i),SubcellAlignmentTargets_PSA(3,2,i),SubcellAlignmentTargets_PSA(1,3,i),SubcellAlignmentTargets_PSA(2,3,i),SubcellAlignmentTargets_PSA(3,3,i));
end;
fprintf(fid,'\n');
fprintf(fid,'10D2: LOCATIONS OF SUBCELL ALIGNMENT TARGETS, IN M1CRS COORDINATE SYSTEMS (meters)\n');
fprintf(fid,'seg#\t   X_M1(sat1)\t   Y_M1(sat1)\t   Z_M1(sat1)\t   X_M1(sat2)\t   Y_M1(sat2)\t   Z_M1(sat2)\t   X_M1(sat3)\t   Y_M1(sat3)\t   Z_M1(sat3)\n');
for j=1:6 % loop on sector number
    for i=1:n_segments 
        k=i+(j-1)*n_segments;
        fprintf(fid,' %c%-2i\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\t%13.9f\n',SectorID(j),i,SubcellAlignmentTargets_M1(1,1,k),SubcellAlignmentTargets_M1(2,1,k),SubcellAlignmentTargets_M1(3,1,k),SubcellAlignmentTargets_M1(1,2,k),SubcellAlignmentTargets_M1(2,2,k),SubcellAlignmentTargets_M1(3,2,k),SubcellAlignmentTargets_M1(1,3,k),SubcellAlignmentTargets_M1(2,3,k),SubcellAlignmentTargets_M1(3,3,k));
    end;
end;
fprintf(fid,' \n');

fprintf(fid,'-----------------------------------------------  SECTION 11: OFF-AXIS PRESCRIPTIONS -------------------------------------------------------\n');
fprintf(fid,'NOTES: * all coefficients are listed in MICRONS\n');
fprintf(fid,'       * all coefficients are calculated based on a nominal segment radius as listed below (in meters)\n');
fprintf(fid,'       * The angle t (in degrees) is the angle measured in the plane XY_PSA from the X_PSA axis to the projection into XY_PSA of the R_M1 line from O_M1 to the projection of O_PSA in the XY_M1 plane; this makes t the rotation from SCRS(k) to PSACRS(k).\n');
fprintf(fid,'       * Zernike coefficients are based on the 7 monomial coefficients listed below (i.e. a20, a22, a31, a33, a40, a42, and a44); all other monomial coefficients are assumed equal to zero\n');
fprintf(fid,'       * Zernike coefficients are identified using the n,m notation, i.e. C20=focus, C22=astigmatism, C31=coma, C33=trefoil, C40=spherical aberration, C42=higher order astigmatism\n');
fprintf(fid,'11A: CYLINDRICAL MONOMIAL COEFFICIENTS (microns)\n');
fprintf(fid,'type#\t   radius (m)\t  t (degrees)\t          a20\t          a22\t          a31\t          a33\t          a40\t          a42\t          a44\t  \n');
for i=1:n_segments    
    fprintf(fid,'%3i\t%13.9f\t%13.9f\t%13.6f\t%13.6f\t%13.6f\t%13.6f\t%13.6f\t%13.6f\t%13.6f\n',i,prescriptionRadius(i), 180*(ToCenterAngle(i)+pi)/pi, 1e6*monomials(i).a20, 1e6*monomials(i).a22, 1e6*monomials(i).a31, 1e6*monomials(i).a33, 1e6*monomials(i).a40, 1e6*monomials(i).a42, 1e6*monomials(i).a44);
end;
fprintf(fid,'11B: BORN & WOLF ZERNIKE COEFFICIENTS (microns)\n');
fprintf(fid,'type#\t   radius (m)\t  t (degrees)\t          C20\t          C22\t          C31\t          C33\t          C40\t          C42\t  \n');
for i=1:n_segments    
    fprintf(fid,'%3i\t%13.9f\t%13.9f\t%13.6f\t%13.6f\t%13.6f\t%13.6f\t%13.6f\t%13.6f\n',i,prescriptionRadius(i), 180*(ToCenterAngle(i)+pi)/pi, 1e6*BWZernike(i).c20, 1e6*BWZernike(i).c22, 1e6*BWZernike(i).c31, 1e6*BWZernike(i).c33, 1e6*BWZernike(i).c40, 1e6*BWZernike(i).c42);
end;
fprintf(fid,'11C: NOLL ZERNIKE COEFFICIENTS (microns)\n');
fprintf(fid,'type#\t   radius (m)\t  t (degrees)\t          C20\t          C22\t          C31\t          C33\t          C40\t          C42\t  \n');
for i=1:n_segments    
    fprintf(fid,'%3i\t%13.9f\t%13.9f\t%13.6f\t%13.6f\t%13.6f\t%13.6f\t%13.6f\t%13.6f\n',i,prescriptionRadius(i), 180*(ToCenterAngle(i)+pi)/pi, 1e6*NollZernike(i).c20, 1e6*NollZernike(i).c22, 1e6*NollZernike(i).c31, 1e6*NollZernike(i).c33, 1e6*NollZernike(i).c40, 1e6*NollZernike(i).c42);
end;
fprintf(fid,'11D: BORN & WOLF ZERNIKE COEFFICIENTS ROTATED INTO THE PSACRS AXES (microns)\n');
fprintf(fid,'type#\t   radius (m)\t          C20\t         C2+2\t         c2-2\t         C3+1\t         C3-1\t         C3+3\t         C3-3\t         C40\t         C4+2\t         C4-2\t \n');
for i=1:n_segments    
    fprintf(fid,'%3i\t%13.9f\t%13.6f\t%13.6f\t%13.6f\t%13.6f\t%13.6f\t%13.6f\t%13.6f\t%13.6f\t%13.6f\t%13.6f\n',i,prescriptionRadius(i), 1e6*PSABWZernike(i).c20, 1e6*PSABWZernike(i).c22, 1e6*PSABWZernike(i).c2m2, 1e6*PSABWZernike(i).c31, 1e6*PSABWZernike(i).c3m1, 1e6*PSABWZernike(i).c33, 1e6*PSABWZernike(i).c3m3, 1e6*PSABWZernike(i).c40, 1e6*PSABWZernike(i).c42, 1e6*PSABWZernike(i).c4m2);
end;
fprintf(fid,' \n');
    

fprintf(fid,'-------------------------------------------  SECTION 12: SEGMENT CLOCKING RELATIVE TO M1 --------------------------------------------------\n');
fprintf(fid,'This section lists - for each segment type - the angle, measured about PSACRS_Z (or, in the PSACRS_XY plane), from the intersection between planes PSACRS_XY and M1CRS_RZ to PSACRS_X\n');
fprintf(fid,'12A: CLOCKING OF SEGMENTS RELATIVE TO M1 RADIAL\n');
fprintf(fid,'type#\t  t (degrees)\n');
for i=1:n_segments    
    fprintf(fid,'%3i\t%13.9f\n',i,180*AngleFromXtoM1RPlane_PSA(i)/pi);
end;


fprintf(fid,' \n');


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
fprintf(fid,'TMT M1 GRAVITY ORIENTATION DATA - Eric Ponslet - %s\n',datestr(now));
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


%################################################################## CANON SECTION ####################################################################

dt=cputime-tstart;
disp(['T=' num2str(dt) '  Computing and outputing ZYZ Euler Angles of the PSA systems, relative to M1 (all angles in Radians)']);

fid=fopen('EulerAngles.txt','wt');
fprintf(fid,'TMT M1 PSA ZYZ EULER ANGLE DATA - Eric Ponslet - %s\n',datestr(now));
fprintf(fid,'Seg    EulerA1(Z)   EulerA2(Y*)  EulerA3(Z**)     Theta         Phi\n');

% loop on all segments in sector A
for i=1:n_segments
    
    R=RPSA_M1(:,:,i);
    
    
    
    % compute ZYZ euler angles that express the PSA system in the M1 system; rotation matrix defining PSA is RPSA_M1
    EulerAnglesPSA_M1(1,i)= atan2( -1.0*R(2,3), -1.0*R(1,3) );
    EulerAnglesPSA_M1(3,i)= atan2( -1.0*R(3,2), R(3,1) );
    EulerAnglesPSA_M1(2,i)= atan2( R(3,2)/sin(EulerAnglesPSA_M1(3,i)), R(3,3) );
    
    CanonPhi(i) = atan(dZdR(i));
    CanonTheta(i) = Tctr_M1(i);
    fprintf(fid,'%3i %13.9f %13.9f %13.9f   %13.9f %13.9f\n',i,EulerAnglesPSA_M1(1,i),EulerAnglesPSA_M1(2,i),EulerAnglesPSA_M1(3,i),CanonTheta(i),CanonPhi(i));
end;

max_rdAlphaGamma = max(abs( (EulerAnglesPSA_M1(1,:) + EulerAnglesPSA_M1(3,:) ) ./ EulerAnglesPSA_M1(1,:) ));
max_rdAlphaTheta = max(abs( (EulerAnglesPSA_M1(1,:) - CanonTheta ) ./ CanonTheta ));
max_rdBetaPhi =    max(abs( (EulerAnglesPSA_M1(2,:) + CanonPhi ) ./ CanonPhi ));

fprintf(fid,'Max absolute relative difference between Alpha and -Gamma = %10.3e\n',max_rdAlphaGamma);
fprintf(fid,'Max absolute relative difference between Alpha and Theta = %10.3e\n',max_rdAlphaTheta);
fprintf(fid,'Max absolute relative difference between Beta and -Phi = %10.3e\n',max_rdBetaPhi);

disp(['         max abs diff between Alpha&-Gamma, Alpha&Theta, and Beta&Phi = ' num2str(max_rdAlphaGamma) '   ' num2str(max_rdAlphaTheta) '   ' num2str(max_rdBetaPhi)]);

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

% 2D plot of a few segments, with position of ESPCRS origin
for i=1:10
    figure('Position',[50+i 50+i scrsz(3)*2/3 scrsz(4)*2/3],'Name',['2D PSA plot of segment ' num2str(i)],'NumberTitle','off');
    hold on;
    text(0,0,num2str(i));
    plot([Vertex_PSA(1,:,i) Vertex_PSA(1,1,i)],[Vertex_PSA(2,:,i) Vertex_PSA(2,1,i)],'-k');
    for k=1:12
        plot(ESPOrigin_PSA(1,i,k),ESPOrigin_PSA(2,i,k),'*r');
    end;
    axis('equal');
end;

% % ----------------------------------------------------------------------------------------------------------------------------------------------------------------
% % 2D plot of array with segment numbers and extreme segments identified
% h=figure('Position',[50 50 scrsz(3)*2/3 scrsz(4)*2/3],'Name','2D M1 Array - Sector A - extreme segments','NumberTitle','off');
% axis([-(a+a/3) cos(30*pi/180)*RmaxBP+a/3 -a/3 RmaxBP+a/3]);
% axis('equal');
% for i=1:n_segments
%     hold on;
%     plot(Center_M1(1,i),Center_M1(2,i),'.k');
%     plot([Vertex_M1(1,:,i) Vertex_M1(1,1,i)],[Vertex_M1(2,:,i) Vertex_M1(2,1,i)],'-k');
%     text(Center_M1(1,i)+0.1,Center_M1(2,i),num2str(i));
% end;
% hold on; i=iMinDiam;  d=0.9; h1=plot([d*Vertex_M1(1,:,i)+(1-d)*Center_M1(1,i) d*Vertex_M1(1,1,i)+(1-d)*Center_M1(1,i)],[[d*Vertex_M1(2,:,i)+(1-d)*Center_M1(2,i) d*Vertex_M1(2,1,i)+(1-d)*Center_M1(2,i)]],'--g','LineWidth',2);
% hold on; i=iMaxDiam;  d=0.9; h2=plot([d*Vertex_M1(1,:,i)+(1-d)*Center_M1(1,i) d*Vertex_M1(1,1,i)+(1-d)*Center_M1(1,i)],[[d*Vertex_M1(2,:,i)+(1-d)*Center_M1(2,i) d*Vertex_M1(2,1,i)+(1-d)*Center_M1(2,i)]],'-g','LineWidth',2);
% hold on; i=iMinArea;  d=0.8; h3=plot([d*Vertex_M1(1,:,i)+(1-d)*Center_M1(1,i) d*Vertex_M1(1,1,i)+(1-d)*Center_M1(1,i)],[[d*Vertex_M1(2,:,i)+(1-d)*Center_M1(2,i) d*Vertex_M1(2,1,i)+(1-d)*Center_M1(2,i)]],'--b','LineWidth',2);
% hold on; i=iMaxArea;  d=0.8; h4=plot([d*Vertex_M1(1,:,i)+(1-d)*Center_M1(1,i) d*Vertex_M1(1,1,i)+(1-d)*Center_M1(1,i)],[[d*Vertex_M1(2,:,i)+(1-d)*Center_M1(2,i) d*Vertex_M1(2,1,i)+(1-d)*Center_M1(2,i)]],'-b','LineWidth',2);
% hold on; i=iMinIrreg; d=0.7; h5=plot([d*Vertex_M1(1,:,i)+(1-d)*Center_M1(1,i) d*Vertex_M1(1,1,i)+(1-d)*Center_M1(1,i)],[[d*Vertex_M1(2,:,i)+(1-d)*Center_M1(2,i) d*Vertex_M1(2,1,i)+(1-d)*Center_M1(2,i)]],'--r','LineWidth',2);
% hold on; i=iMaxIrreg; d=0.7; h6=plot([d*Vertex_M1(1,:,i)+(1-d)*Center_M1(1,i) d*Vertex_M1(1,1,i)+(1-d)*Center_M1(1,i)],[[d*Vertex_M1(2,:,i)+(1-d)*Center_M1(2,i) d*Vertex_M1(2,1,i)+(1-d)*Center_M1(2,i)]],'-r','LineWidth',2);
% legend([h1 h2 h3 h4 h5 h6],'Min Circum. circle','Max Circum. circle','Min Area', 'Max Area','Min Irregularity','Max Irregularity');
% xlabel('X_{M1} (m)');ylabel('Y_{M1} (m)');
% title('Segment Outlines projected into XY_{M1} plane');
% % plot 0/60/120 lines
% for theta=pi/6:2*pi/6:4.01*pi/6
%     hold on;polar([theta theta],[0 RmaxBP],':k');
% end;
% saveas(h,'2D_SectorA_Extremes.fig');
% 
% % ----------------------------------------------------------------------------------------------------------------------------------------------------------------
% % 2D plot of complete array with segment numbers, vertex numbers, and edge sensor numbers
% h=figure('Position',[60 60 scrsz(3)*2/3 scrsz(4)*2/3],'Name','2D M1 Array - All Sectors - with numbering','NumberTitle','off');
% axis([-1*(cos(30*pi/180)*RmaxBP+a/3) (cos(30*pi/180)*RmaxBP+a/3) -1*(RmaxBP+a/3) (RmaxBP+a/3)]);
% axis('equal');
% for i=1:6  % loop on 6 sectors
%     for j=1:n_segments  % loop on all segments
%         iseg=j+(i-1)*n_segments;  % generalized segment number
%         hold on;
%         if mod(i,2)==0  % for every other segment:
%             C=[1 1 1]*0.9;  % set color to gray
%         else
%             C=[1 0.8 0.8];  % set color to pale pink
%         end;
%         patch([Vertex_M1(1,:,iseg) Vertex_M1(1,1,iseg)],[Vertex_M1(2,:,iseg) Vertex_M1(2,1,iseg)],C);  % plot shaded segment with black edge
%         plot(Center_M1(1,iseg),Center_M1(2,iseg),'.k');  % place dot at centers of segments
%         ID=[SectorID(i) num2str(j)];  % segment number string
%         text(Center_M1(1,iseg)+0.1,Center_M1(2,iseg),ID);  % label segments
%     end;
%     if debugParameters % label vertices and edge sensors for segments #1 and 2 only in each sector (labeling all results in unbearably slow plots), unless in debugging mode (much smaller array)
%         jm=n_segments; 
%     else
%         jm=2; 
%     end;
%     for j=1:jm 
%         iseg=j+(i-1)*n_segments;  % generalized segment number
%         d=0.15;
%         plot(Vertex_M1(1,:,iseg),Vertex_M1(2,:,iseg),'.b');  % place dots at vertices
%         text(Vertex_M1(1,:,iseg)-d*(Vertex_M1(1,:,iseg)-Center_M1(1,iseg)),Vertex_M1(2,:,iseg)-d*(Vertex_M1(2,:,iseg)-Center_M1(2,iseg)),num2str([1:6]'),...
%                  'HorizontalAlignment','Center','VerticalAlignment','Middle','Color','b');  % label vertices
%         for k=1:12  % loop on all edge sensors
%             dR=ESOrigin_M1(1:2,iseg,k)-Center_M1(1:2,iseg);  % radial vector from segment center to edge sensor
%             dR=dR/norm(dR);  % normalize to unit length
%             plot(ESOrigin_M1(1,iseg,k)-Gap/2*dR(1),ESOrigin_M1(2,iseg,k)-Gap/2*dR(2),'.g');  % place dots for edges sensor markers (locations of dots are NOT representative)
%             text(ESOrigin_M1(1,iseg,k)-(Gap/2+d*a)*dR(1),ESOrigin_M1(2,iseg,k)-(Gap/2+d*a)*dR(2),num2str(k),'HorizontalAlignment','Center','VerticalAlignment','Middle','Color','g');  % label edge sensors
%         end;
%     end;
% end;
% xlabel('X_{M1} (m)');ylabel('Y_{M1} (m)');
% title('Segment Outlines projected into XY_{M1} plane - with segments (black), vertices (blue), and edge sensors (green) identified');
% % plot 0/60/120 lines
% for theta=pi/6:2*pi/6:4.01*pi/6
%     hold on;polar([theta theta],[0 RmaxBP],':k');
% end;
% saveas(h,'2D_AllSectors_Numbered.fig');
% 
% % ----------------------------------------------------------------------------------------------------------------------------------------------------------------
% % 2D plot of array with BFRH statistics: clocking angle, irregularity, and radius color code and corresponding histograms
% h=figure('Position',[70 70 scrsz(3)*2/3 scrsz(4)*2/3],'Name','Statistics: BFRH Clocking, Irregularity & Radius','NumberTitle','off');
% subplot(231);
%   %axis([-(a+a/3) cos(30*pi/180)*RmaxBP+a/3 -a/3 RmaxBP+a/3]);
%   plot(0,0,'.k');
%   for i=1:n_segments
%       hold on;
%       fill([Vertex_M1(1,:,i) Vertex_M1(1,1,i)],[Vertex_M1(2,:,i) Vertex_M1(2,1,i)],1000*Clocking(i));
%   end;
%   axis('equal');
%   xlabel('X_{M1} (m)');ylabel('Y_{M1} (m)');
%   title(': Dist. of Clocking Angles (mrad)');
%   colorbar('vert');
% subplot(234);
%   hist(Clocking*1000,15);
%   xlabel('Clocking Angle (mrad)');ylabel('count');
%   title('Dist. of Clocking Angles (mrad)');
% subplot(232);
%   %axis([-(a+a/3) cos(30*pi/180)*RmaxBP+a/3 -a/3 RmaxBP+a/3]);
%   plot(0,0,'.k');
%   for i=1:n_segments
%       hold on;
%       fill([Vertex_M1(1,:,i) Vertex_M1(1,1,i)],[Vertex_M1(2,:,i) Vertex_M1(2,1,i)],1000*Irregularity(i));
%   end;
%   axis('equal');
%   xlabel('X_{M1} (m)');ylabel('Y_{M1} (m)');
%   title('Dist. of RMS Irregularity (mm)');
%   colorbar('vert');
% subplot(235);
%   hist(Irregularity*1000,15);
%   xlabel('RMS Irregularity (mm)');ylabel('count');
%   title('Dist. of RMS Irregularity (mm)');
% subplot(233);
%   %axis([-(a+a/3) cos(30*pi/180)*RmaxBP+a/3 -a/3 RmaxBP+a/3]);
%   plot(0,0,'.k');
%   for i=1:n_segments
%       hold on;
%       fill([Vertex_M1(1,:,i) Vertex_M1(1,1,i)],[Vertex_M1(2,:,i) Vertex_M1(2,1,i)],Radius(i));
%   end;
%   axis('equal');
%   xlabel('X_{M1} (m)');ylabel('Y_{M1} (m)');
%   title('Dist. of Best Fit Radius (m)');
%   colorbar('vert');
% subplot(236);
%   hist(Radius,15);
%   xlabel('Best Fit Radius (m)');ylabel('count');
%   title('Dist. of Best Fit Radius (m)');
% saveas(h,'Stats_Clocking_Irregularity_Radius.fig');
% 
% % ----------------------------------------------------------------------------------------------------------------------------------------------------------------
% % 2D plot of array with area and circumscribed radius color code and corresponding histograms
% h=figure('Position',[80 80 scrsz(3)*2/3 scrsz(4)*2/3],'Name','Statistics: Segment Area & Radius of Circumscribed Circle','NumberTitle','off');
% subplot(221);
%   %axis([-(a+a/3) cos(30*pi/180)*RmaxBP+a/3 -a/3 RmaxBP+a/3]);
%   plot(0,0,'.k');
%   for i=1:n_segments
%       hold on;
%       fill([Vertex_M1(1,:,i) Vertex_M1(1,1,i)],[Vertex_M1(2,:,i) Vertex_M1(2,1,i)],HexagonArea(i));
%   end;
%   axis('equal');
%   xlabel('X_{M1} (m)');ylabel('Y_{M1} (m)');
%   title('Dist. of Segment Area ({m^2})');
%   colorbar('vert');
% subplot(222);
%   %axis([-(a+a/3) cos(30*pi/180)*RmaxBP+a/3 -a/3 RmaxBP+a/3]);
%   plot(0,0,'.k');
%   for i=1:n_segments
%       hold on;
%       fill([Vertex_M1(1,:,i) Vertex_M1(1,1,i)],[Vertex_M1(2,:,i) Vertex_M1(2,1,i)],CircumRadii(i));
%   end;
%   axis('equal');
%   xlabel('X_{M1} (m)');ylabel('Y_{M1} (m)');
%   title('Dist. of Radius of Circumscribed Circle (m)');
%   colorbar('vert');
% subplot(223);
%   hist(HexagonArea,15);
%   xlabel('Segment Area ({m^2})');ylabel('count');
%   title('Dist. of Segment Area ({m^2})');
% subplot(224);
%   hist(CircumRadii,15);
%   xlabel('Radius of Circumscribed Circle (m)');ylabel('count');
%   title('Dist. of Radius of Circumscribed Circle (m)');
% saveas(h,'Stats_Area_Circum.fig');
% 
% % ----------------------------------------------------------------------------------------------------------------------------------------------------------------
% % 3D plot of segment outlines and PSA and Edge sensor coordinate systems - sector A
% h=figure('Position',[90 90 scrsz(3)*2/3 scrsz(4)*2/3],'Name','3D Array - Sector A - Gapped and Ungapped, with PSA and ES coordinate systems','NumberTitle','off');
% hold off;
% plot3(Center_M1(1,1:n_segments),Center_M1(2,1:n_segments),Center_M1(3,1:n_segments),'.k');
% axis([-(a/2+a/3) cos(30*pi/180)*RmaxBP+a/3 -a/3 RmaxBP+a/3]);
% axis('equal');
% % plot segment outlines
% for i=1:n_segments
%     hold on;
%     plot3([UngappedVertex_M1(1,:,i) UngappedVertex_M1(1,1,i)],[UngappedVertex_M1(2,:,i) UngappedVertex_M1(2,1,i)],[UngappedVertex_M1(3,:,i) UngappedVertex_M1(3,1,i)],':k');  % plot sector A, ungapped contours
%     plot3([Vertex_M1(1,:,i) Vertex_M1(1,1,i)],[Vertex_M1(2,:,i) Vertex_M1(2,1,i)],[Vertex_M1(3,:,i) Vertex_M1(3,1,i)],'k');  % plot sector A, gapped contours
%     plot3([VertexOpt_M1(1,:,i) VertexOpt_M1(1,1,i)],[VertexOpt_M1(2,:,i) VertexOpt_M1(2,1,i)],[VertexOpt_M1(3,:,i) VertexOpt_M1(3,1,i)],'--k');  % plot sector A, chamfered contours
% end;
% % plot PSA coordinate systems
% PSAaxis=0.3*a;
% for i=1:n_segments
%     for k=1:3  % loop on three axes of PSA coordinate system
%         plot3([Center_M1(1,i) Center_M1(1,i)+PSAaxis*RPSA_M1(1,k,i)],[Center_M1(2,i) Center_M1(2,i)+PSAaxis*RPSA_M1(2,k,i)],[Center_M1(3,i) Center_M1(3,i)+PSAaxis*RPSA_M1(3,k,i)],'b');
%     end;
% 
% end;
% % plot edge sensor coordinate systems
% ESaxis=0.3*a;
% for i=1:n_segments
%     for j=1:12  % loop on all edge sensors of segment 'i'
%         plot3(ESOrigin_M1(1,i,j),ESOrigin_M1(2,i,j),ESOrigin_M1(3,i,j),'.g');  % center of ES coordinate system
%         for k=1:3  % loop on three axes of ES coordinate system
%             plot3([ESOrigin_M1(1,i,j) ESOrigin_M1(1,i,j)+ESaxis*RES_M1(1,k,i,j)],[ESOrigin_M1(2,i,j) ESOrigin_M1(2,i,j)+ESaxis*RES_M1(2,k,i,j)],[ESOrigin_M1(3,i,j) ESOrigin_M1(3,i,j)+ESaxis*RES_M1(3,k,i,j)],'g');
%         end;
%     end;
% end;
% % plot the origins of the edge sensor pocket coordinate systems
% for i=1:n_segments
%     for j=1:12  % loop on all edge sensors of segment 'i'
%         plot3(ESPOrigin_M1(1,i,j),ESPOrigin_M1(2,i,j),ESPOrigin_M1(3,i,j),'.r');  % center of ES coordinate system
%     end;
% end;
% 
% % plot M1CRS Coordinate system
% LM1CRS=0.8*a;
% LM1CRStxt=0.9*a;
% plot3(0,0,0,'.k');  % center of M1CRS
% plot3([0 LM1CRS],[0 0],[0 0],'-k');  % X axis of M1CRS
% text(LM1CRStxt,0,0,'X','HorizontalAlignment','Center','VerticalAlignment','Middle','Color','k');
% plot3([0 0],[0 LM1CRS],[0 0],'-k');  % Y axis of M1CRS
% text(0,LM1CRStxt,0,'Y','HorizontalAlignment','Center','VerticalAlignment','Middle','Color','k');
% plot3([0 0],[0 0],[0 LM1CRS],'-k');  % Z axis of M1CRS
% text(0,0,LM1CRStxt,'Z','HorizontalAlignment','Center','VerticalAlignment','Middle','Color','k');
% 
% xlabel('X_{M1} (m)');ylabel('Y_{M1} (m)');zlabel('Z_{M1} (m)');
% title('Segment Outlines in XYZ_{M1} (black; dotted=ungapped, solid=gapped, dashed=chamfered), PSACRS (blue), and ESCRS (green)');
% saveas(h,'3D_SectorA_PSA_ES.fig');
% 
% 
% % ----------------------------------------------------------------------------------------------------------------------------------------------------------------
% % 3D plot of segment outlines and all fiducials and targets - sector A
% h=figure('Position',[90 90 scrsz(3)*2/3 scrsz(4)*2/3],'Name','3D Array - Sector A - Optical surface fiducials, probe points, and subcell alignment targets','NumberTitle','off');
% hold off;
% plot3(Center_M1(1,1:n_segments),Center_M1(2,1:n_segments),Center_M1(3,1:n_segments),'.k');
% axis([-(a/2+a/3) cos(30*pi/180)*RmaxBP+a/3 -a/3 RmaxBP+a/3]);
% axis('equal');
% % plot segment outlines
% for i=1:n_segments
%     hold on;
%     plot3([Vertex_M1(1,:,i) Vertex_M1(1,1,i)],[Vertex_M1(2,:,i) Vertex_M1(2,1,i)],[Vertex_M1(3,:,i) Vertex_M1(3,1,i)],'k');  % plot sector A, gapped contours
% end;
% % plot the optical surface fiducials
% for i=1:n_segments
%     for j=1:3  % loop on three fiducials
%         plot3(OSFiducials_M1(1,j,i),OSFiducials_M1(2,j,i),OSFiducials_M1(3,j,i),'.r');  % center of ES coordinate system
%     end;
% end;
% % plot the optical surface probe points
% for i=1:n_segments
%     for j=1:3  % loop on three fiducials
%         plot3(OSProbes_M1(1,j,i),OSProbes_M1(2,j,i),OSProbes_M1(3,j,i),'.b');  % center of ES coordinate system
%     end;
% end;
% % plot the subcell alignment targets
% for i=1:n_segments
%     for j=1:3  % loop on three fiducials
%         plot3(SubcellAlignmentTargets_M1(1,j,i),SubcellAlignmentTargets_M1(2,j,i),SubcellAlignmentTargets_M1(3,j,i),'.g');  % center of ES coordinate system
%     end;
% end;
% 
% % plot M1CRS Coordinate system
% LM1CRS=0.8*a;
% LM1CRStxt=0.9*a;
% plot3(0,0,0,'.k');  % center of M1CRS
% plot3([0 LM1CRS],[0 0],[0 0],'-k');  % X axis of M1CRS
% text(LM1CRStxt,0,0,'X','HorizontalAlignment','Center','VerticalAlignment','Middle','Color','k');
% plot3([0 0],[0 LM1CRS],[0 0],'-k');  % Y axis of M1CRS
% text(0,LM1CRStxt,0,'Y','HorizontalAlignment','Center','VerticalAlignment','Middle','Color','k');
% plot3([0 0],[0 0],[0 LM1CRS],'-k');  % Z axis of M1CRS
% text(0,0,LM1CRStxt,'Z','HorizontalAlignment','Center','VerticalAlignment','Middle','Color','k');
% 
% xlabel('X_{M1} (m)');ylabel('Y_{M1} (m)');zlabel('Z_{M1} (m)');
% title('Fiducials and Targets (Optical Surface Fiducials (red), Optical Surface Probe Points (blue), Subcell Alignment Targets (green)');
% saveas(h,'3D_SectorA_Fiducials.fig');



% ----------------------------------------------------------------------------------------------------------------------------------------------------------------
% % 3D plot of entire array with supports and actuators  and supports shown
% h=figure('Position',[100 100 scrsz(3)*2/3 scrsz(4)*2/3],'Name','3D Array - Complete - Gapped, with cell-PSA interfaces','NumberTitle','off');
% hold off
% % % plot segment centers
% % plot3(Center_M1(1,:),Center_M1(2,:),Center_M1(3,:),'.k');
% axis([-1*ceil(cos(30*pi/180)*RmaxBP*1.1) ceil(cos(30*pi/180)*RmaxBP*1.1) -1*ceil(RmaxBP*1.1) ceil(RmaxBP*1.1)]);
% axis equal;
% % ActL=0.4*a;  % length of line segments to represent actuators lines of action
% for j=1:6  % loop on 6 sectors
%     for i=1:n_segments
%         hold on;
%         iseg=i+(j-1)*n_segments;
%         % plot segment outlines
%         if rem(j,2)>0
%             plot3([Vertex_M1(1,:,iseg) Vertex_M1(1,1,iseg)],[Vertex_M1(2,:,iseg) Vertex_M1(2,1,iseg)],[Vertex_M1(3,:,iseg) Vertex_M1(3,1,iseg)],'k');  % plot segment
%         else
%             plot3([Vertex_M1(1,:,iseg) Vertex_M1(1,1,iseg)],[Vertex_M1(2,:,iseg) Vertex_M1(2,1,iseg)],[Vertex_M1(3,:,iseg) Vertex_M1(3,1,iseg)],':k');  % plot segment
%         end
%         
%         hold on;
%         % plot SSA supports
%         CIcenter_M1=1/3*(PSASideAAPCenter_M1(:,1,iseg)+PSASideAAPCenter_M1(:,2,iseg)+PSASideAAPCenter_M1(:,3,iseg)); % center (mean) of cell interface points (for display only)
%         for k=1:3
% %             plot3(PSASideAAPCenter_M1(1,k,iseg),PSASideAAPCenter_M1(2,k,iseg),PSASideAAPCenter_M1(3,k,iseg),'.b');
%             if k==1
%                 style = '-b';
%             else
%                 style = ':b';
%             end;
%             plot3([CIcenter_M1(1) PSASideAAPCenter_M1(1,k,iseg)],[CIcenter_M1(2) PSASideAAPCenter_M1(2,k,iseg)],[CIcenter_M1(3) PSASideAAPCenter_M1(3,k,iseg)],style); % plot lines connecting each of the 3 cell interface points to the mean of those three points (for ease of visualization)
%         end;
%         % plot actuators
% %         plot3(ActCOR_M1(1,iseg),ActCOR_M1(2,iseg),ActCOR_M1(3,iseg),'.b');  % center of rotation
% %         for k=1:3
% %             plot3(ActOrigin_M1(1,k,iseg),ActOrigin_M1(2,k,iseg),ActOrigin_M1(3,k,iseg),'.b');  % ACTCRS origin
% %             plot3([ActOrigin_M1(1,k,iseg) ActOrigin_M1(1,k,iseg)+ActL*ActLOA_M1(1,k,iseg)],...
% %                   [ActOrigin_M1(2,k,iseg) ActOrigin_M1(2,k,iseg)+ActL*ActLOA_M1(2,k,iseg)],...
% %                   [ActOrigin_M1(3,k,iseg) ActOrigin_M1(3,k,iseg)+ActL*ActLOA_M1(3,k,iseg)],'-b');  % lines of action
% %             plot3([ActOrigin_M1(1,k,iseg) ActCOR_M1(1,iseg)],...
% %                   [ActOrigin_M1(2,k,iseg) ActCOR_M1(2,iseg)],...
% %                   [ActOrigin_M1(3,k,iseg) ActCOR_M1(3,iseg)],':b');  % plot lines connecting actuator center of rotation to each of the 3 actuator centers (for ease of visualization)
% %         end;
%         % plot cell top chord
%         for k=1:3
%             kk = k+1;
%             if kk > 3
%                 kk=1;
%             end;
%             plot3([CellNode_M1(1,k,iseg); CellNode_M1(1,kk,iseg)],[CellNode_M1(2,k,iseg); CellNode_M1(2,kk,iseg)],[CellNode_M1(3,k,iseg); CellNode_M1(3,kk,iseg)],'-r');
%             plot3([CellTopChordThirdPoints_M1(1,k,iseg); CellSideAAPCenter_M1(1,k,iseg)],...
%                 [CellTopChordThirdPoints_M1(2,k,iseg); CellSideAAPCenter_M1(2,k,iseg)],...
%                 [CellTopChordThirdPoints_M1(3,k,iseg); CellSideAAPCenter_M1(3,k,iseg)],'-r');
%             plot3(CellSideAAPCenter_M1(1,k,iseg),...
%                 CellSideAAPCenter_M1(2,k,iseg),...
%                 CellSideAAPCenter_M1(3,k,iseg),'.r');
%         end
%         
%     end;
% end;
% % plot M1CRS Coordinate system
% LM1CRS=0.8*a;
% LM1CRStxt=0.9*a;
% plot3(0,0,0,'.k');  % center of M1CRS
% plot3([0 LM1CRS],[0 0],[0 0],'-k');  % X axis of M1CRS
% text(LM1CRStxt,0,0,'X','HorizontalAlignment','Center','VerticalAlignment','Middle','Color','k');
% plot3([0 0],[0 LM1CRS],[0 0],'-k');  % Y axis of M1CRS
% text(0,LM1CRStxt,0,'Y','HorizontalAlignment','Center','VerticalAlignment','Middle','Color','k');
% plot3([0 0],[0 0],[0 LM1CRS],'-k');  % Z axis of M1CRS
% text(0,0,LM1CRStxt,'Z','HorizontalAlignment','Center','VerticalAlignment','Middle','Color','k');
% 
% xlabel('X_{M1} (m)');ylabel('Y_{M1} (m)');zlabel('Z_{M1} (m)');
% title('Segment Outlines in XYZ_{M1} (black; dashed for sectors B,D,E), PSA support frame (blue; dashed for leg // to Y_P_S_A), and cell top chord (red)');
% saveas(h,'3D_AllSectors_Act_Interf.fig');

% ----------------------------------------------------------------------------------------------------------------------------------------------------------------
% % 3D plot of entire array with PSA and M1 Systems Shown
% h=figure('Position',[100 100 scrsz(3)*2/3 scrsz(4)*2/3],'Name','3D Array - Complete - Gapped','NumberTitle','off');
% hold off
% axis([-1*ceil(cos(30*pi/180)*RmaxBP*1.1) ceil(cos(30*pi/180)*RmaxBP*1.1) -1*ceil(RmaxBP*1.1) ceil(RmaxBP*1.1)]);
% axis equal;
% ActL=0.4*a;  % length of line segments to represent actuators lines of action
% for j=1:6  % loop on 6 sectors
%     for i=1:n_segments
%         hold on;
%         iseg=i+(j-1)*n_segments;
%         % plot segment outlines
%         plot3([Vertex_M1(1,:,iseg) Vertex_M1(1,1,iseg)],[Vertex_M1(2,:,iseg) Vertex_M1(2,1,iseg)],[Vertex_M1(3,:,iseg) Vertex_M1(3,1,iseg)],'k');  % plot segment
%         hold on;
%         % plot PSA axes
%         LPSACRS=0.4*a;
%         c=Center_M1(:,iseg);
%         x=RPSA_M1(:,1,iseg)*LPSACRS + c;
%         y=RPSA_M1(:,2,iseg)*LPSACRS + c;
%         z=RPSA_M1(:,3,iseg)*LPSACRS + c;
%         plot3(c(1),c(2),c(3),'.k');  % center of PSACRS
%         plot3([c(1) x(1)],[c(2) x(2)],[c(3) x(3)],'-k');  % X axis of PSACRS
%         plot3([c(1) y(1)],[c(2) y(2)],[c(3) y(3)],'-k');  % Y axis of PSACRS
%         plot3([c(1) z(1)],[c(2) z(2)],[c(3) z(3)],'-k');  % Z axis of PSACRS
%     end;
% end;
% % plot M1CRS Coordinate system
% LM1CRS=0.8*a;
% LM1CRStxt=0.9*a;
% plot3(0,0,0,'.k');  % center of M1CRS
% plot3([0 LM1CRS],[0 0],[0 0],'-k');  % X axis of M1CRS
% %text(LM1CRStxt,0,0,'X','HorizontalAlignment','Center','VerticalAlignment','Middle','Color','k');
% plot3([0 0],[0 LM1CRS],[0 0],'-k');  % Y axis of M1CRS
% %text(0,LM1CRStxt,0,'Y','HorizontalAlignment','Center','VerticalAlignment','Middle','Color','k');
% plot3([0 0],[0 0],[0 LM1CRS],'-k');  % Z axis of M1CRS
% %text(0,0,LM1CRStxt,'Z','HorizontalAlignment','Center','VerticalAlignment','Middle','Color','k');
% 
% xlabel('X_{M1} (m)');ylabel('Y_{M1} (m)');zlabel('Z_{M1} (m)');
% title('Segment Outlines in XYZ_{M1} (black), Actuator Centers & Lines of Action (red), Actuation Center of Rotation (red), and Cell Interface Points (blue)');
% saveas(h,'3D_AllSectors_Act_Interf.fig');
% 
% 
% 
% % ----------------------------------------------------------------------------------------------------------------------------------------------------------------
% % 2D plot of segment outlines in TEMP frame - REDUCED by substracting a regular hexagon of size a_reduced;
% h=figure('Position',[110 110 scrsz(3)*2/3 scrsz(4)*2/3],'Name','Segment Outlines in TEMP frame','NumberTitle','off');
% angles=[0:5]*pi/3;
% a_reduced=(100-50*(Rout/Opt_k)^2)/100*a;
% DVx=a_reduced*cos(angles);
% DVy=a_reduced*sin(angles);
% axis([-1 1 -1 1]*(MaxDiam/2-a_reduced)*1.1);
% axis('equal');
% for i=1:n_segments
%     hold on;
%     plot([Vertex_Temp(1,:,i)-DVx Vertex_Temp(1,1,i)-a_reduced],[Vertex_Temp(2,:,i)-DVy Vertex_Temp(2,1,i)],'-y');
% end;
% hold on; i=iMinDiam;  h1=plot([Vertex_Temp(1,:,i)-DVx Vertex_Temp(1,1,i)-a_reduced],[Vertex_Temp(2,:,i)-DVy Vertex_Temp(2,1,i)],'--g','LineWidth',2);
% hold on; i=iMaxDiam;  h2=plot([Vertex_Temp(1,:,i)-DVx Vertex_Temp(1,1,i)-a_reduced],[Vertex_Temp(2,:,i)-DVy Vertex_Temp(2,1,i)],'-g','LineWidth',2);
% hold on; i=iMinArea;  h3=plot([Vertex_Temp(1,:,i)-DVx Vertex_Temp(1,1,i)-a_reduced],[Vertex_Temp(2,:,i)-DVy Vertex_Temp(2,1,i)],'--b','LineWidth',2);
% hold on; i=iMaxArea;  h4=plot([Vertex_Temp(1,:,i)-DVx Vertex_Temp(1,1,i)-a_reduced],[Vertex_Temp(2,:,i)-DVy Vertex_Temp(2,1,i)],'-b','LineWidth',2);
% hold on; i=iMinIrreg; h5=plot([Vertex_Temp(1,:,i)-DVx Vertex_Temp(1,1,i)-a_reduced],[Vertex_Temp(2,:,i)-DVy Vertex_Temp(2,1,i)],'--r','LineWidth',2);
% hold on; i=iMaxIrreg; h6=plot([Vertex_Temp(1,:,i)-DVx Vertex_Temp(1,1,i)-a_reduced],[Vertex_Temp(2,:,i)-DVy Vertex_Temp(2,1,i)],'-r','LineWidth',2); 
% legend([h1 h2 h3 h4 h5 h6],'Min Circum. circle','Max Circum. circle','Min Area', 'Max Area','Min Irregularity','Max Irregularity');
% xlabel('X_{Temp} (m)');ylabel('Y_{Temp} (m)');
% title('Segment Outlines projected into XY_{Temp} plane (reduced)');
% saveas(h,'Outlines_Temp_top.fig');
% 
% % ----------------------------------------------------------------------------------------------------------------------------------------------------------------
% % 2D plot of segment outlines in PSA frame - REDUCED by substracting a regular hexagon of size a_reduced;
% h=figure('Position',[120 120 scrsz(3)*2/3 scrsz(4)*2/3],'Name','Segment Outlines in PSA frame','NumberTitle','off');
% angles=[0:5]*pi/3;
% DVx=a_reduced*cos(angles);
% DVy=a_reduced*sin(angles);
% axis([-1 1 -1 1]*(MaxDiam/2-a_reduced)*1.1);
% axis('equal');
% for i=1:n_segments
%     hold on;
%     plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-a_reduced],[Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)],'-y');
% end;
% hold on; i=iMinDiam;  h1=plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-a_reduced],[Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)],'--g','LineWidth',2);
% hold on; i=iMaxDiam;  h2=plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-a_reduced],[Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)],'-g','LineWidth',2);
% hold on; i=iMinArea;  h3=plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-a_reduced],[Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)],'--b','LineWidth',2);
% hold on; i=iMaxArea;  h4=plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-a_reduced],[Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)],'-b','LineWidth',2);
% hold on; i=iMinIrreg; h5=plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-a_reduced],[Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)],'--r','LineWidth',2);
% hold on; i=iMaxIrreg; h6=plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-a_reduced],[Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)],'-r','LineWidth',2); 
% legend([h1 h2 h3 h4 h5 h6],'Min Circum. circle','Max Circum. circle','Min Area', 'Max Area','Min Irregularity','Max Irregularity');
% xlabel('X_{PSA} (m)');ylabel('Y_{PSA} (m)');
% title('Segment Outlines projected into XY_{PSA} plane (reduced)');
% saveas(h,'Outlines_PSA_top.fig');
% 
% % ----------------------------------------------------------------------------------------------------------------------------------------------------------------
% % XY and XZ plots of segment outlines in PSA frame - REDUCED by substracting a
% % regular hexagon of size a_reduced;
% h=figure('Position',[130 130 scrsz(3)*2/3 scrsz(4)*2/3],'Name','Segment Profiles in PSA frame','NumberTitle','off');
% angles=[0:5]*pi/3;
% DVx=a_reduced*cos(angles);
% DVy=a_reduced*sin(angles);
% DVz=zeros(size(DVx));
% subplot(211);
%   for i=1:n_segments
%       hold on;
%       plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-DVx(1)],[Vertex_PSA(3,:,i)-DVz Vertex_PSA(3,1,i)],'-y');
%   end;
%   hold on; i=iMinDiam;  h1=plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-DVx(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'--g','LineWidth',2);
%   hold on; i=iMaxDiam;  h2=plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-DVx(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'-g','LineWidth',2);
%   hold on; i=iMinArea;  h3=plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-DVx(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'--b','LineWidth',2);
%   hold on; i=iMaxArea;  h4=plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-DVx(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'-b','LineWidth',2);
%   hold on; i=iMinIrreg; h5=plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-DVx(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'--r','LineWidth',2);
%   hold on; i=iMaxIrreg; h6=plot([Vertex_PSA(1,:,i)-DVx Vertex_PSA(1,1,i)-DVx(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'-r','LineWidth',2); 
%   legend([h1 h2 h3 h4 h5 h6],'Min Circum. circle','Max Circum. circle','Min Area', 'Max Area','Min Irregularity','Max Irregularity');
%   xlabel('X_{PSA} (m)');ylabel('Z_{PSA} (m)');
%   title('Segment Outlines projected into XZ_{PSA} plane (reduced)');
% subplot(212);
%   for i=1:n_segments
%       hold on;
%       plot([Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)-DVy(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'-y');
%   end;
%   hold on; i=iMinDiam;  h1=plot([Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)-DVy(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'--g','LineWidth',2);
%   hold on; i=iMaxDiam;  h2=plot([Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)-DVy(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'-g','LineWidth',2);
%   hold on; i=iMinArea;  h3=plot([Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)-DVy(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'--b','LineWidth',2);
%   hold on; i=iMaxArea;  h4=plot([Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)-DVy(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'-b','LineWidth',2);
%   hold on; i=iMinIrreg; h5=plot([Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)-DVy(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'--r','LineWidth',2);
%   hold on; i=iMaxIrreg; h6=plot([Vertex_PSA(2,:,i)-DVy Vertex_PSA(2,1,i)-DVy(1)],[Vertex_PSA(3,:,i) Vertex_PSA(3,1,i)],'-r','LineWidth',2); 
%   legend([h1 h2 h3 h4 h5 h6],'Min Circum. circle','Max Circum. circle','Min Area', 'Max Area','Min Irregularity','Max Irregularity');
%   xlabel('Y_{PSA} (m)');ylabel('Z_{PSA} (m)');
%   title('Segment Outlines projected into YZ_{PSA} plane (reduced)');
% saveas(h,'Outlines_PSA_side.fig');

dt=cputime-tstart;
disp(['T=' num2str(dt) '  Execution completed']);
