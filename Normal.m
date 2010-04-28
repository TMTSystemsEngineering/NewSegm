function V=Normal(X,Y,k,K);
% function V=Normal(X,Y,k,K);
% Eric Ponslet, March 2010, elp@ericandlucie.com
% This function calculates the components of a vector V (in M1CRS) that is normal to the optical surface (defined by geometric parameters k, and K) 
% at point P(X,Y,Z*) (in M1CRS); 
% * note that Z is implicit such that the point P lies in the optical surface
% 
% INPUTS: 
%   * X (scalar): X coordinate (in M1CRS) of point P
%   * Y (scalar): Y coordinate (in M1CRS) of point P
%   * k (scalar): paraxial radius of curvature of the optical surface
%   * K (scalar): conic constant of the optical surface
% OUTPUTS:
%   * V (3x1, real): M1CRS components of a unit vector normal to the optical surface at point P

R=sqrt(X^2+Y^2);
dZdR=R/sqrt(k^2-(K+1)*R^2);  % slope of optical surface at point P
SlopeNormal=-1/dZdR;  % slope of the NORMAL to the optical surface at point P
RV=-1/sqrt(1+SlopeNormal^2);  % radial component of unit vector normal to the opticval surface at point Pj
ZV=-SlopeNormal/sqrt(1+SlopeNormal^2);  % Z component of unit vector normal to the optical surface at point Pj
V=[RV*X/R; RV*Y/R; ZV]; % normal to optical surface at Pj (positive toward the stars), expressed in cartesian system M1CRS
