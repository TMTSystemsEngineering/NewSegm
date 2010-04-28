function [Q,d]=Intersect(P,V,k,K);
% function [Q,d]=Intersect(P,V,k,K);
% Eric Ponslet, March 2010, elp@ericandlucie.com
% Returns the coordinates of a point Q in the M1 optical surface, along a line defined by vector V, originating from point P
%   * all inputs and outputs are defined in the M1 coordinate system
%   * the optical surface is defined by K and K, according to the equation (in global coordinates) Z= {k-sqrt[k^2-(K+1)*R^2]} / (K+1)
%   * this routine is derived from the more general routine "OpticalLocal.m"
%   * this routine can only process one point at a time (unlike OpticalLocal)
% Technique:
%   By expressing the equation of the optical surface, and that of the line from point P along vector V, and substtituting, one 
%   obtains a second order polynomial for the distance along vector V, from point P to point Q.  
%   That polynomial is solved numerically.
%
% INPUTS:
%   P (3x1, real): coordinates of point P in M1CRS
%   V (3x1, real): components of vector V in M1CRS (does not need to be normalized)
%   k (scalar): paraxial radius of curvature of the optical surface
%   K (scalar): conic constant of the optical surface
% OUTPUTS:
%   Q (3x1, real): coordinates of point Q in M1CRS, the intersection of a straight line from P along V with the optical surface
%   d (scalar): distance along V from point P to point Q

% normalize V
V=V/norm(V);

% calculate coefficients of 2nd order polynomial
a= (K+1)*V(3)^2 + V(1)^2 + V(2)^2;
b= 2*P(1)*V(1) + 2*P(2)*V(2) + 2*P(3)*V(3)*(K+1) - 2*k*V(3);
c= P(1)^2 + P(2)^2 + (K+1)*P(3)^2 - 2*k*P(3);

% calculate root (distance from point P along line V to point Q in the optical surface (Q=P+d*V)
d=(-b-sqrt(b.^2-4.*a.*c))./(2*a);

% calculate point Q
Q=P+d*V;

