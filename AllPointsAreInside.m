function answer = AllPointsAreInside( hexagonRadius, points )
% Checks that a set of 2D points are all inside of a regular hexagon of given
% radius, with two edges parallel to X.

% initialize
answer = 1;

% compute vertices of regular hexagon
vertices = zeros(2,6);
edgeVectors = zeros(2,6);
for i=1:6
    vertices(1,i) = hexagonRadius*cos(pi*(i-1)/3);
    vertices(2,i) = hexagonRadius*sin(pi*(i-1)/3);
end;

% compute edge vectors
for i=1:6
    iNext = i+1;
    if iNext>6
        iNext = 1;
    end;
    edgeVectors(:,i) = vertices(:,iNext)-vertices(:,i);
end;

% test each point against each edge and make sure it is to the left of the edge
for iPoint=1:size(points,2)
    for iEdge=1:6
        vec = points(:,iPoint) - vertices(:,iEdge);
        crossProduct = vec(1)*edgeVectors(2,iEdge) - vec(2)*edgeVectors(1,iEdge);
        if crossProduct > 0
            answer = 0;
        end;
    end;
end;

