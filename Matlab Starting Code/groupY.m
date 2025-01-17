function uy = groupY(u,numPar)
% Group input u by Y (all of the y values are next to each other)
% Assumes original matrix was x rows, y columns and currently grouped by x

% Form the original matrix
uy = reshape(u,numPar.nx,numPar.ny);

uy = transpose(uy);  % y values are on rows

uy = uy(:);



