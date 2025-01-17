function ux = groupX(u,numPar)
% Group input u by X (all of the x values are next to each other)
% Assumes original matrix was x rows, y columns and input is grouped by y

% Form the original matrix
ux = reshape(u,numPar.ny,numPar.nx); % y values are on rows (y grouping)
ux = transpose(ux); % Flip
ux = ux(:);