% example for concatenating two grids.

% in this example, these two sweeps have the same input values except for K
% (fifth dimension)

A = load('202010131551SweepEndpoints.mat');
B = load('202010131611SweepEndpoints.mat');

fields = fieldnames(A);

for i = 1:length(fields)
   fn = fields{i};
   C.(fn) = cat(4,A.(fn),B.(fn));
end

save([date_prefix('yyyymmddHHMM'),'SweepEndpoints.mat'],'-struct','C');