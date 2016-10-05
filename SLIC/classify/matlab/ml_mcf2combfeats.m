function [combfeats,combclass,combcellidx]=ml_mcf2combfeats(features)
%ML_MCF2COMBFEATS Convert MCF to combined feature matrix.
%   COMBFEATS = ML_MCF2COMBFEATS(FEATURES) returns the combined feature
%   matrix from the cell array of feature matrices FEATURES.
%   
%   [COMBFEATS,COMBCLASS,COMBCELLIDX] = ML_MCF2COMBFEATS(...) also returns
%   the combined class labels and cell indices.
%   
%   See also ML_COMBFEATS2MCF

%   AUG-07-2004 Initial write T. Zhao
%   31-OCT-2004 Modified T. Zhao
%       - add comments
%       - change function name tz_combinefeats_mcf --> tz_mcf2combfeats
%   Copyright (c) Murphy Lab, Carnegie Mellon University



nclass=length(features);
combfeats=[];
combclass=[];
combcellidx=[];

for i=1:nclass
    combfeats=[combfeats;features{i}];
    ncell=size(features{i},1);
    combclass=[combclass;zeros(ncell,1)+i];
    combcellidx=[combcellidx;(1:ncell)'];
end
