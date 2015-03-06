% Parameter Estimation and Inverse Problems, 2nd edition, 2011
% by R. Aster, B. Borchers, C. Thurber
% newlist=remove(list,entry);
%
function newlist=remove(list,entry);
newlist=list(find(list-entry*ones(size(list))));
