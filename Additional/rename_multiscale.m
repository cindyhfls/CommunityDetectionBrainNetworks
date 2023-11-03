function newS = rename_multiscale(S)
a = unique(S(~isnan(S) & S>0));
newS = S;
for i = 1:length(a)
    newS(S==a(i)) = i;
end
end