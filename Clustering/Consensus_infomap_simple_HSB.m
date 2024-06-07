function consensusmap = Consensus_infomap_simple_HSB(assignmentsfile,mincol,minsize)
%assignmentsfile = 'rawassn.txt';
% Similar to Jonathan Powers 2011 Neuron paper
% Create initial consensus by accepting all assignments at the mincol threshold and assigning unassigned nodes to their higher threshold assignments
if ~isnumeric(assignmentsfile)
    regularized = importdata(assignmentsfile);
else
    regularized = assignmentsfile;
end
consensusmap = regularized(:,mincol);

unassigned = find(consensusmap<1);
for unassignedindex = unassigned'
    thisassignments = regularized(unassignedindex,mincol:end);
    thisassignments(thisassignments<1) = [];
    if ~isempty(thisassignments)
        consensusmap(unassignedindex) = thisassignments(1);
    end
end


% Reapply size threshold

consensusvals = unique(consensusmap);
for val = consensusvals(:)'
    if nnz(consensusmap==val) < minsize
        toosmallinds = find(consensusmap==val);
        for smallindex = toosmallinds(:)'
            thisassignments = regularized(unassignedindex,mincol:end);
            thisassignments(thisassignments<1) = [];
            thisassignments(thisassignments==val) = [];
            if ~isempty(thisassignments)
                consensusmap(smallindex) = thisassignments(1);
            end
        end
        
    end
end


%%
% dlmwrite([assignmentsfile(1:end-4) '_consensus_mincol',num2str(mincol),'.txt'],consensusmap);

