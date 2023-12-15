function B = get_null_model(A,gamma,nulltype,varargin)
% Jiaxin Cindy Tu 2023.11.28
% This function takes a symmetric input nodal association matrix A and a
% resolution parameter gamma, as well as a nulltype function to state the
% null model type
% modified from single_layer_null_models.m from https://github.com/brain-networks/nonstandard_modularity_maximization.git
% Reference: Esfahlani, Farnaz Zamani, et al. "Modularity maximization as a flexible and generic framework for brain network exploratory analysis." Neuroimage 244 (2021): 118607.
%% Preparation
% check the matrix is symmetric
sz = size(A);
assert(length(sz)==2 & (sz(1)==sz(2)))
n = length(A);
if ~exist('gamma','var')||isempty(gamma)
    gamma = 1;
end

% add more input here for more complex models if needed
parseArgs=inputParser();
addParameter(parseArgs,'d',[])
parse(parseArgs,varargin{:});
d=parseArgs.Results.d;

switch nulltype
    case {'uniform','CPM'}
        P= ~eye(n)*gamma;
        B = A - P; 
        % reference(s):
        %
        %   Bazzi, M., Porter, M. A., Williams, S., McDonald, M., Fenn, D. J., &
        %   Howison, S. D. (2016). Community detection in temporal multilayer
        %   networks, with an application to correlation networks.
        %   Multiscale Modeling & Simulation, 14(1), 1-41.
        %
        %   Traag, V. A., Van Dooren, P., & Nesterov, Y. (2011). Narrow scope for
        %   resolution-limit-free community detection. Physical Review E, 84(1),
        %   016114.
    case {'NG','newman girvan'}
        deg = sum(A,2);
        twom = sum(deg);
        P = (deg*deg')/twom;
        B = A - P;
        % reference:
        %
        %   Newman, M. E., & Girvan, M. (2004). Finding and evaluating community
        %   structure in networks. Physical review E, 69(2), 026113.
    case {'minimum wiring'}
        m = nnz(A)/2;
        edges = find(triu(ones(n),1) > 0);
        di = d(edges);
        [~,ddx] = sort(di,'ascend');
        P = zeros(n);
        w = sort(A(triu(A > 0)),'descend');
        P(edges(ddx(1:m))) = w;
        P = P + P';
        B= sc - P;
    case {'signed_symmetric'}
        A = A.*~eye(n);
        Apos = A .*(A > 0);
        Aneg = -A.*(A < 0);
        kpos = sum(Apos,2);
        twompos = sum(kpos);
        kneg = sum(Aneg,2);
        twomneg = sum(kneg);
        Ppos = (kpos*kpos')/twompos;
        Pneg = (kneg*kneg')/twomneg;
        Bpos = (Apos - Ppos);
        Bneg = (Aneg - Pneg);
        
        B = Bpos/(twompos + twomneg) - Bneg/(twompos + twomneg);
        % reference(s):
        %
        %   Gómez, S., Jensen, P., & Arenas, A. (2009). Analysis of community
        %   structure in networks of correlated data. Physical Review E, 80(1),
        %   016114.
    case {'signed_asymmetric'}
        A = A.*~eye(n);
        Apos = A.*(A > 0);
        Aneg = -A.*(A < 0);
        kpos = sum(Apos,2);
        twompos = sum(kpos);
        kneg = sum(Aneg,2);
        twomneg = sum(kneg);
        Ppos = (kpos*kpos')/twompos;
        Pneg = (kneg*kneg')/twomneg;
        Bpos = (Apos - Ppos);
        Bneg = (Aneg - Pneg);
        B = Bpos/(twompos) - Bneg/(twompos + twomneg);
        % reference(s):
        %
        %   Rubinov, M., & Sporns, O. (2011). Weight-conserving characterization of
        %   complex functional brain networks. Neuroimage, 56(4), 2068-2079.
    case {'binary'}
        avg = mean(nonzeros(A));
        P = (A > 0)*avg;
        B = A - P;
        % reference(s):
        %
        %   Bassett, D. S., Owens, E. T., Porter, M. A., Manning, M. L., & Daniels,
        %   K. E. (2015). Extraction of force-chain network architecture in
        %   granular materials using community detection. Soft Matter, 11(14),
        %   2731-2744.
    otherwise
        error('invalid input for null model type')
end
end