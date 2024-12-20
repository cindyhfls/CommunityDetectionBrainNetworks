function [b,gammas]=eventSamples_HSB(A,n,varargin)
% modified from eventSamples.m from https://github.com/LJeub/HierarchicalConsensus.git
% delete the partition optimization step, only compute gammas

% eventSamples Compute multiresolution ensemble using event sampling.
%
% Syntax
%__________________________________________________________________________
%
%   [gammas]=eventSamples(A,n)
%
%   [gammas]=eventSamples(__,Name,Value)
%
%
% Description
%__________________________________________________________________________
%
%   [gammas]=eventSamples(A,n) computes event sampling ensemble with 'n'
%       partitions.
%
%   [gammas]=eventSamples(__,Name,Value) additionally customizes the
%       behavior of the function by e.g. using a different algorithm to
%       optimize modularity or using a different modularity-like quality
%       function.
%
%
% Input Arguments
%__________________________________________________________________________
%
%   A -- Adjacency matrix of the network
%
%   n -- Number of partitions to be generated
%
%
% Name-Value Pair Arguments
%__________________________________________________________________________
%
% Parameter names can be abbreviated and are not case sensitive.
%
%   'Optimizer' -- Function for finding "optimal" partitions for a
%                  modularity-like quality function with modularity matrix
%                  'B'.
%                  @(B) iterated_genlouvain(B,[],0,1,'moverandw') (default) |
%                  function handle
%
%   'Modularity' -- Function to compute modularity matrix of the form
%                   'B=modularity(A,gamma)' where 'A' is an adjacency
%                   matrix and 'gamma' is a resolution parameter. Note that
%                   the function assumes that 'B=(A+A')/2-gamma*P' for some matrix
%                   'P'.
%                   @modularity (default) | function handle
%
%   'GammaMinSamples' -- Number of partitions to sample at each iteration
%                        when estimating 'gamma_min' (see 'gammaRange' for
%                        more details)
%                        10 (default)| scalar
%
%
% Output Arguments
%__________________________________________________________________________
%
%
%   gammas -- Values of 'gamma' corresponding to each partition in 'S'.
%
% See Also hierarchicalConsensus, exponentialSamples

% Version: 1.1.1
% Date: Thu  8 Mar 2018 15:34:46 CET
% Author: Lucas Jeub
% Email: ljeub@iu.edu

parseArgs=inputParser();
checkFunction=@(x) isa(x,'function_handle');
addParameter(parseArgs,'Optimizer',@(B) iterated_genlouvain(B,[],0,1,...
    'moverandw'),checkFunction);
addParameter(parseArgs,'Modularity',@modularity,checkFunction)
addParameter(parseArgs,'GammaMinSamples',10,@(x) isnumeric(x) && isscalar(x))

parse(parseArgs,varargin{:});

mod_fun=parseArgs.Results.Modularity;
optimizer=parseArgs.Results.Optimizer;
gmin_samples=parseArgs.Results.GammaMinSamples;

P=mod_fun(A,1);
A=(A+A')./2;
P=A-P;
N=length(A);
for i=1:N
    A(i,i)=0;
    P(i,i)=0;
end

[gamma_min, ~]=gammaRange(A,'Modularity',mod_fun,'Samples',gmin_samples);

% get discrete events where interactions change sign
gamma_et=div_0(A,P);
[g_sample,~,ind]=unique(gamma_et);
g_sample=full(g_sample);
PS=sum(sum(P));
AS=sum(sum(A));
Pp=zeros(length(g_sample),1);
Ap=zeros(length(g_sample),1);
for i=1:length(g_sample)
    Pp(i)=sum(P(ind>=i));
    Ap(i)=sum(A(ind>=i));
end

Pp=full(Pp);
Ap=full(Ap);
b_sample=(g_sample.*(PS-Pp)-(AS-Ap))./(g_sample.*(PS-2*Pp)+2*Ap-AS);
b_min=(gamma_min*interp1(g_sample,PS-Pp,gamma_min,'next')-interp1(g_sample,AS-Ap,gamma_min,'next'))/...
    (gamma_min*interp1(g_sample,PS-2*Pp,gamma_min,'next')+interp1(g_sample,2*Ap-AS,gamma_min,'next'));


[b_sample,b_red]=unique(b_sample);
g_sample=g_sample(b_red);
Pp=Pp(b_red);
Ap=Ap(b_red);

% avoid outputting NaN for largest value of gamma due to numerical error
% and 'next' Interpolant not supporting extrapolation
if b_sample(end)<1
    b_sample(end+1)=1;
    Pp(end+1)=Pp(end);
    Ap(end+1)=Ap(end);
    g_sample(end+1)=g_sample(end);
end

% interpolators to handle events
Pminus=griddedInterpolant(b_sample,PS-Pp,'next');
Pplus=griddedInterpolant(b_sample,Pp,'next');
Aminus=griddedInterpolant(b_sample,AS-Ap,'next');
Aplus=griddedInterpolant(b_sample,Ap,'next');

b=linspace(b_min,1,n)';
gammas=(Aminus(b) + b.*(Aplus(b)-Aminus(b)))./...
    ((1-b).*Pminus(b)+b.*Pplus(b));


end

function A=div_0(A,B)
% DIV_0 pointwise division such that 0/0=0
ind=find(A);
A(ind)=A(ind)./B(ind);
end
