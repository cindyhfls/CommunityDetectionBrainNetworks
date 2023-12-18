function D = similarity_measures_HSB(clu, str)
% Jiaxin Cindy Tu 2023.12.14
% Input: 
% clu - matrix (Nnode x Nassignments)
% str - string one of the following methods
% Output:
% D - not an actual distance measure but a square matrix with diagonal of
% zeros
% VI,NMI modified from BCT toolbox partition_distance.m
% Rand Coefficient and Variants from NetworkCommunityToolbox zrand.m
% AMI adapted from https://www.mathworks.com/matlabcentral/fileexchange/33144-the-adjusted-mutual-information
% AMI is useful when the number of clusters is different
% AMI is the AMI_max https://jmlr.csail.mit.edu/papers/volume11/vinh10a/vinh10a.pdf

str = lower(str); % make it case-insensitive by converting to non-capital letters

assert(sum(str==string({'nmi','vi','ami','zrand','arand','rand'})),"unknown similarity measure, please edit similarity_measures_HSB.m to add your measure");

assert(~any(isnan(clu(:)))) % make sure no NaN exist

if any(clu(:)==0)
    warning('The assignments have 0, treating them as valid clusters for now, but you should consider using the cluster without 0')
    clu = clu+1;
end

[n, p] = size(clu);
q = p;
if any(str==string({'nmi','vi','ami'}))
    %% Information theory based
    % MIn = 2MI(X, Y) / [H(X) + H(Y)] Normalized mutual information ([p, q] matrix)
    % VIn = [H(X) + H(Y) - 2MI(X, Y)]/log(n) Normalized variation of information ([p, q] matrix)
    % AMI = (MI-EMI)/(max(H(X),H(Y))-EMI) AMI max
    
    HX = zeros(p, 1);  
    
    d = 10.^ceil(log10(double(1 + max(clu(:)) )));
    for i = 1:p
        Px = nonzeros(accumarray(clu(:, i), 1)) / n;                     % P(x)
        HX(i) = - sum(Px .* log(Px));                                   % H(x)
    end
    HY = HX;
    [VIn,MIn,AMI] = deal(zeros(p, q));
    for i = 1:p
        j_idx = ((i - 1) + 1):q;
        for j = j_idx
            cont = sparse(clu(:, i), clu(:, j), 1); % contingency table
            Pxy = nonzeros(accumarray(d*clu(:, i) +clu(:, j), 1)) / n; 	% P(x,y)
            Hxy = -sum(Pxy .* log(Pxy));                                % H(x,y)
            VIn(i, j) = (2 * Hxy - HX(i) - HY(j)) / log(n);             % VIn
            MI = HX(i) + HY(j) - Hxy; %mutual information not normalized
            MIn(i, j) = 2 * (HX(i) + HY(j) - Hxy) / (HX(i) + HY(j));    % MIn
            EMI = expected_mutual_information(cont);
            AMI(i,j) =(MI-EMI)/(max(HX(i),HY(j))-EMI);
        end
        VIn(j_idx, i) = VIn(i, j_idx);
        MIn(j_idx, i) = MIn(i, j_idx);
        AMI(j_idx, i) = AMI(i,j_idx);
    end
    
    switch str
        case 'nmi'
            D= MIn;
        case 'vi'
            D = VIn;
        case 'ami'
            D = AMI;
    end  
    
end

%% Rand coefficient and variants
if any(str==string({'zrand','arand','rand'}))
    [SR,SAR,zRand]= deal(zeros(p, q));
    for i = 1:p
        j_idx = ((i - 1) + 1):q;
        for j = j_idx
            %Generate contingency table and calculate row/column marginals
            nij = sparse(clu(:, i), clu(:, j), 1);
            ni = sum(nij, 2);
            nj = sum(nij, 1);
            
            %Identify total number of elements, n, numbers of pairs, M, and numbers of
            %classified-same pairs in each partition, M1 and M2.
            n = length(clu(:, i));
            M = n*(n-1)/2;
            M1 = sum(ni.^2-ni)/2;
            M2 = sum(nj.^2-nj)/2;
            
            %Pair counting types:
            a = full(sum(sum(nij.^2-nij)))/2; %same in both
            b = M1-a;                         %same in 1, diff in 2
            c = M2-a;                         %same in 2, diff in 1
            d = M-(a+b+c);                    %diff in both
            
            %Rand and Adjusted Rand indices:
            SR(i, j) = (a+d)/(a+b+c+d);
            meana = M1*M2/M;
            SAR(i, j) = (a-meana)/((M1+M2)/2-meana);
            % PS: The adjusted coefficient is calculated by subtracting the expected
            % value and rescale the result by the difference between the maximum
            % allowed value and the mean value
            
            C1 = 4*sum(ni.^3)-8*(n+1)*M1+n*(n^2-3*n-2);
            C2 = 4*sum(nj.^3)-8*(n+1)*M2+n*(n^2-3*n-2);
            vara = M/16 - (4*M1-2*M)^2*(4*M2-2*M)^2/(256*M^2) + C1*C2/(16*n*(n-1)*(n-2)) + ...
                ((4*M1-2*M)^2-4*C1-4*M)*((4*M2-2*M)^2-4*C2-4*M)/(64*n*(n-1)*(n-2)*(n-3));
            zRand(i,j) = (a-meana)/sqrt(vara);
        end
        SR(j_idx, i) = SR(i, j_idx);
        SAR(j_idx, i) = SAR(i, j_idx);
        zRand(j_idx, i) = zRand(i, j_idx);
    end
    
    switch str
        case 'zrand'
            D = zRand;
        case 'rand'
            D = SR;
        case 'arand'
            D = SAR;
    end
    
end

%% Set diagonals zero so we can get the upper triangle directly with squareform()
for ii = 1:length(D)
    D(ii,ii) = 0;
end

end

function EMI = expected_mutual_information(Cont)
    % Cont is the pre-computed contingency table

    % Calculate the size of the data set
    N = sum(Cont(:));

    % Update the true dimensions
    [R, C] = size(Cont);

    % Calculate the sums of rows and columns
    a = sum(Cont, 2); % Row sums
    b = sum(Cont, 1); % Column sums

    % Initialize EMI
    EMI = 0;

    % Calculate the expected mutual information
    AB = a* b;
    E3 = (AB / N^2) .* log(AB / N^2);
    EPLNP = zeros(R, C);
    LogNij = log([1:min(max(a), max(b))] / N);
    for i = 1:R
        for j = 1:C
            nij = max(1, a(i) + b(j) - N);
            X = sort([nij, N - a(i) - b(j) + nij]);
            if N - b(j) > X(2)
                nom = [a(i)-nij+1:a(i), b(j)-nij+1:b(j), X(2)+1:N-b(j)];
                dem = [N-a(i)+1:N, 1:X(1)];
            else
                nom = [a(i)-nij+1:a(i), b(j)-nij+1:b(j)];
                dem = [N-a(i)+1:N, N-b(j)+1:X(2), 1:X(1)];
            end
            p0 = prod(nom./dem) / N;
            EPLNP(i,j) = nij * LogNij(nij) * p0;
            p1 = p0 * (a(i) - nij) * (b(j) - nij) / (nij + 1) / (N - a(i) - b(j) + nij + 1);
            for nij = max(1, a(i) + b(j) - N) + 1:1:min(a(i), b(j))
                EPLNP(i,j) = EPLNP(i,j) + nij * LogNij(nij) * p1;
                p1 = p1 * (a(i) - nij) * (b(j) - nij) / (nij + 1) / (N - a(i) - b(j) + nij + 1);
            end
        end
    end
    EMI = nansum(nansum(EPLNP - E3));
end

