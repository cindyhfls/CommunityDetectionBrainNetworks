
label_true = repelem(1:5, [30 15 25 20 10]);

mu_mat = [0.4 0.2 0.1 -0.4 -0.4;
        0.2 0.4 0.2 -0.3 -0.4;
        0.1 0.2 0.45 -0.01 -0.2;
        -0.4 -0.3 -0.01 0.5 0.2;
        -0.4 -0.4 -0.2 0.2 0.3];

sigma_mat = repmat(0.1, 5,5);

nROI = 100;
n = 1;
A = NaN(nROI,nROI,n);

rng('default')
for kk= 1:n
    for ii=1:nROI
        for jj=ii:nROI 
           if ii == jj
               A(ii,jj,kk) =1;
               A(jj,ii,kk)=1;
           else
               A(ii,jj,kk) = normrnd( mu_mat(label_true(ii),label_true(jj)) , sigma_mat(label_true(ii),label_true(jj)));
               A(jj,ii,kk) = A(ii, jj, kk);
           end
        end
    end
end

ave_fcmat = mean(A,3);

figure('position',[100 100 800 400]);
subplot(1,2,1);
mat_plot(mu_mat);
subplot(1,2,2);
mat_plot(sigma_mat);

figure;
mat_plot(ave_fcmat);

%%%%%%%%%%%%%%%%%
noise = 0.3*rand(nROI)-0.15;
noise(eye(nROI)==1) = 0;

ave_fcmat_wnoise = ave_fcmat + noise;
figure;
mat_plot(ave_fcmat_wnoise);
title('With noise');
%%%%%%%%%%%%%%%%%


