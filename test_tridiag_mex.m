% test tridiag mex

if ~(exist('col','file') == 2)
    run('~/Documents/mai_code/mai_setup.m');
end
addpath('~/Documents/mai_code/ADMM_tridiag/');
addpath('~/Documents/mai_code/pthread_tutor/');
addpath('~/Documents/mai_code/util/');

% d = [1 2 1]';
% 
% a = [1 2]';
% b = [3 4 5]';
% c = [6 7]';

% d = [1 2]';
% 
% a = [1]';
% b = [3 4]';
% c = [6]';

rng(0);
N = 10;
M = 5;
scale = 10;
d = scale*randn(N,M);
%d = d + 1i*scale*rand(N,M);
a = scale*rand(N-1,1)-scale/2;
b = scale*rand(N,1)-scale/2;
c = scale*rand(N-1,1)-scale/2;


% d = col(1:N);
% %d = d + 1i*col(1:N);
% a = col(1:N-1) - 2;
% c = flipud(a);
% b = col(N:-1:1);

d = single(d);
a = single(a);
b = single(b);
c = single(c);


T = diag(a,-1) + diag(b) + diag(c,1);
x0 = T\d;

tic
x1 = apply_tridiag_inv(a, b, c, d);
toc 
x1_real = apply_tridiag_inv(a, b, c, real(d));

% mex tridiag_inv_mex_nopar.c
% try
%     tic
%     x2 = tridiag_inv_mex_nopar(a, b, c, d);
%     toc
% catch
%     display('failed, prob seg fault');
% end
% norm(x1-x2)

%%x
% mex -O CFLAGS="\$CFLAGS -std=c99" -I./def/ tridiag_inv_mex.c
% 
% try
%     tic
%     x3 = tridiag_inv_mex(single(a), single(b), single(c), single(d));
%     toc
% catch
%     display('failed, prob seg fault');
% end

mex -O CFLAGS="\$CFLAGS -std=c99" -I./def/ tridiag_inv_mex_varnthread.c

ncores = int16(2);
try
    tic
    x3 = tridiag_inv_mex_varnthreads(single(a), single(b), single(c), ...
        single(d), ncores);
    toc
catch
    display('failed, prob seg fault');
end

%norm(x1-x3)
norm(col(x0)-col(x3))
% norm(x1_real-x3)
return
%% plot times as func of N

Ns = 100:100:10000;
Ns = 10:10:1000;

for ii = 1:length(Ns)
	N = Ns(ii);
	M = N;
	scale = 10;
	d = 3 + scale*randn(N,M);
	d = d + 1i*scale*rand(N,M);
	a = scale*rand(N-1,1)-scale/2;
	b = scale*rand(N,1)-scale/2;
	c = scale*rand(N-1,1)-scale/2;
	tic
	x1 = ir_apply_tridiag_inv(a, b, c, d);
	ml_toc(ii) = toc;     
	tic
    	x3 = tridiag_inv_mex(a, b, c, d);
   	pthr_toc(ii) = toc;

end

figure; plot(Ns, ml_toc); hold on; plot(Ns, pthr_toc,'r')
xlabel('size of tridiag matrix') 
ylabel('computation time')
ylabel('computation time (sec)')
legend('matlab', 'pthreaded MEX')
title('Tridiagonal Solver Compute Time')


