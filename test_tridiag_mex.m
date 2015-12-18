% test tridiag mex
run('~/Documents/mai_code/mai_setup.m');
addpath('~/Documents/mai_code/ADMM_tridiag/');
addpath('~/Documents/mai_code/pthread_tutor/');

d = [1 2 1]';

a = [1 2]';
b = [3 4 5]';
c = [6 7]';

% d = [1 2]';
% 
% a = [1]';
% b = [3 4]';
% c = [6]';
rng(0);
N = 200;
M = 200;
d = 10*rand(N,M);
a = 10*rand(N-1,1);
b = 10*rand(N,1);
c = 10*rand(N-1,1);

tic
x1 = apply_tridiag_inv(a, b, c, d);
toc 
mex tridiag_inv_mex_nopar.c
try
    tic
    [x2, new_c, new_d] = tridiag_inv_mex_nopar(a, b, c, d);
    toc
catch
    display('failed, prob seg fault');
end
norm(x1-x2)