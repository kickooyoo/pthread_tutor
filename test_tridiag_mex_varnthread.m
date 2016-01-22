%% test tridiag mex
%% check value against \ and ir_apply_tridiag

rng(0);
N = 5; % 256
M = 4; % 256
% hyperthreading means up to 160 on mpel8
scale = 10;
d = scale*randn(N,M);
d = d + 1i*scale*rand(N,M);
a = scale*rand(N-1,1)-scale/2;
b = scale*rand(N,1)-scale/2;
c = scale*rand(N-1,1)-scale/2;

d = single(d);
a = single(a);
b = single(b);
c = single(c);

T = diag(a,-1) + diag(b) + diag(c,1);
x0 = T\d;

x1 = ir_apply_tridiag_inv(a, b, c, d);

% compile as needed
mex -O CFLAGS="\$CFLAGS -std=c99" -I./def/ tridiag_inv_mex_varnthread2.c

ncores = int16(jf('ncore'));

try
    x2 = tridiag_inv_mex_varnthread2(a, b, c, d, ncores);
catch
    display('tridiag_inv_mex_varnthread.c failed');
end

% [col(x0) col(x2)]

equivs(x0, x2)
equivs(x1, x2)