%% test tridiag mex
%% check value against \ and ir_apply_tridiag

N = 50; % 256
M = 30; % 256
scale = 10;
ncores = int16(2);%int16(jf('ncore'));
nrep = 200;
for jj = 1:nrep
    rng(jj);
    % hyperthreading means up to 160 on mpel
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
    % mex -O CFLAGS="\$CFLAGS -std=c99" -I./def/ tridiag_inv_mex_varnthread2.c
    
    try
%         x2 = tridiag_inv_mex_varnthread2(a, b, c, d, ncores);
        x2 = tridiag_inv_mex2(a, b, c, d, ncores);
    catch
        display('tridiag_inv_mex_varnthread2.c failed');
    end
%     err(jj) = norm(x0-x2)/N;
end
% if any(err > 1e-3)
%     display('bad err');
%     keyboard;
% end

% [col(x0) col(x2)]

% norm(x1-x2)
% equivs(x0, x2)
% equivs(x1, x2)
return;
%% timing test
nrep = 4;
ncores = int16(1:jf('ncore'));
ncores = int16(2);
for ii = 1:ncores
        ncore = int16(jf('ncore'));
        for jj = 1:nrep
                tic
                T = diag(a,-1) + diag(b) + diag(c,1);
                x0 = T\d;
                bs_toc(ii,jj) = toc;
        end
        for jj = 1:nrep
                tic
                x1 = ir_apply_tridiag_inv(a, b, c, d);
                ir_toc(ii,jj) = toc;
        end
        for jj = 1:nrep
                tic
                x2 = tridiag_inv_mex2(a, b, c, d, ii);
                pth_toc(ii,jj) = toc;
        end
end

figure; plot(1:ncores, mean(bs_toc,2));
hold on; plot(1:ncores, mean(ir_toc,2),'r');
hold on; plot(1:ncores, mean(pth_toc,2),'g');

return;
%% bad inputs

% mixed row and col vectors for a, b, c OK
x3 = tridiag_inv_mex_varnthread2(a, b', c, d, ncores);

% bad sizes
x3 = tridiag_inv_mex_varnthread2(a(3:end), b, c, d, ncores);
x3 = tridiag_inv_mex_varnthread2(a, b, [c; 1; -1], d, ncores);
x3 = tridiag_inv_mex_varnthread2(a, b, c, scale*rand(N+1, M), ncores);

% bad types
x3 = tridiag_inv_mex_varnthread2(double(a), b, c, d, ncores);
x3 = tridiag_inv_mex_varnthread2(a, b, c, double(d), ncores);
x3 = tridiag_inv_mex_varnthread2(a, int16(b), c, d, ncores);
