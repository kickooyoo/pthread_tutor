%% test tridiag mex
%% check value against \ and ir_apply_tridiag

% compile as needed
mex -O CFLAGS="\$CFLAGS -std=c99 -DMmex" -I./def/ tridiag_inv_mex_noni.c
    
if 1
N = 50; % 256
M = 50; % 256
scale = 10;
sparse_i = cat(1, col(2:N*M), col(1:N*M), col(1:N*M-1));
sparse_j = cat(1, col(1:N*M-1), col(1:N*M), col(2:N*M));

ncores = int32(4);%int16(jf('ncore'));
nrep = 1;

for jj = 1:nrep
    rng(jj);
    % hyperthreading means up to 160 on mpel
    d = scale*randn(N,M);
    d = d + 1i*scale*rand(N,M);
    a = scale*rand(N-1,M) - scale;%-scale/2;
    b = scale*rand(N,M) + scale*N;
    c = scale*rand(N-1,M) - scale;%-scale/2;
    
    d = single(d);
    a = single(a);
    b = single(b);
    c = single(c);
    
    a_pad = col(cat(1, a, zeros(1, M)));
    a_long = a_pad(1:end-1);
    c_pad = col(cat(1, c, zeros(1, M)));
    c_long = c_pad(1:end-1);    
    T = sparse(sparse_i, sparse_j, double(cat(1, a_long, b(:), c_long)));
    if (condest(T) > 1e15)
            keyboard;
    end
%     Tf = diag(a_long, -1) + diag(b(:)) + diag(c_long, 1);
    x0 = T\double(d(:));
%     x0b = Tf\d(:);

    x1 = ir_apply_tridiag_inv(a_long, b(:), c_long, d(:));
    
    try
        x2 = tridiag_inv_mex_noni(a, b, c, d, ncores);
    catch
        display('tridiag_inv_mex_varnthread.c failed');
    end
    
    norm(x0-x2)/(N*M)
    norm(x1 - x2)/(N*M)
    err(jj) = norm(x0-x2)/(N*M);
end
% if any(err > 1e-3)
%     display('bad err');
%     keyboard;
% end

% [col(x0) col(x2)]

% norm(x1-x2)
% equivs(x0, x2)
% equivs(x1, x2)
% return;

end
%% timing test
nrep = 16;
warmup = 4;
%ncores = int16(jf('ncore'));
ncores = int32([4]);% 8 16]);
for ii = 1:ncores
    for jj = 1:nrep
        if (jj > warmup) tic; end
        T = sparse(sparse_i, sparse_j, double(cat(1, a_long, b(:), c_long)));
        x0 = T\double(d(:));
        if (jj > warmup)
            bs_toc(ii,jj-warmup) = toc;
        end
    end
    if 1
        for jj = 1:nrep
            if (jj > warmup) tic; end
            x1 = ir_apply_tridiag_inv(a_long, b(:), c_long, d(:));
            if (jj > warmup)
                ir_toc(ii,jj-warmup) = toc;
            end
        end
    end
    for jj = 1:nrep
        if (jj > warmup) tic; end
        x2 = tridiag_inv_mex_noni(a, b, c, d, ii);
        if (jj > warmup)
            pth_toc(ii,jj-warmup) = toc;
        end
    end
end

figure; plot(1:ncores, mean(bs_toc,2));
hold on; plot(1:ncores, mean(ir_toc,2),'r');
hold on; plot(1:ncores, mean(pth_toc,2),'g');
legend('backslash', 'ir apply', 'pthread');
 return;
%% bad inputs

% mixed row and col vectors for a, b, c OK
x3 = tridiag_inv_mex_noni(a, b', c, d, ncores);

% bad sizes
x3 = tridiag_inv_mex_noni(a(3:end), b, c, d, ncores);
x3 = tridiag_inv_mex_noni(a, b, [c; 1; -1], d, ncores);
x3 = tridiag_inv_mex_noni(a, b, c, scale*rand(N+1, M), ncores);

% bad types
x3 = tridiag_inv_mex_noni(double(a), b, c, d, ncores);
x3 = tridiag_inv_mex_noni(a, b, c, double(d), ncores);
x3 = tridiag_inv_mex_noni(a, int16(b), c, d, ncores);
