% test tridiag_inv_mex

n = 10;
n = 5;
%sub = ones(n-1,1);
sub = [11 2 3 4]';
sup = sub;
%diags = 3*ones(n,1);
diags = [1.1 2.1 3.1 4.1 5.1]';
%rhs = 2*ones(n,1); % make complex later
rhs = [1.2 2.2 3.2 4.2 5.2]';

% just choose one block, so block_size = n
output = tridiag_inv_mex(double(sub),double(diags),double(sup),double(rhs),int32(n));
% how to make mex file accept horizontal vectors?


