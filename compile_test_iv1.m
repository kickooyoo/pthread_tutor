% test simple mex on iv1

% gcc -v in term to see I have 4.8
% mex -setup in ml to see I need gcc 4.7

% sudo apt-get install g++.4.7
% sudo apt-get install gcc.4.7

% changes to /usr/local/MATLAB/R2015a/bin/mexopts.sh on gcc and g++
% http://www.walkingrandomly.com/?p=1959

% for C99, had to remove -ansi flag from mexopts.sh on line 59 and 75
% apparently ansi means -std=c89, hence the conflict
mex -O CFLAGS="\$CFLAGS -std=c99" hello_mex.c

%mex -O CFLAGS="\$CFLAGS -Wp,-lang-c-c++-comments" tridiag_inv_mex.c

mex timestwo.c
