function Dir1Shift = Dir1ShiftCalc(nu,T,N,M)
kPrimeVec = transpose(0:N-1);                                                                                      % Vector of k'
Dir1 = sin(pi*(nu*N*T+kPrimeVec))./sin(pi/N*(nu*N*T+kPrimeVec))...
    .*exp(1i*pi*(nu*N*T+kPrimeVec)*(N-1)/N);                                                                       % First Dirichlet kernel function
Dir1(isnan(Dir1)) = N;
Dir1Shift = zeros(N*M,N);
for jj = 0:N-1
    Dir1Shift(:,jj+1) = repelem(circshift(Dir1,jj),M);
end
end