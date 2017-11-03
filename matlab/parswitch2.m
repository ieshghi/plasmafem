function D = parswitch2(epsilon,delta,kappa,A,gamma)
    M = [1 (1+epsilon)^2 (1+epsilon)^4 0;
    1 (1-epsilon)^2 (1-epsilon)^4 0;...
    1 (1-delta*epsilon)^2 ...
    (1-delta*epsilon)^4-4*(1-delta*epsilon)^2*kappa^2*epsilon^2 ...
    kappa*epsilon;
    1 (1-delta*epsilon)^2 ...
    (1-delta*epsilon)^4-4*(1-delta*epsilon)^2*gamma^2*kappa^2*epsilon^2 ...
    -gamma*kappa*epsilon];

%Matrix entries for the particular solutions

    B=-[(1+epsilon)^4/8+A*(1/2*(1+epsilon)^2*log(1+epsilon)-(1+epsilon)^4/8);
    (1-epsilon)^4/8+A*(1/2*(1-epsilon)^2*log(1-epsilon)-(1-epsilon)^4/8);
    (1-delta*epsilon)^4/8+A*(1/2*(1-delta*epsilon)^2*log(1-delta*epsilon)...
    -(1-delta*epsilon)^4/8);
    (1-delta*epsilon)^4/8+A*(1/2*(1-delta*epsilon)^2*log(1-delta*epsilon)...
    -(1-delta*epsilon)^4/8)];

% Compute coefficients d_1, d_2, d_3, and d_4

    D=M\B;

end
