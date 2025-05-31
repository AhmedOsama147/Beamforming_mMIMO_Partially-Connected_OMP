function x = randLaplacian(m, n, mu, var)
    % Generate Laplacian noise
    u = rand(m, n)-0.5;
    b = sqrt(var/2);
    x = mu - b*sign(u).*log(1- 2*(1-exp(-pi/b))*abs(u));
end