function [Cext, Csca, Cabs] = MieScattering(lambda, R, n_m, npart)

m = npart / n_m;
k = 2 * pi * n_m * R / lambda;
v = k * R;
w = m * v;
N = round(2+v+4*v^(1 / 3));

j = (1:N); % finite upper bond for the sum

% Ricatti-Bessel functions and their derivatives
psi_j = @(j, x) sqrt(pi*x/2) * besselj(j+0.5, x);
deriv_psi = @(j, x) psi_j(j-1, x) - (j / x) * psi_j(j, x);

kzi = @(j, x) sqrt(pi*x/2) * (besselj(j+0.5, x) + 1i * bessely(j+0.5, x));
deriv_kzi = @(j, x) kzi(j-1, x) - (j / x) * kzi(j, x);

% Definition of the coefficients
global a b;
a = @(j,m,w,v) (m * psi_j(j, w) * deriv_psi(j, v) - psi_j(j, v) * deriv_psi(j, w)) / ...
    (m * psi_j(j, w) * deriv_kzi(j, v) - psi_j(j, v) * deriv_kzi(j, w));
b = @(j,m,w,v) (psi_j(j, w) * deriv_psi(j, v) - m * psi_j(j, v) * deriv_psi(j, w)) / ...
    (psi_j(j, w) * deriv_kzi(j, v) - m * deriv_psi(j, w) * kzi(j, v));

% Definition of the coeffients of extinction, scattering and absorption respectively
Cext = (2 * pi / k^2) * sum_sig_ext(j,m,w,v);
Csca = (2 * pi / k^2) * sum_sig_sca(j,m,w,v);
Cabs = Cext - Csca;

end
