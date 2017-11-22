function [Cext,Csca,Cabs]=MieScattering(lambda,R,n_m,npart)

m=npart/n_m;
k=2*pi*n_m/lambda;
v=k*R;
w=m*v;
N=round(2+v+4*v^(1/3));

j=(1:N); % finite upper bond for the sum

% Ricatti-Bessel functions and their derivatives
fpsi = @(j,x) sqrt(pi*x/2)*besselj(j+1/2,x);
deriv_psi = @(j,x) fpsi(j-1,x) - (j/x) * fpsi(j,x);

kzi = @(j,x) sqrt(pi*x/2)*(besselj(j+1/2,x)+1i*bessely(j+1/2,x));
deriv_kzi = @(j,x) kzi(j-1,x) - (j/x) * kzi(j,x);

% Definition of the coefficients
global a b;
a = @(j) (m*fpsi(j,w)*deriv_psi(j,v)-fpsi(j,v)*deriv_psi(j,w))/...
    (m*fpsi(j,w)*deriv_kzi(j,v)-fpsi(j,v)*deriv_kzi(j,w));
b = @(j) (fpsi(j,w)*deriv_psi(j,v)-m*fpsi(j,v)*deriv_psi(j,w))/...
    (fpsi(j,w)*deriv_kzi(j,v)-m*deriv_psi(j,w)*kzi(j,v));

% Definition of the coeffients of extinction, scattering and absorption respectively
Cext=(2*pi/k^2)*sum_sig_ext(j);
Csca=(2*pi/k^2)*sum_sig_sca(j);
Cabs=Cext-Csca;

end
