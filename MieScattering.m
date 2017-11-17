function [Cext,Csca,Cabs]=MieScattering(lambda,R,n_m,npart)

m=npart/n_m;
k=2*pi*n_m/lambda;
v=k*R;
w=m*v;
N=round(2+v+4*v^(1/3));

% computation
j=(1:N);

%Implement here a procedure to calculate extinction, scattering and
%absorption cross-sections

%Cext=...;
%Csca=...;
%Cabs=...;

end
