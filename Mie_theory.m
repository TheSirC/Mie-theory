% SCATTERING BY A SPHERICAL NANOPARTICLE USING MIE THEORY
clear all;

lambda = 1e-9; %light wavelength in m 
R=50e-9; %particle radius in m
A=pi*R^2;  %particle geometrical cross-section

n_m=1.5; %refractive index of the medium

% Load experimental data of the particle refractive index
% 2 files are available: for gold and silver
data = dlmread('./données/ag_palik.txt');
wlgths = data(:,1);
n = data(:,2);
k = data(:,3);

lambda1 = 30e-9; % starting wavelength in m
lambda2 = 800e-9; % end wavelength in m
step = 1e-9; 
nbrLambdas = int16((lambda2-lambda1)/step+1);

Qabs = zeros(nbrLambdas,1); %absorption efficiency
Qsca = zeros(nbrLambdas,1); %scattering efficiency
Qext = zeros(nbrLambdas,1); %extinction efficiency
lambdas = zeros(nbrLambdas,1);

count = 1; 
for lambda=lambda1:step:lambda2
    % Interpolate the experimental refractive index of the particle...
    npart = interp1q(wlgths,n,lambda) + 1i*interp1q(wlgths,k,lambda);
    
    [e s a] = MieScattering(lambda,R,n_m,npart);
    Qsca(count,1) = s/A;
    Qext(count,1) = e/A;
    Qabs(count,1) = a/A;
    lambdas(count,1) = lambda*1e9;
    count = count + 1  
end

figure(1); 
plot(lambdas,Qabs,lambdas,Qsca,lambdas,Qext);
xlabel('Wavelength (nm)');
ylabel('Efficiency');
legend('Qabs','Qsca','Qext');