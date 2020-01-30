function [H1, H2, H3, H4, A1, A2, A3, A4] = Scanlan(fs)

Vs= 1./fs;
k = pi*fs; % et donc K = 2*k;

if k==0; k=1e-6; end

J0 = besselj(0,k);
J1 = besselj(1,k);
Y0 = bessely(0,k);
Y1 = bessely(1,k);

denom = (J1+Y0).^2 + (Y1-J0).^2;
k_=k;
while denom==0
    k_=k_*0.95;
    J0 = besselj(0,k_);
    J1 = besselj(1,k_);
    Y0 = bessely(0,k_);
    Y1 = bessely(1,k_);
    denom = (J1+Y0).^2 + (Y1-J0).^2;
end

F = (J1.*J1+J1.*Y0 + Y1.*Y1-Y1.*J0) ./ ((J1+Y0).^2 + (Y1-J0).^2);
G = - (J1.*J0 + Y1.*Y0) ./ ((J1+Y0).^2 + (Y1-J0).^2);

H1 = -Vs .* F;
H2 =  Vs/4 .* (1+F+2/pi*Vs.*G);
H3 = Vs.^2 /2/pi .* (F-pi/2./Vs.*G);
H4 = pi/2*(1+2/pi*Vs.*G);

A1 = -Vs/4.*F;
A2 = -Vs/16 .* (1-F-2/pi*Vs.*G);
A3 = Vs.^2 /8/pi .* (F-pi/2./Vs.*G);
A4 = Vs/4.*G;