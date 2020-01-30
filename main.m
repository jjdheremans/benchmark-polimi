% main.

clear all 
close all
clc


%% Input parameters
% Geometry/Material
ml = 22740;     % [kg/m]
Jl = 2.47e+06;  % [kgm2/m]

B = 31;         % [m]
xi = 3e-3;      % [-]
fz = 0.1;       % [Hz]
ft = 0.278;     % [Hz]

rho = 1.22 ;    % [kg/m3]
U_mat = [15];%,30,45,60,70] ; % [m/s]
Iw = 0.05 ;     % [-]
Lw = 20 ;       % [m]
typspec = 1 ;   % VK Spectrum

% Time/Frequency discretization
tmax = 600 ;    % [s]
nvalues = 2^13 ; 

%% Definition of wind PSD
sigmaw = @(U) Iw * U ;    % [m/s]
if typspec == 1 
    Sw = @(f,U) sigmaw(U)^2./f.*( 4*f*Lw/U ) .* ( 1 + 755.2 * ( f* Lw / U).^2) ./ ( 1 + 283.2 * ( f* Lw / U).^2 ).^(11/6) ;
else
    fprintf('PSD type not implemented yet\n')
end

% Bessel functions
dt = tmax / nvalues;                    % [s]
fmax = 1 / dt;                          % [Hz]
fny = fmax / 2;                         % [Hz]

Vs = @(f) 1./f ;
J0 = @(f) besselj(0, f*pi);
J1 = @(f) besselj(1, f*pi);
Y0 = @(f) bessely(0, f*pi);
Y1 = @(f) bessely(1, f*pi);
F =  @(f) (J1(f) .* (J1(f) + Y0(f)) + Y1(f) .* (Y1(f) - J0(f))) ./ ((J1(f) + Y0(f)).^2 + (Y1(f) - J0(f)).^2);
G =  @(f) -(J1(f) .* J0(f) + Y1(f) .* Y0(f)) ./ ((J1(f) + Y0(f)).^2 + (Y1(f) - J0(f)).^2);

% Mass/Stiffness/Damping system matrices
M = diag([ml, Jl]) ;
C = 4*pi*xi*M*diag([fz, ft]) ;
K = 4*pi^2*M*diag([fz, ft].^2) ;

% determination of flutter derivatives 
h = cell(4, 1);
h{1} = @(f) 2 * pi * F(f);
h{2} = @(f) -2 * pi * (1/4 + F(f) / 4 + Vs(f) .* G(f) / (2 * pi));
h{3} = @(f) 2 * pi * (F(f) - pi * G(f) ./ (2 * Vs(f)));
h{4} = @(f) 1 + 2 / pi * Vs(f) .* G(f);

a = cell(4, 1);
a{1} = @(f) pi / 2 * F(f);
a{2} = @(f) pi / 2 * (1/4 - F(f) / 4 - Vs(f) .* G(f) / (2 * pi));
a{3} = @(f) pi / 2 * (F(f) - pi * G(f) ./ (2 * Vs(f)));
a{4} = @(f) Vs(f) .* G(f) / (2 * pi);

% Self-excited forces
Fse1_fun = @(f,U) - 1/2 * rho * U^2 * B * [h{1}(f)/U h{2}(f)*B/U ; a{1}(f)*B/U a{2}(f)*B^2/U ] ;
Fse2_fun = @(f,U)   1/2 * rho * U^2 * B * [2*pi^3*h{4}(f)./(B*Vs(f)^2) h{3}(f) ; 2*pi^3*a{4}(f)./(Vs(f)^2)  a{3}(f)*B ] ;

% Buffetting force
A = @(f) 2 ./ (7 * f).^2 .* (7 * f - 1 + exp(-7 * f));
X = @(f,U) ( 0.5 * rho * U * B ) * [2*pi; pi/2*B].*A(f) ;
Sq = @(f,U) X(f,U) * Sw(f,U) * X(f,U)';


%% LOOP ON AVG Wind speeds
fs = linspace(0.001, 3, nvalues); % [Hz]

for looper = 1 : length(U_mat )
    
    uliege = U_mat(looper) ;

    %% Modal Analysis
    [phi, ev] = eig(K, M);
    phi = phi ./ max(phi,[],1) ;
    if ( ~prod(diag(ev) - (2*pi*([fz ft].')).^2 < 1e-03 ) ) % Are the EV well recovered?
        fprintf('Assertion failed: the initial eigen values haven''t been recovered\n')
        return
    end
    Ms = phi' * M * phi;
    Cs = phi' * C * phi;
    Ks = phi' * K * phi;
    HM1 = @(f) (-Ms * (2 * pi * f).^2 + 1i * 2 * pi * f * (Cs - Fse1_fun(f, uliege)) + (Ks - Fse2_fun(f, uliege)));
    Sx = zeros(2, 2, nvalues);
    for i = 1 : length( fs )
        freq = fs(i) ;
        H(:,:,i) = inv( HM1(freq) );
        Sx(:, :, i) = H(:,:,i) * Sq(freq, uliege) * conj(transpose(H(:,:,i))) ;
    end
    Sx = real( Sx ) ;

    %% Graphical Representation
    fig = figure ;
    strings = cell(size(Sx,1)*size(Sx,2),1) ;
    for i = 1 : size( Sx, 1) 
        for j = 1 : size( Sx, 2) 
            semilogy( fs, squeeze(Sx(i,j,:)),'LineWidth',1.5 )
            hold on 
            index = sub2ind([size(Sx,1), size(Sx,2)],i,j) ;
            strings{index} = sprintf('Sx: i=%d, j=%d',i,j) ;
        end
    end
    legend(strings, 'interpreter','latex')
    grid on 
    xlabel('Frequences [Hz]', 'FontSize',12)
    ylabel('PSD [m^2/s^3]', 'FontSize',12)
    ttitle = sprintf('Wind Speed: %d [m/s]', uliege ) ;
    title( ttitle )
    
    figure 
    fplot(F,'DisplayName','F(fs)'); hold on
    fplot(G,'DisplayName','G(fs)')
    legend show
    grid on

    for kk = 1:length(fs)
        for k = 1:length(a)
            aa(k,kk) = a{k}(fs(kk));
            hh(k,kk) = h{k}(fs(kk));
        end
    end

    figure 
    subplot(1,2,1)
    for jj = 1:size(hh,1)
        plot(Vs(fs),hh(jj,:),'DisplayName',['h_',num2str(jj)]); hold on
    end
    title('Flutter derivatives h_i')
    xlabel('V*')
    xlim([0 50])
    legend show
    grid on
    subplot(1,2,2)
    for ii = 1:size(aa,1)
        plot(Vs(fs),aa(ii,:),'DisplayName',['a_',num2str(ii)]); hold on
    end
    title('Flutter derivatives a_i')
    xlabel('V*')
    xlim([0 50])
    legend show
    grid on
    
    figure 
    fplot(A,[0 10]); hold on
    title('Admittance A(f*)')
    grid on

end
