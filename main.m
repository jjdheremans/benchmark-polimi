% main.

clear all 
close all
clc

%% Graph Ouput:
plotPSD = false ;   % plot of PSDs
plotAdm = false ;   % plot of admittance function
plotfde = false ;   % plot of fractional derivatives
plotFEG = false ;   % plot of theodorsen functions F and G
plotfrq = true  ;   % plot of eigen frequencies evolution with resp. to U

%% Input parameters
% Geometry/Material
ml = 22740;     % [kg/m]
Jl = 2.47e+06;  % [kgm2/m]

B = 31;         % [m]
xi = 0.003;     % [-]
fz = 0.1;       % [Hz]
ft = 0.278;     % [Hz]

rho = 1.22 ;                    % [kg/m3]
U_mat = linspace(15,100,50) ;% [m/s]
Iw = 0.05 ;                     % [-]
Lw = 20 ;                       % [m]
typspec = 1 ;                   % VK Spectrum

% Time/Frequency discretization
tmax = 600 ;    % [s]
nvalues = 2^10 ; 



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
F =  @(f)  (J1(f) .* (J1(f) + Y0(f)) + Y1(f) .* (Y1(f) - J0(f))) ./ ((J1(f) + Y0(f)).^2 + (Y1(f) - J0(f)).^2);
G =  @(f) -(J1(f) .* J0(f) + Y1(f) .* Y0(f)) ./ ((J1(f) + Y0(f)).^2 + (Y1(f) - J0(f)).^2);

% Mass/Stiffness/Damping system matrices
M = diag([ml, Jl]) ;
C = 4*pi*xi*M.*diag([fz, ft]) ;
K = 4*pi^2*M.*diag([fz, ft].^2) ;

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

% Equivalent Stiffness/Damping/Mass Matrix
Ceq_fun = @(f,U) ( C - Fse1_fun(f,U) ) ;
Keq_fun = @(f,U) ( K - Fse2_fun(f,U) ) ;
Meq = M ;
Meqm1 = Meq\eye(2) ;

% FRF Matrix
HM1 = @(f,U) (-M * (2 * pi * f).^2 + 1i * 2 * pi * f * Ceq_fun(f,U) + Keq_fun(f,U)) ;

%% LOOP ON AVG Wind speeds
tol = 1e-4 ;


if plotPSD
    fig1 = figure ;  hold on ;
end

freqSto = zeros( 2, length(U_mat ) ) ;
xiSto   = zeros( 2, length(U_mat) ) ;

for looper = 1 : length(U_mat )
    
    
    uliege = U_mat(looper) ;
    
    % Evolution of eigen frequencies with respect to avg wind speeds U
    counter = 0 ;
    for mode = 1 : 2
        counter = 0 ;
        if mode == 1 && looper ==1
            ev_im1(mode) = 2*pi*fz ;
        elseif mode == 2 && looper ==1
            ev_im1(mode) = 2*pi*ft ; 
        end
        freq = ev_im1(mode) / ( 2*pi );
        fs = freq*B/uliege ;
        Keq = Keq_fun( fs, uliege ) ;
        Ceq = Ceq_fun( fs, uliege ) ;
        cdt = 0 ; counter = 0 ;
        while cdt ~= 1 || counter > 10
            [em,evc,~] = polyeig(Keq,1i*Ceq,-Meq) ;
            [~,index] = min ( abs ( evc - ev_im1(mode) ) ) ; 
            ev = evc(index) ;
            freq =  (ev) / ( 2 * pi ) ;
            fs = freq*B/uliege ;
            Keq = Keq_fun( fs, uliege ) ;
            Ceq = Ceq_fun( fs, uliege ) ;
            abstol = tol * mean( abs( Meq(mode,:)*ev^2 + Ceq(mode,:)*ev + Keq(mode,:) ) ) ; 
            cdt = sum( abs( ( -Meq*ev^2+1i*Ceq*ev+Keq ) * em(:,index) ) < abstol ) / 2 ;
            ev_im1(mode) = ev ;
            counter = counter + 1 ;
        end
        if counter > 10 
            warning('Iterative procedure didn''t converge')
        end
        freqSto(mode,looper) = abs(freq) ;
        xiSto(mode,looper) = imag( ev ) / abs( ev ) ;
    end
   
    % PSD
    fs_vec = linspace(0.001, 3, 100); % [Hz]    
    for i = 1 : length( fs_vec ) 
        fs = fs_vec(i) ;
        Keq = Keq_fun( fs, uliege ) ;
        Ceq = Ceq_fun( fs, uliege ) ;
        [phi,e] = polyeig(Keq,1i*Ceq,-Meq) ; 
        [~,ind] = min( abs( e - 2*pi*freqSto(:,looper).' ) ) ;
        phi = phi( : , ind(1,:) ) ;
        H(:,:) = HM1(fs, uliege)\eye(2);
        Sx.mod(:, :, i) = H * Sq(fs, uliege) * H' ;
        Sx.nod(:, :, i) =  phi * Sx.mod(:, :, i) * phi.' ;
        HH(:,:,i) = H ;
        Sx.mod = real( Sx.mod ) ;
        Sx.nod = real( Sx.nod ) ;
    end
    stdvalues.nod1(looper,:) = sqrt( trapz( fs_vec, squeeze(Sx.mod(1,1,:)) ).' ) ;
    stdvalues.nod2(looper,:) = sqrt( trapz( fs_vec, squeeze(Sx.mod(2,2,:)) ).' ) ;
    
    %% Graphical Representation
    if plotPSD
        % Plot of Nodal PSDs
        figure(fig1) ;
        strings = cell(size(Sx.mod,1)*size(Sx.mod,2),1) ;
        for i = 1 : size( Sx.mod, 1) 
            for j = 1 : size( Sx.mod, 2) 
                semilogy( fs, squeeze(Sx.mod(i,j,:)),'LineWidth',1.5 )
                hold on 
                index = sub2ind([size(Sx.mod,1), size(Sx.mod,2)],i,j) ;
                strings{index} = sprintf('Sx: i=%d, j=%d',i,j) ;
            end
        end
    end

end

figure
plot (fs_vec, abs(squeeze(HH(1,1,:))) ) 
%% Graphical post-processing

% Annotation of fig1
if plotPSD
    figure(fig1) ; grid on
    % legend(strings, 'interpreter','latex')
    xlabel('Frequences [Hz]', 'FontSize',12, 'interpreter','latex')
    ylabel('PSD [m$^2$/s$^3$]', 'FontSize',12, 'interpreter','latex')
end
    

% Annotation of fig2
if plotfrq
    fig2 = figure ; grid on ; hold on
    plot( U_mat, freqSto(1,:), '*-' )
    plot( U_mat, freqSto(2,:), '*-' )  
    legend({'Mode 1','Mode 2'}, 'interpreter','latex','FontSize',12)
    xlabel('Avg. Wind Speed $U$ [m/s]', 'FontSize',12, 'interpreter','latex')
    ylabel('Eigen Frequencies $f_i$ [Hz]', 'FontSize',12, 'interpreter','latex')
    ylim([0 0.30])
end

% Plot of Theodorsen functions F and G
if plotFEG
    figure 
    fplot(F); hold on
    fplot(G);
    xlabel('$f^*$ [Hz]','interpreter','latex','FontSize',12)
    ylabel('$F(f^*)$ and $G(f^*)$','interpreter','latex','FontSize',12)
    legend({'$F(f^*)$','$G(f^*)$'},'interpreter','latex','FontSize',12)
    grid on
    for kk = 1:length(fs)
        for k = 1:length(a)
            aa(k,kk) = a{k}(fs(kk));
            hh(k,kk) = h{k}(fs(kk));
        end
    end
end


% Plot of flutter derivatives h_i and a_i
if plotfde
    figure 
    subplot(2,1,1)
    for jj = 1:size(hh,1)
        plot(Vs(fs),hh(jj,:),'DisplayName',['h_',num2str(jj)]); hold on
    end
    title('Flutter derivatives $h_i$', 'interpreter','latex','FontSize',12)
    xlabel('V* [Hz$^{-1}$]', 'interpreter','latex','FontSize',12)
    xlim([0 50])
    legend show
    grid on
    subplot(2,1,2)
    for ii = 1:size(aa,1)
        plot(Vs(fs),aa(ii,:),'DisplayName',['a_',num2str(ii)]); hold on
    end
    title('Flutter derivatives $a_i$', 'interpreter','latex','FontSize',12)
    xlabel('V* [Hz$^{-1}$]', 'interpreter','latex','FontSize',12)
    xlim([0 50])
    legend show
    grid on 
end


% Plot of admittance function
if plotAdm
    figure 
    hold on
    fplot(A,[0 10]); 
    ylabel('Admittance A(f*)', 'interpreter','latex','FontSize',12)
    xlabel('$f^*$ [Hz]', 'interpreter','latex','FontSize',12)
    grid on
end







    
