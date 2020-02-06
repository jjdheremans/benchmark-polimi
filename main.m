% main.

clear all 
close all
clc

%% Graph Ouput:

NbPlots = 10 ;       % Number of U to display

plotPSD = true  ;   % plot of PSDs
plotFRF = true  ;   % plot of FRFs
plotSTD = true  ;   % plot of Standard Deviation
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

rho = 1.22 ;                            % [kg/m3]
nUvalues = 50 ;                         % Nb of values for avg. speed U
U_mat = linspace(15,100,nUvalues) ;     % [m/s]
Iw = 0.05 ;                             % [-]
Lw = 20 ;                               % [m]
typspec = 1 ;                           % VK Spectrum

% Time/Frequency discretization
tmax = 600 ;    % [s]
nvalues = 2^10 ; 



%% Definition of wind PSD
sigmaw = @(U) Iw * U ;    % [m/s]
if typspec == 1 
    Sw = @(f,U) sigmaw(U)^2./f.*( 4*f*Lw/U ) .* ( 1 + 755.2 * ( f* Lw / U).^2) ./ ( ( 1 + 283.2 * ( f* Lw / U).^2 ).^(11/6) ) ;
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
Sq = @(f,U) X(f*B/U,U) * Sw(f,U) * X(f*B/U,U)';

% Equivalent Stiffness/Damping/Mass Matrix
Ceq_fun = @(f,U) ( C - Fse1_fun(f,U) ) ;
Keq_fun = @(f,U) ( K - Fse2_fun(f,U) ) ;
Meq = M ;
Meqm1 = Meq\eye(2) ;

% FRF Matrix
HM1 = @(f,U) (-M * (2 * pi * f).^2 + 1i * 2 * pi * f * Ceq_fun(f*B/U,U) + Keq_fun(f*B/U,U)) ;

%% LOOP ON AVG Wind speeds
tol = 1e-4 ;

if plotFRF 
    fig3 = figure ;         
    subplot(1,2,1)
    hold on
    grid on 
    subplot(1,2,2)
    hold on 
    grid on 
end

if plotPSD 
    fig1 = figure ;         
end

varvalues = zeros(2,2,length(U_mat) ) ;
freqSto = zeros( 2, length(U_mat ) ) ;
xiSto   = zeros( 2, length(U_mat) ) ;
ev_im1 = zeros(  2, 1 ) ;
strings = cell(NbPlots,1) ;
plotcounter = 1 ; indexPlot = round( linspace(1,nUvalues,NbPlots) ) ;

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
    f_vec = linspace(0.0001, 5, 500); % [Hz]    
    for i = 1 : length( f_vec ) 
        f = f_vec(i) ;
        fs = f * B/uliege ;
        Keq = Keq_fun( fs, uliege ) ;
        Ceq = Ceq_fun( fs, uliege ) ;
        H(:,:) = HM1( f, uliege)\eye(2);
        Sx.mod(:, :, i) = H * Sq(f, uliege) * H' ;
        HH(:,:,i) = H ;
        Sx.mod = real( Sx.mod ) ;
    end
    
    for ii = 1 : 2 
        for jj = 1 : 2
            varvalues(ii,jj,looper) = trapz( f_vec, squeeze(Sx.mod(ii,jj,:)) ) ;
        end 
    end
    
    %% Graphical Representation
    
    if looper == indexPlot( plotcounter ) 

        strings{plotcounter} = sprintf('U=%3.2f [m/s]',uliege) ;

        % Plot of Nodal PSDs
        if plotPSD
            figure(fig1) ;
            subplot(1,2,1)
            semilogy( f_vec, squeeze(Sx.mod(1,1,:)),'LineWidth',1.5 )
            hold on
            grid on 
            subplot(1,2,2)
            semilogy( f_vec, squeeze(Sx.mod(2,2,:)),'LineWidth',1.5 )
            hold on
            grid on 
        end

        % Plot of FRF
        if plotFRF
            figure(fig3)
            subplot(1,2,1)
            plot (f_vec,abs(squeeze( HH(1,1,:)))) 
            subplot(1,2,2)
            plot (f_vec,abs(squeeze( HH(2,2,:))))
        end

        plotcounter = plotcounter + 1 ;
        
    end

end


%% Graphical post-processing

% Annotation of fig PSDs
if plotPSD
    figure(fig1) ; grid on
    subplot(1,2,1)
    legend(strings, 'interpreter','latex')
    title('Modal PSD $S_{1,1}(\omega)$','interpreter','latex')
    xlabel('Frequences [Hz]', 'FontSize',12, 'interpreter','latex')
    ylabel('PSD [m$^2$/s$^3$]', 'FontSize',12, 'interpreter','latex')
    subplot(1,2,2)
    title('Modal PSD $S_{2,2}(\omega)$','interpreter','latex')
    legend(strings, 'interpreter','latex')
    xlabel('Frequences [Hz]', 'FontSize',12, 'interpreter','latex')
    ylabel('PSD [m$^2$/s$^3$]', 'FontSize',12, 'interpreter','latex')
end

% Annotation of fig FRFs
if plotFRF
    figure(fig3) ; grid on
    subplot(1,2,1)
    legend(strings, 'interpreter','latex')
    xlabel('Frequences [Hz]', 'FontSize',12, 'interpreter','latex')
    ylabel('FRF $H(\omega)$ [s$^2$.kg$^{-1}$]', 'FontSize',12, 'interpreter','latex')
    subplot(1,2,2)
    legend(strings, 'interpreter','latex')
    xlabel('Frequences [Hz]', 'FontSize',12, 'interpreter','latex')
    ylabel('FRF $H(\omega)$ [s$^2$.kg$^{-1}$]', 'FontSize',12, 'interpreter','latex')
end

% Annotation of fig2
if plotfrq
    fig2 = figure ;
    grid on ; hold on
    plot( U_mat, freqSto(1,:), '*-' )
    plot( U_mat, freqSto(2,:), '*-' )  
    plot( U_mat, xiSto(1,:), 'o-' )
    plot( U_mat, xiSto(2,:), 'o-' )
    legend({'$\omega_1$','$\omega_2$','$\xi_1$','$\xi_2$'}, 'interpreter','latex','FontSize',12)
    xlabel('Avg. Wind Speed $U$ [m/s]', 'FontSize',12, 'interpreter','latex')
    ylabel('Eigen Frequencies $f_i$ [Hz]', 'FontSize',12, 'interpreter','latex')
    ylim([0 0.30])
end

% Plot of standard deviation
if plotSTD
   figure ; 
   hold on ;
   grid on ;
   for ii = 1 : 2
       for jj = 1 : 2
           if ii == jj && ii == 1
                plot( U_mat, squeeze(sqrt(varvalues(ii,jj,:))), 'LineWidth',1.5)
           elseif ii == jj && ii == 2
                plot( U_mat, squeeze(sqrt(varvalues(ii,jj,:)))*B/2, 'LineWidth',1.5)
           else
               plot( U_mat, squeeze(varvalues(ii,jj,:)), 'LineWidth',1.5)
           end
       end
   end
   legend({'$\sigma_{x_{zz}}$', '$\sigma_{x_{zt}}$', '$\sigma_{x_{tz}}$', '$\sigma_{x_{tt}}*B/2$'},'interpreter','latex','location','nw')
   xlabel('Frequencies [Hz]','interpreter','latex')
   ylabel('Standard deviation - Co-variance','interpreter','latex')
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
    fs = linspace( 0, 5, 200) ;
    for kk = 1:length(fs)
        for k = 1:length(a)
            aa(k,kk) = a{k}(fs(kk));
            hh(k,kk) = h{k}(fs(kk));
        end
    end
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







    
