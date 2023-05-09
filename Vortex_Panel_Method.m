%% Patrick Kennedy - Term Project
clear, clc, close all; clf;
set(0, 'DefaultFigureWindowStyle', 'docked');
should_plot = false;
airfoil = '2412';
V = 1;
m = 120;

%% Model NACA 2412 body
figure(1);
title(sprintf('NACA 2412 (%d panels)',m));
xlabel('x/c');
ylabel('y/c');
hold on; grid on; axis equal;
xlim([0 1]);

% read airfoil
naca2412 = naca4digit(airfoil, m);

% plot mcl and body
plot(naca2412.x, naca2412.mcl, 'b--'); % lower
plot(naca2412.XY(1,:), naca2412.XY(2,:), 'b-'); % body
scatter(naca2412.XY(1,:), naca2412.XY(2,:), 5); % vertecies
legend('MCL', 'Body', 'Points');

%% iterate accross AoAs
alphas = linspace(-4, 12, 10); %deg
output = cell(size(alphas));
newfig = 1;
for i = 1:length(alphas)
    a = alphas(i); % deg
    newfig = newfig + 1;
    if ~should_plot
        newfig = -1;
    end
    output{i} = vortext_panel_method(naca2412, a, newfig);
    fprintf([ ...
        'AoA %.2f\n\t' ...
        'Cl = %f\n\t' ...
        'Cd = %f\n\t' ...
        'Cmle = %f\n\t' ...
        'xcp = %f\n\t' ...
        'cmc4 = %f\n\n'], a, output{i}.cl, output{i}.cd, output{i}.cmle, ...
                        output{i}.xcp, output{i}.cmc4);
end

%% plot Cl vs a

plot_legend = ["",""];

Cl = zeros(size(output));
Cmc4 = zeros(size(output));
for i = 1:length(output)
    Cl(i) = output{i}.cl;
    Cmc4(i) = output{i}.cmc4;
end

newfig = newfig + 1;
if ~should_plot
    newfig = 2;
end

figure(newfig);
grid on; hold on;
title(sprintf('Flow Characteristics vs AoA (%d panels)', m));
xlabel('AoA (deg)');
plot(alphas, Cl, 'b-');
coeffs = polyfit(alphas, Cl, 1);
fplot(@(a) coeffs(1)*a + coeffs(2), [-4 15], 'r.');

if coeffs(2) >= 0
    plot_legend(end+1) = 'Cl';
    plot_legend(end+1) = sprintf('%.2fa + %.2f',coeffs);
else
    plot_legend(end+1) = 'Cl';
    plot_legend(end+1) = sprintf('%.2fa %.2f',coeffs);
end

% get slope
ao = coeffs(1);
a_Cl0 = -coeffs(2) / coeffs(1);
Cl_a0 = coeffs(2);
fprintf('a0 = %f\n', ao);
fprintf('a_Cl0 = %f\n', a_Cl0);
fprintf('Cl_a0 = %f\n', Cl_a0);

%% plot Cmc/4 vs a
plot(alphas, Cmc4, 'g-');
coeffs = polyfit(alphas, Cmc4, 1);
fplot(@(a) coeffs(1)*a + coeffs(2), [-4 15], '.');

if coeffs(2) >= 0
    plot_legend(end+1:end+2) = ["Cm_c_/_4", sprintf('%.2fa + %.2f',coeffs)];
else
    plot_legend(end+1:end+2) = ["Cm_c_/_4", sprintf('%.2fa %.2f',coeffs)];
end

% get slope
dcmc4_da = coeffs(1);

%% get aerodynamic center
xac = 0.25 - dcmc4_da/ao;
fprintf('x_ac = %f\n', xac);

%% plot xcp vs AoA
Xcp = zeros(size(output));
for i = 1: length(output)
   Xcp(i) = output{i}.xcp;
end

plot(alphas, Xcp);
legend([plot_legend(3:end) 'X_c_p']);

%% plot Cmac
Cmac = -Cl.*(Xcp - xac);
newfig = newfig + 1;
figure(newfig);
grid on; hold on;
plot(Cl, Cmac);
title(sprintf('Cm_a_c vs Cl (%d panels)',m));
xlabel('Cl');
ylabel('Cm_a_c');
coeffs = polyfit(Cl, Cmac, 0);
cmac = coeffs;
fprintf('Cmac fit = %f\n', cmac);

%% VORTEX PANEL METHOD ====================================================
function outs = vortext_panel_method(naca4, aoa_deg, fig)
%% parse parameters

alpha = deg2rad(aoa_deg);

if ~exist('fig','var')
    fig = -1;
end

num_panels = naca4.n;
Xbody = naca4.XY(1,:);
Ybody = naca4.XY(2,:);
num_points = num_panels + 1;

Ybody(1,1) = 0;
Ybody(end,1) = 0;

%% define geometry

X_mid = zeros(num_panels, 1); % x mid point
Y_mid = zeros(num_panels, 1); % y mid point
Sj = zeros(num_panels, 1); % length
theta = zeros(num_panels, 1); % panel angle

for i = 1:num_panels
    X_mid(i) = 0.5*( Xbody(i) + Xbody(i+1) );
    Y_mid(i) = 0.5*( Ybody(i) + Ybody(i+1) );
    Sj(i) = sqrt( (Xbody(i+1)-Xbody(i))^2 + (Ybody(i+1)-Ybody(i))^2 );
    theta(i,1) = atan2( (Ybody(i+1)-Ybody(i)) , (Xbody(i+1)-Xbody(i)) );
end

%% find coeffs

% init
Cn1 = zeros(num_panels);
Cn2 = zeros(num_panels);
Ct1 = zeros(num_panels);
Ct2 = zeros(num_panels);

for i = 1:num_panels
    for j = 1:num_panels

        if (i == j) % identity line

            Cn1(i,j) = -1;
            Cn2(i,j) = 1 ;
            Ct1(i,j) = pi/2;
            Ct2(i,j) = pi/2;

        else

            % current state
            sj = Sj(j);
            xi = X_mid(i);
            Xj = Xbody(j);
            thetai = theta(i);
            thetaj = theta(j);
            yi = Y_mid(i);
            Yj = Ybody(j);

            A = -(xi-Xj)*(cos(thetaj))-(yi-Yj)*(sin(thetaj));
            B = (xi-Xj)^2+(yi-Yj)^2;
            C = sin(thetai-thetaj);
            D = cos(thetai-thetaj);
            E = (xi-Xj)*sin(thetaj)-(yi-Yj)*cos(thetaj);
            F = log(1+(sj^2+2*A*sj)/B);
            G = atan2((E*sj),(B+A*sj));
            P = ((xi-Xj)*sin(thetai-2*thetaj))+((yi-Yj)*cos(thetai-2*thetaj));
            Q = ((xi-Xj)*cos(thetai-2*thetaj))-((yi-Yj)*sin(thetai-2*thetaj));

            Cn1(i,j) = 0.5*D*F+C*G-Cn2(i,j);
            Cn2(i,j) = D+((0.5*Q*F)/Sj(j))-((A*C+D*E)*(G/Sj(j)));
            Ct1(i,j) = 0.5*C*F-D*G-Ct2(i,j);
            Ct2(i,j) = C+((0.5*P*F)/Sj(j))+((A*D-C*E)*(G/Sj(j)));

        end
    end
end
clear sj xi Xj thetai thetaj yi Yj;

%% determine An, RHS, and At 

% create An and RHS
An = zeros(num_points);
RHS = zeros(num_panels, 1);
for i = 1:num_points-1
    An(i,1) = Cn1(i,1);
    An(i,num_points) = Cn2(1,num_points-1);
    RHS(i) = sin(theta(i) - alpha);
    for j = 2:num_points-1
        An(i,j) = Cn1(i,j) + Cn2(i,j-1);
    end
end
An(num_points,1) = 1;
An(num_points,num_points) = 1;
for j = 2:num_points-1
    An(num_points,j) = 0;
end
RHS(num_points) = 0;

% create At
At = zeros(num_panels, num_points);
for i = 1:num_points-1
    At(i,1) = Ct1(i,1);
    At(i,num_points) = Ct2(i,num_points-1);
    for j = 2:num_points-1
        At(i,j) = Ct1(i,j) + Ct2(i,j-1);
    end
end

%% determine Cp
Gamma = An\RHS;

% calc velocity across each panel
V = zeros(1,num_panels);
for i = 1:num_panels
    % do summation
    Vsum = 0;
    for j = 1:num_points
        Vsum = Vsum + At(i,j)*Gamma(j);
    end

    % set velocity
    V(i) = cos(theta(i)-alpha) + Vsum;
end

% calc coeff of pressure across each panel
Cp = zeros(1,num_panels);
for i = 1:num_panels
    Cp(i) = 1 - (V(i))^2;
end

%% numerical integration to find outputs

Cp_l = flip(Cp(1:num_panels/2));
X_l = flip(X_mid(1:num_panels/2));
Cp_u = Cp(num_panels/2+1:end);
X_u = X_mid(num_panels/2+1:end);

% normal coeff
cn = 0;
for i = 1:length(Cp_l)
    dx = naca4.x(i) - naca4.x(i+1);
    cn = cn + (Cp_l(i) - Cp_u(i))*dx;
end

% axial coeff
ca = 0;
for i = 1:length(Cp_l)
    dx = naca4.x(i) - naca4.x(i+1);
    dy_l = Ybody(i) - Ybody(i+1);
    dy_u = Ybody(num_panels/2+i+1) - Ybody(num_panels/2+i);
    ca = ca + (Cp_u(i)*dy_u/dx - Cp_u(i)*dy_l/dx)*dx;
end

% calculate cl and cd
ClCd = [cos(alpha) -sin(alpha) ; sin(alpha) cos(alpha)] * [cn ; ca];
cl = ClCd(1);
cd = ClCd(2);

% coeff moment about LE
cmle = 0;
for i = 1:length(Cp_l)-1
    
    xil = Xbody(i);
    xi1l = Xbody(i+1);
    dxl = xil - xi1l;
    xiu = Xbody(num_points - i+1);
    xi1u = Xbody(num_points - i);
    dxu = xiu - xi1u;
  
    cmle = cmle + ...
        ((Cp_u(i) + Cp_u(i+1))/2)*((xiu+xi1u)/2)*dxu - ...
        ((Cp_l(i) + Cp_l(i+1))/2)*((xil+xi1l)/2)*dxl;
end

% xcp
xcp = -cmle / cl;

% cmc/4
cmc4 = -cl*(xcp-0.25);

% generate outs struct
outs = struct();
outs.cl = cl;
outs.cd = cd;
outs.cmle = cmle;
outs.xcp = xcp;
outs.cmc4 = cmc4;

%% plot
if fig ~= -1 % if should plot
    figure(fig);
    hold on;
    plot(naca4.XY(1,:), naca4.XY(2,:), 'b-'); % body
    plot(X_u,Cp_u, 'g');
    plot(X_l,Cp_l, 'r');
    plot(X_u,Cp_u+Cp_l);
    xline(xcp);
    set(gca,'Ydir','reverse');
    xlabel('x/c');
    ylabel('Coefficient of Pressure');
    grid on;
    title(sprintf('NACA %s @ AoA = %d', naca4.name, aoa_deg));
    legend('Body','Upper','Lower','Cp_u+Cp_l');
end

end

%% NACA 4-DIGIT ===========================================================

% naca4(mpxx) returns a characterized airfoild using A&V mean cord line
% equations. This function will parse 4 digit numbers and return a struct
% of the airfoil characteristics
% 
% digits = string: format 'mpxx' that represents NACA 4 digits
% n = int: number of panels for XY model, will only produce even number of
% panels, any odd number will be round down with floor()
function naca4 = naca4digit(digits, n)
naca4 = struct(); % init struct
naca4.name = digits;
naca4.n = n;
m = str2double(digits(1)) / 100; % parse m
p = str2double(digits(2)) / 10; % parse p
xx = str2double(digits(3:4)) / 100; % parse xx

% add demensions to airfoil
naca4.m = m;
naca4.p = p;
naca4.xx = xx;
naca4.Rle = 1.1019 * xx^2; % nose radius

% Model MCL
syms x_
yt = xx*(1.4845.*sqrt(x_) - 0.63.*x_ - 1.758.*x_.^2 ...
    + 1.4215.*x_.^3 - 0.5075.*x_.^4); % thickness distrubution
dyc_dx = diff(yt,x_);

% convert yt and dyc/dt to function handel
yt = matlabFunction(yt);
dyc_dx = matlabFunction(dyc_dx);

% add to airfoil
naca4.yt = yt;
naca4.dyt_dx = dyc_dx;

% mcl equations for fwd and aft
mcl_fwd = (2*m/p).*x_ - (m/p^2) .* x_.^2;
mcl_aft = (-m/(p^2 -2*p +1))*(2*p-1 - 2*p.*x_ + x_.^2);

% account for singularity
if p == 0
    mcl_fwd = 0;
end

% find derivate of y_mcl(x)
dmcl_fwd = diff(mcl_fwd, x_);
dmcl_aft = diff(mcl_aft, x_);

% convert syms to function handel
mcl_fwd = matlabFunction(mcl_fwd,'Vars',x_);
mcl_aft = matlabFunction(mcl_aft,'Vars', x_);
dmcl_fwd = matlabFunction(dmcl_fwd,'Vars',x_);
dmcl_aft = matlabFunction(dmcl_aft,'Vars',x_);

% conditional f(x) for mcl and dyc/dx
mcl = @(x) (x < p) .* mcl_fwd(x) ... % fwd
        + ~(x < p) .* mcl_aft(x); % aft

dyc_dx = @(x) (x < p) .* dmcl_fwd(x) ... % fwd
           + ~(x < p) .* dmcl_aft(x); % aft

% iterate from theta = 0 to 180 and collect x locations on r=0.5 circle
thetas = linspace(0, 180, floor(n/2)+1);
x = zeros(1,length(thetas));
r = 0.5;
for i = 1:length(thetas)
    theta = thetas(i);
    x(i) = r + r*cosd(theta);
end

% format x coords for upper and lower surface
xu = flip(x);
xl = x;

% trim x coords from LE and TE
xu = xu(2:end-1);
xl = xl(2:end-1);

% define body coords
Xu = @(x) x - yt(x).*sin(atan(dyc_dx(x)));
Xl = @(x) x + yt(x).*sin(atan(dyc_dx(x)));
Yu = @(x) mcl(x) + yt(x).*cos(atan(dyc_dx(x)));
Yl = @(x) mcl(x) - yt(x).*cos(atan(dyc_dx(x)));

% parse into vec of coords
% add (1, 0) in front and back, and (0, 0) in between upper and lower
XY = [1 Xl(xl) 0 Xu(xu) 1 ; % TE -> lower -> LE -> upper -> TE
      0 Yl(xl) 0 Yu(xu) 0]; % TE -> lower -> LE -> upper -> TE

% add values to airfoil
naca4.x = x;
naca4.mcl = mcl(x);
naca4.XY = XY;

end
