%% Assignment 5, First Order Upwind Differencing Scheme

% Setting Variables
n = 5;
L = 1;
gamma = 1;
v_in = 1;
phi1 = 1;
phi0 = 0;
rho = 1;

analyticalSoln = [1, 0.9387, 0.7963, 0.6244, 0.4100, 0.1505, 0]; %from Versteeg Eg 5.2

[plotX, numericalSolnUDS] = FOUpwindScheme(n, L, gamma, v_in, phi1, phi0, rho);

disp(['Numerical solution of First Order Upwind differencing scheme =', num2str(numericalSolnUDS)])

% Calculate the error between numerical and analytical solutions
error = numericalSolnUDS - analyticalSoln(1:end);
figure(1);
plot(plotX, numericalSolnUDS, "go--", LineWidth=1)
hold on;
plot(plotX, analyticalSoln(1:end), "r--", LineWidth=1)
xlabel('Length');
ylabel('Flux');
title('Flux Distribution');
grid on;

figure (2);
plot(plotX, error, "b.-", LineWidth=1)
xlabel('Length');
ylabel('Error');
title('Error between Numerical and Analytical Solutions');
grid on;
hold off;

figure(3);
[plotX, numericalSolnUDSR] = FOUpwindScheme(n, L, gamma, -v_in, 0, 1, rho);
disp(['Numerical solution of central differencing scheme (Reversed) =', num2str(numericalSolnUDS)])
plot(plotX,numericalSolnUDSR,"bo--",LineWidth=1)
hold on;
xlabel('Length');
ylabel('Flux');
title('Flux Distribution and error');
grid on;
legend('Numerical Solution', 'Analytical Solution','Error','Reversed UDS')

legend("Position", [0.68182,0.58719,0.27857,0.15595])
%%
function[x,phi_N]=FOUpwindScheme(n,L,gamma,v_in,phi1,phi0,rho)
%------FIRST ORDER UPWIND DIFFERENCING SCHEME------
dx = L/n; % cell spacing

x(1) = 0; % generate x spacing for plotting
x(2) = dx/2;
for i=3:n+1
 x(i) = x(i-1) + dx;
end
x(n+2) = L;

dxw = zeros(1,n+2); 
dxe = zeros(1,n+2); 
phi_N=zeros(1,n+2);
Fw = zeros(1,n+2); 
Fe = zeros(1,n+2); 
Dw = zeros(1,n+2); 
De = zeros(1,n+2);

% Calculate neighbor distances and convection terms
for i=2:n+1
  dxw(i) = (x(i) - x(i-1));  % Distance to west neighbor
  dxe(i) = (x(i+1) - x(i));  % Distance to east neighbor
  Fe(i) = rho*v_in;           % East face convection
  Fw(i) = rho*v_in;           % West face convection (same inlet velocity)
  De(i) = gamma/dxe(i);       % East face diffusion
  Dw(i) = gamma/dxw(i);       % West face diffusion
end

% Initialize coefficients for tridiagonal matrix system
aw = zeros(1,n+2);      % West face coefficient
ae = zeros(1,n+2);      % East face coefficient
ap = zeros(1,n+2);      % Diagonal coefficient
d = zeros(1,n+2);       % Right-hand side of the system
a = zeros(1,n+2);       % Lower diagonal coefficient
b = zeros(1,n+2);       % Upper diagonal coefficient
c = zeros(1,n+2);       % Source terms
dnew = zeros(1,n+2);    % Modified diagonal for forward substitution
cnew = zeros(1,n+2);    % Modified source term for forward substitution

% Calculate coefficients for interior nodes (2nd to n-th)
for i=2:n+1
 aw(i) = Dw(i) + max(0,Fw(i)); ae(i) = De(i) + max(0,-Fe(i)); 
 ap(i) = aw(i) + ae(i) + Fe(i) - Fw(i);
 a(i) = -ae(i); b(i) = -aw(i); d(i) = ap(i);
end

d(1) = 1; 
b(1) = 0; 
a(1) = 0; 
c(1) = phi1;
d(n+2) = 1; 
b(n+2) = 0; 
a(n+2) = 0; 
c(n+2) = phi0;
dnew(1) = d(1);
cnew(1) = c(1);

for i=2:n+2
 dnew(i) = d(i) - b(i)*a(i-1)/dnew(i-1); cnew(i) = c(i) - b(i)*cnew(i-1)/dnew(i-1);
end

phi_N(n+2) = cnew(n+2)/dnew(n+2);
for i = n+1:-1:1
 phi_N(i) = (cnew(i) - a(i)*phi_N(i+1))/dnew(i);
 hold on
end

end