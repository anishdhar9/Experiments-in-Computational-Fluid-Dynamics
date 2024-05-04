%% Assignment 5, Central Differencing Scheme
% 

% Setting Variables
n = 5;
L = 1;
gamma = 1;
v_in = 1;
phi1 = 1;
phi0 = 0;
rho = 1;

analyticalSoln = [1, 0.9387, 0.7963, 0.6244, 0.4100, 0.1505, 0];

[plotX, numericalSolnCDS] = CentralDifferencingScheme(n, L, gamma, v_in, phi1, phi0, rho);

disp(['Numerical solution of central differencing scheme =', num2str(numericalSolnCDS)])

% Calculate the error between numerical and analytical solutions
error = numericalSolnCDS - analyticalSoln(1:end);
figure;
plot(plotX, numericalSolnCDS, "go--", LineWidth=1)
hold on;
plot(plotX, analyticalSoln(1:end), "r--", LineWidth=1)
xlabel('Length');
ylabel('Flux');
title('Flux Distribution');
legend('Numerical Solution', 'Analytical Solution')
grid on;

figure;
plot(plotX, error, "b.-", LineWidth=1)
xlabel('Length');
ylabel('Error');
title('Error between Numerical and Analytical Solutions');
grid on;

hold off;

[plotX, numericalSolnCDS] = CentralDifferencingScheme(n, L, gamma, -v_in, 0, 1, rho);
disp(['Numerical solution of central differencing scheme (Reversed) =', num2str(numericalSolnCDS)])
disp(['Numerical solution of central differencing scheme (Reversed) =',num2str(numericalSolnCDS)])
plot(plotX,numericalSolnCDS,"bo--",LineWidth=1)
xlabel('Length');
ylabel('Flux');
title('Flux Distribution');
grid on;
hold off;
%%
function[x,phi_N]=CentralDifferencingScheme(n,L,gamma,v_in,phi1,phi0,rho)
%------CENTRAL DIFFERENCING SCHEME------
dx = L/n; % cell spacing

x(1) = 0; % generate x spacing for plotting
x(2) = dx/2;
for i=3:n+1
 x(i) = x(i-1) + dx;
end
x(n+2) = L;

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
for i=3:n
  aw(i) = Dw(i) + 0.5*Fw(i);     % West face coefficient (diffusion + convection)
  ae(i) = De(i) - 0.5*Fw(i);     % East face coefficient (diffusion - convection)
end

% Adjust coefficients for boundary nodes (1st and n+1th)
aw(2) = Dw(2) + Fw(2);          % West face coefficient at left boundary (full convection)
ae(2) = De(2) - 0.5*Fe(2);      % East face coefficient at left boundary
aw(n+1) = Dw(n+1) + 0.5*Fw(n+1); % West face coefficient at right boundary
ae(n+1) = De(n+1) - Fe(n+1);     % East face coefficient at right boundary (full convection)

% Build tridiagonal coefficient matrix system
for i=2:n+1
  ap(i) = aw(i) + ae(i) + Fe(i) - Fw(i); % Diagonal coefficient (sum of all terms)
  a(i) = -ae(i);                         % Lower diagonal coefficient
  b(i) = -aw(i);                         % Upper diagonal coefficient
  d(i) = ap(i);                           % Right-hand side
end

% Apply boundary conditions
d(1) = 1; b(1) = 0; a(1) = 0; c(1) = phi1; % Left boundary (fixed phi)
d(n+2) = 1;
b(n+2) = 0; 
a(n+2) = 0; 
c(n+2) = phi0;
dnew(1) = d(1); cnew(1) = c(1);

for i=2:n+2
 dnew(i) = d(i) - b(i)*a(i-1)/dnew(i-1); cnew(i) = c(i) - b(i)*cnew(i-1)/dnew(i-1);
end

phi_N(n+2) = cnew(n+2)/dnew(n+2);
for i = n+1:-1:1
 phi_N(i) = (cnew(i) - a(i)*phi_N(i+1))/dnew(i);
end

end