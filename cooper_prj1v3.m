clear; clc;

d = (1/12)*[2 2.5 1.5]; %diameter (ft)
L = [20 10 30]; %length (ft)
e = [0.00085 0.00085 0.00085]; %roughness factor (ft)
visc = .0000205; %viscosity (lbf*s/ft^2)
rho = 1.94; %density (slug/ft^3)
A = (pi/4)*d.^2; %area (ft^2)
Q = 750 * 231/(12^3)/60; %volumetric flow rate (ft^3/s)

syms x %solution vector x = [v1 v2 v3 Pressuredrop]

% Equations 1, 2, and 3 are derived from the darcy-weisbach equation, with
% x(1) representing v1, x(2) representing v2, x(3) representing v3, and
% x(4) representing delta P, the pressure drop. All equations are set equal
% to zero, so the pressure drop x(4) is brought over to the same side as
% all of the other terms. Equation 4 comes from the relation Q1 + Q2 + Q3 =
% Q, where Q is brought over to the left side such that V1*A1 + V2*A2 +
% V3*A3 - Q = 0. See derivation of these equations in the pdf in onedrive.

eq1 = @(x) (32*visc*x(1)*L(1)) / (d(1)^2) - x(4);
eq2 = @(x) (32*visc*x(2)*L(2)) / (d(2)^2) - x(4);
eq3 = @(x) (32*visc*x(3)*L(3)) / (d(3)^2) - x(4);
eq4 = @(x) x(1)*A(1) + x(2)*A(2) + x(3)*A(3) - Q;

F = @(x) [eq1(x); eq2(x); eq3(x); eq4(x)];

val = fsolve(F,[-31 14 422 1]); 
% this solution vector val = [v1, v2, v3, Pressure drop]
% with velocities in units (ft/s) and pressure drop in units (lbf/ft^2 or
% psf) converges to specific values no matter the initial guess.


v = val(1:3); deltaP = val(4);
% storing the solved values for velocity and pressure drop

Qfin = [0 0 0];
for i = 1:3
    Qfin(i) = v(i)*A(i)*60*(12^3)/231;
end
% taking the velocity results and plugging in to find Q1, Q2, Q3, and then
% converting from ft^3/s to gal/min
% (ft^3/s)*(60s/1min)*(12^3in^3/1ft^3)*(1gal/231in^3) = gal/min

deltaP = deltaP/144; %converting from psf to psi

% NOTE: here the pressure drop is wildly incorrect BECAUSE the reynolds
% number for this scenario ends up being around 400000, indicating that our
% assumption of f = 64/Re is incorrect. 

eq1b = @(x) ((log10( (e(1)/d(1))/3.7 + 5.74/(rho*x(1)*d(1)/visc)^0.9 ))^-2) ...
    * (L(1)/d(1)) * (x(1)^2)/8 * rho - x(4); 
eq2b = @(x) ((log10( (e(2)/d(2))/3.7 + 5.74/(rho*x(2)*d(2)/visc)^0.9 ))^-2) ...
    * (L(2)/d(2)) * (x(2)^2)/8 * rho - x(4); 
eq3b = @(x) ((log10( (e(3)/d(3))/3.7 + 5.74/(rho*x(3)*d(3)/visc)^0.9 ))^-2) ...
    * (L(3)/d(3)) * (x(3)^2)/8 * rho - x(4); 
% these equations come from the approximate friction factor expression for
% Re > 2300

F2 = @(x) [eq1b(x); eq2b(x); eq3b(x); eq4(x)];

% the "function evaluation limit" for fsolve is maxed out when solving this
% system, so the intermediate value of each fsolve is plugged back into the
% command. this loops until a solution is reached.
temp = fsolve(F2, [1 1 1 1], optimoptions('fsolve','Display','iter'));

it = 1; %iteration counter

% loop checks if F2(x) is getting close to zero. If this requires more than
% 500 iterations of fsolve, the script aborts.
while norm(F2(temp)) > 1e-5
    temp = fsolve(F2, temp);
    it = it +1;
    if it>500
        clc;
        disp('Convergence failed, try another initial guess.')
        return
    end
end

clc;


% results are displayed to the screen in correct units
disp('Assuming Laminar flow: ')
disp(['The flow through pipe 1 is: ' num2str(Qfin(1)) ' gal/min'])
disp(['The flow through pipe 2 is: ' num2str(Qfin(2)) ' gal/min'])
disp(['The flow through pipe 3 is: ' num2str(Qfin(3)) ' gal/min'])
disp(['The pressure drop across the system is: ' num2str(deltaP) ' psi'])

disp('--------------------------------------------------------------------')

disp('Assuming turbulent flow: ')
disp(['The flow through pipe 1 is: ' num2str(temp(1)*A(1)*60*(12^3)/231) ' gal/min'])
disp(['The flow through pipe 2 is: ' num2str(temp(2)*A(2)*60*(12^3)/231) ' gal/min'])
disp(['The flow through pipe 3 is: ' num2str(temp(3)*A(3)*60*(12^3)/231) ' gal/min'])
disp(['The pressure drop across the system is: ' num2str(temp(4)/144) ' psi'])

