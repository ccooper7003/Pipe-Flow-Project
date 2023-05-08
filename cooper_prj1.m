clear; clc;

d = (1/12)*[2 2 2]; %diameter (ft)
L = [10 10 10]; %length (ft)
visc = .0000205; %viscosity (lbf*s/ft^2)
rho = 1.94*32.2; %density (lbm/ft^3)
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

% eq1 = @(x) (32*visc*x(1)*L(1)) / (d(1)^2) - x(4);
% eq2 = @(x) (32*visc*x(2)*L(2)) / (d(2)^2) - x(4);
% eq3 = @(x) (32*visc*x(3)*L(3)) / (d(3)^2) - x(4);
% eq4 = @(x) x(1)*A(1) + x(2)*A(2) + x(3)*A(3) - Q;
% 
% F = @(x) [eq1(x); eq2(x); eq3(x); eq4(x)];
% 
% val = fsolve(F,[1 14 422 1]) 
% this solution vector val = [v1, v2, v3, Pressure drop]
% with velocities in units (ft/s) and pressure drop in units (lbf/ft^2 or
% psf) converges to specific values no matter the initial guess.

% 
% FAILED ATTEMPT AT SOLVING FOR FRICTION FACTORS
eq1 = @(x) ((32.2*2048*visc^2*L(1)) / (x(1)*(d(1)^3)*rho))^-1 - x(4)^-1;
eq2 = @(x) ((32.2*2048*visc^2*L(2)) / (x(2)*(d(2)^3)*rho))^-1 - x(4)^-1;
eq3 = @(x) ((32.2*2048*visc^2*L(3)) / (x(3)*(d(3)^3)*rho))^-1 - x(4)^-1;
eq4 = @(x) A(1)/(x(1)*d(1)) + A(2)/(x(2)*d(2)) + A(3)/(x(3)*d(3)) - (rho*Q)/(32.2*64*visc);

F = @(x)  [(32.2*2048*visc^2*L(1)) / (x(4)*(d(1)^3)*rho) - x(1); ...
    (32.2*2048*visc^2*L(2)) / (x(4)*(d(2)^3)*rho) - x(2); ... 
    (32.2*2048*visc^2*L(3)) / (x(4)*(d(3)^3)*rho) - x(3); ... 
    A(1)/(x(1)*d(1)) + A(2)/(x(2)*d(2)) + A(3)/(x(3)*d(3)) - (rho*Q)/(32.2*64*visc)];

F([1 1 1 1])

