clear; clc;

% The function takes inputs in the form pipewizard([x1 x2 x3 x4], Q), where
% the vector x is an INITIAL GUESS for the problem. x(1:3) are volumetric
% flow rates (gal/min) and x(4) is a pressure drop (psi). Q, the total flow
% rate through the network, is given and should be inputted in gal/min.

% NOTE: the length and diameter parameters currently inside the function
% are for problem #2.


Q = 750;
x = [250 102 300 10];



% pipewizard([10 10 10 10],750);



% function [x] = pipewizard(x,Q)

d = (1/12)*[2 2 2]; %diameter (ft)
L = [10 10 10]; %length (ft)
e = [0.00085 0.00085 0.00085]; %roughness factor (ft)
visc = 2.05e-5; %viscosity (lbf*s/ft^2)
rho = 1.94; %density (slug/ft^3) kept in slugs to avoid 32.2 ft/s^2 factor
A = (pi/4)*d.^2; %area (ft^2)
Q = Q / 60/(12^3)*231; %volumetric flow rate (ft^3/s)
% INPUT GIVEN PARAMETERS HERE ^^^

% the input to the function is the initial guess with x(1:3) as volumetric
% flow rates (gal/min) and x(4) as a pressure drop (pound-force/square inch
% -- psi)

x = [x(1)/A(1)*.00222801 x(2)/A(2)/60/(12^3)*231 ...
    x(3)/A(3)/60/(12^3)*231 x(4)*144];
% converting the initial guess from flow rate in gallon/min to velocity in
% ft/s and pressure drop from psf to psi. Our system is written in terms of
% velocities v1, v2, v3 of the respective pipes and the pressure drop.
% These are our four unknowns. 

eq1 = @(x) myf1(x(1),e(1),d(1),rho,visc) * (L(1)/d(1))* (x(1)^2)/2 - x(4)/rho;
eq2 = @(x) myf1(x(2),e(2),d(2),rho,visc) * (L(2)/d(2))* (x(2)^2)/2 - x(4)/rho;
eq3 = @(x) myf1(x(3),e(3),d(3),rho,visc) * (L(3)/d(3))* (x(3)^2)/2 - x(4)/rho;
eq4 = @(x) x(1)*A(1) + x(2)*A(2) + x(3)*A(3) - Q;
F = @(x) [eq1(x); eq2(x); eq3(x); eq4(x)];
% functions are defined to initialize the following loop. note the 'norm'
% command is the 'magnitude' (sqrt(F1^2 + F2^2 + F3^2 + F4^2)) of the F
% vector. If the norm of F is close to zero, the system is solved. These
% functions are repeatedly redefined during the loop (the values of x, and
% subsequently the value of all the friction factors, are changing after
% each fsolve iteration).

% Equations 1, 2, and 3 are derived from the darcy-weisbach equation, with
% x(1) representing v1, x(2) representing v2, x(3) representing v3, and
% x(4) representing delta P, the pressure drop. All equations are set equal
% to zero, so the pressure drop x(4) is brought over to the same side as
% all of the other terms. Equation 4 comes from the relation Q1 + Q2 + Q3 =
% Q, where Q is brought over to the left side such that V1*A1 + V2*A2 +
% V3*A3 - Q = 0. See derivation of these equations in the pdf in onedrive.

x = abs(fsolve(F,x,optimoptions('fsolve','Display','iter','MaxFunctionEvaluations',1000)));

disp(['The flow through pipe 1 is: ' num2str(x(1)*A(1)*60*(12^3)/231) ' gal/min'])
disp(['The flow through pipe 2 is: ' num2str(x(2)*A(2)*60*(12^3)/231) ' gal/min'])
disp(['The flow through pipe 3 is: ' num2str(x(3)*A(3)*60*(12^3)/231) ' gal/min'])
disp(['The pressure drop across the system is: ' num2str(x(4)) ' psi'])
% results are displayed to the command window in appropriate units

% myf1(x(1),e(1),d(1),rho,visc)

% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fric = myf1(vel, e, d, rho, visc)
% this function allows the caluclation of the friction factor given the
% velocity through the pipe (m/s), the roughness factor (ft), the diameter
% (ft), the density (slug/ft^3) and viscosity (lbf*s/ft^2)

fricPipe = @(fric) (fric)^(-1/2) + ....
    2*log10( (e/d)/3.7 +2.51/(abs(vel)*rho*d/visc * sqrt(fric)) );
% function in the form fricPipe(fric) = 0 which allows us to solve the
% Colebrook equation numerically

Re = abs(vel)*rho*d/visc;
% Reynold's number is calculated so that the correct friction equation can
% be selected using the if/else statement below.

frApprox = ((log10( (e/d)/3.7 + 5.74/(Re)^0.9 ))^-2);
% the approximate solution to the Colebrook equation is calculated as used
% as an initial guess for its solution.

if Re <= 2300
    fric = 64/Re;
else
    fric = fsolve(fricPipe,frApprox,optimoptions('fsolve','Display','none'));
    % If the Reynolds number is greater than 2300, the Colebrook equation
    % is solved using the function fricPipe(fric) = 0 and the initial guess
    % fApprox
end

end
