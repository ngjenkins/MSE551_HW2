%% MSE 551 HW2

%% Question 1
% Define Conversion Parameters for LJ Units (Input Actual Values For
% Different Atom Types) should you decide to use them
mass = 1; sigma = 1; epsilon = 1;

% Define Stepping Parameters for Interatomic Separation Distance r
r_i = 0.9; r_f = 4.0; N = 25;
r_step = (r_f-r_i)/N;
r = (r_i:r_step:r_f);

% Include Imported Energy Values From Summary Data Files as E_LJ
PATH = mfilename('fullpath');
DATA = struct2cell(importdata([PATH(1:end-4),'\HW2_1_SUMMARY.txt'],' '));
E_LJ = DATA{1,1}(:,3);

figure
plot(r,E_LJ)
title('Q1: Energy vs Separation [LJ]')
xlabel('$r/\sigma$','interpreter','latex')
ylabel('$E/\epsilon$','interpreter','latex')

%% Question 2
% Relevant units of type metal are distance=Ang, energy=eV,

% List minimized energies and corresponding lattice constants
% for each crystal structure
E_sc = -3.0627
a_sc = 2.687
E_bcc = -6.6032
a_bcc = 3.219
E_fcc = -13.273
a_fcc = 4.050

% Given that we know Al takes an FCC structure the results are not
% surprising.

%% Question 3
% Relevant units of type metal are distance=Ang, energy=eV, volume=Ang^3

% Define Stepping Parameters for lattice parameter a
a_i = 3.0; a_f = 5.0; N = 25;
a_step = (a_f-a_i)/N;
a = (a_i:a_step:a_f);
V = a.^3;

% Include Imported Energy Values From Summary Data Files as E_Al
PATH = mfilename('fullpath');
DATA = struct2cell(importdata([PATH(1:end-4),'\HW2_3_SUMMARY.txt'],' '));
E_Al = DATA{1,1}(:,3)';
FIT = polyfit(V,E_Al,10); x_fit = linspace(min(V),max(V),100);
P = polyval(FIT,x_fit);

figure
plot(V,E_Al,'o',x_fit,P)
title('Q1: Energy vs UC Volume [Al]')
xlabel('V [$\AA^3$]','interpreter','latex')
ylabel('E [$eV$]','interpreter','latex')

% Calculate first and second order differences as ~derivatives, scaling for
% step size, then evaluating at global maximum to find curvature at minimum
[M,I] = min(P);
V0 = x_fit(I);
S = diff(P)/(x_fit(2)-x_fit(1));
D = diff(P,2)/(x_fit(2)-x_fit(1));
dE2d2v = D(I);
% Use above parameters to ~bulk modulus
B = V0*dE2d2v;
% Scale to typical units
% 1 eV = 1.619e-19 J
% 1 ang = 1e-10 m
B = B*(1.619e-19)*(1e30)*1e-9; % GPa

disp(['The bulk modulus is ~',char(string(B)),' GPa, which is on the proper order of magnitude and in the neighborhood of potential values for Al'])

%% Question 4

% My project will focus on an ideal graphene oxide system. The goal is to
% simulate the vibrational behavior of chains of common GO functional
% groups (OH, COOH, maybe bridged -O- groups). The project will focus on
% two key requirements for the program: the ability to generate any number
% of arbitrary chains of functional groups in various periodicities, and
% the ability to analyze the vibrational phase (and amplitude) of each of 
% these groups for different meta-structures. Ideally, the position of
% funtional groups should be able to be procedurally generated from 2-D
% unit cell structures (both lattice and basis).







