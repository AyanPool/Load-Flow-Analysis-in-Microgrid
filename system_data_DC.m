function [nDG,nPQ,Rv,R_line,PL0,r,LF,LT,nbus,nline] = system_data_DC(MG_system)

switch(MG_system)
    case 1

Pbase = 5000;
Vbase = 380;

%% System parameters of six-bus DC Microgrid

nDG = 3;
nPQ = 3;

% Control Parameters 
R_v = [0.3; 0.4; 0.2];

% Line Parameters 
R_line = [ 1          4      0.20; 
           2          5      0.40; 
           3          6      0.30;
           4          5      0.25;
           5          6      0.20;];

LF = R_line(:,1); 
LT = R_line(:,2);

% Load Parameters 
PLdata = [0; 0; 0; 4500; 5700; 6350];
PL0 = PLdata/Pbase;

r_base = Vbase^2/Pbase;
r = R_line(:,3) / r_base;
Rv = R_v/r_base;
nbus = 6;
nline = length(R_line(:,1));
end