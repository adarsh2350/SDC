%% Data sharing
function [n, K, T, D, T_avg, opr_chrgs, base_fare, reg_chrg, reg_lim, add_chrg, tax, profit_factor, dis_factor, DS, CS] = SDC_input


n = 4;   % no. of zones

K = 2;   % no. of types of cabs

T = 600;  % Max. time cab can be driven throughout the day

D = 2;   % Dry run

T_avg = [20 35];  % avg. time for one ride

opr_chrgs = [26000 33000]; % operational charges for a cab
base_fare = [25 35]; % base fare
reg_chrg = [9 12]; % regular distance charges at d=0
reg_lim = 15; % max. distance to implement add_chrg
add_chrg = [12 16]; % additional distance charges
tax = 18; % taxes 

profit_factor = 1.2;  % minimum profit on operational charges
dis_factor = 0.8; % minimum discount from current pricing 
DS = 0.8;  % Driver's share on profit
CS = 0.2;  % Company's share on profit

