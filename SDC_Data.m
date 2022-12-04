%% Data sharing
function [dist, time, t_freq, cst_C] = SDC_Data


%% City is divided into n zones

% Distance between two zones (in km)

dist = [0 2.4 5.7 13.4;
        2.4 0 2.3 17.7;
        6.1 2.2 0 6.4
        13.4 17.7 6.4 0];

% time taken to travel between two zones (in min)

time = [0 7 14 33;
        7 0 6 35;
        18 6 0 22;
        33 35 22 0];

% frequency of rides booked between two zones (on a particular day)

t_freq = zeros(4,4,2);

t_freq(:,:,1) = [0 150 180 300;
                 210 0 120 130;
                 200 140 0 130;
                 200 160 120 0];

t_freq(:,:,2) = [0 100 100 120;
                120 0 50 60;
                90 30 0 30;
                180 50 40 0];

% cost of current rides

cst_C = zeros(4,4,2);

cst_C(:,:,1) = [0 90 160 350;
                100 0 100 450;
                180 90 0 160;
                400 425 150 0];

cst_C(:,:,2) = [0 140 225 400;
                150 0 140 520;
                250 140 0 220;
                450 500 210 0];
