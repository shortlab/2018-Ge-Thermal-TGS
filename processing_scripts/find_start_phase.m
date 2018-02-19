function [start_time]=find_start_phase(time,amplitude,num_half_periods,grat)
%Function to find the start time for a given fixed null-point start
%   time: time vector for TGS data
%   amplitude: vector of recorded TGS amplitudes
%   num_half_periods: is the null-point at which you decide to start. Best-choice based on analysis is for null-point 2
%   grat: grating spacing at which measurement was conducted

%%There are four possibly different topgies for the inital TGS trace based on what the acoustic phase looks like at the
%% zero time registered by the amplitude grating measurement. Based on performance, a robost search method for arbitrary 
%% fixed null-point starts has not been developed. Rather, manual determinations have been coded for the first four null-point
%% for each of the four morphologies. You cannot select num_half_periods>4. The morphology is detected by the locations of
%% initial maxima and minima in the trace and the null-point start time selected based on that morphology determination.

%%contact cdennett (at) mit.edu for a more complete description of each of the possible cases here.

if num_half_periods>4
    display('Cannot start this far out, defaulting to zero start')
    start_time=time(1);
end

if grat>6
    [~,pos_loc]=findpeaks(amplitude(20:end),time(20:end),'NPeaks',5);
    [~,neg_loc]=findpeaks(-1*amplitude(20:end),time(20:end),'NPeaks',5);
    start_point=0.5e-9;
else
    [~,pos_loc]=findpeaks(amplitude(2:end),time(2:end),'NPeaks',5);
    [~,neg_loc]=findpeaks(-1*amplitude(2:end),time(2:end),'NPeaks',5);
    start_point=0;
end


if neg_loc(1)<pos_loc(1)
    check_length=pos_loc(1)-neg_loc(1);
    if neg_loc(1)-check_length/2<start_point %Case #3
        if num_half_periods==1
            start_time=neg_loc(1)+1/2*(pos_loc(1)-neg_loc(1));
        elseif num_half_periods==2
            start_time=pos_loc(1)+1/2*(neg_loc(2)-pos_loc(1));
        elseif num_half_periods==3
            start_time=neg_loc(2)+1/2*(pos_loc(2)-neg_loc(2));
        elseif num_half_periods==4
            start_time=pos_loc(2)+1/2*(neg_loc(3)-pos_loc(2));
        end
    else
        if num_half_periods==1
            start_time=neg_loc(1)-1/2*(pos_loc(1)-neg_loc(1));
        elseif num_half_periods==2
            start_time=neg_loc(1)+1/2*(pos_loc(1)-neg_loc(1));
        elseif num_half_periods==3
            start_time=pos_loc(1)+1/2*(neg_loc(2)-pos_loc(1));
        elseif num_half_periods==4
            start_time=neg_loc(2)+1/2*(pos_loc(2)-neg_loc(2));
        end
    end
else
    check_length=neg_loc(1)-pos_loc(1);
    if pos_loc(1)-check_length/2<start_point %Case #1
        if num_half_periods==1
            start_time=pos_loc(1)+1/2*(neg_loc(1)-pos_loc(1));
        elseif num_half_periods==2
            start_time=neg_loc(1)+1/2*(pos_loc(2)-neg_loc(1));
        elseif num_half_periods==3
            start_time=pos_loc(2)+1/2*(neg_loc(2)-pos_loc(2));
        elseif num_half_periods==4
            start_time=neg_loc(2)+1/2*(pos_loc(3)-neg_loc(2));
        end
    else
        if num_half_periods==1
            start_time=pos_loc(1)-1/2*(neg_loc(1)-pos_loc(1));
        elseif num_half_periods==2
            start_time=pos_loc(1)+1/2*(neg_loc(1)-pos_loc(1));
        elseif num_half_periods==3
            start_time=neg_loc(1)+1/2*(pos_loc(2)-neg_loc(1));
        elseif num_half_periods==4
            start_time=pos_loc(2)+1/2*(neg_loc(1)-pos_loc(2));
        end
    end
end

end