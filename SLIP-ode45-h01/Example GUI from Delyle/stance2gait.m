function [ F ] = stance2gait( F_S, t0, tf,t)
%stance2gait Converts force data from stance time reference frame to gait
%            reference frame. Will also interpolate the data so that all
%            the arrays are through the same time steps. For use with 
%            Farell_2014_analysis.m
% F_S = Array where the first column is time normalized to stance period 
%        (t = 0 is touchdown, t = 1 is liftoff). Time should be in 
%        ASCENDING order!
% t0 = touchdown time normalized by gait period
% tf = liftoff time normalized by gait period
% t = the time array to interpolate over.

F_T = F_S;

% test if first column not in ascending order
test = sort(F_T(:,1));
test = abs(F_T(:,1) - test);
if any(test>0)
    error(['First column of input array not in ascending order!',char(10),...
            'You should do something about that'])
end

% delete any values where ts < 0 or ts > 1
[i,~] = find(F_T(:,1)<0 | F_T(:,1)>1);
F_T(i,:) = [];

% Make the first row [0,0]
if F_T(1,1) == 0;
    F_T(1,2) = 0;
else
    F_T = [[0,0];F_T];
end

%Make the last value [1,0]
if F_T(end,1) == 1;
    F_T(1,2) = 0;
else
    F_T = [F_T;[1,0]];
end


if tf < t0
    tf = 1 + tf;
end

F_T(:,1) = F_T(:,1)*(tf-t0) + t0;
F_T(F_T(:,1) > 1,1) = F_T(F_T(:,1) > 1,1) - 1;
F_T = wrtsort(F_T);
% These next two lines important to avoid interpolation errors at ends.
F_T(end+1,:) = [1+F_T(1,1),F_T(1,2)]; 
F_T = [[F_T(end-1,1)-1,F_T(end-1,2)]; F_T];
F = interp1(F_T(:,1),F_T(:,2),t,'Pchip');
F = F';
F(isnan(F)) = 0;

end

