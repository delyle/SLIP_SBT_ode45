function i_closest = closestapproach(ti,vi,xi,yi,T2)

% find initial position

x0 = xi(1); y0 = yi(1);

% create takeoff vector after apex

i_flight = ti >= T2(end) & vi <= 0;

% find distance between flight trajectory and initial position

d2 = (xi(i_flight) - x0).^2 + (yi(i_flight) - y0).^2;

[~,i_closest_flight] = min(d2);
i_closest = i_closest_flight + find(i_flight,1,'first');


