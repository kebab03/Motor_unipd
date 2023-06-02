function [ ga, gb, gc ] = calc_dq2abc( gd, gq, theta )

theta = theta*pi/180;

ga = gd*cos(theta)          - gq*sin(theta);
gb = gd*cos(theta - 2*pi/3) - gq*sin(theta - 2*pi/3);
gc = gd*cos(theta - 4*pi/3) - gq*sin(theta - 4*pi/3);

end

