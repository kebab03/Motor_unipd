function [ gd, gq ] = calc_abc2dq( ga, gb, gc, theta )

theta = theta*pi/180;

gd =  2/3*( ga*cos(theta) + gb*cos(theta - 2*pi/3) + gc*cos(theta - 4*pi/3) );
gq = -2/3*( ga*sin(theta) + gb*sin(theta - 2*pi/3) + gc*sin(theta - 4*pi/3) );

end
