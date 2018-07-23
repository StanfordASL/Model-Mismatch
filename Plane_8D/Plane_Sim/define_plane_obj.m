function [ plane, plane_box ] = define_plane_obj()

alpha_0 = 5*pi/180;
lambda = 0.1;

%% Vertices

body_V =...
    [10., -0.5, -1.;
    10., 0.5, -1.;
    10., -0.5, 0.;  
    10., 0.5, 0.;
    5., -0.5, 1.;  
    5., 0.5, 1.;
    -10., -0.5, 1.;  
    -10., 0.5, 1.;
    -10., -0.5, 0.; 
    -10., 0.5, 0.;
    -7.5, -0.5, -1.; 
    -7.5, 0.5, -1.]*lambda;

wing_V = [ 3., 0., 3*sin(alpha_0)+.1;
         3., 0., 3*sin(alpha_0)-.1;
         0., -10., .1;
         0., -10., -.1;
         0., 10., .1;
         0., 10., -.1;
        -2., -10., -2*sin(alpha_0)+.1;
        -2., -10., -2*sin(alpha_0)-.1;
        -2., 10., -2*sin(alpha_0)+.1;
        -2., 10., -2*sin(alpha_0)-.1]*lambda;
    
stab_V = [ -8., 0., .5; 
           -8., 0., .1;
           -9., -4., .5;  
           -9., -4., .1;
           -9., 4., .5;  
           -9., 4., .1;
           -10., -4., .5;
           -10., -4., .1;
           -10., 4., .5;
           -10., 4., .1]*lambda;

tail_V =  [-7., -.2, 1.;  
           -6., .2, 1.;
            -10., -.2, 1.;
            -10., .2, 1.;
            -10., -.2, 4.;
            -10., .2, 4.;
            -9., -.2, 4.;
            -9., .2, 4.]*lambda;
        
        
box_V = [-1, -1,-0.1;
             1, -1, -0.1;
             1, 1, -0.1;
             -1,1,-0.1;
             -1, -1,0.4;
             1, -1, 0.4;
             1, 1, 0.4;
             -1,1,0.4];             
        
%% Define convex hull 

body = convhull(body_V(:,1),body_V(:,2),body_V(:,3));
wing = convhull(wing_V(:,1),wing_V(:,2),wing_V(:,3));
stab = convhull(stab_V(:,1),stab_V(:,2),stab_V(:,3));
tail = convhull(tail_V(:,1),tail_V(:,2),tail_V(:,3));

box = convhull(box_V(:,1),box_V(:,2),box_V(:,3));

%% Plane obj

plane.body = struct('hull',body,'V',body_V,'nV',size(body_V,1));
plane.wing = struct('hull',wing,'V',wing_V,'nV',size(wing_V,1));
plane.stab = struct('hull',stab,'V',stab_V,'nV',size(stab_V,1));
plane.tail = struct('hull',tail,'V',tail_V,'nV',size(tail_V,1));

plane.box = struct('hull',box,'V',box_V,'nV',size(box_V,1));

%% Plot

plot_plane(plane,zeros(3,1),zeros(3,1),1);
                
end
            
