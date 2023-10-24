%these are the global variables used in functions below.
global D_O 
global beam_angle
global beam_dia
global ptv_rad 
global ptv_center
global oar_a
global oar_b
global oar_c
global oar_center
global head_a
global head_b
global head_c
global head_center
%initializing some global variables.
D_O = 1;
beam_angle = 30;
beam_dia = 30;
ptv_rad = 15;
ptv_center = [30,0,15];
oar_a = 15;
oar_b = 25;
oar_c = 15;
oar_center = [0,30,45];
head_a = 80;
head_b = 100;
head_c = 80;
head_center = [0,0,0];
global struct_arr
global daf_table
global rdf_table
%Draw_3D_Scene()

daf_table = Compute_Depth_Dose();

rdf_table = Compute_Radial_Dose();

struct_arr = Compute_Beam_Directions();

Compute_Skin_Entry_Points();

Compute_Beam_Safety_Flags();
%these lines are for testing or running the functions.
%radial_distance = Compute_Radial_Distance([45,0,15],18);
%Depth = Compute_Depth_from_Skin([45,0,15], 18);
%dose = Compute_Point_Dose_from_Beam([30,0,15], 18)
%dose = Compute_Point_Dose_from_All_Beams([30,0,15])
%Compute_Surface_Dose_PTV();
%Compute_Surface_Dose_OAR();
%Compute_Volume_Dose_PTV();
%Compute_Volume_Dose_OAR();
%this function is for drawing 3D scene of head, OAR, PTV and the isocenter.
function Draw_3D_Scene()
    %These are variables used in the function.
    global ptv_rad
    global ptv_center
    global oar_center
    global oar_a
    global oar_b
    global oar_c
    global head_center
    global head_a
    global head_b
    global head_c
    
    %this is for creating the structure of the head.
    [X_head, Y_head, Z_head] = ellipsoid(head_center,head_center,head_center,head_a,head_b,head_c);
    %this is for creating the structure of PTV.
    [X,Y,Z] = sphere;
    X_ptv = X * ptv_rad + ptv_center(1);
    Y_ptv = Y * ptv_rad + ptv_center(2);
    Z_ptv = Z * ptv_rad + ptv_center(3);
    %this is for creating the structure of OAR.
    [X_oar, Y_oar, Z_oar] = ellipsoid(oar_center(1),oar_center(2),oar_center(3),oar_a,oar_b,oar_c);
    
    %the isocenter equals to the center of ptv.
    isocenter = ptv_center;
    
    %this is for plotting the head.
    head_plot = surf(X_head,Y_head,Z_head,'FaceAlpha',0.1,'FaceColor','green');
    hold on
    %this is for plotting axis.
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    %this is for plotting the structure of PTV
    ptv_plot = surf(X_ptv,Y_ptv,Z_ptv,'FaceAlpha',0.5,'FaceColor','yellow');
    hold on 
    %this is for plotting the structure of OAR
    oar_plot = surf(X_oar,Y_oar,Z_oar,'FaceAlpha',0.5,'FaceColor','red');
    hold on
    %this is for plotting isocenter
    isocenter_plot = plot3(isocenter(1),isocenter(2),isocenter(3),'.','markersize',10);

end
%this function is for creating depth dose function table
%the table will be created based on the given graph for depth dose.
function depth_dose_function_table = Compute_Depth_Dose()
    %the table will contain depth and dose in two columns
    depth_dose_function_table = [];
    %the depth starts from 0 and ends at 200
    ini_depth = 0;
    final_depth = 200;
    while ini_depth <= final_depth
        %computing dose using DAF function.
        dose = DAF(ini_depth);
        data = [ini_depth, dose];
        depth_dose_function_table = [depth_dose_function_table; data];
        ini_depth = ini_depth + 1;
    end
end

%this function is for creating radial dose function table
%the table will store radial distance and dose in two columns.
function radial_dose_function_table = Compute_Radial_Dose()
    radial_dose_function_table = [];
    ini_distance = -23;
    final_distance = 23;
    while ini_distance <= final_distance
        %computing dose using RDF function
        dose = RDF(ini_distance);
        data = [ini_distance,dose];
        radial_dose_function_table = [radial_dose_function_table; data];
        ini_distance = ini_distance + 1;
    end
end

%this function is for computing dose depending on depth.
function dose = DAF(d)
%this is created based on a given graph for depth dose
    if d <= 20
        dose = 0.5 + 0.025 * d;
    else
        dose = 1 - 0.005 * (d-20);
    end
end
%this function is for computing dose depending on radial distance.
function dose = RDF(d)
%this is created based on a given graph for radial distance dose.
    if d > 22.5
        dose = 0;
    elseif d > 7.5
        dose = - 1/15 * (d - 22.5);
    elseif d > -7.5
        dose = 1;
    elseif d > -22.5
        dose = 1/15 * (d + 22.5);
    else
        dose = 0;
    end
end

%this function is for comping beam direction vectors.
function data_set = Compute_Beam_Directions()
    
    data_set = [];
    %first the horizontal vectors are created separated by 30 degrees.
    arc = [];
    for i = 30:30:150
        x = 30*sind(i);
        y = 30*cosd(i);
        z = 0;
        temp_vec = [x, y, z];
        temp_vec = temp_vec/norm(temp_vec);
        %latitude and longitude are computed.
        lat = atan2(z,sqrt(x*x+y*y)) / pi * 180;
        long = atan2(y,x) / pi * 180;
        data = [temp_vec, long, lat];
        
        data_set = [data_set; data];
        arc = [arc;temp_vec];
        
    end
    %this for loop is for computing direction vectors after rotating the
    %horizontal cluster of vector about y axis by 30 degrees every time. 
    [m n] = size(arc);
    for i = -30:-30:-180
        %computing rotation matrix
        [rot_mat,~] = rotation_axis('y',i); 
        for j = 1:m
            %new vector created
            temp_vec = transpose(rot_mat * transpose(arc(j,:)));
            %computing longitude and latitude.
            x = temp_vec(1);
            y = temp_vec(2);
            z = temp_vec(3);
            lat = atan2(z,sqrt(x*x+y*y)) / pi * 180;
            long = atan2(y,x) / pi * 180;
            data = [temp_vec, long, lat];
            data_set = [data_set; data];
        end    
    end
    %this is for storing data for vectors that are parallel to y-axis
    temp_vec = [0,1,0];
    x = temp_vec(1);
    y = temp_vec(2);
    z = temp_vec(3);
    %computing longitude and latitude.
    lat = atan2(z,sqrt(x*x+y*y)) / pi * 180;
    long = atan2(y,x) / pi * 180;
    data = [temp_vec, long, lat];
    data_set = [data_set; data];
    
    temp_vec = [0,-1,0];
    x = temp_vec(1);
    y = temp_vec(2);
    z = temp_vec(3);
    %computing longitude and latitude.
    lat = atan2(z,sqrt(x*x+y*y)) / pi * 180;
    long = atan2(y,x) / pi * 180;
    data = [temp_vec, long, lat];
    data_set = [data_set; data];

    %this is for plotting all the direction vectors.
    %{
    [m n] = size(data_set);
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    hold on
    for i = 1:m
        data = data_set(i,:);
        dir = data(1:3);
        plot3([0;dir(1)],[0;dir(2)],[0;dir(3)])
        hold on
    end
    %}
end


%{
This function is for computing homogeneous rotation matrix depending on the
axis and an angle of rotation.
%}
function [Mat3, Mat4] = rotation_axis(axis, angle)
    angle = angle / 180 * pi;
    %this case is where rotation about X axis is happening
    if (axis == 'x') || (axis == 'X')
        Mat3 = [[1,0,0];[0,cos(angle),-sin(angle)];[0,sin(angle),cos(angle)]];
        Mat4 = [[1,0,0,0];[0,cos(angle),-sin(angle),0];[0,sin(angle),cos(angle),0];[0,0,0,1]];
    %this case is where rotation about Y axis is happening
    elseif (axis == 'y') || (axis == 'Y')
        Mat3 = [[cos(angle),0,sin(angle)];[0,1,0];[-sin(angle),0,cos(angle)]];
        Mat4 = [[cos(angle),0,sin(angle),0];[0,1,0,0];[-sin(angle),0,cos(angle),0];[0,0,0,1]];
    %this case is where rotation about Z axis is happening
    elseif (axis == 'z') || (axis == 'Z')
        Mat3 = [[cos(angle),-sin(angle),0];[sin(angle),cos(angle),0];[0,0,1]];
        Mat4 = [[cos(angle),-sin(angle),0,0];[sin(angle),cos(angle),0,0];[0,0,1,0];[0,0,0,1]];
    else
        disp("Enter the appropriate axis");
    end
end

%this function is for computing skin entry points.
function Compute_Skin_Entry_Points()
    global head_center
    global head_a
    global head_b
    global head_c
    global ptv_center
    global struct_arr
    [m n] = size(struct_arr);
    new_arr = [];
    for i = 1:m
        %direction vector is computed
        data = struct_arr(i,:);
        dir = data(1:3);
        
        %computing skin entry point
        P = [30,0,15];
        [n,inter] = intersect_line_and_elipsoid(P,dir,80,100,80);
        
        %the correct intersecting skin entry point is found.
        intersection = inter(2,:);
        %computing depth from each skin entry point
        depth = round(sqrt(dot((P-intersection),(P-intersection))));
        data = [data, intersection, depth];
        new_arr = [new_arr; data];
    end
    struct_arr = new_arr;

    %this is for plotting direction vectors, head, and entry points.
    [m n] = size(struct_arr);
    %this is for constructing head
    [X_head, Y_head, Z_head] = ellipsoid(0,0,0,80,100,80);
    head_plot = surf(X_head,Y_head,Z_head,'FaceAlpha',0.1,'FaceColor','green');
    hold on
    
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    hold on
    %this is for plotting direction vectors.
    for i = 1:m
        data = struct_arr(i,:);
        %this is the direction vector
        dir = data(1:3);
        plot3([ptv_center(1);ptv_center(1)+150*dir(1)],[ptv_center(2);ptv_center(2)+150*dir(2)],[ptv_center(3);ptv_center(3)+150*dir(3)])
        hold on
    end

    %this is for plotting skin entry points.
    for i = 1:m
        data = struct_arr(i,:);
        %this is the skin entry point
        point = data(6:8);
        plot3(point(1),point(2),point(3),'.','markersize',10);
        hold on
    end
    
end
%{
this function is for computing the intersection of a line and a elipsoid
and the number of the intersections
%}
function [n_intersection, intersection] = intersect_line_and_elipsoid(P,v,a,b,c)
    syms t;
    %{
    this is the equation of a elipsoid where the coordinate of the point 
    from the line is substituted.
    %}
    eqn = (P(1) + t * v(1))^2 / a^2 + (P(2) + t * v(2))^2 / b^2 + (P(3) + t * v(3))^2 / c^2 == 1;
    %{
    solving for t is for computing the exact coordinate of the
    intersections
    %}
    sltn = solve(eqn,t,'Real',true);
    intersection = [P(1) + sltn * v(1), P(2) + sltn * v(2), P(3) + sltn * v(3)];
    [m n] = size(intersection);
    n_intersection = m;
    %this is the case where there is no intersections.
    if n_intersection == 0
        %an empty matrix is returned for having no intersections.
        intersection = [];

    %this is the case where the 2 solutions are the same.
    %However, this case would be considered as having one intersection.
    elseif isequal(intersection(1,:),intersection(2,:))
        intersection = intersection(1,:);
        n_intersection = 1;
    end
end


%this function is for computing beam safety of each beam
%this the function will store 'True' value if it is safe and store 'False'
%if the function is hitting OAR.
function Compute_Beam_Safety_Flags()
    global struct_arr;
    global oar_center;
    global oar_a;
    global oar_b;
    global oar_c;
    [m n] = size(struct_arr);
    new_arr = [];
    syms x y z;
    unsafe = [];
    %this is for computing safety of each beam
    for i = 1:m
        data = struct_arr(i,:);
        dir = data(1:3);
        eqn1 = x*x / 225 + (y-30)*(y-30) / 625 + (z-45)*(z-45) / 225 == 1;
        eqn2 = distance_of_point_from_line([30,0,15],dir,[x, y, z]) == 15;
        %vpasolve is used for solving system of equations.
        sltn = vpasolve([eqn1,eqn2],[x, y, z]);
        %this is the case where the beam hits OAR.
        if size(sltn.x) >= 1
            safe = 'False';
            unsafe = [unsafe; eqn2];
        %this is the case where the beam does not hit OAR.
        else
            safe = 'True';
        end
        data = [data,safe];
        new_arr = [new_arr; data];
    end
    struct_arr = new_arr;
    
    %This is for plotting OAR and unsafe beams.
    %{
    [X_oar,Y_oar,Z_oar] = ellipsoid(oar_center(1),oar_center(2),oar_center(3),oar_a,oar_b,oar_c);
    surf(X_oar,Y_oar,Z_oar,'FaceAlpha',0.5,'FaceColor','red');
    hold on
    [p q] = size(unsafe);
    %this is for plotting each unsafe beams
    for i = 1:p
        fimplicit3(unsafe(i,:))
        hold on
    end
    %}
    
end

%this function is for computing radial distance
function dist = Compute_Radial_Distance(P,ind)
    global struct_arr;
    data = struct_arr(ind,:);
    dir = data(1:3);
    %a function from HW1 is used for computing the distance
    dist = round(distance_of_point_from_line([30,0,15],dir,P));
end

%this function is for computing depth from the skin
%a point is treated as if it was located on a plane with the normal
%vector that is equal to the beam direction vector.
function Depth = Compute_Depth_from_Skin(Point, ind)
    global struct_arr;
    data = struct_arr(ind,:);
    %direction vector
    dir = data(1:3);
    %skin entry point
    skin = data(6:8);
    %the intersect between a line and a plane is computed.
    depth_point = intersect_line_and_plane(Point,dir,skin,dir);
    %the distance between the intersect and the skin entry point
    Depth = round(sqrt(dot((skin-depth_point),(skin-depth_point))));
end

%this is for computing dose from a beam.
function dose = Compute_Point_Dose_from_Beam(Point, ind)
    global daf_table
    global rdf_table
    global D_O
    
    depth = Compute_Depth_from_Skin(Point, ind);
    rad = Compute_Radial_Distance(Point, ind);
    %if the radius is smaller than 23, it will look up the rdf table.
    if rad <23
        rdf = rdf_table(rad + 23,2);
    %if the radius is not smaller than 23, the rdf will be 0
    else
        rdf = 0;
    end
    %looking up daf table
    daf = daf_table(depth+1,2);
    
    dose = D_O * rdf * daf;
    
end

%this is for computing point dose from all beams
function dose = Compute_Point_Dose_from_All_Beams(Point)
    global struct_arr
    [m n] = size(struct_arr);
    
    dose = 0;
    %the dose will be sum of all dose from each beam
    for i = 1:m
        if struct_arr(i,10) == 'True'
            dose = dose + Compute_Point_Dose_from_Beam(Point,i);
        end
    end
end

%this is for computing surface dose of PTV
function Compute_Surface_Dose_PTV()
    global D_O
    global struct_arr
    global ptv_rad
    global ptv_center
    %this is for creating the structure of PTV
    [p q] = size(struct_arr);
    [X, Y, Z] = sphere;
    Xptv = X * ptv_rad + ptv_center(1);
    Yptv = Y * ptv_rad + ptv_center(2);
    Zptv = Z * ptv_rad + ptv_center(3);

    ptv_hot_point = [];
    %the hottest dose will start from the smallest value 
    ptv_hot_dose = 0;
    ptv_cold_point = [];
    %the coldest dose will start from the largest value.
    ptv_cold_dose = D_O * p;
    [m n] = size(Xptv);
    %this is for computing all the dose of points on the surface.
    for i = 1:m
        for j = 1:n
            %this is a point on the surface
            P = [Xptv(i,j),Yptv(i,j),Zptv(i,j)];
            
            dose = Compute_Point_Dose_from_All_Beams(P);
            %the hottest dose will increase after each trial of comparison
            if dose > ptv_hot_dose
                ptv_hot_point = P;
                ptv_hot_dose = dose;
            end
            %the coldest dose will decrease after each trial of comparison
            if dose < ptv_cold_dose
                ptv_cold_point = P;
                ptv_cold_dose = dose;
            end
        end
    end
    %print the hottest dose, the coldest dose and locations
    disp("Hottest PTV surface dose");
    disp(ptv_hot_dose);
    disp("Hottest PTV surface location");
    disp(ptv_hot_point);
    disp("Coldest PTV surface dose");
    disp(ptv_cold_dose);
    disp("Coldest PTV surface location");
    disp(ptv_cold_point);
    
    %this is for plotting ptv and the location of dose using different
    %colours
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    ptv_plot = surf(Xptv,Yptv,Zptv,'FaceAlpha',0.5,'FaceColor','yellow');
    hold on
    plot3(ptv_hot_point(1),ptv_hot_point(2),ptv_hot_point(3),'.','markersize',20,'Color','magenta');
    hold on
    plot3(ptv_cold_point(1),ptv_cold_point(2),ptv_cold_point(3),'.','markersize',20,'Color','green');
    hold on

end

%this is for computing surface dose of OAR
function Compute_Surface_Dose_OAR()
    global D_O
    global struct_arr
    global oar_a
    global oar_b
    global oar_c
    global oar_center
    [p q] = size(struct_arr);
    %this is for creating structure of an ellipsoid.
    [Xoar, Yoar, Zoar] = ellipsoid(oar_center(1),oar_center(2),oar_center(3),oar_a,oar_b,oar_c);
    
    oar_hot_point = [];
    %the hottest dose starting from the smallest value
    oar_hot_dose = 0;
    oar_cold_point = [];
    %the coldest dose starting from the coldest value.
    oar_cold_dose = D_O * p;
    [m n] = size(Xoar);
    %this is for computing dose for each point on the surface
    for i = 1:m
        for j = 1:n
            %this is the point on the surface
            P = [Xoar(i,j),Yoar(i,j),Zoar(i,j)];
            dose = Compute_Point_Dose_from_All_Beams(P);
            %the hottest dose will increase
            if dose > oar_hot_dose
                oar_hot_point = P;
                oar_hot_dose = dose;
            end
            %the coldest dose will decrease
            if dose < oar_cold_dose
                oar_cold_point = P;
                oar_cold_dose = dose;
            end
        end
    end
    %printing the hottest and the coldest dose and the locations
    disp("Hottest OAR surface dose");
    disp(oar_hot_dose);
    disp("Hottest OAR surface location");
    disp(oar_hot_point);
    disp("Coldest OAR surface dose");
    disp(oar_cold_dose);
    disp("Coldest OAR surface location");
    disp(oar_cold_point);
    
    %this is for plotting OAR and the location of the hottest and the
    %coldest dose using different colours.
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    oar_plot = surf(Xoar,Yoar,Zoar,'FaceAlpha',0.1,'FaceColor','red');
    hold on
    plot3(oar_hot_point(1),oar_hot_point(2),oar_hot_point(3),'.','markersize',20,'Color','yellow');
    hold on
    plot3(oar_cold_point(1),oar_cold_point(2),oar_cold_point(3),'.','markersize',20,'Color','green');
    hold on

end
%this function is for computing DVH for PTV
function Compute_Volume_Dose_PTV()
    global ptv_rad
    global ptv_center
    
    %this is for compuing dose for each point inside the structure.
    dose_list = [];
    for i = 1:ptv_rad
        for j = 1:ptv_rad
            for k = 1:ptv_rad
                %it will be moving in 8 different direction from center
                P1 = ptv_center + [i*2,j*2,k*2];
                P2 = ptv_center + [-i*2,j*2,k*2];
                P3 = ptv_center + [i*2,-j*2,k*2];
                P4 = ptv_center + [-i*2,-j*2,k*2];
                P5 = ptv_center + [i*2,j*2,-k*2];
                P6 = ptv_center + [-i*2,j*2,-k*2];
                P7 = ptv_center + [i*2,-j*2,-k*2];
                P8 = ptv_center + [-i*2,-j*2,-k*2];
                P = [P1;P2;P3;P4;P5;P6;P7;P8];
                for l = 1:8
                    point = P(l,:);
                    %make sure that the point is inside the structure.
                    if sqrt(dot((point - ptv_center),(point - ptv_center))) <= ptv_rad
                        %dose inside the structure
                        dose = Compute_Point_Dose_from_All_Beams(point);
                        
                        dose_list = [dose_list, dose];
                    end
                end
                P = [];
            end
        end
    end
    %this is for compuing dose for each point inside the structure.
    for i = 1:ptv_rad
        %this will be moving in 6 different directions.
        P1 = ptv_center + [i*2,0,0];
        P2 = ptv_center + [-i*2,0,0];
        P3 = ptv_center + [0,i*2,0];
        P4 = ptv_center + [0,-i*2,0];
        P5 = ptv_center + [0,0,i*2];
        P6 = ptv_center + [0,0,-i*2];
        P = [P1;P2;P3;P4;P5;P6];
        for l = 1:6
            point = P(l,:);
            %make sure that the point is inside the structure.
            if sqrt(dot((point - ptv_center),(point - ptv_center))) <= ptv_rad
                dose = Compute_Point_Dose_from_All_Beams(point);
                dose_list = [dose_list, dose];
            end
        end
        P = [];
        
    end
    %this is for computing dose of the center.
    dose = Compute_Point_Dose_from_All_Beams(ptv_center);
    dose_list = [dose_list, dose];
    %creating histogram
    histo = histogram(dose_list,10)
end

%this function is for computing dose for point inside OAR
function Compute_Volume_Dose_OAR()
    global oar_a
    global oar_b
    global oar_c
    global oar_center
    count = 0;
    dose_list = [];
    %this is for computing dose for each point inside the structure
    for i = 1:oar_a
        for j = 1:oar_b
            for k = 1:oar_c
                %it will move in 8 different direction start from center
                P1 = oar_center + [i*2,j*2,k*2];
                P2 = oar_center + [-i*2,j*2,k*2];
                P3 = oar_center + [i*2,-j*2,k*2];
                P4 = oar_center + [-i*2,-j*2,k*2];
                P5 = oar_center + [i*2,j*2,-k*2];
                P6 = oar_center + [-i*2,j*2,-k*2];
                P7 = oar_center + [i*2,-j*2,-k*2];
                P8 = oar_center + [-i*2,-j*2,-k*2];
                P = [P1;P2;P3;P4;P5;P6;P7;P8];
                for l = 1:8
                    point = P(l,:);
                    %make sure the point is inside the structure
                    if (point(1)-oar_center(1))^2/(oar_a*oar_a) + (point(2)-oar_center(2))^2/(oar_b*oar_b) + (point(3)-oar_center(3))^2/(oar_c*oar_c)<= 1
                        dose = Compute_Point_Dose_from_All_Beams(point)
                        count = count +1
                        
                        dose_list = [dose_list, dose];
                    end
                end
                P = [];
            end
        end
    end
    %this is for computing dose for each point inside the structure
    for i = 1:oar_a
        %it will move in 2 direction from center
        P1 = oar_center + [i*2,0,0];
        P2 = oar_center + [-i*2,0,0];
        P = [P1;P2];
        for l = 1:2
            point = P(l,:);
            %make sure the point is inside the structure
            if (point(1)-oar_center(1))^2/(oar_a*oar_a) + (point(2)-oar_center(2))^2/(oar_b*oar_b) + (point(3)-oar_center(3))^2/(oar_c*oar_c)<= 1
                dose = Compute_Point_Dose_from_All_Beams(point);
                dose_list = [dose_list, dose];
                
            end
        end
        P = [];    
    end
    %this is for computing dose for each point inside the structure
    for i = 1:oar_b
        %it will move in 2 direction from center
        P1 = oar_center + [0,i*2,0];
        P2 = oar_center + [0,-i*2,0];
        P = [P1;P2];
        for l = 1:2
            point = P(l,:);
            %make sure the point is inside the structure
            if (point(1)-oar_center(1))^2/(oar_a*oar_a) + (point(2)-oar_center(2))^2/(oar_b*oar_b) + (point(3)-oar_center(3))^2/(oar_c*oar_c)<= 1
                dose = Compute_Point_Dose_from_All_Beams(point);
                dose_list = [dose_list, dose];
            end
        end
        P = [];    
    end
    %this is for computing dose for each point inside the structure
    for i = 1:oar_c
        %it will move in 2 direction from center
        P1 = oar_center + [0,0,i*2];
        P2 = oar_center + [0,0,-i*2];
        P = [P1;P2];
        for l = 1:2
            point = P(l,:);
            %make sure the point is inside the structure
            if (point(1)-oar_center(1))^2/(oar_a*oar_a) + (point(2)-oar_center(2))^2/(oar_b*oar_b) + (point(3)-oar_center(3))^2/(oar_c*oar_c)<= 1
                dose = Compute_Point_Dose_from_All_Beams(point)
                dose_list = [dose_list, dose];
            end
        end
        P = [];    
    end
    %computing dose for the center
    dose = Compute_Point_Dose_from_All_Beams(oar_center);
    dose_list = [dose_list, dose];
    %creating histogram
    histo = histogram(dose_list,10)

end

function dist = distance_of_point_from_line(P,v,A)
    %This is the distance between the point P from the line and a point A.
    PA = sqrt(dot(P-A,P-A));
    %v will be normalized
    v = v / norm(v);
    %this compute the length of projection of PA onto vector v.
    proj = abs(dot(v,(P-A)));
    %the distance between between the line and the point will be computed 
    %using the length of PA and the projection.
    dist = sqrt((PA+proj)*(PA-proj));
end
%{
this function is used for computing the intersect 
between a line and a plane
%}
function intersection = intersect_line_and_plane(A,n,P,v)
    %{
    The dot product of normal vector of a plane and direction vector of a
    line is computed. 
    %}
    if (dot(n,v) == 0)
        %this is the case where the line is on the plane.
        if (dot((A-P),n) == 0)
            syms t;
            disp("The line is on the plane")
            %t is any scalar value. 
            %This would imply that there are infinite solutions.
            intersection = P + t*v;

        %this is the case where the line and plane are parallel.
        else
            disp("The line and a plane are parallel")
            %an empty matrix is returned.
            intersection = [];
        end
    else
        %each vector is normalized
        n = n/norm(n);
        v = v/norm(v);
        AP = A-P;
        %t is a coefficient of the direction vector of a line.
        t = dot(AP,n) / dot(v,n);
        %computing intersection using t value.
        intersection = P+v*t;
    end
end


