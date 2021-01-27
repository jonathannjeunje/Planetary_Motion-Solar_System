%% ODE Solving describing the motion of the bodies in the solar system
function [] = SolarSystem_Sim()
    %% GLOBAL VARIABLE DEFINITION *****************************************
    global nb; % The number of considered bodies.
    global ib; % The channel-indice (in the 3D matrix) of the planet. 
    global is; % The indice of the sun.
    global nv; % The Number of Initial Values per bodies
    global ip_xyz; % The indices of the xyz-coordinates respectively of the bodies' position and position.
    global iv_xyz; % The indices of the xyz-coordinates respectively of the bodies' position and velocity.
    global days_in_a_yr % Average number of seconds in an Earth year.
    global names_b; % List of bodies' name
    global colors_b; % List of bodies' colors
    global scales_b; % List of bodies' scales gotten from the log of their mean radius (in km)
    global inner_planets; % List of planets closer to the sun
    global outer_planets; % List of planets further from the sun
    global interacting_bodies; %
    
    global G; % The Universal gravitational constant (in m^3 kg^-1 s^-2)
    global Mb; % The Mass of the Bodies (in kg)
    
    global reset_env % Resets the environement
    global time_of_sim_in_yrs; % The time of simulation (in Earth years)
    global show_plot; % Show plot.
    global show_animation; % Show Animation.
    global sun_fixed; % Use the sun as reference origin for the system.
    global animation_speed; % The speed of simulation.
    global animation_bodies; % List of bodies to animate.
    global view_from_top; % View from top.
    global grid_on; % Grid on.
    global h; % The Step size for ODE Solver.
    global import_file_name; % The file's name containing previously computed data.
    global export_data; % If computed data should be exported or not.
    global generate_video; % Enables video of the simulation to be generated.
    
    %% INITIALISATION: Core settings **************************************
    nb = 10; ib = 1; is = 1; nv = 6; ip_xyz = [1 2 3]; iv_xyz = [4 5 6]; 
    days_in_a_yr = 365.25; inner_planets = 2:5; outer_planets = 6:10;
    names_b = {'Sun','Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune','Pluto'};
    colors_b = {[1 0.5 0],[0.5 0.5 0.5],[1 0.9 0],[0.3 0.6 0.8],[0.6 0.2 0.4],[0.6 0 0.3],[1 1 0],[0.3 0.8 0.8],[0.1 0.7 0.8],'c'};
    scales_b = log([6.95500E+08 2.44000E+06 6.05180E+06 6.37101E+06 3.38990E+06 6.99110E+07 5.82320E+07 2.53620E+07 2.46240E+07 1.19500E+06]); 
    
    %% INITIALISATION: User settings **************************************
    reset_env = 'yes'; view_from_top = 'yes'; grid_on = 'yes'; 
    sun_fixed = 'no'; show_plot = 'yes'; show_animation = 'yes';
    export_data = 'yes'; generate_video = 'no';
    ode_solver = 'import'; % ode45, ode113, ode652, ode15s, import
%     import_file_name = '2018.05.06.212048-165yrs-by1-ode45.dat';
%     import_file_name = '2018.04.20.094756-165yrs-by1-ode113.dat';
    import_file_name = '2018.04.20.123252-165yrs-by1-ode652.dat';
    time_of_sim_in_yrs = 165; % The period of revolution of Neptune:165yrs Pluto:248yrs
    h = 1; animation_speed = 30;
    interacting_bodies = [1:9];
    animation_bodies = interacting_bodies([1:9]); % Select the bodies to annimate.
    
    %% INITIALISATION: Reset Environment **********************************
    if strcmp(reset_env,'yes')
        clc; % Clears Command line.
        delete(findall(0,'Type','figure')); % Clears Previous figures.
    end
    
    %% INITIALISATION: Solar system constants and Initial Values. *********
    G = 6.67E-11; % (in m^3 kg^-1 s^-2)
    G = ((6.67E-11)*(24*3600)^2)/(1000^3); % (in km^3 kg^-1 day^-2)
    Mb = [1.98854E+30,3.30200E+23,4.86850E+24,5.97219E+24,6.41850E+23,...
            1.89813E+27,5.68319E+26,8.68103E+25,1.02410E+26,1.30700E+22];
    
    bp_INIT = [1.81899E+08	9.83630E+08     -1.58778E+07
            -5.67576E+10	-2.73592E+10	2.89173E+09
            4.28480E+10     1.00073E+11     -1.11872E+09
            -1.43778E+11	-4.00067E+10	-1.38875E+07
            -1.14746E+11	-1.96294E+11	-1.32908E+09
            -5.66899E+11	-5.77495E+11	1.50755E+10
            8.20513E+10     -1.50241E+12	2.28565E+10
            2.62506E+12     1.40273E+12     -2.87982E+10
            4.30300E+12     -1.24223E+12	-7.35857E+10
            1.65554E+12     -4.73503E+12	2.77962E+10
            ]./1000; % (in Km)

    bv_INIT = [-1.12474E+01	7.54876E+00     2.68723E-01
            1.16497E+04     -4.14793E+04	-4.45952E+03
            -3.22930E+04	1.36960E+04     2.05091E+03
            7.65151E+03     -2.87514E+04	2.08354E+00
            2.18369E+04     -1.01132E+04    -7.47957E+02
            9.16793E+03     -8.53244E+03	-1.69767E+02
            9.11312E+03     4.96372E+02     -3.71643E+02
            -3.25937E+03	5.68878E+03     6.32569E+01
            1.47132E+03     5.25363E+03     -1.42701E+02
            5.24541E+03     6.38510E+02     -1.60709E+03
            ].*(24*3600/1000); % (in km/day)

    %% INITIALISATION: Initial values column vector construction **********
    y_INIT = zeros(nv*nb,1); % Initialisation
    for ib = interacting_bodies
        y_INIT(ip_xyz+(ib-1)*nv) = bp_INIT(ib,:);
        y_INIT(iv_xyz+(ib-1)*nv) = bv_INIT(ib,:);
    end
    
    %% ODE SOLVER: Run Selected ode solver ********************************
    tspan = (0:h:time_of_sim_in_yrs*days_in_a_yr);
    switch ode_solver
        case 'ode45'
            [T, Y] = ode45(@SolarSystem_ODE, tspan, y_INIT);
        case 'ode113'
            [T, Y] = ode113(@SolarSystem_ODE, tspan, y_INIT);
        case 'ode652'
            [T, Y] = ode652(@SolarSystem_ODE, tspan, y_INIT);
        case 'ode15s'
            [T, Y] = ode15s(@SolarSystem_ODE, tspan, y_INIT);
        case 'import'
            [T, Y] = Import_Data(import_file_name);
    end

    %% Export Data: Save computed data to file ****************************
    if strcmp(export_data,'yes') && not(strcmp(ode_solver,'import'))
        Export_Data(T, Y);
    end
    
    %% PLOT: Distance of bodies from the reference ************************
    if strcmp(show_plot,'yes')
        SolarSystem_Plot(T,Y);
    end
    
    %% ANIMATION: Animation of the Solar System ***************************
    if strcmp(show_animation,'yes')
        Solarsystem_Animation(T,Y);
    end
    
    %% FUNCTION: SolarSystem_ODE ******************************************
    function [output] = SolarSystem_ODE(t,y)
        output = zeros(nv*nb,1); % Initialisation
        
        %% Distance between bodies.
        r = zeros(nb,nb); % Initialisation
        for i = interacting_bodies
            for j = interacting_bodies
                r(i,j) = sqrt(sum((y(ip_xyz+(i-1)*nv)-y(ip_xyz+(j-1)*nv)).^2)); % Distance between body i and j
            end
        end
        
        %% Body's motion
        for ib = interacting_bodies
            output(ip_xyz+(ib-1)*nv) = y(iv_xyz+(ib-1)*nv); % v_x, v_y, v_z 
            iib = setdiff(interacting_bodies,ib); % Indices of the interacting bodies wrt ib.
            output(iv_xyz(1)+(ib-1)*nv) = G*sum(Mb(iib).*(y(ip_xyz(1)+(iib-1)*nv)-y(ip_xyz(1)+(ib-1)*nv))'./(r(ib,iib).^3));
            output(iv_xyz(2)+(ib-1)*nv) = G*sum(Mb(iib).*(y(ip_xyz(2)+(iib-1)*nv)-y(ip_xyz(2)+(ib-1)*nv))'./(r(ib,iib).^3));
            output(iv_xyz(3)+(ib-1)*nv) = G*sum(Mb(iib).*(y(ip_xyz(3)+(iib-1)*nv)-y(ip_xyz(3)+(ib-1)*nv))'./(r(ib,iib).^3)); 
        end
        
        Solver_Progress(t)
    end
    
    %% FUNCTION: SolarSystem_Plot *****************************************
    function [] = SolarSystem_Plot(T,Y)
        T = T./days_in_a_yr; 
        
        fig_Sun = figure('name','The Sun');
        for ib = is
            if strcmp(sun_fixed,'yes')
                plot(T,sqrt(sum((Y(:,ip_xyz+(ib-1)*nv)-Y(:,ip_xyz+(is-1)*nv)).^2,2)),'Color',colors_b{ib});hold on;
            else
                plot(T,sqrt(sum((Y(:,ip_xyz+(ib-1)*nv)).^2,2)),'Color',colors_b{ib});hold on;
            end
        end
        xlabel('Time (in Earth years)');
        ylabel('The distance of the Sun to the reference (in m)');
        legend(names_b{is},'Location','northeast');
        xlim([0 time_of_sim_in_yrs]); grid on;
        
        fig_IP = figure('name','Inner Planets');
        for ib = inner_planets
            if strcmp(sun_fixed,'yes')
                plot(T,sqrt(sum((Y(:,ip_xyz+(ib-1)*nv)-Y(:,ip_xyz+(is-1)*nv)).^2,2)),'Color',colors_b{ib});hold on;
            else
                plot(T,sqrt(sum((Y(:,ip_xyz+(ib-1)*nv)).^2,2)),'Color',colors_b{ib});hold on;
            end
        end
        xlabel('Time (in Earth years)');
        ylabel('The distance of the planet to the reference (in m)');
        legend(names_b{inner_planets},'Location','northeast');
        xlim([0 time_of_sim_in_yrs]); grid on;

        fig_OP = figure('name','Outer planets');
        for ib = outer_planets
            if strcmp(sun_fixed,'yes')
                plot(T,sqrt(sum((Y(:,ip_xyz+(ib-1)*nv)-Y(:,ip_xyz+(is-1)*nv)).^2,2)),'Color',colors_b{ib});hold on;
            else
                plot(T,sqrt(sum((Y(:,ip_xyz+(ib-1)*nv)).^2,2)),'Color',colors_b{ib});hold on;
            end
        end
        xlabel('Time (in Earth years)');
        ylabel('The distsnce of the planet to the reference (in m)');
        legend(names_b{outer_planets},'Location','northeast');
        xlim([0 time_of_sim_in_yrs]); grid on;
    end

    %% FUNCTION: SolarSystem_Animation ************************************
    function [] = Solarsystem_Animation(T,Y)
        N = length(T);
        
        if  strcmp(generate_video,'yes') 
            name = strcat(datestr(now,'yyyy.mm.dd.HHMMSS'),...
                  '-',ode_solver,...
                  '.avi');
            v = VideoWriter(name);
            open(v);
        end
         
        scrsz = get(0,'ScreenSize');
        fig_ANIM = figure('position', [0.05*scrsz(3) 0.05*scrsz(4) 0.75*scrsz(3) 0.85*scrsz(4)]);
        
        set(gcf,'Name','Solar System');

        % Plot the overall orbits of the planets
        for ib = animation_bodies
            if strcmp(sun_fixed,'yes')
                plots_b(ib) = plot3(Y(:,ip_xyz(1)+(ib-1)*nv)-Y(:,ip_xyz(1)+(is-1)*nv),...
                                    Y(:,ip_xyz(2)+(ib-1)*nv)-Y(:,ip_xyz(2)+(is-1)*nv),...
                                    Y(:,ip_xyz(3)+(ib-1)*nv)-Y(:,ip_xyz(3)+(is-1)*nv),...
                                    'linewidth',0.5,'Color',colors_b{ib}...
                                );hold on;
            else
                plots_b(ib) = plot3(Y(:,ip_xyz(1)+(ib-1)*nv),...
                                    Y(:,ip_xyz(2)+(ib-1)*nv),...
                                    Y(:,ip_xyz(3)+(ib-1)*nv),...
                                    'linewidth',0.5,'Color',colors_b{ib}...
                                );hold on;
            end
        end
        legend(names_b{animation_bodies},'Location','northeast');
        
        if  strcmp(grid_on,'yes') 
            grid on; 
        end
        if  strcmp(view_from_top,'yes') 
            axis equal; 
        end 
        
        xlabel('x (in m)');
        ylabel('y (in m)');
        zlabel('z (in m)');
        
        % Loop over all the orbits
        for i = 1:animation_speed:N
            for ib = animation_bodies
                if strcmp(sun_fixed,'yes')
%                     plots_b(ib) = plot3(Y(1:i,ip_xyz(1)+(ib-1)*nv)-Y(:,ip_xyz(1)+(is-1)*nv),...
%                                         Y(1:i,ip_xyz(2)+(ib-1)*nv)-Y(:,ip_xyz(2)+(is-1)*nv),...
%                                         Y(1:i,ip_xyz(3)+(ib-1)*nv)-Y(:,ip_xyz(3)+(is-1)*nv),...
%                                         'linewidth',0.5,'Color',colors_b{ib}...
%                                     );hold on;
                    plots_b(ib) = plot3(Y(i,ip_xyz(1)+(ib-1)*nv)-Y(i,ip_xyz(1)+(is-1)*nv),...
                                        Y(i,ip_xyz(2)+(ib-1)*nv)-Y(i,ip_xyz(2)+(is-1)*nv),...
                                        Y(i,ip_xyz(3)+(ib-1)*nv)-Y(i,ip_xyz(3)+(is-1)*nv),...
                                        '.','markersize',scales_b(ib),'Color',colors_b{ib}...
                                      ); hold on;
                else
%                     plots_b(ib) = plot3(Y(1:i,ip_xyz(1)+(ib-1)*nv),...
%                                         Y(1:i,ip_xyz(2)+(ib-1)*nv),...
%                                         Y(1:i,ip_xyz(3)+(ib-1)*nv),...
%                                         'linewidth',0.5,'Color',colors_b{ib}...
%                                     );hold on;
                    plots_b(ib) = plot3(Y(i,ip_xyz(1)+(ib-1)*nv),...
                                        Y(i,ip_xyz(2)+(ib-1)*nv),...
                                        Y(i,ip_xyz(3)+(ib-1)*nv),...
                                        '.','markersize',scales_b(ib),'Color',colors_b{ib}...
                                      ); hold on;
                end
            end
            frame = getframe(gcf);
            
            if  strcmp(generate_video,'yes') 
                writeVideo(v,frame);
            end
            
            if i < N
                delete(plots_b);
            end
        end
        
        if  strcmp(generate_video,'yes') 
            close(v);    
        end
                    
    end
       
    %% FUNCTION: ode652 ***************************************************
    function [T,Y] = ode652(ODE,TSPAN,Y_INIT)
        % Definition of the RK-Method parameters
        A_matrix = [0 0 0 0; .5 0 0 0; 0 .5 0 0; 0 0 1 0];
        b_weights = [1/6 1/3 1/3 1/6];
        c_nodes = [0 .5 .5 1];
        
        nu = length(c_nodes); %Det of RK-Method's # of stages
        
        STEP = TSPAN(2)-TSPAN(1); % Def of the step size
        T = TSPAN'; % Def of the time span.
        N = length(T); %Det of the # of steps.
        M = length(Y_INIT); % The number of values for a given step.
        Y = zeros(M,N); % Initialisation of the vector of Y.
        Y(:,1) = Y_INIT; % Def of the Initial value.

        %% RK-Method
        Xi = zeros(M,nu);%Initialisation of the Xi's.
        for n = 1:N-1

            So = zeros(M,1); %Init of the outter sum for each Y.
            for j = 1:nu

                Si = zeros(M,1); %init of the inner sum for each Xi.
                for i = 1:j-1
                    Si = Si + A_matrix(j,i)*ODE(T(n)+c_nodes(i)*STEP,Xi(:,i)); %Det of the inner sum.
                end;

                Xi(:,j) = Y(:,n)+STEP*Si; %Determination of Xi's

                So = So + b_weights(j)*ODE(T(n)+c_nodes(j)*STEP,Xi(:,j)); %Det of the outter sum.

            end;

            Y(:,n+1) = Y(:,n)+STEP*So; %Det of the Approximation of the next y by RK-method
        end;

        Y = Y';
    end
    
    %% FUNCTION: Import_Data **********************************************
    function [T, Y] = Import_Data(import_file_name)
        % Import
        data = csvread(import_file_name);
        
        % Set variables
        T = data(:,1)';
        Y = data(:,2:end);
    end
    
    %% FUNCTION: Export_Data **********************************************
    function [] = Export_Data(T, Y)
        % Build file name
        name = strcat(datestr(now,'yyyy.mm.dd.HHMMSS'),...
                      '-',num2str(time_of_sim_in_yrs),'yrs',...
                      '-by', num2str(h),...
                      '-',ode_solver,...
                      '.dat');
        % Build data
        data = horzcat(T,Y);
        
        % Export data
        csvwrite(name,data);
    end
    
    %% FUNCTION: Solver_Progress ******************************************
    function [] = Solver_Progress(t)
        progress = (t/(time_of_sim_in_yrs*days_in_a_yr))*100;
        progress = sprintf('%.2f',progress);
        msg = ['Solar System Simulation for ',...
                num2str(time_of_sim_in_yrs),' yrs -- ',...
                progress,'%'];
        
%         fprintf(repmat('\b',1,line_size)); 
%         fprintf(msg);
        disp(msg);
    end
    
end