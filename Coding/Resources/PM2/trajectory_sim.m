%this function is to calculate the motion of the planets based on their
%initial velocity, as well as the gravitational influence of other planets
%solar_system.mat includes all of the planets, moons of earth, jupiter &
%saturn, as well as halleys comet, Voyager 1 and 2.
% Stephen Walker 2009 (stephen.walker@student.uts.edu.au)
clear;
load solar_system.mat;
%loads initial positions of all of the bodies being modeled 
start_time = clock;
% to determine the runtime of this code - for this version, 1100 days
% at res 120 delta_t 60 takes 2hrs 57min on 2.53 core 2 duo
G = 6.67e-11;
%gravitational constant
day_count =  365 * 10;

res = 120;
delta_t = 60;
%delta_t is the resolution of step size between calculations.
%res is how frequently this data is written to the binary file. 
%res shouldnt be smaller than delta_t - waste of resources for no benefit


scenario_runtime = day_count * 86400 /delta_t;
[total_planets,zzzzz] = size(planets_au);
%scenario executions determined by the delta_t and the number of days of
%runtime requested
start_time = clock;
AU=149597870691;
KM = 1000;
for p =1:total_planets;
    planets(p,1) = planets_au(p,1)* KM;
    planets(p,2) = planets_au(p,2)* KM;
    planets(p,3) = planets_au(p,3)* KM;

    planets(p,4) = planets_au(p,4)* KM;
    planets(p,5) = planets_au(p,5)* KM;
    planets(p,6) = planets_au(p,6)* KM;
    planets(p,7) = planets_au(p,7); % mass import
    


%planets will soon include flags to indicate body type and current status
% eg, powered flight, within atmosphere,at boundary conditions, eceeded
% boundary conditions at high velocity
% for now, this is the main class. setting the import ephemerides will
% define missions etc. this class will call functions for drag, launch
% guidance etc as required. 
% eventually, this will be a function of an
% independant main class which will implement missions calling this and
% other modules

 end
%this is the import function for the ephemerides. The KM & AU values are
%used to format the imported ephemerides if needed.


fid_b = fopen('test.bin','w')

%open the binary file to record the calculated planet positions.
%this was more efficient than growing an array within Matlab - it tended to
%slow down on longer simulations
for t = 1:scenario_runtime;

    %LOOP A start
    %LOOP A selects all elements to be modeled - one at a time
    for a = 1 : total_planets ;
        
        %delta_t = planets(a,11);
        %this is future work to try and make the code more scaleable

        pos_a = [planets(a,1),planets(a,2),planets(a,3)];
        vel_a = [planets(a,4),planets(a,5),planets(a,6)];
        mass_a = planets(a,7);
   %initial positions , mass and velocity of planet/body a      
 net_grav_a = [0,0,0];
 
        % LOOP B start
        %Loop B compares all of the planets/ bodies (except a)with planet/body a
        %the stepping is wierd (a+1) as it is only necessary to calculate
        %the grav vector from a to b. once calculated, it is easy to
        %reverse the calculated vector (multiply by -1) to get the b to a
        %vector. This halved the work of B loop and halved the execution
        %time

        for b = (a+1):total_planets; %the a+1 is in to create diagonal


            mass_b = planets(b,7);
            pos_b = [planets(b,1),planets(b,2),planets(b,3)];
            %reads initial positions and mass of body b

            diff =  pos_b -pos_a; 
            %with both pos values, a vector between a and b is determined
            pos_diff.mag = sqrt((diff(1)^2) + (diff(2)^2) + (diff(3)^2));
            %the magnitude of this vector is calculated
            
            unit_diff = diff * (1 / pos_diff.mag);
            %The vector between a and b is turned into a unit vector
            pos_diff.grav = ((mass_a * mass_b * G) / (pos_diff.mag^2));
            %the gravity magnitude is calculated 
            grav = unit_diff * pos_diff.grav;
            % the gravity magnitude is applied to the unit vector giving a
            % gravity vector between a and b


            grav_array.x(a,b) = grav(1)* -1;
            grav_array.y(a,b) = grav(2)* -1;
            grav_array.z(a,b) = grav(3)* -1;
            %the gravity vector values are loaded into an array (virtually
            %3d, one array per axis - this only does HALF of the array
            grav_array.x(b,a) = grav(1);
            grav_array.y(b,a) = grav(2);
            grav_array.z(b,a) = grav(3);
            % this does the other HALF of the array

         


        end;%loop_b_end


        net_grav_a(1) = sum(grav_array.x(:,a));
        net_grav_a(2) = sum(grav_array.y(:,a));
        net_grav_a(3) = sum(grav_array.z(:,a));
        %the total gravitational force on body a 

        mt = delta_t/mass_a;
        
        net_grav_a = net_grav_a * mt;
        %this turns the net grav value into an acceleration value, by removing mass, 
        %and factoring in the delta_t
        new_vel_a = vel_a + net_grav_a;
        %the new velocity is the old velocity added to the delta_velocity
        new_pos_a = pos_a + (new_vel_a * delta_t);
        %new position is the old position plus the delta position.
        planets(a,1:3) = new_pos_a;
        
        planets(a,4:6) = new_vel_a;
        %The calculated position and velocity are written to  planets 


    end
    %LOOP A END

    if rem(t,res) == 0;
        % every "res" seconds the position values for all planets are
        % sampled 
        disp((t /scenario_runtime)*100);
        for a = 1:total_planets
            planet_out(((3*a) - 2)) = planets(a,1); %x position values
            planet_out(((3*a) - 1)) = planets(a,2); %y position values
            planet_out((3*a)) = planets(a,3);       %z position values
%         %the sampled values arewritten to planet_out.
        end
       fwrite(fid_b,planet_out(:),'*double');  
       %planet out is written to a binary file
       

    end
  
end %END OF IMPLEMENTATION LOOP
%time step is finished
fclose(fid_b);
% once all timesteps are done, binary file is closed


stop_time = clock;
run_time = stop_time - start_time
save run_time_duration.mat run_time
%the execution time is stored in a .mat file