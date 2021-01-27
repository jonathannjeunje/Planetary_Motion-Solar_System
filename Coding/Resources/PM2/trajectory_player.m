%this code is to play back the solar system simulated data
% Stephen Walker 2009 (stephen.walker@student.uts.edu.au)
clear

fid=fopen('test.bin','r');
[fname,permission,machineformat,encoding]=fopen(fid)
% open the binary file for reading
 planet_out_2=reshape(fread(fid,'*double'),75,[])';
 fclose(fid);
%read planet values out of file to array
whos planet_out*
%  return

planet_out = planet_out_2;
%unnecesscary, but it allowed two different pieces of code to be easily joinedspliced together
[zzzz,count] = size(planet_out)
%gets dimensions of planet data array
 figure(3)
  clf
 i = 1:1:zzzz/3;
 %define the range of data being displayed
hold on
colours=lines(zzzz);
%setup unique colour for each modeled body
j=0;
try
for inc = 3:3:(count);
    j=j+1;

 plot3 (planet_out(i,(inc-2)),planet_out(i,(inc-1)), planet_out(i,inc),'-','color',colours(j,:));
 % a 3d plot of the planetary data - it is hard to tell if it is 3d without
 % zooming in , as planets are all in very similar planes. Pluto & Halleys
 % comet give it away
end
%catch me
    % do nothing...
end
axis equal

 % axis([-4.3948    5.0089   -3.5114    3.9054]*1e11)
%pause(10)
%end