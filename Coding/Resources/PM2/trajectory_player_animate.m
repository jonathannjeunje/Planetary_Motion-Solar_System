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
hold on
set(gcf,'Color',[1 1 1]*0.2)
set(gca,'Color',[1 1 1]*0.2)


i = 1:10:zzzz;
%i = 1:10:min(43200,zzzz);
inc = 3:3:(count);

a = plot3(planet_out(i,(inc-2)),planet_out(i,(inc-1)), planet_out(i,inc));
axis auto
view(3)
axis manual
grid on
axis equal


for l = 1:100:zzzz%min(43200,zzzz)
    k = 1:100:l;
    for j = 1:length(a);
        inc = j*3;
        set(a(j),'XData',planet_out(k,(inc-2)),'YData',planet_out(k,(inc-1)),'ZData',planet_out(k,inc));
    end
    drawnow;
    
end
