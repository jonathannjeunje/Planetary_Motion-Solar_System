
a = 0; b = 10; N = 100; h = (b-a)/N;
g = 9.8; % m/sec
l = 1; % m
%F = @(t,y)[y(2);(-g/l)*sin(y(1))];
[T,Y] = ode45(@pendulumMotion,a:h:b,[pi/4 0]);
plot(T,Y(:,1),'k');


for i=1:N
    figure(1);
    x=(l*cos(Y(i,1)-pi/2));
    y=(l*sin(Y(i,1)-pi/2));
    plot(x,y,'.','markersize',40); 
    axis([-2 2 -1.5 0]);
    line([0 x], [0 y],'Linewidth',2);
    Yr=getframe;
end

movie(Yr)

