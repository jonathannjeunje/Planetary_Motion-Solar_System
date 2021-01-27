
%Math 652 HW6 RK-Method Implementation due Friday
%% 
clear; %Clear memory.
close all; %Close previous figures.
%clc, clf;

%% Variable declaration
a = 0; b = 1; %Definition of the interval [a,b]
f = @(t,y)(-y); %Definition of the function f(t,y) = y'(t) = -y
m = 4; %Def of the max(k), k = 0,1,2,3,...,m

%Definition of the RK-Method parameters
A_matrix = [0 0 0 0; .5 0 0 0; 0 .5 0 0; 0 0 1 0];
b_weights = [1/6 1/3 1/3 1/6];
c_nodes = [0 .5 .5 1];

nu = length(c_nodes); %Det of RK-Method's # of stages

e_max = zeros(1,m); %Initialisation of the vector of error.
e_ratio = zeros(1,m-1); %Initialisation of the vector of error ratios

%% Numerical approximation for various k's
for k = 0:m

    h =.1/2^k; %Def of the Step-size.
    t = a:h:b; %Det of the vector of steps.
    N = length(t); %Det of the # of steps.
    y = zeros(1,N); %Initialisation of the vector of y's.
    y(1) = exp(-t(1)); %Def of the Initial value.
   
    %% RK-Method
    xi = zeros(1,nu);%Initialisation of the Xi's.
    for n = 1:N-1
        
        So = 0; %Init of the outter sum for each y.
        for j = 1:nu
            
            Si = 0; %init of the inner sum for each Xi.
            for i = 1:j-1
                Si = Si + A_matrix(j,i)*f(t(n)+c_nodes(i)*h,xi(i)); %Det of the inner sum.
            end;
            
            xi(j) = y(n)+h*Si; %Determination of Xi's
            
            So = So + b_weights(j)*f(t(n)+c_nodes(j)*h,xi(j)); %Det of the outter sum.
            
        end;
        
        y(n+1) = y(n)+h*So; %Det of the Approximation of the next y by RK-method
    end;
    
    %% Det of the max error.
    e_max(k+1) = norm(y-exp(-t),inf);
      
    %% Plots
    if mod((k+1),4)==1
        figure('units','normalized','outerposition',[0 0 1 1]);
        %figure('rend','painters','pos',[250 50 900 600]);
        subplot(2,2,1);
    elseif mod((k+1),4)==0
        subplot(2,2,4);
    else
        subplot(2,2,mod((k+1),4));
    end;
    plot(t,y,'*r'); hold on;
    plot(t,exp(-t),'b');
    legend('Approxition','Exact','Location','northeast');
    ylabel('y');xlabel('t');grid on;
    title(strcat('Approximation vs exact solution with k=', num2str(k),', h=', num2str(h)));
    
end;

%% Claim on the order of convergence: The order of convergence of our 
%implementation of RK-method is 4.

%Let's verify our claim: Since we ran our method with the step size
%h=.1/2^k
%and recorded the errors. We need to compute the ratio [e(k)/e(k+1)], verify it
%is ~(2^4) i.e. ~16 and conclude.

for k=1:m-1
    e_ratio(k)=e_max(k)/e_max(k+1); %Det of the error ratios.
end;

%% Plots
figure('units','normalized','outerposition',[0 0 1 1]); %Instantiating another figure window.
subplot(1,2,1); %Def our first subplot.
bar(0:m,e_max);ylabel('e\_max');xlabel('k');grid on; %Error plot.
text(0:m,e_max,num2str(e_max'),'vert','bottom','horiz','center');
title('Errors bar chart');
subplot(1,2,2); %Def our 2nd subplot.
plot(1:m-1,e_ratio,'*--'); ylim([min(e_ratio)-3 max(e_ratio)+3]);;ylabel('e\_ratio');xlabel('k'); grid on; %Error ratio plot.
%text(1:m-1,e_ratio,strcat('\leftarrow ', num2str(e_ratio')));
text(1:m-1,e_ratio,num2str(e_ratio'),'vert','bottom','horiz','left');
title('Error ratios plot');

%Conclusion: Hence, by interpreting the error ratio plot, the order of
%convergence of our RK-method is 4.
