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