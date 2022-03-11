
time_step = 0.000005;        % Sim resolution
pos_step = 0.01;            

total_pos_steps = 100;      % Total number of position steps
simulation_time = 25;       % Run the animation for this long (sec)

alpha = 0.2;                % Angle of incline (rad)
viscocity = 320;            % Viscocity nu (m^2/s), value taken from notes 11.8
g = 9.8;                    % Gravitational acc (m/s^2)

% Useful quantities:
height = pos_step * total_pos_steps;
beta = (viscocity * time_step) / (pos_step^2);

f0 = ones(total_pos_steps,1);                   % Initial solution
u_SS = ones(total_pos_steps,1);                 % Steady State solution

x = linspace(0.0, height, total_pos_steps );    % Plot x-axis

for i = 1:total_pos_steps

    f0(i,1) = 0.02 * sin( 0.004*(i-1)^2 );                                 % Some strange initial condition
    u_SS(i,1) = - g * sin(alpha) * x(i) * ( 0.5*x(i) - 1.0) / viscocity;    % Known result from class
end

% Construct the A matrix
A = (1 + 2*beta) * diag(ones(1,total_pos_steps));
A = A - beta * diag(ones(1,total_pos_steps - 1), 1);
A = A - beta * diag(ones(1,total_pos_steps - 1), -1);

% No slip boundary condition on f
A(1,1) = 1.0;
for i = 2:total_pos_steps
    
    A(1,i) = 0.0;
end

% No stress boundary condition on f
A(total_pos_steps, total_pos_steps) = 1.0 + beta;
A(total_pos_steps, total_pos_steps - 1) = -beta;
for i = 1:(total_pos_steps - 2)
    
    A(total_pos_steps, i) = 0.0;
end

A_inv = inv(A);
f_new = f0;
time = 0.0;

while time < simulation_time % main loop
    
    u = f_new + u_SS; % How velovity is obtained
    
    cla;
    plot(x, u(:,1) );
    hold on;
    plot(x, u_SS(:,1) ); % show Steady State plot

    % Make the plot look pretty:
    grid on;
    ylim([-0.01 0.03]);
    xlabel('x (m)');
    ylabel('velocity u(x) (m/s)');
    legend({'Time dep. solution','Steady State'},'Location','southeast');

    pause(0.3); 

    f_new = A_inv * f_new; % Update solution
    time = time + 0.3;     % Update time
end
