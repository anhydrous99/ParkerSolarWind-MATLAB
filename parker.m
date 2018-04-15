k = 1.38064852*10^(-23);  % Boltzmann Constant (m^2kgs^-2K^-1)
T = 1.3*10^6;             % Typical Coronal Temprature (K)
M_p = 1.6726219*10^(-27); % Mass of a Proton (kg)
U = 2*k*T/M_p;            % Ion Thermal Velocity Squared
G = 6.67408*10^(-11);     % Gravitational Constant (m^3kg^-1s^-2)
M_sun = 1.989*10^30;      % Mass of sun (kg)
r_corona = 5360500000;    % radius from center of sun to
                          %   outer edge of corona (m)
r_scorona = 2100000;      % start of corona
AU = 149597870700;        % 1 Astronomical Unit (m)
v_corona = 400000;        % Typical velocity of solar wind
                          %   at the corona
a = sqrt(k*T/M_p);        % used With Maxwell Distribution of Ion Velocities
% Plot Maxwell Distribution where T is 1.3*10^6
fplot(@(v) sqrt(2./pi).*(v.^2/a.^3).*exp(-v.^2./(2.*a.^2)), [0 500000]);
line([4*10^5 4*10^5],[0 (5.8)*10^(-6)], 'Color','red');
title('Maxwell-Boltzmann Velocity Distribution of Ions in Corona');
xlabel('Velocity (m/s)');
ylabel('Count');

v_es_c = sqrt(2*G*M_sun/r_scorona);  % escape velocity at start of Corona
v_es_a = sqrt(2*G*M_sun/r_corona); % escape velocity at end of Corona
fprintf('Escape Velocity at Start of Corona: %8.2e m/s\n', v_es_c);
fprintf('Escape Velocity at End of Corona  : %8.2e m/s\n', v_es_a);

odefun = @(t,y) (2*U/t-G*M_sun/t^2)/(y-U/y);
[r,v] = ode45(odefun, [r_corona AU], v_corona);
figure(2);
plot(r, v, 'r-o');
title('Solution of Parker Diff. Eq. with ODE45');
xlabel('Distance (m)');
ylabel('Velocity (m/s)');
clear k T M_p U G M_sun r_corona AU v_corona odefun a