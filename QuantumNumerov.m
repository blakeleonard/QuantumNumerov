% Calculates the Radial Wave Function for a given potential using the Numerov Method.

% Blake Leonard 2011

% Physics Department
% Washington University in St. Louis


clear;

% Variable Declaration and initial conditions

h = 0.001;                          % step size

ediv = 1000;                        % energy step size

mass = 2.402028*10^-24;             % Units of Kelvin / c^2 ( c in angstroms )  

A = 0.0823454302;                   % 2*mass / h_bar^2       1 / (K*ang^2)   
                     
step_outside = 9000;

maxstep = 10000;

estep_max = 1000;

theta_step_max = 1000;

theta_step = ( 2*pi ) / theta_step_max;

energy_step1 = 100;

energy_step2 = 500;


% Calculate V(r) & Plot

for istep = 1:maxstep

		r(istep) = 2 + h * (istep);

		V(istep) = 40.88*( (2.556/r(istep))^12 - (2.556/r(istep))^6);     % Kelvin
		
end

axis([2.5,10]);

xlabel("r(femtometers)");

ylabel("Potential(MeV)");

plot(r, V);

disp('');

z = input('Hit Enter to Advance');

disp('');


% Loop through Energy Values

for estep = 1:estep_max

	U_0(1) = 0;               U_1(1) = 0;
	U_0(2) = 0.0001;          U_1(2) = 0.0001;      % initial condition  

	E(estep) = (estep)*ediv;

	
	%  Calculate functions f_0(r) (L=0 case) and f_1(r) (L=1 case) 

	for istep = 1:maxstep
	
		f_0(istep) = A * ( E(estep) - V(istep) );
	
		f_1(istep) = f_0(istep) - ( 2 / (r(istep))^2);
	
	end


	% Numerov Method

	for istep =3:maxstep

		U_0(istep) = ( ( ( 2 - ( (5/6)*h^2 ) * f_0(istep - 1 ) ) * U_0(istep-1) ) - ( 1 + ((h^2)/12) * f_0(istep-2) ) * U_0(istep-2) ) / ( 1 + ((h^2)/12)*f_0(istep) );
	
		U_1(istep) = ( ( ( 2 - ( (5/6)*h^2 ) * f_1(istep - 1 ) ) * U_1(istep-1) ) - ( 1 + ((h^2)/12) * f_1(istep-2) ) * U_1(istep-2) ) / ( 1 + ((h^2)/12)*f_1(istep) );
		
		R_0(istep) = U_0(istep) / r(istep);
		
		R_1(istep) = U_1(istep) / r(istep);

	end
	
	if (estep == 1)
	
		axis([2.5,10]);

		xlabel("r");

		ylabel("Radial Wave function l = 0");

		plot(r, R_0);

		disp('');

		z = input('Hit Enter to Advance');

		disp('');
	
	
		axis([2.5,10]);

		xlabel("r");

		ylabel("Radial Wave function l = 1");

		plot(r, R_1);

		disp('');

		z = input('Hit Enter to Advance');

		disp('');
		
	end


	% Calculate Avg Derivatives at point outside potential range

	R_0prime = ( ( ( R_0(step_outside) - R_0(step_outside-1) ) / h ) + ( ( R_0(step_outside+1) - R_0(step_outside) ) / h ) ) / 2;     

	R_1prime = ( ( ( R_1(step_outside) - R_1(step_outside-1) ) / h ) + ( ( R_1(step_outside+1) - R_1(step_outside) ) / h ) ) / 2; 


	% Calculate beta

	beta_0 = ( ( step_outside * h ) / ( R_0(step_outside) ) ) * R_0prime;

	beta_1 = ( ( step_outside * h ) / ( R_1(step_outside) ) ) * R_1prime;


	% Calculate Bessel and Neumann functions

	k = sqrt(2*mass*E(estep));

	j_0 = ( sin( k*step_outside*h) )/ (k*step_outside*h);

	j_0prime = ( ( cos(k*step_outside*h) ) / (k*step_outside*h)  ) - ( ( sin(k*step_outside*h) ) / ( (k*step_outside*h)^2 ) );

	j_1 = ( ( sin(k*step_outside*h) ) / ( (k*step_outside*h)^2 )  ) - ( ( cos(k*step_outside*h) ) / (k*step_outside*h) );

	j_1prime = ( ( 2 * cos(k*step_outside*h) ) / ( (k*step_outside*h)^2 )  ) - ( ( 2 * sin(k*step_outside*h) ) / ( (k*step_outside*h)^3 ) ) + ( sin( k*step_outside*h) )/ (k*step_outside*h);


	n_0 = - ( cos( k*step_outside*h) )/ (k*step_outside*h);

	n_0prime = ( ( sin(k*step_outside*h) ) / (k*step_outside*h)  ) + ( ( cos(k*step_outside*h) ) / ( (k*step_outside*h)^2 ) );

	n_1 = - ( ( cos(k*step_outside*h) ) / ( (k*step_outside*h)^2 )  ) - ( ( sin(k*step_outside*h) ) / (k*step_outside*h) );

	n_1prime = ( ( 2 * sin(k*step_outside*h) ) / ( (k*step_outside*h)^2 )  ) + ( ( 2 * cos(k*step_outside*h) ) / ( (k*step_outside*h)^3 ) ) - ( cos( k*step_outside*h) )/ (k*step_outside*h);


	% Calculate delta

	delta_0(estep) = atan( ( ( (k*step_outside*h) * j_0prime ) - ( beta_0 * j_0 ) ) / ( ( (k*step_outside*h) * n_0prime ) - ( beta_0 * n_0 ) ) );

	delta_1(estep) = atan( ( ( (k*step_outside*h) * j_1prime ) - ( beta_1 * j_1 ) ) / ( ( (k*step_outside*h) * n_1prime ) - ( beta_1 * n_1 ) ) );


	% Calculate Total Cross Section
	
	sigma_tot(estep) = ( ( 4*pi ) / ( 2*mass*E(estep) ) ) * ( ( sin(delta_0(estep)) )^2 + 3 * ( sin(delta_1(estep)) )^2 );	
	
end


% Calculate Diff Cross Section for a couple diff energies

for istep = 1:theta_step_max

	theta(istep) = theta_step*(istep-1);

	sigma_diff1(istep) = ( ( 1 ) / ( 2*mass*E(energy_step1) ) ) * ( ( sin(delta_0(energy_step1)) )^2 + 9 * ( ( sin(delta_1(energy_step1)) )^2 ) * ( ( cos(theta(istep)) )^2 ) + 6 * sin(delta_0(energy_step1)) * sin(delta_1(energy_step1)) * cos(delta_0(energy_step1) - delta_1(energy_step1) ) * cos(theta(istep)) );

	sigma_diff2(istep) = ( ( 1 ) / ( 2*mass*E(energy_step2) ) ) * ( ( sin(delta_0(energy_step2)) )^2 + 9 * ( ( sin(delta_1(energy_step2)) )^2 ) * ( ( cos(theta(istep)) )^2 ) + 6 * sin(delta_0(energy_step2)) * sin(delta_1(energy_step2)) * cos(delta_0(energy_step2) - delta_1(energy_step2) ) * cos(theta(istep)) );

end


% Plot 

axis("auto");

xlabel("Energy(MeV)");

ylabel("Delta_0 Phase Shift");

plot(E, delta_0);

disp('');

z = input('Hit Enter to Advance');

disp('');


axis("auto");

xlabel("Energy(MeV)");

ylabel("Delta_1 Phase Shift");

plot(E, delta_1);

disp('');

z = input('Hit Enter to Advance');

disp('');


axis("auto");

xlabel("Energy(MeV)");

ylabel("Total Cross Section");

plot(E, sigma_tot);

disp('');

z = input('Hit Enter to Advance');

disp('');


axis("auto");

xlabel("Theta(radians)");

ylabel("Diff. Cross Section at Energy 1");

plot(theta, sigma_diff1);

disp('');

z = input('Hit Enter to Advance');

disp('');


axis("auto");

xlabel("Theta(radians)");

ylabel("Diff. Cross Section at Energy 2");

plot(theta, sigma_diff2);

disp('');

z = input('Hit Enter to Advance');

disp('');