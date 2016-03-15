#! /usr/bin/octave -q
#La línea de arriba señala que es un script que se le pasará a este programa. Recuerda que el archivo tiene que ser ejecutable. 
#La opción -q es para evitar el mensajillo de saludo al principio.
## Copyright (C) 2016 Gonzalo Rodríguez Prieto <gonzalo.rprieto@uclm.es>
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.
###########################################################################################
#
#
#  OCTAVE SCRIPT TO OBTAIN THE SHOCK WAVE VELOCITY FROM A SHADOW STREAK IMAGE.
#  IN TWO POSSIBILITIES: BY DATA SMOOTHING AND BY POLYNOMIA ADJUSTMENT.
#    Made by Gonzalo Rodríguez Prieto
#              Version 0.85
#
#
#########################################################
#
#   It uses the functions:
#        polyadj
#        deri
#	 	 display_rounded_matrix
#  They must be in the same directory.
#
###########################################################################################


#Charge the fitting package:
pkg load optim
#Charge the smooth package:
pkg load data-smoothing

tic; #Control total computing time.

#Empty memory from user created variables:
clear;

###
# Parameters:
###
mm = 1e3; %meters
us = 1e6; %seconds


###########################################################################################
#### RADIAL PART: FINDING THE RADIAL EXPANSION, AND FIND VELOCITY AND ACCELERATION.
###########################################################################################



###
# Reading the file with the shock wave radial data:
###
###
# Data extracted by Mario Barbaglia, and stored in *.txt files with this format:
# radius(mm) TAB time(µs)
# Its name indicates the experimental conditions.
###
arg_list = argv (); %Aquí está el comando que se usó, arg_list{1} y todo lo que le pongas detrás.

printf("Argumentos: \n");
for i = 1:nargin %nargin es el número de argumentos, incluidos el comando de entrada.
  printf (" %s", arg_list{i});
endfor
printf ("\n");

#String with the file name:
filename = "500µm-15kV-All.txt"; %Testing purposes.
%filename = arg_list{1}

[file, msg] = fopen(filename, "r");
if (file == -1) 
   error ("shadow script: Unable to open file name: %s, %s",filename, msg); 
endif;

r_rows =  "%f %f"; %2 rows to read, but three things, including the TAB in between columns.
if feof(file)==0 #Read the file until it finds the EndOfFile character.
   data = textscan(file, r_rows);
endif;
fclose(file);

#Transforming the data in numbers:
t = cell2mat(data(:,1)); #Time vector (in µs)
  t = t(isfinite(t)); %Removing Nans and other garbage from the row.
  t = sort(t); %Rearranging t from 0 to maximum for some cases...
r = cell2mat(data(:,2)); #Radius vector (in mm)
  r = r(isfinite(r)); %Removing Nans and other garbage from the row.
  r = sort(r);%Rearranging r from 0 to maximum for some cases...



###
# Smoothing the data and extracting from there velocity and acceleration:
###

trad = linspace(t(1),t(end),100)'; %Equispaciated vector. Necessary for smoothing
dev = std(r); %Standard deviation. Best for not over smoothing.
#Radial data smoothing:
r_smooth = regdatasmooth(t,r, "xhat", trad, "stdev", dev, "relative" ); % "xhat": points for x values; "stdev": value of admited standard deviation;
% "relative":to use relative differences;.

#Smooth velocity (In mm/µs)
paso = abs(trad(1)-trad(2)); %For numerical derivative, length between points
v_smooth = deri(r_smooth,paso); %Numerical derivative  of shock wave expansion.
  v_smooth(v_smooth==0) = []; %Removing garbage zeroes.

#Smooth acceleration (In mm/µs²)
acc_smooth = deri(v_smooth',paso);  %Numerical derivative  of shock wave expansion velocity.
 acc_smooth(acc_smooth==0) = []; %Removing garbage zeroes.
 acc_smooth = acc_smooth'; %Make it a column.
 v_smooth = v_smooth'; %Make it a column.



###
# r² fitting to autosimilar model:
###
#Making the r² data and calculating errors:
[r_sq_pol,r_sq_str] = polyfit(t,r.^2,1);
r_sq_data = polyval(r_sq_pol,t);
	R = r_sq_str.R; %The "R", whatever this it...
	d = ( R' * R)\eye(2); %The covariance matrix
	d = diag(d)'; %Vector of the diagonal elements of something...
	MSE = (r_sq_str.normr^2)/ r_sq_str.df; %Variance of the residuals.
	r_sq_err = sqrt(MSE*d); %Standard errors of the polynomial values
%Velocity: (derivated from r = sqrt(linear appr. r²)) (in m/µs)
vel_r_sq = ( 0.5*r_sq_pol(1) ) ./ sqrt(polyval(r_sq_pol,t)  );
%Acceleration (Derived from r = sqrt(linear approx. r²) (in mm/µs²):
accel_r_sq = - (0.25*r_sq_pol(1)^2 ./ polyval(r_sq_pol,t).^1.5 );



###
# Polynomial fitting(It does not work very well):
###
#Polynomial fit to the data:
[r_pol, r_struc] = polyadj(t,r,"order",3,"max_order",5);
r_data = polyval(r_pol,t); %Radius data.

#Polynomial derivative of radial data, velocity (In mm/µs):
v_pol = polyder(r_pol);
v_data = polyval(v_pol,t); %Velocity data.

#Polynomial acceleration: (In mm/µs):
acc_pol = polyder(v_pol);
acc_data = polyval(acc_pol,t); %Acceleration data.



###
# Cheking graphs:
###
figure("visible","off");
#Smoothing and experimental data:
 graphics_toolkit ("gnuplot"); %To save LATEX symbols properly.
plot(t,r,"*k",trad,r_smooth,"-g"); %Show the radial shock wave expansion and the smoothing.
title(filename,"interpreter","tex"); %Graph title: Filename of data.
xlabel('t(\mus)',"interpreter","tex");%Graph labels
ylabel("r(mm)");
graphname = horzcat(filename(1:index(filename,".","first")),"rad_smooth.jpg"); %Graph filename
print(graphname); %Save a file with the graph
close; %Close the graph window
#Radial acceleration:
figure("visible","off");
plot(trad(1:size(acc_smooth)),acc_smooth,"-r", t, acc_data, "-b"); %Show the radial shock wave acceleration smoothing and polynomia.
title(filename); %Graph title: Filename of data.
xlabel('t(\mus)',"interpreter","tex");%Graph labels
ylabel('accel (mm/\mus^2)',"interpreter","tex");
legend('Smoothing', 'Pol. interpolation');
graphname = horzcat(filename(1:index(filename,".","first")),"accel_smooth.jpg"); %Graph filename
print(graphname); %Save a file with the graph
close; %Close the graph window

graphics_toolkit ("fltk"); %To come back to normal.



###
# Saving the results:
###
#Saving acceleration:
#Output file name:
name = horzcat(filename(1:index(filename,".","first")),"accel.dat"); %Adding the right sufix to the filename variable.
output = fopen(name,"w"); %Opening the file.
#First line (Veusz system, acepted as garbage line by QtiPlot):
fdisp(output,"descriptor `t(µs)`	`acc_smooth(mm/µs²)`	`t_pol(µs)`	`acc_pol(mm/µs²)`");
redond = [3 3 3 3]; %Saved precision.
#Matching vectors size:
t2 = t;
t2(end:size(acc_smooth)) = 0.0;
acc_data(end:size(acc_smooth)) = 0.0; 
#Matrix to save:
acceleration = [trad(1:size(acc_smooth)), acc_smooth, t2, acc_data]; %Puting the vectors in columns.
display_rounded_matrix(acceleration, redond, output); 
fclose(output); %Closing the file.

#Saving velocity:
#Output file name:
name = horzcat(filename(1:index(filename,".","first")),"vel.dat"); %Adding the right sufix to the filename variable.
output = fopen(name,"w"); %Opening the file.
#First line (Veusz system, acepted as garbage line by QtiPlot):
fdisp(output,"descriptor `t(µs)`	`vel_smooth(mm/µs)`	`t_pol(µs)`	`vel_pol(mm/µs)`");
redond = [3 3 3 3]; %Saved precision  for all the files!!!
#Matching vectors size:
t3 = t;
t3(end:size(v_smooth)) = 0.0;
v_data(end:size(v_smooth)) = 0.0;
velocity = [trad(1:size(v_smooth)), v_smooth, t3, v_data]; %Puting the vectors in columns.
display_rounded_matrix(velocity, redond, output); 
fclose(output); %Closing the file.


#Saving r² data:
#Output file name:
name = horzcat(filename(1:index(filename,".","first")),"r_2.dat"); %Adding the right sufix to the shot name chosen from the filename variable.
output = fopen(name,"w"); %Opening the file.
#First line ((Veusz system, acepted as garbage line by QtiPlot):
fdisp(output,"descriptor `t_r²(µs)`  `radius²(mm²)`");
redond = [3 3];
rad2 = [t, r_sq_data]; %Puting the vectors in columns.
display_rounded_matrix(rad2, redond, output); 
fclose(output); %Closing the file.
#Saving r² polynomial values:
#Output file name:
name = horzcat(filename(1:index(filename,".","first")),"r_2_parameters.txt"); %Adding the right sufix to the shot name chosen from the filename variable.
output = fopen(name,"w"); %Opening the file.
fdisp(output,filename(1:index(filename,".","first")-1) ); %Experimental conditions.
fdisp(output,'');
fdisp(output, 'Pol. values (Term 1, term 0)');
fdisp(output, r_sq_pol);
fdisp(output, 'Err. values (Term 1, term 0)');
fdisp(output, r_sq_err);
Ener = ( ( r_sq_pol(1)* us)/mm^2 ) * 1.27; %Last parameter is air density (kg/m³) and mm and us are conversion factor to SI units. Joules.
t_0 = abs(r_sq_pol(2)/r_sq_pol(1)) / us;  %Initial time for self-similar behaviour. Seconds
En = [Ener, t_0]; %To store the data.
fdisp(output,'Energy and t_0 values (J and sec.)');
fdisp(output, En);

fclose(output);



###
# Total processing time
###

timing = toc; 

disp("Script shadow execution time:")
disp(timing)
disp(" seconds")  


#That's...that's all folks!!!
