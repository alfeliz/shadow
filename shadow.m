#! /usr/bin/octave -q
#La línea de arriba señala que es un script que se le pasará a este programa. Recuerda que el archivo tiene que ser ejecutable. 
#La opción -q es para evitar el mensajillo de saludo al principio.


#Charge the fitting package:
pkg load optim
#Charge the smooth package:
pkg load data-smoothing



###########################################################################################
#
#
#  OCTAVE SCRIPT TO OBTAIN THE SHOCK WAVE VELOCITY FROM A SHADOW STREAK IMAGE.
#  IN TWO POSSIBILITIES: BY DATA SMOOTHING AND BY POLYNOMIA ADJUSTMENT.
#    Made by Gonzalo Rodríguez Prieto
#       (gonzalo#rprieto AT uclm#es)
#       (Mail: change "#" by "." and "AT" by "@"
#              Version 0.50
#
#
#########################################################
#
#   It uses the functions:
#        polyadj
#        deri
#	 	 display_rounded_matrix (THIS NOT NOW)
#  They must be in the same directory.
#
###########################################################################################

tic; #Control total computing time.

#Empty memory from user created variables:
clear;


###########################################################################################
#### RADIAL PART: FINDING THE RADIAL EXPANSION, AND FIND THE VELOCITY BY DERIVATION:
###########################################################################################

###
# Reading the file with the shock wave radial data:
###
###
# Data extracted by Mario Barbaglia, and stored in *.txt files with this format:
# radius(mm) TAB time(µs)
# Its name indicates the experimental conditions.
###


arg_list = argv (); %Aquí está el comando que se usó, arg_list{1} y todo o que le pongas detrás.

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

#Making a matrix with the radial data and sorting it:
rad = [t,r];
rad = sortrows(rad);


#Smoothing of radial data to velocity and acceleration:
trad = linspace(t(1),t(end),100)'; %Equispaciated vector. Necessary for smoothing
dev = max(abs(r))*0.01; %Use of 1% of maximum radial expasion as "standard deviation" of data.
#Radial data smoothing:
[r_smooth,lambda] = regdatasmooth(t,r, "xhat", trad, "stdev", dev);

#Smooth velocity (In mm/µs)
paso = abs(trad(1)-trad(2)); %For numerical derivative, length between points
v_smooth = deri(r_smooth,paso); %Numerical derivative  of shock wave expansion.
  v_smooth(v_smooth==0) = []; %Removing garbage zeroes.

#Smooth acceleration (In mm/µs²)
a_smooth = deri(v_smooth',paso);  %Numerical derivative  of shock wave expansion velocity.
 a_smooth(a_smooth==0) = []; %Removing garbage zeroes.

#Making the r² data:
r_sq_pol = polyfit(t,r.^2,1);
r_sq_data = polyval(r_sq_pol,t);
%Velocity: (derivated from r = sqrt(linear appr. r²)) (in m/µs)
vel_r_sq = ( 0.5*r_sq_pol(1) ) ./ sqrt(polyval(r_sq_pol,rad(:,1))  );
%Acceleration (Derived from r = sqrt(linear approx. r²) (in mm/µs²):
accel_r_sq = - (0.25*r_sq_pol(1)^2 ./ polyval(r_sq_pol,rad(:,1)).^1.5 );

#Polynomial fit to the data:
[r_pol, r_struc] = polyadj(rad(:,1),rad(:,2));
r_data = polyval(r_pol,rad(:,1)); %Radius data.

#Polynomial derivative of radial data, velocity (In mm/µs):
v_pol = polyder(r_pol);
v_data = polyval(v_pol,rad(:,1)); %Velocity data.

#Polynomial acceleration: (In mm/µs):
acc_pol = polyder(v_pol);
acc_data = polyval(acc_pol,rad(:,1)); %Acceleration data.

%###
%# Showing graph with radial smoothing and radius data:
%###
 %graphics_toolkit ("gnuplot"); %To save LATEX symbols properly.
%plot(rad(:,1),rad(:,2),"*k",trad,r_smooth,"-g"); %Show the radial shock wave expansion and the smoothing.
%title(filename,"interpreter","tex"); %Graph title: Filename of data.
%xlabel('t(\mus)',"interpreter","tex");%Graph labels
%ylabel("r(mm)");

%graphname = horzcat(filename(1:index(filename,".","first")),"rad_smooth.jpg"); %Graph filename
%print(graphname); %Save a file with the graph
%close; %Close the graph window


%###
%# Showing graph with radial acceleration:
%###
%plot(trad(1:size(a_smooth,2)),a_smooth,"-r"); %Show the radial shock wave expansion and the smoothing.
%title(filename); %Graph title: Filename of data.
%xlabel('t(\mus)',"interpreter","tex");%Graph labels
%ylabel('accel (mm/\mus^2)',"interpreter","tex");

%graphname = horzcat(filename(1:index(filename,".","first")),"accel_smooth.jpg"); %Graph filename
%print(graphname); %Save a file with the graph
%close; %Close the graph window

%graphics_toolkit ("fltk"); %To come back to normal.

###
# Saving the results:
###

%#Saving acceleration from smooth data:
%#Output file name:
%name = horzcat(filename(1:index(filename,".","first")),"accel.dat"); %Adding the right sufix to the filename variable.
%output = fopen(name,"w"); %Opening the file.
%#First line (Veusz system, acepted as garbage line by QtiPlot):
%fdisp(output,"descriptor `t(µs)`  `acceleration(mm/µs²)`");
%redond = [4 6]; %Saved precision  for all the files!!!
%acceleration = [trad(1:size(a_smooth,2))'; a_smooth]'; %Puting the vectors in columns. Notice the ' to make the columns right!!!
%display_rounded_matrix(acceleration, redond, output); 
%fclose(output); %Closing the file.

%#Saving velocity from smooth data:
%#Output file name:
%name = horzcat(filename(1:index(filename,".","first")),"vel.dat"); %Adding the right sufix to the filename variable.
%output = fopen(name,"w"); %Opening the file.
%#First line (Veusz system, acepted as garbage line by QtiPlot):
%fdisp(output,"descriptor `t(µs)`  `velocity(mm/µs)`");
%redond = [4 6]; %Saved precision  for all the files!!!
%velocity = [trad(1:size(v_smooth,2))'; v_smooth]'; %Puting the vectors in columns. Notice the ' to make the columns right!!!
%display_rounded_matrix(velocity, redond, output); 
%fclose(output); %Closing the file.

%#Saving radius:
%#Output file name:
%name = horzcat(filename(1:index(filename,".","first")),"rad_smooth.dat"); %Adding the right sufix to the shot name chosen from the filename variable.
%output = fopen(name,"w"); %Opening the file.
%#First line ((Veusz system, acepted as garbage line by QtiPlot):
%fdisp(output,"descriptor `t(µs)`  `radius(mm)`");
%radius = [trad'; r_smooth']'; %Puting the vectors in columns. Notice the ' to make the columns right!!!
%display_rounded_matrix(radius, redond, output); 
%fclose(output); %Closing the file.

%#Saving radius "al cuadrado":
%#Output file name:
%name = horzcat(filename(1:index(filename,".","first")),"rad_2.dat"); %Adding the right sufix to the shot name chosen from the filename variable.
%output = fopen(name,"w"); %Opening the file.
%#First line ((Veusz system, acepted as garbage line by QtiPlot):
%fdisp(output,"descriptor `t(µs)`  `radius²(mm²)`");
%rad2 = [trad'; r_sq']'; %Puting the vectors in columns. Notice the ' to make the columns right!!!
%display_rounded_matrix(rad2, redond, output); 
%fclose(output); %Closing the file.



###
# Total processing time
###

timing = toc; 

disp("Script shadow execution time:")
disp(timing)
disp(" seconds")  


#That's...that's all folks!!!
