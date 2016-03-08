# shadow

*shadow.m* is a Octave script that can be used from the command line (It has been tested in Linux).
	
Its task is to make numericaly the first and second derivative of some vector with two columns that have not equal separation on the horizontal axis. 

To do that, it approximates first the initial data to a large vector with equal separation, and there performs the first derivative and the derivative of the derivative, so the second.

It works with just the name of the file to be used and produced these files:

* Two JPG graphs of which one shows the original data and its approximation, and the other the "aceleration" of the first vector.
* Three text files withthe sufix .dat that have the smoothed data, the value of its square and the "acceleration" data.

An example of usage of the script is:

		./shadow.m testfile.txt
With these assumptions:

1. Octave is installed and able to run on the computer.
2. Octave packages that are called in the first lines of the script are also there.
3. In the same folder, it is also placed the file to work with and the octave files "display_rounded_matrix.m" and "deri.m"
