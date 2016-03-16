## Copyright (C) 2013 Gonzalo Rodríguez Prieto
##
## This file is part of Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {} polyadj (@var{x},@var{fx},options)
##   Adjust some points to the best polynomium among some.
##
## @code{padj = polyadj (@var{x},@var{fx},options)} returns best polynomium that adjusts to the function
## @var{fx} of @var{x} values. Options: "order", ord; minimum order of polynomia to adjust.
## "max_order", max_order; and maximum order of polynomia.
## By default, @var{ord}=1 and @var{max_order} = 10.
##
## @code{[padj, padjstruc, st_er, t] = polyadj (@var{x},@var{fx}, @var{order}, @var{max_order}))} returns 
## polynomium in @var{padj}, structure in @var{padjstruc}, standard errors of polynomium values in 
## @var{st_er} and t value in @var{t}.
## @end deftypefn

## Author: Gonzalo Rodríguez Prieto (gonzalo.rprietoATuclm.es)
## Created: Nov 2013
## Modified: Jun 2015
## Modified: Mar 2016



function [padj, padjstruc, st_er, t] = polyadj(x,fx, varargin)

#########
# Control error part
#########
if (nargin < 2) #If too few parameters.
 print_usage;
endif;

if (nargin > 6) #If too many parameters.
 print_usage;
endif;

if (isscalar(x) ==1) #Works when a vector is passed, not a escalar or other type.
 error("polyadj: Variable x must be a vector.");
endif;

if (isscalar(fx) ==1) #Works when a vector is passed, not a escalar or other type.
 error("polyadj: Variable fx must be a vector.");
endif;

if ( length(x) != length(fx) ) #Length of both vectors are not the same.
 error("polyadj: Lengths of vectors x and fx must coincide.");
endif;

#Defaults:
ord = 1;
max_order = 10;

  ## parse options for order and precision
  idx = [];
  if ( nargin > 2)
    for i = 1:nargin-2
      arg = varargin{i};
      if ischar(arg)
        switch arg
          case "order"
            ord = varargin{i+1};
            idx = [idx,i,i+1];
          case "max_order"
            max_order = varargin{i+1};
            idx = [idx,i,i+1];
        endswitch
      endif
    endfor
  endif
  varargin(idx) = [];
  %~ options = varargin;




#########
# Adjusting loop part
#########

n = ord;
t_vector = [];
t_tot = 666;

#Fitting loop :
while ( (n<max_order) && (t_tot>0) ) %n is polynomial order.
	#Fitting:
	[padj, padjstruc] = polyfit(x,fx, n);

	#Calculation of the statistical errors of the fit:
		%From http://www.facstaff.bucknell.edu/maneval/help211/fitting.html
	R = padjstruc.R; %The "R", whatever this it...
	d = ( R' * R)\eye(n+1); %The covariance matrix
	d = diag(d)'; %ROW vector of the diagonal elements.
	MSE = (padjstruc.normr^2)/ padjstruc.df; %Variance of the residuals.
	st_er = sqrt(MSE*d); %Standard errors
	t = abs(padj)./st_er; %Observed t-values. The bigger, the better.
	t_tot = sum(t);
	t_vector = [t_vector,t_tot];
	n = n +1;
endwhile;

[tmax,pos] = max(t_vector); %Better adjustment to the data when t is maximum.
best_fit = ord + (pos-1);

[padj, padjstruc] = polyfit(x,fx, best_fit );

#Calculation of the statistical errors of the fit:
	%From http://www.facstaff.bucknell.edu/maneval/help211/fitting.html
R = padjstruc.R; %The "R", whatever this it...
d = ( R' * R)\eye(best_fit+1); %The covariance matrix
d = diag(d)'; %ROW vector of the diagonal elements.
MSE = (padjstruc.normr^2)/ padjstruc.df; %Variance of the residuals.
st_er = sqrt(MSE*d); %Standard errors
t = abs(padj)./st_er; %Observed t-value. The bigger, the better.

endfunction;

#That's...that's all, folks! 
