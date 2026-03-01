We have encountered two bugs in the lifting line code that can create problems in certain situations.
There are two easy ways to deal with these bugs.

--------------------------------------------------------------
Bug 1: If the initial reference speed "uoo" in main_simple is set to a low speed, fsolve may not converge.

Solution: The easiest way to solve this is to either keep the speed at a high value or to decrease the step size in "valpha" (input vector to plainpolar) i.e. increase the number of points in the linspace function.

If this does not work, it is possible to remove the stop condition for when fsolve did not converged. This is done in the "pointsolve" function by removing (or commenting) the last if statement. No gurantee what happens with the results by doing this, but atleast the code runs regardless of velocity and number of vortices. Keep an eye on the residuals manually.

--------------------------------------------------------------
Bug 2: If the variable "delta" defined in build_h15() (or any other build function) is not a vector the program will not run. This is a problem if you only have xfoil data for one flap angle.

Solution: The easy fix here is to define a "dummy" element in the delta vector for a second flap angle (which you don't have data for) and then use the same xfoil data for both flap angles. 

For h15 it could look like this:

function sps = build_h15()
	% sps = build_h15()
	%
	% Example of how to use assemble_afgrid with computed XFoil polar files.
	% This  particular call constructs the interpolation tables for the H15
	% airfoil used in the wing lab model.

    	dir = '../xfoil/';
    	fnames = cell(3, 2);
    	reynolds = [3 7 15]' * 1e5;
    	delta = [0 1]' * pi/180;
    
    	fnames{1,1} = [dir 'h15re3e5.txt'];
    	fnames{2,1} = [dir 'h15re7e5.txt'];
    	fnames{3,1} = [dir 'h15re15e5.txt'];
    
    	fnames{1,2} = [dir 'h15re3e5.txt'];
    	fnames{2,2} = [dir 'h15re7e5.txt'];
    	fnames{3,2} = [dir 'h15re15e5.txt'];


    	sps = assemble_afgrid(fnames, reynolds, delta);
    
end
--------------------------------------------------------------