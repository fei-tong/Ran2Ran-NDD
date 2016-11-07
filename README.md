# Ran2Ran-NDD

This project is to get the Ran2Ran distance distribution within an arbitrary polygon,
    based on the line kinematic measure approach and the probabilistic sum method.
Node density is considered.

The project contains two folders: 
* /single_triangle for calculating the Ran2Ran distance distribution within a single triangle, and
* /any_two_triangles for calculating the Ran2Ran distance distribution between any two triangles.
 
Two sections will be introduced below: Sec1: Approach; Sec2: Simualtion
## Sec1: Approach:
   The genergal function is 
   * -> function [ d_array, cdf_array ] = f_rand2rand_arbitrary_polygon( triangle_cell,varargin )</br>
        % This function is to get the numerical Ran2Ran distance distribution </br>
        %   within an arbitrary polygon.</br>
        %</br>
        % Input:</br>
        %   triangle_cell: Triangulated triangules of an arbitrary polygon. Each</br>
        %       cell element contains a triangle, which is like</br>
        %       [x1 y1;x2 y2;x3 y3], where [xi yi] is a vertex of the triangle.</br>
        %</br>
        %   varargin contains the following one optional argument:</br>
        %       -> density array: (n-2)*1, containing "node densities" of each</br>
        %          triangle.</br>
        % Output:</br>
        %   [ d_array, cdf_array ]: Ran2Ran distance distribution within the polygon</br>

    This general function will call the following to functions to get the results associated with the triangulated triangles. 
    Then with the probabilistic sum method, the result assiciated with the polygon can be obtained:
   * -> in folder /single_triangle</br>
       function [d_array,pdd_pdf,pdd_cdf] = f_formula_pdd_pdf_triangle(x,y,d_step)%(X,Y,Z,a,b,c)</br>
        % This function is to get the Ran2Ran distance distribution within a</br>
        % single triangle, including both pdf and cdf.</br>
        % Input:</br>
        %   (x,y): the 3 vertexes of the triangle.</br>
        %   d_step: A large number which determines the accuracy of the results.</br>
        % Output:</br>
        %   numerical results of pdf: [d_array,pdd_pdf] and cdf: [d_array,pdd_cdf]</br>
   * -> in folder /any_two_triangles</br>
       function [ d_array, pdf_array, cdf_array ] = f_rand2rand_between_any_2_triangles( t1,t2,d_step )</br>
        % This function is to get the Ran2Ran distance distribution between any</br>
        % two triangles.</br>
        %</br>
        % Input:</br>
        %   t1: triangle 1, t1 is like [x1 y1;x2 y2;x3 y3], where [xi yi] is a</br>
        %       vertex of t1.</br>
        %   t2: triangle 2, t2 is similar to t1;</br>
        %   d_step: A large number which determines the accuracy of the results.</br>
        %</br>
        % Output:</br>
        %   [d_array, cdf_array]: Numerical result of the Ran2Ran distance</br>
        %                         distribution between 2 triangles</br>
        
       The above function: f_rand2rand_between_any_2_triangles will call the following function:
   * -> function [ L ] = f_theta_func_between_any_2_triangles( theta,t1,t2 )</br>
        % This function is to get the lengths of 3 segments which are generated due</br>
        % to the intersection between a line with angle theta (with regard to</br>
        % x-axis) and the two triangles. </br>
        % The function will be called by the following function:</br>
        % function [ d_array, pdf_array, cdf_array ] = f_rand2rand_between_any_2_triangles( t1,t2,d_step )</br>
        %</br>
        % Input:</br>
        %   theta: the angle of an arbitrary line passed the origin.</br>
        %   t1: triangle 1, t1 is like [x1 y1;x2 y2;x3 y3], where [xi yi] is a</br>
        %       vertex of t1.</br>
        %   t2: triangle 2, t2 is similar to t1.</br>
        % Output:</br>
        %   L: n*3, Lengths of three segments. n is determined by the step distance</br>
        %       between two adjacent parallel lines with angle of theta.</br>
        

## Sec2: Simulation:
   * -> function [d_array,pdd_cdf,sim_density] = f_sim_rand2rand_arbitrary_polygon( triangle_cell,varargin)</br>
        % Simulation: to get the Ran2Ran distance distribution within an</br>
        %             arbitrary polygon triangulated into several triangules,</br>
        %             stored in triangle_cell</br>
        % Input: </br>
        %   triangle_cell: Triangulated triangules of an arbitrary polygon.Each</br>
        %       cell element contains a triangle, which is like</br>
        %       [x1 y1;x2 y2;x3 y3], where [xi yi] is a vertex of the triangle. </br>
        %</br>
        %   varargin contains the following one optional argument:</br>
        %       -> density array: (n)*1, containing "node densities" of each</br>
        %          triangle (n: the number of triangles).</br>
        % Output:</br>
        %   [ d_array, pdd_cdf, sim_density ]: Ran2Ran distance distribution within the</br>
        %                      polygon, with the simulated density of each triangle.</br>
