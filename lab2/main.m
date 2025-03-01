% File Name:        main.m
% -------------------------------------------------------------------------
% 
% Purpose of file:  Contains the main function
%
% -------------------------------------------------------------------------
%
% Remarks:          You will need to write your own functions to...
%                   1)  Read the satellite coordinates that you generate
%						from the broadcast ephemeris
%					2) 	Read the igs ephemeris files using the readSP3()
%						function
%					3) 	Generate sky plot of satellites based on the
%					    satellite  and the pillar coordiantes
%                   4)  Calculte the error of different kinds of ephemeris
%                       with precise orbits (igswwwwd.sp3) as reference

%  You can modify any part of this code to make it more convenient for you
%  to use. This just provides a general procedure.
%-------------------------------------------------------------------------%


% Main function...

% clear memory, figures, command window, output numbers in standard notation
clear all;
close all;
clc;
format long g;

%% User Inputs
prn = 11;

preciseFile = 'igswwwwd.sp3';
brdc24File = 'youroutput24h.txt';
brdcbestFile = 'youroutputbest.txt';

almanacFile = 'Almanac-SatPos-12.txt';
urPredFile = 'iguwwwwd_00.sp3';
urObsFile = 'iguwwwwd+1_00.sp3';
rapidFile = 'igrwwwwd.sp3';

%% Read Files
refsat=readSP3(preciseFile);   %%Read the sp3 file
%% Sky plot of satellkite
   %calculate the vector from receiver to satellites in xyz
   %transform the vector to ENU
   %calculate the azimuth and elevation
   %plot the sky plot
   %calculate the DOP
%% Parse out PRN 11
idx=find(refsat(:,1)==prn);
PRN11=refsat(idx,:);
%% Plots

%% Statistics
