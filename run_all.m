% ----- Runs the whole project ------
%  Results and figures go into results 

projectFolder = fileparts(mfilename('fullpath'));
cd(projectFolder);
addpath(projectFolder);
% ...

part1_2d5          % 2D model, numerical error analysis    (Table 2, Figure 3)
part2_3d_verify14   % checks of the 3D solver               (Table 2, Figure 2)
part3_3d_sweep13    % iso vs aniso, kappa/theta sweep       (Table 3, Figures 4-5)
part4_optimise15    % best contact and voltage              (Table 4, Figures 1 and 6)
