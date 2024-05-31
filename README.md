This program is the source code for simulating and fitting the recurrent equation model to experimental data.

System requirements
  MATLAB v2020b or v2022a
  
Installation guide
  These source codes are MATLAB scripts which do not need to be installed. 

Demo and Instructions for use
  The functions of the recurrent equation model are in the folder “matlab_code\model_functions”. The time needed to simulate each set of parameters for the recurrent equation model is ~0.02 s.
  
  The original data is in the file “data_and_parameters\Imaging_24hr_data.xlsx”.
  The code for fitting the experimental data is in the folder “matlab_code\model_fitting”. This code reads the original data and outputs the optimized data for the recurrent equation model. The time needed for parameter optimization is ~2 min. 
  The parameters from curve fitting are the files “Dx_steady_state_nonlinear_3_27-Mar-2023_2modules.mat” and “Dx_steady_state_nonlinear_3_27-Mar-2023_3modules.mat” in folder “data_and_parameters”.
  
  The functions used for reading and writing data and figures are in the folder “matlab_code\output_functions”.
  The example codes for drawing figures in the paper are in the folder “matlab_code\figure_code_examples”.

License
  
  Copyright (C) <2024>  <Junjie Luo, Cheng Huang, Mark J. Schnitzer>
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <https://www.gnu.org/licenses/>.
  
