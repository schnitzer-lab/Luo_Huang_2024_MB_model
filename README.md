This program is the source code for simulating and fitting the recurrent equation model to experimental data.

System requirements

    MATLAB v2020b or v2022a
  
Installation guide

    These source codes are MATLAB scripts which do not need to be installed. 

Demo and Instructions for use
  
    The functions of the recurrent equation model are in the folder “matlab_code\model_functions”. The time needed to simulate each set of parameters for the recurrent equation model is ~0.02 s. To generate the confidence interval (CI) in figures by simulating 10000 sets of randomly sampled parameters takes ~4 min. 
    
    The original data is in the file “data_and_parameters\Imaging_24hr_data.xlsx”.
    
    The code for fitting the experimental data is in the folder “matlab_code\model_fitting”. The code "fit_nonlinear_models.m" reads the original data file “Imaging_24hr_data.xlsx” and outputs the optimized data for the recurrent equation model into the folder “data_and_parameters”. The time needed for parameter optimization is ~2 min. 
    
    The parameters from curve fitting are the files “Dx_steady_state_nonlinear_3_27-Mar-2023_2modules.mat” and “Dx_steady_state_nonlinear_3_27-Mar-2023_3modules.mat” in folder “data_and_parameters”. “Dx_steady_state_nonlinear_3_27-Mar-2023_2modules.mat” stores the parameters of the 2-module model. “Dx_steady_state_nonlinear_3_27-Mar-2023_3modules.mat” stores the parameters of the 3-module model. 
    
    The functions used for reading and writing data and figures are in the folder “matlab_code\output_functions”.
    
    The example codes for drawing figures in the paper are in the folder “matlab_code\figure_code_examples”. 
 
        'plotPanel_5c_S10j_R2_CI.m'  plots Figure 5c and Extended data Figure 10j.
        
        'plotPanel_5d.m'             plots Figure 5d.
        
        'plotPanel_5e_R2.m'          plots Figure 5e.
        
        'plotPanel_5f_R2_CI.m'       plots Figure 5f.
        
        'plotPanel_5g_R2_CI.m'       plots Figure 5g and Extended data Figure 10l.
        
        'plotPanel_5h_R2_CI.m'       plots Figure 5h.
        
        'plotPanel_5i_R2.m'          plots Figure 5i, Extended data Figure 10m and 10h.
        
        'plotPanel_5j_R2_CI.m'       plots Figure 5j.
        
        'plotPanel_5k_R2_CI.m'       plots Figure 5k.
        
        'supplemental_figure_1.m'    plots Supplemental Figure 1 in Appendix.

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
    
