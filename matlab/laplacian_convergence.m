%convergenza Laplaciano 

clc
clear all
close all

  hh = [1/20; 1/40; 1/80; 1/160];
  file1 = fopen('laplacianBasic_errorsL2','r');
  file2 = fopen('laplacianBasic_errorsH1','r');

for i=1:4
  test1 = fscanf(file1, '%[qwertyuioplkjhgfdssazxcvbnm01]');
  test2 = fscanf(file2, '%[qwertyuioplkjhgfdssazxcvbnm01]');
  ERR_L2=fscanf(file1,'%f');
  ERR_H1=fscanf(file2,'%f');
  
  
  figure()
  loglog(hh,ERR_L2,'r','LineWidth',2)
  hold on
  grid on
  loglog(hh,hh.^2,'r--','LineWidth',2)
  axis tight
  legend('||u-u_h||_{L^2(\Omega)}','h^2')
  title(test1)
  
  figure()
  loglog(hh,ERR_H1,'b','LineWidth',2)
  hold on
  grid on
  loglog(hh,hh,'b--','LineWidth',2)
  axis tight
  legend('||u-u_h||_{H^1(\Omega)}','h')
  title(test2)
end