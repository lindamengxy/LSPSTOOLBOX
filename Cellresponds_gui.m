function [ output_args ] =Cellresponds()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%average Maps based on age.
eventwindow=50;
direct_t=8;

filepath='/Volumes/disk1/kawens/Norm/eventwindow50';
cd(filepath);
a=dir('*.mat')
Cellfile=size(a);
C=cell(Cellfile);
 foname='/Volumes/disk1/kawens/Norm/eventwindow50/singletraces';
 
for fn=1:1:Cellfile
    fn
    filename=fullfile(filepath,a(fn).name);
    C{fn}=load(filename);
    
    
    cells=C{fn}.cells;
   
    
    
    
    
      close all
       
       multipletraces_gui(cells,foname);
    
end 
end
  
  
  
  