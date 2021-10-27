function [fun,opt] = cviconfig(cvi)
%CVICONFIG CVI configuration object.
%   CVICONFIG(CVI) creates a CVI configuration object. CVI is a string 
%   representing the clustering criterion to be used.
%
%   CVICONFIG returns an error if the CVI does not exist in the toolbox.
%
%   [FUN, OPT] = CVICONFIG(CVI) returns a function handle FUN corresponding 
%   to the parameter CVI and string value OPT indicating if the CVI should
%   be maximied ('mx') or minimized ('mn').
%
% ------------------------------------------------------------------------
%   Version 1.0 (Matlab R2020b Unix)
%   Copyright (c) 2021, A. Jose-Garcia and W. Gomez-Flores
% ------------------------------------------------------------------------
% Type of cluster vality index
if strcmpi(cvi,'db')
    % Davies-Bouldin index
    fun = @dbindex;     
    opt = 'mn';
elseif strcmpi(cvi,'db2')
    % Davies-Bouldin index enhanced
    fun = @db2index;    
    opt = 'mn';    
elseif strcmpi(cvi,'dunn')
    % Dunn index
    fun = @dunnindex;   
    opt = 'mx';
elseif strcmpi(cvi,'gd31')
    % Dunn index variant 3,1
    fun = @gd31index;
    opt = 'mx';
elseif strcmpi(cvi,'gd41')
    % Dunn index variant 4,1
    fun = @gd41index; 
    opt = 'mx';
elseif strcmpi(cvi,'gd51')
    % Dunn index variant 5,1
    fun = @gd51index; 
    opt = 'mx';
elseif strcmpi(cvi,'gd33')
    % Dunn index variant 3,3
    fun = @gd33index; 
    opt = 'mx';
elseif strcmpi(cvi,'gd43')
    % Dunn index variant 4,3
    fun = @gd43index; 
    opt = 'mx';
elseif strcmpi(cvi,'gd53')
    % Dunn index variant 5,3
    fun = @gd53index; 
    opt = 'mx';
elseif strcmpi(cvi,'ch')
    % Calinski-Harabasz index
    fun = @chindex; 
    opt = 'mx';
elseif strcmpi(cvi,'pbm')
    % I or PBM index
    fun = @pbmindex;  
    opt = 'mx';
elseif strcmpi(cvi,'xb')
    % Xie-Beni index
    fun = @xbindex;
    opt = 'mn';
elseif strcmpi(cvi,'sdbw')
    % S_Dbw validity index
    fun = @sdbwindex;
    opt = 'mn';
elseif strcmpi(cvi,'sf')
    % Score function index
    fun = @sfindex;
    opt = 'mx';
elseif strcmpi(cvi,'sil')
    % Silhouette index
    fun = @silindex;
    opt = 'mx';
elseif strcmpi(cvi,'cs')
    % CS index
    fun = @csindex;
    opt = 'mn';    
elseif strcmpi(cvi,'cop')
    % COP index
    fun = @copindex;
    opt = 'mn';
elseif strcmpi(cvi,'sv')
    % SV index
    fun = @svindex;
    opt = 'mx';    
elseif strcmpi(cvi,'sym') 
    % Symmetry index
    fun = @symindex;
    opt = 'mx';
elseif strcmpi(cvi,'sdb')
    % Sym-Davies-Bouldin index
    fun = @symdbindex; 
    opt = 'mn';
elseif strcmpi(cvi,'sdunn')
    % Sym-Dunn index
    fun = @symdunnindex; 
    opt = 'mx';
elseif strcmpi(cvi,'cind')
    % C index
    fun = @cindex;
    opt = 'mn';
elseif strcmpi(cvi,'cvdd')
    % cvdd index
    fun = @cvddindex;
    opt = 'mx';
elseif strcmpi(cvi,'cvnn')
    % cvnn index
    fun = @cvnnindex;
    opt = 'mn';    
elseif strcmpi(cvi,'dbcv')
    % dbcv index
    fun = @dbcvindex;
    opt = 'mx';
elseif strcmpi(cvi,'lccv')
    % lccv index
    fun = @lccvindex;
    opt = 'mx';  
elseif strcmpi(cvi,'ssdd')
    % ssdd index
    fun = @ssddindex;
    opt = 'mn';
elseif strcmpi(cvi,'wb')
    % wb index
    fun = @wbindex;
    opt = 'mn';
else
    error('Unknown cluster validity index');
end