function kchic_coreg(inDir, threshMR,threshCT)

%% Preliminaries
% Set input directory 
if nargin<1
    inDir = spm_select(1,'dir');
    threshMR = 0.1;
    threshCT = 2000;
end

% Set SPM home directory (for tissure probability maps)
myspm = what('spm12');
tpmPath = [myspm.path filesep 'tpm'];
addpath(tpmPath);

% Select data
preT1Files = spm_select('FPListRec',inDir,'pre_T1w');
pstT1Files = spm_select('FPListRec',inDir,'post_T1w');
CTFiles    = spm_select('FPListRec',inDir,'CT');

%% Coregistration
% Set number of subjects for loop
nSub = size(preT1Files,1);

% Fire up SPM
spm('defaults','fmri');
spm_jobman('initcfg');

% Loop
for si = 1:nSub
    % Define things for current subject
    preT1Sub   = preT1Files(si,:);
    pstT1Sub   = pstT1Files(si,:);
    CTSub      = CTFiles(si,:);
    [path,nam] = spm_fileparts(preT1Sub);
    outnam = nam(1:10);
    
    % Define batch structure
    matlabbatch = [];
    
    % Segmentation of pre-op T1
    matlabbatch{1}.spm.spatial.preproc.channel.vols = {preT1Sub};
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[tpmPath filesep 'TPM.nii,1']};
    matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[tpmPath filesep 'TPM.nii,2']};
    matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[tpmPath filesep 'TPM.nii,3']};
    matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[tpmPath filesep 'TPM.nii,4']};
    matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[tpmPath filesep 'TPM.nii,5']};
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[tpmPath filesep 'TPM.nii,6']};
    matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
    
    % Coregistration of  post-op T1 to (bias-corrected) pre-op T1
    matlabbatch{2}.spm.spatial.coreg.estwrite.ref(1) = cfg_dep('Segment: Bias Corrected (1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
    matlabbatch{2}.spm.spatial.coreg.estwrite.source = {pstT1Sub};
    matlabbatch{2}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
    
    % Coregistration of implantation CT to (bias-corrected) pre-op T1
    matlabbatch{3}.spm.spatial.coreg.estwrite.ref(1) = cfg_dep('Segment: Bias Corrected (1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
    matlabbatch{3}.spm.spatial.coreg.estwrite.source = {CTSub};
    matlabbatch{3}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
    
    % Brain extraction
    matlabbatch{4}.spm.util.imcalc.input(1) = cfg_dep('Segment: Bias Corrected (1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
    matlabbatch{4}.spm.util.imcalc.input(2) = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
    matlabbatch{4}.spm.util.imcalc.input(3) = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
    matlabbatch{4}.spm.util.imcalc.input(4) = cfg_dep('Segment: c3 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
    matlabbatch{4}.spm.util.imcalc.output = [outnam '_brain.nii'];
    matlabbatch{4}.spm.util.imcalc.outdir = {path};
    matlabbatch{4}.spm.util.imcalc.expression = 'i1.*((i2+i3+i4)>0.5)';
    matlabbatch{4}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{4}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{4}.spm.util.imcalc.options.mask = 0;
    matlabbatch{4}.spm.util.imcalc.options.interp = 1;
    matlabbatch{4}.spm.util.imcalc.options.dtype = 4;
    
    % Mask coregistered CT with brain mask and intensity threshold
    matlabbatch{5}.spm.util.imcalc.input(1) = cfg_dep('Coregister: Estimate & Reslice: Resliced Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rfiles'));
    matlabbatch{5}.spm.util.imcalc.input(2) = cfg_dep('Image Calculator: ImCalc Computed Image: brain', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{5}.spm.util.imcalc.output = [outnam '_electrodes.nii'];
    matlabbatch{5}.spm.util.imcalc.outdir = {path};
    matlabbatch{5}.spm.util.imcalc.expression = ['(i1.*(i2>' num2str(threshMR) '))>' num2str(threshCT)];
    matlabbatch{5}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{5}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{5}.spm.util.imcalc.options.mask = 0;
    matlabbatch{5}.spm.util.imcalc.options.interp = 1;
    matlabbatch{5}.spm.util.imcalc.options.dtype = 4;
    
    % Run
    spm_jobman('run', matlabbbatch);   
end
end