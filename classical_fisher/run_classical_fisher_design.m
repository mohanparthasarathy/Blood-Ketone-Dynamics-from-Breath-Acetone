function D=run_classical_fisher_design(dataFile,outputDir,userOpts)
%RUN_CLASSICAL_FISHER_DESIGN Classical FIM and sparse measurement design.
% Sensitivities use log-physical parameters. Wavelength design includes
% spectral baseline/interferent nuisance columns; no quantum optimization.

here=fileparts(mfilename('fullpath')); root=fileparts(here);
if nargin<1 || isempty(dataFile)
    dataFile=fullfile(root,'synthetic_data','outputs','cohort','synthetic_cohort.mat');
    if ~isfile(dataFile); dataFile=fullfile(root,'matlab','outputs','cohort','synthetic_cohort.mat'); end
end
if nargin<2 || isempty(outputDir); outputDir=fullfile(here,'outputs'); end
if nargin<3; userOpts=struct; end
o=defaults(userOpts); if ~exist(outputDir,'dir'); mkdir(outputDir); end
addpath(fullfile(root,'synthetic_data'),'-begin'); addpath(fullfile(root,'matlab'),'-begin');
S=load(dataFile); id=o.patient_id; idx=find(S.meta.patient_id==id);
if isempty(idx); error('Patient %d not found.',id); end
[~,ord]=sort(S.meta.breath_time_h(idx)); idx=idx(ord); M=S.meta(idx,:);
scenario=string(M.scenario(1)); theta=truth_theta(M);
names=["log_r_pre","log_r_post","log_tau_c","log_k_r", ...
    "log_Gamma_A","log_Pe","log_Da_aw"];

% Candidate measurement grid.
phaseIdx=unique(round(linspace(1,size(S.spectra,1),min(o.phase_candidates,size(S.spectra,1)))));
waveIdx=unique(round(linspace(1,numel(S.wn),min(o.wavelength_candidates,numel(S.wn)))));
wn=double(S.wn(waveIdx)); tq=double(S.time_s(phaseIdx,idx));
if max(abs(tq-tq(:,1)),[],'all')<1e-8; phaseTime=tq(:,1); else; phaseTime=mean(tq,2); end
breathTimes=double(M.breath_time_h);

[Y0,B0,S0]=simulate(theta,M,S.base,wn,phaseTime,scenario,o.fast_nx);
nP=numel(theta); J=zeros([size(Y0) nP]); JB=zeros(numel(B0),nP); JS=JB;
for p=1:nP
    tp=theta; tm=theta; tp(p)=theta(p)*exp(o.fd_step); tm(p)=theta(p)*exp(-o.fd_step);
    [yp,bp,sp]=simulate(tp,M,S.base,wn,phaseTime,scenario,o.fast_nx);
    [ym,bm,sm]=simulate(tm,M,S.base,wn,phaseTime,scenario,o.fast_nx);
    J(:,:,:,p)=(yp-ym)/(2*o.fd_step); JB(:,p)=(bp-bm)/(2*o.fd_step);
    JS(:,p)=(sp-sm)/(2*o.fd_step);
end

% Add normalized spectral nuisance parameters: offset, linear drift, broad
% water-like background, and an overlapping VOC. These create a meaningful
% wavelength discrimination problem while leaving the acetone simulator intact.
x=(wn-mean(wn))/max(eps,max(wn)-min(wn)); scale=o.nuisance_scale_abs;
shapes=[ones(numel(wn),1),x(:), ...
    exp(-.5*((wn(:)-1640)/95).^2), ...
    exp(-.5*((wn(:)-1460)/42).^2)];
shapes=shapes./max(abs(shapes),[],1)*scale;
nN=size(shapes,2); Jall=zeros([size(Y0) nP+nN]); Jall(:,:,:,1:nP)=J;
for k=1:nN
    Jall(:,:,:,nP+k)=repmat(reshape(shapes(:,k),1,numel(wn),1), ...
        size(Y0,1),1,size(Y0,3));
end
allNames=[names,"offset","linear_baseline","water_background","VOC_overlap"];
sigma=max(o.noise_floor_abs,median(double(M.noise_sd_abs)));

% Contributions by wavelength, metabolic breath, and within-breath phase.
Fw=block_information(Jall,2,sigma); Ft=block_information(Jall,1,sigma);
Fs=block_information(Jall,3,sigma);
Ffull=sum(Fw,3); prior=diag(1./o.prior_sd_log(:).^2);
Freg=Ffull+prior;
[diagTable,corrMatrix,eigenTable]=fim_diagnostics(Ffull,Freg,allNames,nP);

targetRow=find(abs(breathTimes-o.target_time_h)==min(abs(breathTimes-o.target_time_h)),1);
G=zeros(2,nP+nN); G(1,1:nP)=JB(targetRow,:)/o.target_B_scale_mM;
G(2,1:nP)=JS(targetRow,:)/o.target_slope_scale_mM_h;

subsetTable=[]; selected=struct;
for m=o.wavelength_subset_sizes
    if m>numel(wn); continue; end
    sets=struct;
    sets.greedy_D=greedy(Fw,m,prior,[],"D");
    sets.target_BHB=greedy(Fw,m,prior,G,"target");
    [~,peakOrder]=sort(squeeze(sum(Y0.^2,[1 3])),'descend'); sets.peak=peakOrder(1:m).';
    sets.uniform=unique(round(linspace(1,numel(wn),m)));
    rng(o.random_seed+m,'twister'); sets.random=sort(randperm(numel(wn),m));
    methods=fieldnames(sets);
    for q=1:numel(methods)
        method=methods{q}; sel=sets.(method); F=prior+sum(Fw(:,:,sel),3);
        [ld,tv]=criteria(F,G); list=strjoin(compose('%.1f',wn(sel)),',');
        row=table(m,string(method),string(list),ld,tv,cond(F), ...
            'VariableNames',{'N_wavelengths','method','wavenumbers_cm1', ...
            'regularized_logdet','target_variance_score','condition_number'});
        subsetTable=append(subsetTable,row); selected.(sprintf('%s_m%d',method,m))=wn(sel);
    end
end

% Greedy targeted selection of breath times and breath-phase samples.
timeTable=selection_path(Ft,prior,G,breathTimes,"breath_time_h",o.time_subset_sizes);
phaseTable=selection_path(Fs,prior,G,phaseTime,"phase_time_s",o.phase_subset_sizes);
schur=nuisance_adjusted(Ffull,1:5,6:(nP+nN),o.schur_ridge);

D=struct('patient_id',id,'scenario',scenario,'theta',theta,'parameter_names',allNames, ...
    'noise_sd_abs',sigma,'wavenumber_cm1',wn,'phase_time_s',phaseTime, ...
    'breath_time_h',breathTimes,'F_full',Ffull,'F_regularized',Freg, ...
    'F_metabolic_adjusted',schur,'diagnostics',diagTable, ...
    'correlation_matrix',corrMatrix,'eigenvalues',eigenTable, ...
    'wavelength_designs',subsetTable,'time_design',timeTable, ...
    'phase_design',phaseTable,'selected_wavelengths',selected, ...
    'target_time_h',breathTimes(targetRow),'options',o);
writetable(diagTable,fullfile(outputDir,'fim_diagnostics.csv'));
writetable(array2table(corrMatrix,'VariableNames',cellstr(allNames), ...
    'RowNames',cellstr(allNames)),fullfile(outputDir,'parameter_correlations.csv'), ...
    'WriteRowNames',true);
writetable(eigenTable,fullfile(outputDir,'fim_eigenvalues.csv'));
writetable(subsetTable,fullfile(outputDir,'wavelength_designs.csv'));
writetable(timeTable,fullfile(outputDir,'breath_time_design.csv'));
writetable(phaseTable,fullfile(outputDir,'breath_phase_design.csv'));
save(fullfile(outputDir,'classical_fisher_design.mat'),'D','-v7.3');
make_plots(D,outputDir); disp(diagTable); disp(subsetTable); disp(timeTable); disp(phaseTable)
end

function [Y,B,dB]=simulate(th,M,base,wn,tq,scenario,nx)
p=base; p.fast.nx=nx; p.slow.tau_c=th(3); p.slow.k_r=th(4);
p.slow.Gamma_A=th(5); p.fast.Pe=th(6); p.fast.Da_aw=th(7);
if scenario=="developing"; p.slow.r0=th(1); p.slow.delta_r=max(0,th(2)-th(1));
elseif scenario=="resolving"; p.slow.r0=th(2); p.slow.delta_r=max(0,th(1)-th(2));
else; p.slow.r0=th(1); p.slow.delta_r=0; end
p.scales.A_star_mM=p.slow.kappa_QA*p.scales.K_star_mM/p.slow.Gamma_A;
p.scales.C_star_ppm=p.scales.A_star_mM/p.scales.lambda_ba*24450;
slow=solve_slow_metabolism(p,char(scenario)); t=double(M.breath_time_h);
B=interp1(slow.time_h,slow.B_mM,t,'pchip');
b=interp1(slow.time_h,slow.b,t,'pchip'); q=interp1(slow.time_h,slow.q,t,'pchip');
r=interp1(slow.time_h,slow.r,t,'pchip');
dB=p.scales.K_star_mM/p.scales.T_m_h*(p.slow.rho_B*r- ...
    (1+p.slow.kappa_BQ)*b+p.slow.kappa_QB*q);
a=interp1(slow.time_h,slow.a,t,'pchip'); Y=zeros(numel(t),numel(wn),numel(tq));
for k=1:numel(t)
    z=solve_fast_breath(a(k),p,M.U_scale(k)); C=interp1(z.time_s,z.C_exh_ppm,tq,'pchip');
    [wf,A]=synthetic_acetone_spectrum(C,p);
    spec=interp1(wf,A.',wn,'pchip').'; % phase x wavelength
    Y(k,:,:)=permute(spec,[3 2 1]);
end
end

function th=truth_theta(M)
s=string(M.scenario(1)); r0=M.r0(1); dr=M.delta_r(1);
if s=="developing"; rp=r0; rq=r0+dr; elseif s=="resolving"; rp=r0+dr; rq=r0;
else; rp=r0; rq=r0; end
th=[rp rq M.tau_c(1) M.k_r(1) M.Gamma_A(1) M.Pe(1) M.Da_aw(1)];
end

function F=block_information(J,blockDim,sigma)
sz=size(J); np=sz(4); nb=sz(blockDim); F=zeros(np,np,nb);
for k=1:nb
    if blockDim==1; A=reshape(J(k,:,:,:),[],np);
    elseif blockDim==2; A=reshape(J(:,k,:,:),[],np);
    else; A=reshape(J(:,:,k,:),[],np); end
    F(:,:,k)=(A.'*A)/sigma^2;
end
end

function sel=greedy(Fblocks,m,prior,G,mode)
n=size(Fblocks,3); sel=zeros(1,m); available=true(1,n); F=prior;
for k=1:m
    score=-inf(1,n);
    for j=find(available)
        C=F+Fblocks(:,:,j);
        if mode=="D"; score(j)=stable_logdet(C);
        else; score(j)=-trace(G*(C\G.')); end
    end
    [~,sel(k)]=max(score); available(sel(k))=false; F=F+Fblocks(:,:,sel(k));
end
sel=sort(sel);
end
function T=selection_path(Fb,prior,G,values,name,sizes)
T=[];
for m=sizes
    if m>numel(values); continue; end
    sel=greedy(Fb,m,prior,G,"target"); F=prior+sum(Fb(:,:,sel),3); [ld,tv]=criteria(F,G);
    row=table(m,string(strjoin(compose('%.4g',values(sel)),',')),ld,tv, ...
        'VariableNames',["N_selected",name,"regularized_logdet","target_variance_score"]);
    T=append(T,row);
end
end
function [ld,tv]=criteria(F,G); ld=stable_logdet(F); tv=trace(G*(F\G.')); end
function x=stable_logdet(F); R=chol((F+F.')/2+1e-12*eye(size(F))); x=2*sum(log(diag(R))); end
function [T,corr,eigenTable]=fim_diagnostics(F,Freg,names,nPhysical)
e=eig((F+F.')/2); tol=max(e)*max(size(F))*eps(max(e)); rankRaw=sum(e>tol);
cov=inv(Freg); se=sqrt(max(0,diag(cov))); corr=cov./max(eps,se*se.');
T=table(names.',diag(F),diag(Freg),se, ...
    [rankRaw;nan(numel(names)-1,1)],[cond(Freg);nan(numel(names)-1,1)], ...
    'VariableNames',{'parameter','raw_information','regularized_information', ...
    'CRLB_scaled_SE','raw_FIM_rank','regularized_condition_number'});
T.Properties.UserData.correlation=corr; T.Properties.UserData.nPhysical=nPhysical;
eigenTable=table((1:numel(e)).',sort(e,'descend'), ...
    'VariableNames',{'index','raw_FIM_eigenvalue'});
end
function S=nuisance_adjusted(F,interest,nuisance,ridge)
A=F(interest,interest); B=F(interest,nuisance); C=F(nuisance,nuisance);
S=A-B*((C+ridge*eye(size(C)))\B.'); S=(S+S.')/2;
end
function T=append(T,A); if isempty(T); T=A; else; T=[T;A]; end; end

function make_plots(D,out)
f=figure('Color','w','Position',[100 100 1100 430]); tiledlayout(1,2,'TileSpacing','compact');
nexttile; imagesc(log10(abs(D.F_regularized)+eps)); colorbar; axis square
xticks(1:numel(D.parameter_names)); yticks(1:numel(D.parameter_names));
xticklabels(D.parameter_names); yticklabels(D.parameter_names); xtickangle(45); title('log_{10}|F|');
nexttile; T=D.wavelength_designs; methods=unique(T.method); hold on
for k=1:numel(methods); Q=T(T.method==methods(k),:); plot(Q.N_wavelengths,Q.target_variance_score,'-o','DisplayName',methods(k)); end
set(gca,'YScale','log'); grid on; xlabel('Wavelength count'); ylabel('Target variance score'); legend('Location','best');
exportgraphics(f,fullfile(out,'fisher_design_summary.png'),'Resolution',200); close(f)
end

function o=defaults(o)
d=struct('patient_id',14,'fd_step',.01,'fast_nx',61,'phase_candidates',21, ...
 'wavelength_candidates',61,'noise_floor_abs',2e-6,'nuisance_scale_abs',1e-4, ...
 'wavelength_subset_sizes',[3 5 8 10],'time_subset_sizes',[2 3 5 8], ...
 'phase_subset_sizes',[3 5 10 15],'target_time_h',4.25, ...
 'target_B_scale_mM',1,'target_slope_scale_mM_h',.5,'schur_ridge',1e-8, ...
 'prior_sd_log',[2 2 1 1 1 1 1 10 10 10 10],'random_seed',9217);
f=fieldnames(d); for k=1:numel(f); if ~isfield(o,f{k}); o.(f{k})=d.(f{k}); end; end
end
