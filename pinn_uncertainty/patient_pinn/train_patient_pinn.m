function fit=train_patient_pinn(tObs,Cobs,gain,base,fitMask,o,trueTheta,scenario)
%TRAIN_PATIENT_PINN Multistart, curriculum-trained patient-specific PINN.
scenario=string(scenario); fitIdx=find(fitMask); valMask=false(size(fitMask));
if numel(fitIdx)>=4; valMask(fitIdx(end-1))=true; end
trainMask=fitMask & ~valMask;
gammaPrior=NaN;
if string(o.gamma_mode)=="prior"
    sigmaLog=sqrt(log(1+o.gamma_prior_cv^2));
    rng(o.random_seed+round(1e4*trueTheta(8)),'twister');
    gammaPrior=trueTheta(8)*exp(sigmaLog*randn-.5*sigmaLog^2);
end

candidates=cell(o.num_restarts,1); scores=inf(o.num_restarts,1);
for restart=1:o.num_restarts
    fprintf('  restart %d/%d\n',restart,o.num_restarts);
    candidates{restart}=initialize_candidate(restart);
    candidates{restart}=adam_stage(candidates{restart},trainMask,o.max_iterations,NaN,"train");
    scores(restart)=score_candidate(candidates{restart},valMask,trainMask);
    candidates{restart}.selection_score=scores(restart);
    fprintf('    selection score %.3e; %d Adam iterations\n',scores(restart), ...
        candidates{restart}.iterations_completed);
end
[~,best]=min(scores); fit=candidates{best}; fit.best_restart=best;
if o.refit_iterations>0
    fit=adam_stage(fit,fitMask,o.refit_iterations,1,"refit");
end
if o.use_lbfgs && exist('lbfgsupdate','file')==2 && o.lbfgs_iterations>0
    adamFit=fit;
    try
        fit=lbfgs_stage(fit,fitMask);
    catch ME
        fit=adamFit;
        fit.lbfgs_error=string(ME.message);
        warning('L-BFGS refinement skipped: %s',ME.message);
    end
elseif o.use_lbfgs && o.verbose
    warning('lbfgsupdate unavailable; returning converged Adam solution.');
end
fit.final_loss=fit.history(end,1);

    function cand=initialize_candidate(restart)
        rng(o.random_seed+1009*restart+round(100*trueTheta(4)),'twister');
        layers=featureInputLayer(1,Normalization="none",Name="time");
        for kk=1:o.hidden_layers
            layers=[layers;fullyConnectedLayer(o.hidden_width,Name="fc"+kk); ...
                tanhLayer(Name="tanh"+kk)]; %#ok<AGROW>
        end
        net=dlnetwork([layers;fullyConnectedLayer(3,Name="states")]);
        lb=[.05 .05 .5 .5 .3]; ub=[3.5 3.5 7 5 4.5];
        trend=(Cobs(end)-Cobs(1))/max(.1,abs(Cobs(1)));
        if scenario=="baseline"
            rr=.12+.18*rand; start=[rr rr 3 2 1.5];
        elseif trend>.10
            start=[.12+.18*rand 1.8+rand 2.5+rand 1.4+1.4*rand 1.5];
        else
            start=[1.8+rand .12+.18*rand 2.5+rand 1.4+1.4*rand 1.5];
        end
        if string(o.gamma_mode)=="oracle"; start(5)=trueTheta(8);
        elseif string(o.gamma_mode)=="population"; start(5)=o.gamma_population;
        elseif string(o.gamma_mode)=="prior"; start(5)=gammaPrior; end
        raw=dlarray(single(logit((start-lb)./(ub-lb))));
        cand=struct('net',net,'raw_parameters',raw,'history',zeros(0,5), ...
            'lb',lb,'ub',ub,'iterations_completed',0,'selection_score',Inf);
        cand=decode(cand);
    end

    function cand=adam_stage(cand,mask,nIter,fixedRamp,label)
        net=cand.net; raw=cand.raw_parameters;
        avgNet=[]; avgSqNet=[]; avgRaw=[]; avgSqRaw=[];
        tauC=dlarray(single(linspace(0,o.t_max_h,o.num_collocation)),'CB');
        tauO=dlarray(single(tObs(mask).'),'CB'); Cfit=dlarray(single(Cobs(mask).'),'CB');
        Gfit=dlarray(single(gain(mask).'),'CB'); bestMetric=Inf; stale=0;
        for it=1:nIter
            if isnan(fixedRamp)
                progress=it/max(1,nIter);
                ramp=min(1,max(0,(progress-o.curriculum_data_fraction))/ ...
                    max(eps,o.curriculum_ramp_fraction));
            else; ramp=fixedRamp; end
            w=weights(ramp);
            [loss,gNet,gRaw,parts,physical]=dlfeval(@pinn_model_loss,net,raw, ...
                tauC,tauO,Cfit,Gfit,base,o,cand.lb,cand.ub,trueTheta, ...
                gammaPrior,scenario,w);
            lr=o.learning_rate; if it>.7*nIter; lr=.25*lr; end
            [net,avgNet,avgSqNet]=adamupdate(net,gNet,avgNet,avgSqNet,it,lr);
            [raw,avgRaw,avgSqRaw]=adamupdate(raw,gRaw,avgRaw,avgSqRaw,it,lr);
            row=[double(gather(extractdata(loss))) double(gather(extractdata(parts)))];
            cand.history(end+1,:)=row; %#ok<AGROW>
            if it==1 || mod(it,o.print_every)==0
                fprintf('    %s %5d loss %.3e data %.3e physics %.3e Gamma %.3f\n', ...
                    label,it,row(1),row(2),row(3),double(gather(extractdata(physical(5)))));
            end
            if it>=o.min_iterations && mod(it,o.early_stop_check_every)==0
                if row(1)<bestMetric*(1-o.early_stop_rel_tol); bestMetric=row(1); stale=0;
                else; stale=stale+o.early_stop_check_every; end
                if stale>=o.early_stop_patience
                    if o.verbose; fprintf('    convergence plateau detected\n'); end
                    break
                end
            end
        end
        cand.net=net; cand.raw_parameters=raw;
        cand.iterations_completed=cand.iterations_completed+it; cand=decode(cand);
    end

    function cand=lbfgs_stage(cand,mask)
        tauC=dlarray(single(linspace(0,o.t_max_h,o.num_collocation)),'CB');
        tauO=dlarray(single(tObs(mask).'),'CB'); Cfit=dlarray(single(Cobs(mask).'),'CB');
        Gfit=dlarray(single(gain(mask).'),'CB'); rawFixed=cand.raw_parameters; w=weights(1);
        state=lbfgsState; net=cand.net;
        lossFcn=@(n) dlfeval(@network_loss,n,rawFixed,tauC,tauO,Cfit,Gfit, ...
            base,o,cand.lb,cand.ub,trueTheta,gammaPrior,scenario,w);
        for jj=1:o.lbfgs_iterations
            [net,state]=lbfgsupdate(net,lossFcn,state);
            if o.verbose && mod(jj,o.lbfgs_print_every)==0
                fprintf('    L-BFGS %4d loss %.3e\n',jj,double(state.Loss));
            end
            if state.GradientsNorm<o.lbfgs_gradient_tolerance || ...
                    state.StepNorm<o.lbfgs_step_tolerance || ...
                    state.LineSearchStatus=="failed"; break; end
        end
        cand.net=net; cand.lbfgs_iterations_completed=jj;
        [L,~,~,parts]=dlfeval(@pinn_model_loss,net,rawFixed,tauC,tauO,Cfit,Gfit, ...
            base,o,cand.lb,cand.ub,trueTheta,gammaPrior,scenario,w);
        cand.history(end+1,:)=[double(gather(extractdata(L))) double(gather(extractdata(parts)))];
    end

    function s=score_candidate(cand,vMask,trMask)
        if any(vMask); mm=vMask; else; mm=trMask; end
        state=evaluate_pinn(cand.net,tObs(mm),o.t_max_h); aa=state(3,:).';
        cs=base.slow.kappa_QA*base.scales.K_star_mM/cand.Gamma_A/ ...
            base.scales.lambda_ba*24450;
        pred=gain(mm).*cs.*aa; scale=max(.25,max(Cobs(trMask))-min(Cobs(trMask)));
        s=mean(((pred-Cobs(mm))/scale).^2)+o.selection_physics_weight*cand.history(end,3);
    end

    function w=weights(ramp)
        w=struct('data',o.lambda_data,'physics',o.lambda_physics*ramp, ...
            'ic',o.lambda_ic,'positive',o.lambda_positive, ...
            'gamma',o.lambda_gamma,'forcing',o.lambda_forcing*ramp);
    end
    function cand=decode(cand)
        z=double(gather(extractdata(cand.raw_parameters)));
        p=cand.lb+(cand.ub-cand.lb)./(1+exp(-z));
        if scenario=="baseline"; p(2)=p(1); end
        if string(o.gamma_mode)=="oracle"; p(5)=trueTheta(8);
        elseif string(o.gamma_mode)=="population"; p(5)=o.gamma_population; end
        cand.r_pre=p(1); cand.r_post=p(2); cand.tau_c=p(3); cand.k_r=p(4);
        cand.Gamma_A=p(5); cand.Gamma_prior=gammaPrior;
    end
end

function [loss,gradients]=network_loss(net,raw,varargin)
[loss,gradients,~]=pinn_model_loss(net,raw,varargin{:});
end
function z=logit(p); p=max(1e-6,min(1-1e-6,p)); z=log(p./(1-p)); end
