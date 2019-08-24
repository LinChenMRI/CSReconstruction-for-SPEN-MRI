function x = fnlCg(x0,params)
%-----------------------------------------------------------------------
%
% res = fnlCg(x0,params)
%
% implementation of a L1 penalized non linear conjugate gradient reconstruction
%
% The function solves the following problem:
%
% given k-space measurments y, and a fourier operator F the function 
% finds the image x that minimizes:
%
% Phi(x) = ||F*x - y||^2 + lambda1*|W*x|_1 + lambda2*TV(x) 
%
%
% the optimization method used is non linear conjugate gradient with fast&cheap backtracking
% line-search.
% 
% (c) Michael Lustig 2007
%-------------------------------------------------------------------------
x = x0;


% line search parameters
maxlsiter = params.lineSearchItnlim ;
gradToll = params.gradToll ;
alpha = params.lineSearchAlpha;     beta = params.lineSearchBeta;
t0 = params.lineSearchT0;
k = 0;
t = 1;

% copmute g0  = grad(Phi(x))

g0 = wGradient(x,params);

dx = -g0;


% iterations
while(1)

% backtracking line-search

	% pre-calculate values, such that it would be cheap to compute the objective
	% many times for efficient line-search
	[FTtx, FTtdx, Dtx, Dtdx,XFMtx,XFMtdx] = preobjective(x, dx, params);
	f0 = objective(FTtx, FTtdx, Dtx, Dtdx,XFMtx,XFMtdx, 0, params);
	t = t0;
        [f1, ERRobj, RMSerr]  =  objective(FTtx, FTtdx, Dtx, Dtdx,XFMtx,XFMtdx, t, params);
	
	lsiter = 0;

	while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:)))^2 && (lsiter<maxlsiter)
		lsiter = lsiter + 1;
		t = t * beta;
		[f1, ERRobj, RMSerr]  =  objective(FTtx, FTtdx, Dtx, Dtdx,XFMtx,XFMtdx, t, params);
	end

	if lsiter == maxlsiter
		disp('Reached max line search,.... not so good... might have a bug in operators. exiting... ');
		return;
	end

	% control the number of line searches by adapting the initial step search
	if lsiter > 2
		t0 = t0 * beta;
	end 
	
	if lsiter<1
		t0 = t0 / beta;
	end

	x = (x + t*dx);

	%--------- uncomment for debug purposes ------------------------	
% 	disp(sprintf('%d   , obj: %f, RMS: %f, L-S: %d', k,f1,RMSerr,lsiter));

	%---------------------------------------------------------------
	
    %conjugate gradient calculation
    
	g1 = wGradient(x,params);
	bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
	g0 = g1;
	dx =  - g1 + bk* dx;
	k = k + 1;
	
	%TODO: need to "think" of a "better" stopping criteria ;-)
	if (k > params.Itnlim) || (norm(dx(:)) < gradToll) 
		break;
	end

end


return;


function [FTtx, FTtdx, Dtx, Dtdx,XFMtx,XFMtdx] = preobjective(x, dx, params)

% precalculates transforms to make line search cheap

FTtx = params.FT*x;
FTtdx = params.FT*dx;

if params.TVWeight
    Dtx = params.TV*x;
    Dtdx = params.TV*dx;
else
    Dtx = 0;
    Dtdx = 0;
end

if params.xfmWeight
    XFMtx = params.XFM*x;
    XFMtdx = params.XFM*dx;
else
    XFMtx = 0;
    XFMtdx = 0;
end


function [res, obj, RMS] = objective(FTtx, FTtdx, Dtx, Dtdx,XFMtx,XFMtdx, t, params)
%calculated the objective function

p = params.pNorm;

obj = FTtx + t*FTtdx - params.data;

obj = obj(:)'*obj(:);

if params.TVWeight
    w = Dtx(:) + t*Dtdx(:);
    TV = (w.*conj(w)+params.l1Smooth).^(p/2); 
else
    TV = 0;
end

if params.xfmWeight
   w = XFMtx(:) + t*XFMtdx(:); 
   XFM = (w.*conj(w)+params.l1Smooth).^(p/2);
else
    XFM=0;
end



TV = sum(TV.*params.TVWeight(:));
XFM = sum(XFM.*params.xfmWeight(:));
RMS = sqrt(obj/sum(abs(params.data(:))>0));

res = obj + (TV) + (XFM) ;

function grad = wGradient(x,params)

gradXFM = 0;
gradTV = 0;

gradObj = gOBJ(x,params);
if params.xfmWeight
gradXFM = gXFM(x,params);
end
if params.TVWeight
gradTV = gTV(x,params);
end

grad = (gradObj +  params.xfmWeight.*gradXFM + params.TVWeight.*gradTV);



function gradObj = gOBJ(x,params)
% computes the gradient of the data consistency

gradObj = params.FT'*(params.FT*x - params.data);

% gradObj = (1-params.BW).*gradObj;    %  加入权重    修改了TV算子

gradObj = 2*gradObj ; 


function grad = gXFM(x,params)
% compute gradient of the L1 transform operator

p = params.pNorm;

Dx = params.XFM*x;

G = p*Dx.*(Dx.*conj(Dx) + params.l1Smooth).^(p/2-1);

grad = params.XFM'*G;

% grad = params.BW.*grad;   % 加入权重


function grad = gTV(x,params)
% compute gradient of TV operator

p = params.pNorm;

Dx = params.TV*x;
% 
% Dx(:,:,1) = params.BW.*Dx(:,:,1);   %  加入权重    修改了TV算子
% Dx(:,:,2) = params.BW.*Dx(:,:,2);

G = p*Dx.*(Dx.*conj(Dx) + params.l1Smooth).^(p/2-1);

grad = params.TV'*G;
