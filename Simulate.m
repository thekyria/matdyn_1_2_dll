function [angles, speeds, eq_tr, ed_tr, efd, PM, voltages, ...
          stepsize, errest, time, success] = ...
  Simulate(baseMVAIn, busIn, genIn, branchIn, areasIn, genCostIn, ...
           SSqlim, SSdc, SSalg, SStol, SSmax_it, ...
           genDyn, excDyn, govDyn, ...
           TDfreq, TDstep, TDstoptime, TDmethod, TDtol, TDstepMin, TDstepMax, ...
           event, buschange, linechange)
%SIMULATE Runs a dynamic simulation
%  [ANGLES, SPEEDS, EQ_TR, ED_TR, EFD, PM, VOLTAGES, STEPSIZE, ERREST,
%   TIME, SUCCESS] = ...
%    rundyn(BASEMVAIN, BUSIN, GENIN, BRANCHIN, AREASIN, GENCOSTIN, ...
%    SSQLIM, SSDC, SSALG, SSTOL, SSMAX_IT, ...
%    GENDYN, EXCDYN, GOVDYN, ...
%    TDFREQ, TDSTEP, TDSTOPTIME, TDMETHOD, TDTOM, TDSTEPMIN, TDSTEPMAX, ...
%    EVENT, BUSCHANGE, LINECHANGE)
% 
%  INPUTS
%    baseMVA, bus, gen, branch, areas, gencost = as per loadflow()
%    SSqlim = enforce gen reactive power limits at expense of |V|
%                 0 - do NOT enforce limits
%                 1 - enforce limits, simultaneous bus type conversion
%                 2 - enforce limits, one-at-a-time bus type conversion
%    SSdc = DC modeling for power flow & OPF
%                 0 - use AC formulation & corresponding algorithm opts
%                 1 - use DC formulation, ignore AC algorithm options
%    SSalg = AC power flow algorithm
%                 1 - Newton's method            
%                 2 - Fast-Decoupled (XB version)
%                 3 - Fast-Decoupled (BX version)
%                 4 - Gauss-Seidel               
%    SStol = termination tolerance on per unit P & Q mismatch (eg. 1e-8)
%    SSmax_it = maximum number of iteration
%                 10 - for Newton's method
%                 30 - for Fast-Decoupled
%                 1000 - for Gauss-Seidel
%    TDmethod = integration method to use
%                 1 - modified Euler
%                 2 - Runge-Kutta
%                 3 - Runge-Kutta Fehlberg
%                 4 - Higham and Hall
%                 5 - Modified Euler with interface error control
%    TDtol = Specify the tolerance of the error. This argument is only used
%                 for the Runge-Kutta Fehlberg and Higham and Hall methods
%    TDminStepSize = sets the step size; for adaptive step size algorithms
%                 (Runge-Kutta Fehlberg & Higham and Hall methods) this is 
%                 the minimum step size
%    TDmaxStepSize = sets the maximum step size; used only in adaptive step
%                 size algorithms (Runge-Kutta Fehlberg & Higham and Hall)
% 
%  OUTPUTS
%    angles = generator angles
%    speeds = generator speeds
%    eq_tr = q component of transient voltage behind reactance
%    ed_tr = d component of transient voltage behind reactance
%    efd = Excitation voltage
%    PM = mechanical power
%    voltages = bus voltages
%    stepsize = step size integration method
%    errest = estimation of integration error
%    time = time points
%    success = 0: failure / 1: success

%% Constants
STEADYSTATETOL = 1e-6;

%% define named indices into bus, gen, branch matrices
% bus related indexes
PD = 3;
QD = 4;
VM = 8;
VA = 9;
% generator related indexes
GEN_BUS = 1;
GEN_STATUS = 8;

%% initialize output arguments
angles = 0;
speeds = 0;
eq_tr = 0;
ed_tr = 0;
efd = 0;
PM = 0;
voltages = 0;
stepsize = 0;
errest = 0;
time = 0;
success = 0;
global freq

%% input argument validation
if nargin < 24
	return;
end

% debug print
display(baseMVAIn); display(busIn); display(genIn); display(branchIn);
display(areasIn); display(genCostIn);
%-
display(SSqlim); display(SSdc); display(SSalg); display(SStol);
display(SSmax_it);
%-
display(genDyn); display(excDyn); display(govDyn);
%-
display(TDfreq); display(TDstep); display(TDstoptime); display(TDmethod);
display(TDtol); display(TDstepMin); display(TDstepMax);
%-
display(event); display(buschange); display(linechange);

%% Load all data
freq = TDfreq;
step = TDstep;
stoptime = TDstoptime;

% ----- Ext2Int indexing -----
[i2e, busIn, genIn, branchIn, areasIn] = ...
  ext2int(busIn, genIn, branchIn, areasIn);

% ----- Load generator data ----- 
Pgen0 = genDyn;
genmodel = Pgen0(:,1);
% Define generator models
d = (1:length(genmodel))';
type1 = d(genmodel==1);
type2 = d(genmodel==2);
% Check transient saliency
xd_tr = Pgen0(type2,8);
xq_tr = Pgen0(type2,9);
if sum(xd_tr~=xq_tr)>=1
  Pgen0(type2,9) = Pgen0(type2,8);
end
excmodel = Pgen0(:,2);
govmodel = Pgen0(:,3);

% ----- Load exciter data -----
Pexc0 = excDyn;
for i = 1:length(excDyn(1,:))
    Pexc0(excDyn(:,1),i) = excDyn(:,i);
end

% ----- Load governor data -----
Pgov0 = govDyn;
for i = 1:length(govDyn(1,:))
    Pgov0(govDyn(:,1),i) = govDyn(:,i);
end

% ----- Load event data -----
% buschange and linechange matrices are padded to size = size(event,1)
type1 = buschange;                 % input arg 'buschange'
type2 = linechange;                % input arg 'linechange'
buschange=zeros(size(event,1),4);  % init internal variable 'buschange'
linechange=zeros(size(event,1),4); % init internal variable 'linechange'
i1=1; i2=1;                        % counters
for i=1:size(event,1)
  if event(i,2)==1
    buschange(i,:) = type1(i1,:);
    i1=i1+1;     
  elseif event(i,2)==2
    linechange(i,:) = type2(i2,:);
    i2=i2+1;
  end
end

%% Initialization: Power Flow 
% Run power flow
[baseMVAIn, busIn, gen, branchIn, success] = ...
  loadflow(baseMVAIn, busIn, genIn, branchIn, areasIn, genCostIn, ...
           SSqlim, SSdc, SSalg, SStol, SSmax_it);
if ~success
  % Loadflow failed
  fprintf('Simulate(): Error: Power flow did not converge. Exiting...\n')
  return;
end

% Bus voltages in complex form
U0=busIn(:,VM).*(cos(busIn(:,VA)*pi/180) + 1i.*sin(busIn(:,VA)*pi/180));
U00=U0;

% Get generator info
on = gen(:, GEN_STATUS) > 0;     % which generators are on?
gbus = gen(on, GEN_BUS);         % what buses are they at?
ngen = length(gbus);             % number of generators

nbus = length(U0);               % number of buses

%% Construct augmented Ybus 
Pl=busIn(:,PD)./baseMVAIn;                  % load power
Ql=busIn(:,QD)./baseMVAIn;

xd_tr = zeros(ngen,1);
xd_tr(genmodel==2) = Pgen0(genmodel==2,8);  % 4th order mdl: xd_tr column 8
xd_tr(genmodel==1) = Pgen0(genmodel==1,7);  % classical mdl: xd_tr column 7

[Ly, Uy, Py] = AugYbus(baseMVAIn, busIn, branchIn, xd_tr, gbus, Pl, Ql, U0);


%% Calculate Initial machine state
% Xgen - array of generator variables: ngen x 4
% columns
% 1          2                      3         4 
% delta[rad] omega(absolute)[rad/s] Eq_tr[pu] Ed_tr[pu]
[Efd0, Xgen0] = GeneratorInit(Pgen0, U0(gbus), gen, baseMVAIn, genmodel);

omega0 = Xgen0(:,2);

[Id0,Iq0,Pe0] = MachineCurrents(Xgen0, Pgen0, U0(gbus), genmodel);
Vgen0 = [Id0, Iq0, Pe0];


%% Exciter initial conditions
Vexc0 = abs(U0(gbus));
[Xexc0,Pexc0] = ExciterInit(Efd0, Pexc0, Vexc0, excmodel);


%% Governor initial conditions
Pm0 = Pe0;
[Xgov0, Pgov0] = GovernorInit(Pm0, Pgov0, omega0, govmodel);
Vgov0 = omega0;


%% Check Steady-state
Fexc0 = Exciter(Xexc0, Pexc0, Vexc0, excmodel);
Fgov0 = Governor(Xgov0, Pgov0, Vgov0, govmodel);
Fgen0 = Generator(Xgen0, Xexc0, Xgov0, Pgen0, Vgen0, genmodel);

% Check Generator Steady-state
if sum(sum(abs(Fgen0))) > STEADYSTATETOL
  fprintf('> Error: Generator not in steady-state\n> Exiting...\n')
  return;
end
% Check Exciter Steady-state
if sum(sum(abs(Fexc0))) > STEADYSTATETOL
	fprintf('> Error: Exciter not in steady-state\n> Exiting...\n')
	return;
end
% Check Governor Steady-state
if sum(sum(abs(Fgov0))) > STEADYSTATETOL
  fprintf('> Error: Governor not in steady-state\n> Exiting...\n')
  return;
end

%% Initialization of main stability loop
t = -0.02; % simulate 0.02s without applying events
err = 0;
failed = 0;
eulerfailed = 0;

if TDmethod==3 || TDmethod==4 % 3: Runge-Kutta Fehlberg, 4: Higham and Hall
  step = TDstepMin;
end

ev=1;
eventhappened = false;
i=0;

%% Allocate memory for variables
chunk = 5000; % chunk of memory allocation

time = zeros(chunk,1);
time(1,:) = t;
errest = zeros(chunk,1);
errest(1,:) = err;
stepsize = zeros(chunk,1);
stepsize(1,:) = step;

% System variables
voltages = zeros(chunk, length(U0)); voltages(1,:) = U0.';

% Generator
angles = zeros(chunk,ngen); angles(1,:) = Xgen0(:,1).*180./pi;
speeds = zeros(chunk,ngen); speeds(1,:) = Xgen0(:,2)./(2.*pi.*freq);
eq_tr = zeros(chunk,ngen); eq_tr(1,:) = Xgen0(:,3);
ed_tr = zeros(chunk,ngen); ed_tr(1,:) = Xgen0(:,4);

% Exciter and governor
efd = zeros(chunk,ngen); efd(1,:) = Efd0(:,1);  
PM = zeros(chunk,ngen); PM(1,:) = Pm0(:,1);


%% Main stability loop
while t < stoptime + step

  %% Output    
  i=i+1;
%   fprintf('Simulate(): %6.2f%% completed\n', t/stoptime*100)

  %% Numerical Method
  switch TDmethod
    case 1 % Modified Euler
        [Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, ...
        U0, t, newstepsize] = ...
          ModifiedEuler(t, Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, ...
            Xgov0, Pgov0, Vgov0, Ly, Uy, Py, gbus, genmodel, excmodel, ...
            govmodel, step);
          
    case 2 % Runge-Kutta
        [Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, ...
        U0, t, newstepsize] = ...
          RungeKutta(t, Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, ...
            Xgov0, Pgov0, Vgov0, Ly, Uy, Py, gbus, genmodel, excmodel, ...
            govmodel, step);
          
    case 3 % Runge-Kutta Fehlberg
        [Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, ...
        U0, err, failed, t, newstepsize] = ...
          RungeKuttaFehlberg(t, Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, ...
            Vexc0, Xgov0, Pgov0, Vgov0, U0, Ly, Uy, Py, gbus, genmodel, ...
            excmodel, govmodel, TDtol, TDstepMax, step);
          
    case 4 % Higham and Hall
        [Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, ...
        U0, err, failed, t, newstepsize] = ...
          RungeKuttaHighamHall(t, Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, ...
          Vexc0, Xgov0, Pgov0, Vgov0, U0, Ly, Uy, Py, gbus, genmodel, ...
          excmodel, govmodel, TDtol, TDstepMax, step);                 
    case 5 % Modified Euler with interface error control
        [Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, ...
        U0, t, eulerfailed, newstepsize] = ...
          ModifiedEuler2(t, Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, ...
          Xgov0, Pgov0, Vgov0, Ly, Uy, Py, gbus, genmodel, excmodel, ...
          govmodel, step);
  end

  if eulerfailed
    fprintf('Simulate(): No solution found. ')
    fprintf('Try lowering tol or increasing max_it count in ModifEuler2. ')
    fprintf('Exiting... \n')
    return;
  end        
    
  if failed
    t = t-step;
  end

  % End exactly at stop time
  if t + newstepsize > stoptime
    newstepsize = stoptime - t;
  elseif step < TDstepMin
    fprintf('Simulate(): Error: No solution found with minimum step size')
    fprintf('Exiting... \n')
    return;
  end;
    
    
  %% Allocate new memory chunk if matrices are full
  if i>size(time,1)
    angles = [angles;zeros(chunk,ngen)];
    speeds = [speeds;zeros(chunk,ngen)];
    eq_tr = [eq_tr;zeros(chunk,ngen)];
    ed_tr = [ed_tr;zeros(chunk,ngen)];
    efd = [efd; zeros(chunk,ngen)];
    PM = [PM; zeros(chunk,ngen)];
    voltages = [voltages; zeros(chunk,length(U0))];
    stepsize = [stepsize; zeros(chunk,1)];
    errest = [errest; zeros(chunk,1)];
    time = [time; zeros(chunk,1)];
  end

    
  %% Save values
  stepsize(i,:) = step.';
  errest(i,:) = err.';
  time(i,:) = t;

  voltages(i,:) = U0.';

  % exc
  efd(i,:) = Xexc0(:,1).*(genmodel>1); % Set Efd to zero when using
                                       % classical generator model  

  % gov
  PM(i,:) = Xgov0(:,1);
   
  % gen
  angles(i,:) = Xgen0(:,1).*180./pi;
  speeds(i,:) = Xgen0(:,2)./(2.*pi.*freq);
  eq_tr(i,:) = Xgen0(:,3);
  ed_tr(i,:) = Xgen0(:,4);   
    
    
  %% Adapt step size if event will occur in next step
  if    ~isempty(event)     ...
     && ev <= size(event,1) ...
     && (TDmethod == 3 || TDmethod == 4)
    if t + newstepsize >= event(ev,1)
      if event(ev,1) - t < newstepsize
        newstepsize = event(ev,1) - t;
      end
    end
  end 

    
  %% Check for events
  if ~isempty(event) && ev <= size(event,1)    

    for k=ev:size(event,1)                    % cycle through all events ..   
      if abs(t-event(ev,1))>10*eps ||  ev > size(event,1) %. that happen @t
        break;
      else
        eventhappened = true;
      end

      switch event(ev,2)
        case 1
          busIn(buschange(ev,2),buschange(ev,3)) = buschange(ev,4);
        case 2
          branchIn(linechange(ev,2),linechange(ev,3)) = linechange(ev,4);
      end
      ev=ev+1;
    end          

    if eventhappened
      % Refactorise
      [Ly, Uy, Py] = AugYbus(baseMVAIn, busIn, branchIn, xd_tr, gbus, ...
                             busIn(:,PD)./baseMVAIn, ...
                             busIn(:,QD)./baseMVAIn, U00);
      U0 = SolveNetwork(Xgen0, Pgen0, Ly, Uy, Py, gbus, genmodel);            

      [Id0,Iq0,Pe0] = MachineCurrents(Xgen0, Pgen0, U0(gbus), genmodel);
      Vgen0 = [Id0,Iq0,Pe0];
      Vexc0 = abs(U0(gbus));

      % decrease stepsize after event occured
      if TDmethod==3 || TDmethod==4
        newstepsize = TDstepMin;
      end

      i=i+1; % if event occurs, save values at t- and t+

      %% Save values
      stepsize(i,:) = step.';
      errest(i,:) = err.';
      time(i,:) = t;

      voltages(i,:) = U0.';

      % exc
      efd(i,:) = Xexc0(:,1).*(genmodel>1); % Set Efd to zero when using
                                           % classical generator model  

      % gov
      PM(i,:) = Xgov0(:,1);

      % gen
      angles(i,:) = Xgen0(:,1).*180./pi;
      speeds(i,:) = Xgen0(:,2)./(2.*pi.*freq);
      eq_tr(i,:) = Xgen0(:,3);
      ed_tr(i,:) = Xgen0(:,4);

      eventhappened = false;
    end
  end

    
  %% Advance time    
  step = newstepsize;
  t = t + step;
    
end % end of main stability loop

%% Save only the first i elements
angles = angles(1:i,:);
speeds = speeds(1:i,:);
eq_tr = eq_tr(1:i,:);
ed_tr = ed_tr(1:i,:);

efd = efd(1:i,:);
PM = PM(1:i,:);

voltages = voltages(1:i,:);

stepsize = stepsize(1:i,:);
errest = errest(1:i,:);
time = time(1:i,:);

return;
