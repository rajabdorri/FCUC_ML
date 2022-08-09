
***************************************************************
*** PARAMETERS
***************************************************************
*$set season spring
*$set day d1
*$set last_day d0
$include Input_Data_for_UC_Lapalma.gms


***************************************************************
*** VARIABLES
***************************************************************

variable obj objective function variable

variable c_aux(t) auxilliary variable

variable co(t,i) operation cost in each time period

positive variable g(t,i) generator outputs

positive variable g_lin(t,i,b) generator block outputs

binary variable suc(t,i,j) start up cost

binary variable x(t,i) binary variable equal to 1 if generator is producing, and 0 otherwise

binary variable y(t,i) binary variable equal to 1 if generator is start-up, and 0 otherwise

binary variable z(t,i) binary variable equal to 1 if generator is shut-down, and 0 otherwise

positive variable curt(t) wind curtailment

positive variable ll(t) unserved load

positive variable inertia_AO(t,i) inertia

positive variable k_AO(t,i) k

positive variable r_AO(t,i) r

parameter k0,kh,kk,kp,kr;
*$include %method%_coefficients.txt


*logistic regression
*k0=-7.53970158;
*kh=-0.63260783;
*kk=0.10089387;
*kp=-4.71837262;
*kr=0.86408313;

*SVM, c=1
*k0=-75.96200017;
*kh=-4.47285656;
*kk=0.91201327;
*kp=-61.19083931;
*kr=11.72095723;

*SVM, c=0.1
*k0=-10.29703624;
*kh=-0.64044481;
*kk=0.12205558;
*kp=-7.39730728;
*kr=1.33332431;

*SVM, c=10
*k0=-74.9710448;
*kh=-4.38458316;
*kk=0.89111734;
*kp=-59.61915712;
*kr=11.4237979;

*Perceptron
*k0=-3.8716;
*kh=-2.499;
*kk=0.5470;
*kp=-42.4038;
*kr=5.4978;

***************************************************************
*** EQUATION DECLARATION
***************************************************************

equations

cost objective function
cost_aux(t) auxilliary equation
bin_set1(t,i) setting start-up binary variables
bin_set10(t,i) setting start-up binary variables
bin_set2(t,i) setting start-up binary variables
gen_sum(t,i) summing the generation of blocks per generator
gen_min(t,i) genertor minimum output
cost_sum(t,i) generation cost summation
block_output(t,i,b) limiting the output of each generator block
min_updown_1(t,i) minimum updown time constraint 1
min_updown_2(t,i) minimum updown time constraint 2
min_updown_3(t,i) minimum updown time constraint 3
ramp_limit_min(t,i) ramp-down limit
ramp_limit_max(t,i) ramp-up limit
ramp_limit_min_1(i) ramp-down limit for the first time period
ramp_limit_max_1(i) ramp-up limit for the first time period
start_up_cost1(t,i,j) stairwise linear cost function - equation 1
start_up_cost2(t,i) stairwise linear cost function - equation 2
power_balance(t) power balance for each bus
curt_cons(t) wind curtailment constraint
ML_const(t,i) reserve requirement
inertia_eq(t,i)
k_eq(t,i)
r_eq(t,i)
qss(t,i)
rocof(t,i)
;

***************************************************************
*** SETTINGS
***************************************************************
scalar max_rocof /2.5/;
scalar qss_max /0.5/;
*needed for running twice through the same set in a single equation
alias (t,tt);
alias (i,ii);

***************************************************************
*** EQUATIONS
***************************************************************

cost..
         obj =e= sum(t,c_aux(t));

cost_aux(t)..
         c_aux(t) =e= sum(i,co(t,i))+curt(t)*ws_penalty
*+sum(s,ll(t,s))*voll;
         ;

bin_set1(t,i)$(ord(t) gt 1)..
         y(t,i) - z(t,i) =e= x(t,i) - x(t-1,i);

bin_set10(t,i)$(ord(t) = 1)..
         y(t,i) - z(t,i) =e= x(t,i) - onoff_t0(i);

bin_set2(t,i)..
         y(t,i) + z(t,i) =l= 1;

cost_sum(t,i)..
         co(t,i) =e= a(i)*x(t,i) + sum(b,g_lin(t,i,b)*k(i,b)) + sum(j,suc_sw(i,j)*suc(t,i,j));

gen_sum(t,i)..
         g(t,i) =e= sum(b,g_lin(t,i,b));

gen_min(t,i)..
         g(t,i) =g= g_min(i)*x(t,i);

block_output(t,i,b)..
         g_lin(t,i,b) =l= g_max(i,b)*x(t,i);

min_updown_1(t,i)$(L_up_min(i)+L_down_min(i) gt 0 and ord(t) le L_up_min(i)+L_down_min(i))..
         x(t,i) =e= onoff_t0(i);

min_updown_2(t,i)..
         sum(tt$(ord(tt) ge ord(t)-g_up(i)+1 and ord(tt) le ord(t)),y(tt,i)) =l= x(t,i);

min_updown_3(t,i)..
         sum(tt$(ord(tt) ge ord(t)-g_down(i)+1 and ord(tt) le ord(t)),z(tt,i)) =l= 1-x(t,i);

ramp_limit_min(t,i)$(ord(t) gt 1)..
         -ramp_down(i) =l= g(t,i) - g(t-1,i);

ramp_limit_max(t,i)$(ord(t) gt 1)..
         ramp_up(i) =g= g(t,i) - g(t-1,i);

ramp_limit_min_1(i)..
         -ramp_down(i) =l= g('t1',i) - g_0(i);

ramp_limit_max_1(i)..
         ramp_up(i) =g= g('t1',i) - g_0(i);

start_up_cost1(t,i,j)..
         suc(t,i,j) =l= sum(tt$(ord(tt) lt ord(t) and ord(tt) ge suc_sl(i,j) and ord(tt) le suc_sl(i,j+1)-1),z(t-ord(j),i))+
         1$(ord(j) lt card(j) and count_off_init(i)+ord(t)-1 ge suc_sl(i,j) and count_off_init(i)+ord(t)-1 lt suc_sl(i,j+1))+
         1$(ord(j) = card(j) and count_off_init(i)+ord(t)-1 ge suc_sl(i,j));

start_up_cost2(t,i)..
         sum(j,suc(t,i,j)) =e= y(t,i);

power_balance(t)..
         sum(i,g(t,i)) + RES_det(t,"%day%","%season%")-curt(t) =e= d(t,"%day%","%season%")
*-ll(t,s);
         ;

curt_cons(t)..
         curt(t) =l= RES_det(t,"%day%","%season%");
         
inertia_eq(t,i)..
        inertia_AO(t,i)=e=sum(ii$(ord(ii) <> ord(i)),inertia(ii)*mbase(ii)*x(t,ii));
        
k_eq(t,i)..
        k_AO(t,i)=e=sum(ii$(ord(ii) <> ord(i)),kg(ii)*x(t,ii)*mbase(ii));
        
r_eq(t,i)..
        r_AO(t,i)=e=sum(ii$(ord(ii) <> ord(i)),sum(b,g_max(ii,b))*x(t,ii)-g(t,ii));
        
qss(t,i)..
        r_AO(t,i)=g=g(t,i)-0.01*d(t,"%day%","%season%")*qss_max;
        
rocof(t,i)..
        2*inertia_AO(t,i)*max_rocof/50=g=g(t,i);
         
ML_const(t,i).. k0 + kh*inertia_AO(t,i)
                        + kk*k_AO(t,i)
                        + kp*g(t,i)
                        + kr*r_AO(t,i)
*                        + 100*(1-x(t,i))
                        =g=0; 

***************************************************************
*** SOLVE
***************************************************************

model ep /
cost,
cost_aux,
bin_set1,
bin_set10,
bin_set2,
gen_sum,
gen_min,
cost_sum,
block_output,
min_updown_1,
min_updown_2,
min_updown_3,
ramp_limit_min,
ramp_limit_max,
ramp_limit_min_1,
ramp_limit_max_1,
start_up_cost1,
start_up_cost2,
power_balance,
curt_cons,
*ML_const,
inertia_eq,
k_eq,
r_eq,
qss,
rocof
/;

option reslim = 1200;
option Savepoint=1;
option optcr=0.01;
ep.optfile = 1;

ep.limrow =0;
ep.limcol =0;

file opt cplex option file /cplex.opt/;
put opt;
put 'threads 4'/;
put 'miptrace _UC_deterministic_4_04_0.csv'/;
putclose;


solve ep using mip minimizing obj;

parameter overall_suc, output_cost,time_elapsed,generation(t,i),online_reserve(t,i),online_inertia(t,i),online_k(t,i);

overall_suc=sum((t,i,j), suc_sw(i,j)*suc.l(t,i,j));
output_cost=sum((t,i), a(i)*x.l(t,i) + sum(b,g_lin.l(t,i,b)*k(i,b)));
time_elapsed = ep.etSolver;
generation(t,i)=g.l(t,i)+eps;
online_reserve(t,i)=sum(b,g_max(i,b))*x.l(t,i)-g.l(t,i)+eps;
online_inertia(t,i)=inertia_AO.l(t,i);
online_k(t,i)=k_AO.l(t,i);



Execute_Unload 'DM_ML_%method%_%season%_%day%.gdx'
*, obj, generation, online_reserve,online_inertia,xone,xtwo, lambda_p, lambda_xone,lambda_xtwo,c_aux, c, g,x, y, z, curt, overall_suc,output_cost,time_elapsed;
*$ontext
execute 'gdxxrw.exe DM_ML_%method%_%season%_%day%.gdx par=generation rng=gen!a1:zz28000'
execute 'gdxxrw.exe DM_ML_%method%_%season%_%day%.gdx par=online_reserve rng=reserve!a1:zz28000'
execute 'gdxxrw.exe DM_ML_%method%_%season%_%day%.gdx par=online_inertia rng=inertia!a1:zz28000'
execute 'gdxxrw.exe DM_ML_%method%_%season%_%day%.gdx par=online_k rng=k!a1:zz28000'