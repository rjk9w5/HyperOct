% Main driver for hypersonics Project. 
clear all
close all
clc
addpath('./Subroutines')
warning('off')
fpr = @(T) sutherlands(T, 1.716*10^-5, 110.6,273.1)*1005.4/sutherlands(T,0.0241, 194,273.1);
fmu = @(T) sutherlands(T, 1.716*10^-5, 110.6,273.1)*1005.4;
fk = @(T) sutherlands(T,0.0241, 194,273.1);
flr = @(Pr) Pr^(1/2);

deg2rad = pi/180;
% simple wedge assumption
defl = 18*deg2rad;
t_run = [538.103; 538.179; 538.734; 538.805];
alpha = [-5.012;5.540;-5.081;4.617]*deg2rad;
Pr = 0.715;
c = 0.0:0.01:0.3512;

FreeStream
fileID = fopen('cp_results.txt','w');
fileID2 = fopen('hf_ss_results.txt','w');
for i = 1:4
    fsst = state_init('Mach', M(i), 'Velocity', u(i), 'Density', r2(i), 'Temperature', T2(i), 'Enthalpy', 1005.4*T(i), 'Pressure', p2(i));
    % Newtonian
    cp_newtonian = 2*sin(defl+alpha(i))^2;

    % CPG Assumption
    [ds_cpg, B_cpg] = cpg_obs(fsst, defl+alpha(i), 1.4);
    cp_cpg = (ds_cpg.p-fsst.p)/q_inf(i);
    % Equillibrium Assumption
    [~, ~, ~, ~, ds_equil, B_equil] = equil_obs(fsst, defl+alpha(i));
    cp_equil = (ds_equil.p-fsst.p)/q_inf(i);

    (cp_newtonian - cp_equil)/cp_equil*100
    (cp_cpg - cp_equil)/cp_equil*100
    
    fprintf(fileID,'Newtonian\t\tCPG\t\t\tEquillibrium\n')
    fprintf(fileID,'%8.7f\t\t%8.7f\t%8.7f\n',cp_newtonian, cp_cpg, cp_equil)
    
    % (ds_equil.V - ds_cpg.V)/ds_equil.V
    % (ds_equil.M - ds_cpg.M)/ds_equil.M
    % (ds_equil.T - ds_cpg.T)/ds_equil.T
    % (ds_equil.p - ds_cpg.p)/ds_equil.p
    % (ds_equil.r - ds_cpg.r)/ds_equil.r
    % disp('\n')
    % Laminar Case
    est = ds_cpg;
    r0 = flr(Pr);
    Taw = Taw_function(est, r0);
    Tw = 339.5;
    % Reference Temperature State
        Tref = est.T + 0.5 * (Tw - est.T) + 0.22 * (Taw - est.T);
        refst = state_init('Temperature', Tref, 'Pressure', est.p, 'Density', est.p/(287*Tref));

        Reref = refst.r*est.V*c/sutherlands(refst.T, 1.716*10^-5, 110.6,273.1);
        m = fmu(refst.T)
        k = fk(refst.T)
        Prref = fpr(refst.T);

        ltwall = (.5*refst.r*est.V^2*0.664)./sqrt(Reref)/(1*10^3);
        lqdotwall = refst.r * est.V * est.cp * (Taw - Tw) * 0.332 ./ sqrt(Reref) * Prref^(-2/3)/(1*10^7);
        % size(refst.r * est.V * est.cp * (Taw - Taw) * 0.332 ./ sqrt(Reref) * Prref^(-2/3)/(1*10^7))

    lqdotwall(1) = qdot_stag(fsst, est, Tw, .001)/1e7;
    lqdotwall(1);
    ltwall(1) = 0;
    fprintf(fileID2,'Results for %2.0f\n',i)
    fprintf(fileID2,'0\t\t25\t\t50\t\t75\t\t100\n')
    n = max(size(c));
    id1 = 1;
    id2 = floor(n*.25);
    id3 = floor(n*.5);
    id4 = floor(n*.75);
    id5 = n;
    fprintf(fileID2,'%6.5f\t\t%6.5f\t\t%6.5f\t\t%6.5f\t\t%6.5f\n',lqdotwall(id1),lqdotwall(id2),lqdotwall(id3),lqdotwall(id4),lqdotwall(id5))
    fprintf(fileID2,'%6.5f\t\t%6.5f\t\t%6.5f\t\t%6.5f\t\t%6.5f\n\n',ltwall(id1),ltwall(id2),ltwall(id3),ltwall(id4),ltwall(id5))

    time = num2str(t_run(i));
	figure(1)
	    subplot(2,2,i)
	    plot(c,lqdotwall)
	    title(['Laminar Heat Flux at wall @ t = ', time])
	    xlabel('Chord [m]')
	    ylabel('q_{w} [kW/cm^2]')

    figure(2)
        subplot(2,2,i)
        plot(c,ltwall)
        title(['Laminar Shear Stress at wall @ t = ', time])
        xlabel('Chord [m]')
        ylabel('\tau_{w} [kPa]')
end %for
fclose(fileID)
fclose(fileID2)
rmpath('./Subroutines')
