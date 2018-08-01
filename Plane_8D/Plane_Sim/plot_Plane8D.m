%% Check Plot

%Dubins4D check
% figure(1)
% plot3(Xp(:,1),Xp(:,2),Xp(:,3),'b-','linewidth',2.0);
% grid on

% figure()
% plot(t_ref(1:end-1),Rp(1:end-1,:)-R,'linewidth',2);
% grid on
% xlabel('Time [s]')
% title('r_p vs r')

figure()
plot(t_ref(1:T_steps-1),dVdt(:,1),'r-','linewidth',2); hold on
plot(t_ref(1:T_steps-1),dVdt(:,2),'b-','linewidth',2);
grid on
xlabel('Time [s]')
title('dV/dt: true vs approx');

%% Trajectory plot

%Trajectories
figure(fig_FMT)
grid on
hold on
plot3(X(:,1),X(:,2),X(:,3),'k-','linewidth',2.0);
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
for i = 1:(0.5/dT):length(t_ref)
    proj_Ellipse([1:3],E_sol,1,X_ref(i,1:3)',30,'b',X_ref(i,4));
end

%% Rest of the plots

%States
figure()
title('Relative states');

subplot(2,2,1)
plot(t_ref(1:end-1),R(:,1:3),'r-','linewidth',2);
xlabel('time [s]');
ylabel('pos');

subplot(2,2,2)
plot(t_ref(1:end-1),R(:,5),'r-','linewidth',2);
xlabel('time [s]');
ylabel('v');

subplot(2,2,3)
plot(t_ref(1:end-1),R(:,[4;6;7;8]),'linewidth',2);
xlabel('time [s]');
ylabel('angles');
legend('\psi_r','\gamma','\phi','\alpha');

%Control
figure()
subplot(2,2,1)
hold on
plot(t_ref(1:end-1),U(:,1),'r-','linewidth',2);
line([t_ref(1),t_ref(end)],[-u_lim(1),-u_lim(1)],'color','k','linewidth',2);
line([t_ref(1),t_ref(end)],[u_lim(1),u_lim(1)],'color','k','linewidth',2);
grid on; xlabel('time [s]'); ylabel('a');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

subplot(2,2,2)
hold on
plot(t_ref(1:end-1),U(:,2),'r-','linewidth',2);
line([t_ref(1),t_ref(end)],[-u_lim(2),-u_lim(2)],'color','k','linewidth',2);
line([t_ref(1),t_ref(end)],[u_lim(2),u_lim(2)],'color','k','linewidth',2);
grid on; xlabel('time [s]'); ylabel('$\dot{\phi}$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

subplot(2,2,3)
hold on
plot(t_ref(1:end-1),U(:,3),'r-','linewidth',2);
line([t_ref(1),t_ref(end)],[-u_lim(3),-u_lim(3)],'color','k','linewidth',2);
line([t_ref(1),t_ref(end)],[u_lim(3),u_lim(3)],'color','k','linewidth',2);
grid on; xlabel('time [s]'); ylabel('$\dot{\alpha}$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%V,rho
figure()
plot(t_ref(1:end-1),V_bnds,'linewidth',2);
line([0,t_ref(end)],[rho_sol,rho_sol],'color','k','linewidth',2);
grid on
xlabel('Time [s]');
ylabel('V(r)');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
set(gcf,'Color','w');

%% Create video

keyboard;
create_Plane_vid();