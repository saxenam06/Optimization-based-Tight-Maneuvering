function Draw_MPC_Goal_stop_Mv_obs_param(t,xx,xx1,ux,x_ref,N,ref_obs,ref_traj,CG,lf,lr,w)

set(0,'DefaultAxesFontName', 'Times NewRoman')
set(0,'DefaultAxesFontSize', 12)

line_width=1.5;
fontsize_labels=14;

x_r_1=[];
y_r_1=[];

x_obs_r_1=[];
y_obs_r_1=[];
x_obs_r_2=[];
y_obs_r_2=[];
x_obs_r_3=[];
y_obs_r_3=[];


figure(500)

set(gcf, 'PaperPositionMode', 'auto')
set(gcf, 'Color', 'w')
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1.75 0.30]);

tic
cut=0;
F=size(xx,2)-cut;
P=fix(F/2);
r=diff(fix(linspace(0,F,P+1)));
X=1:size(xx,2)-cut;
C=mat2cell(X,1,r)

for k=1:P
    h_t=10; w_t=5;
    
    plot(x_ref(1),x_ref(2),'*b');
    hold on;
    x1=xx(1, fix(median(C{k})),1); y1=xx(2,fix(median(C{k})),1); th1=xx(3,fix(median(C{k})),1); 
    x_r_1=[x_r_1; x1];
    y_r_1= [y_r_1; y1];
    plot(x_r_1, y_r_1, '-g','linewidth', line_width); hold on
    plot(x1,y1,'--g');
    
    R=[cos(th1) sin(th1) 0;-sin(th1) cos(th1) 0;0  0 1]';
    T=[1 0 x1;0 1 y1;0 0 1];
    pr1=T*R*[lf;w/2;1];
    pr2=T*R*[-lf;w/2;1];
    pr3=T*R*[-lf;-w/2;1];
    pr4=T*R*[lf;-w/2;1];
    
    fill([pr1(1) pr2(1) pr3(1) pr4(1)], [pr1(2) pr2(2) pr3(2) pr4(2)],'g');
    
     x_obs1=mean(xx(9,C{k},1)); y_obs1=mean(xx(10,C{k},1));
     x_obs_r_1=[x_obs_r_1 x_obs1];
     y_obs_r_1 = [ y_obs_r_1 y_obs1];
     plot(x_obs_r_1, y_obs_r_1, '-b','linewidth', line_width); hold on
     plot(x_obs1,y_obs1,'--b');
     
     x_obs1_rct=[x_obs1+lf, x_obs1+lf x_obs1-lr x_obs1-lr];
     y_obs1_rct=[y_obs1+(w/2), y_obs1-(w/2) y_obs1-(w/2) y_obs1+(w/2)];
     fill(x_obs1_rct, y_obs1_rct,'b');
     
     x_obs2=mean(xx(11,C{k},1)); y_obs2=mean(xx(12,C{k},1));
     x_obs_r_2=[x_obs_r_2 x_obs2];
     y_obs_r_2 = [ y_obs_r_2 y_obs2];
     plot(x_obs_r_2, y_obs_r_2, '-b','linewidth', line_width); hold on
     plot(x_obs2,y_obs2,'--b');
     
     x_obs2_rct=[x_obs2+lf, x_obs2+lf x_obs2-lr x_obs2-lr];
     y_obs2_rct=[y_obs2+(w/2), y_obs2-(w/2) y_obs2-(w/2) y_obs2+(w/2)];
     fill(x_obs2_rct, y_obs2_rct,'b');
     
     x_obs3=mean(xx(13,C{k},1)); y_obs3=mean(xx(14,C{k},1));
     x_obs_r_3=[x_obs_r_3 x_obs3];
     y_obs_r_3 = [ y_obs_r_3 y_obs3];
     plot(x_obs_r_3, y_obs_r_3, '-b','linewidth', line_width); hold on
     plot(x_obs3,y_obs3,'--b');
     
     x_obs3_rct=[x_obs3+lf, x_obs3+lf x_obs3-lr x_obs3-lr];
     y_obs3_rct=[y_obs3+(w/2), y_obs3-(w/2) y_obs3-(w/2) y_obs3+(w/2)];
     fill(x_obs3_rct, y_obs3_rct,'b');
     
     
     if k<size(xx,2)
         plot(xx1(2:N+1,1,fix(median(C{k}))), xx1(2:N+1,2,fix(median(C{k}))),'r--*')
         
     end
     
     plot(ref_obs{1}(:,1),ref_obs{1}(:,2),'-k','linewidth',1.25);
     plot(ref_obs{2}(:,1),ref_obs{2}(:,2),'-k','linewidth',1.25);
     
     plot(ref_traj(1,:),ref_traj(2,:),'.m');
     plot(ref_traj(1,:),ref_traj(2,:)+1.75,'--c');          
     
    hold off
    
    ylabel('$y$-position (m)','interpreter','latex','FontSize',fontsize_labels)
    xlabel('$x$-position (m)','interpreter','latex','FontSize',fontsize_labels)
    axis([0 130 -4 4])
    pause(0.1)
    box on;
    grid on
    %aviobj = addframe(aviobj,gcf);
    drawnow
    % for video generation
    FIG(k) = getframe(gcf); % to get the current frame
end
toc

video = VideoWriter('Scenario_Tight_Evasive_Maneuvering.mp4','MPEG-4');
video.FrameRate = 5;  % (frames per second) this number depends on the sampling time and the number of frames you have
open(video)
writeVideo(video,FIG)
close (video)

figure(1)
subplot(3,3,2)
plot(xx(1,1:1:end-1-cut),atan(((lr./(lf+lr))*tan(xx(8,1:1:end-1-cut))))*180/pi,'-k','linewidth',1.5);
ylabel('beta angle (deg)')
xlabel('Longitudinal Distance (m)')
grid on
subplot(3,3,3)
plot(xx(1,1:1:end-1-cut),xx(3,1:1:end-1-cut)*180/pi,'-k','linewidth',1.5);
ylabel('Yaw angle (deg)')
xlabel('Longitudinal Distance (m)')
grid on
subplot(3,3,1)
plot(xx(1,1:1:end-1-cut),ux(2,1:1:end-1-cut)*180/pi,'-k','linewidth',1.5);
ylabel('Steering rate (deg/s)')
xlabel('Longitudinal Distance (m)')
grid on
subplot(3,3,4)
plot(xx(1,1:1:end-1-cut),xx(4,1:1:end-1-cut).*ux(3,2:1:end-cut)*N,'-k','linewidth',1.5);
ylabel('Lidar Limit (m)')
xlabel('Longitudinal Distance (m)')
grid on
subplot(3,3,5)
plot(xx(1,1:1:end-1-cut),xx(4,1:1:end-1-cut),'-k','linewidth',1.5);
ylabel('Longitudinal Speed (m/s)')
xlabel('Longitudinal Distance (m)')
grid on
subplot(3,3,6)
plot(xx(1,1:1:end-1-cut),ux(3,1:1:end-1-cut),'-k','linewidth',1.5);
ylabel('Time Step (s)')
xlabel('Longitudinal Distance (m)')
grid on
subplot(3,3,7)
plot(xx(1,1:1:end-1-cut),xx(8,1:1:end-1-cut)*180/pi,'-k','linewidth',1.5);
ylabel('Steering angle (deg)')
xlabel('Longitudinal Distance (m)')
grid on
subplot(3,3,8)
plot(xx(1,1:1:end-1-cut),xx(7,1:1:end-1-cut),'-k','linewidth',1.5);
ylabel('Acceleration (m/s^2)')
xlabel('Longitudinal Distance (m)')
grid on
subplot(3,3,9)
plot(xx(1,1:1:end-1-cut),ux(1,1:1:end-1-cut),'-k','linewidth',1.5);
ylabel('Jerk (m/s^3)')
xlabel('Longitudinal Distance (m)')
grid on
print -dpng 'Results_States_Controls'





    
    
    
