clear all 
close all
nlist=[30 50 100 300 500 1000 3000 5000 10000 30000 50000 100000];
p=20;
numave=20;

for i=1:size(nlist,2)
  n=nlist(i);
  file=['out/objS_p',num2str(p),'n',num2str(n),'m10r110r280numave',num2str(numave),'.mat']
  load(file)
  file=['out/ctime_p',num2str(p),'n',num2str(n),'m10r110r280numave',num2str(numave),'.mat']
  load(file)
  file=['out/objG_p',num2str(p),'n',num2str(n),'m10r110r280numave',num2str(numave),'.mat']
  load(file)
  file=['out/objK_p',num2str(p),'n',num2str(n),'m10r110r280numave',num2str(numave),'.mat']
  load(file)

  if(isempty(objS.random) == 0)
    rand.obj(i)=mean(objS.random,2);
    rand.obj2(i)=mean(objG.random,2);
    rand.obj3(i)=mean(objK.random,2);
    rand.ctime(i)=mean(ctime.random,2);
    nrand=i;
  end
  if(isempty(objS.AG) == 0)
    AG.obj(i)=mean(objS.AG,2);
    AG.obj2(i)=mean(objG.AG,2);
    AG.obj3(i)=mean(objK.AG,2);
    AG.ctime(i)=mean(ctime.AG,2);
    nAG=i;
  end
  if(isempty(objS.GSDP) == 0)
    GSDP.obj(i)=mean(objS.GSDP,2);
    GSDP.obj2(i)=mean(objG.GSDP,2);
    GSDP.obj3(i)=mean(objK.GSDP,2);
    GSDP.ctime(i)=mean(ctime.GSDP,2);
    nGSDP=i;
  end
  if(isempty(objS.WSDPKF) == 0)
    WSDPKF.obj(i)=mean(objS.WSDPKF,2);
    WSDPKF.obj2(i)=mean(objG.WSDPKF,2);
    WSDPKF.obj3(i)=mean(objK.WSDPKF,2);
    WSDPKF.ctime(i)=mean(ctime.WSDPKF,2);
    nWSDPKF=i;
  end
  if(isempty(objS.CRconvexKF) == 0)
    CRconvexKF.obj(i)=mean(objS.CRconvexKF,2);
    CRconvexKF.obj2(i)=mean(objG.CRconvexKF,2);
    CRconvexKF.obj3(i)=mean(objK.CRconvexKF,2);
    CRconvexKF.ctime(i)=mean(ctime.CRconvexKF,2);
    nCRconvexKF=i;
  end
  if(isempty(objS.convexKF) == 0)
    convexKF.obj(i)=mean(objS.convexKF,2);
    convexKF.obj2(i)=mean(objG.convexKF,2);
    convexKF.obj3(i)=mean(objK.convexKF,2);
    convexKF.ctime(i)=mean(ctime.convexKF,2);
    nconvexKF=i;
  end
  if(isempty(objS.GKF) == 0)
    GKF.obj(i)=mean(objS.GKF,2);
    GKF.obj2(i)=mean(objG.GKF,2);
    GKF.obj3(i)=mean(objK.GKF,2);
    GKF.ctime(i)=mean(ctime.GKF,2);
    nGKF=i;
  end
end

ref(:)=GKF.obj(:);
figure(1)
semilogx(nlist(1:nrand),rand.obj(1:nrand)./ref(1:nrand),'+-','LineWidth',2,'MarkerSize',10,'DisplayName','Rand');
hold on
semilogx(nlist(1:nAG),AG.obj(1:nAG)./ref(1:nAG),'o-','LineWidth',2,'MarkerSize',10,'DisplayName','Greedy (nondynam.)');
semilogx(nlist(1:nGSDP),GSDP.obj(1:nGSDP)./ref(1:nGSDP),'<-','LineWidth',2,'MarkerSize',10,'DisplayName','SDP (Gform)');
semilogx(nlist(1:nWSDPKF),WSDPKF.obj(1:nWSDPKF)./ref(1:nWSDPKF),'>-','LineWidth',2,'MarkerSize',10,'DisplayName','WSDP (Wform)');
semilogx(nlist(1:nCRconvexKF),CRconvexKF.obj(1:nCRconvexKF)./ref(1:nCRconvexKF),'^-','LineWidth',2,'MarkerSize',10,'DisplayName','CR Approx. Convex');
semilogx(nlist(1:nconvexKF),convexKF.obj(1:nconvexKF)./ref(1:nconvexKF),'v-','LineWidth',2,'MarkerSize',10,'DisplayName','Approx Convex');
semilogx(nlist(1:nGKF),GKF.obj(1:nGKF)./ref(1:nGKF),'square-','LineWidth',2,'MarkerSize',10,'DisplayName','Greedy');
set( gca, 'FontName','Times New Roman','FontSize',16 );
xlabel('{\it n}')
ylabel('{\it f}_{SS}/{\it f}_{Greedy}')
grid on
legend
saveas(gcf,'n-fss.fig')

ref(:)=GKF.obj2(:);
figure(2)
semilogx(nlist(1:nrand),rand.obj(1:nrand)./ref(1:nrand),'+-','LineWidth',2,'MarkerSize',10,'DisplayName','Rand');
hold on
semilogx(nlist(1:nAG),AG.obj2(1:nAG)./ref(1:nAG),'o-','LineWidth',2,'MarkerSize',10,'DisplayName','Greedy (nondynam.)');
semilogx(nlist(1:nGSDP),GSDP.obj2(1:nGSDP)./ref(1:nGSDP),'<-','LineWidth',2,'MarkerSize',10,'DisplayName','SDP (Gform)');
semilogx(nlist(1:nWSDPKF),WSDPKF.obj2(1:nWSDPKF)./ref(1:nWSDPKF),'>-','LineWidth',2,'MarkerSize',10,'DisplayName','SDP (Wform)');
semilogx(nlist(1:nCRconvexKF),CRconvexKF.obj2(1:nCRconvexKF)./ref(1:nCRconvexKF),'^-','LineWidth',2,'MarkerSize',10,'DisplayName','CR Approx. Convex');
semilogx(nlist(1:nconvexKF),convexKF.obj2(1:nconvexKF)./ref(1:nconvexKF),'v-','LineWidth',2,'MarkerSize',10,'DisplayName','Approx. Convex');
semilogx(nlist(1:nGKF),GKF.obj2(1:nGKF)./ref(1:nGKF),'square-','LineWidth',2,'MarkerSize',10,'DisplayName','Greedy');
set( gca, 'FontName','Times New Roman','FontSize',16 );
xlabel('{\it n}')
ylabel('{\it f}_{OG}/{\it f}_{Greedy}')
grid on
legend
saveas(gcf,'n-fog.fig')

ref(:)=GKF.obj3(:);
figure(3)
semilogx(nlist(1:nrand),rand.obj(1:nrand)./ref(1:nrand),'+-','LineWidth',2,'MarkerSize',10,'DisplayName','Random sel.');
hold on
semilogx(nlist(1:nAG),AG.obj3(1:nAG)./ref(1:nAG),'o-','LineWidth',2,'MarkerSize',10,'DisplayName','Greedy (nondynam.)');
semilogx(nlist(1:nGSDP),GSDP.obj3(1:nGSDP)./ref(1:nGSDP),'<-','LineWidth',2,'MarkerSize',10,'DisplayName','SDP (Gform)');
semilogx(nlist(1:nWSDPKF),WSDPKF.obj3(1:nWSDPKF)./ref(1:nWSDPKF),'>-','LineWidth',2,'MarkerSize',10,'DisplayName','SDP (Wform)');
semilogx(nlist(1:nCRconvexKF),CRconvexKF.obj3(1:nCRconvexKF)./ref(1:nCRconvexKF),'^-','LineWidth',2,'MarkerSize',10,'DisplayName','CR Approx. Convex');
semilogx(nlist(1:nconvexKF),convexKF.obj3(1:nconvexKF)./ref(1:nconvexKF),'v-','LineWidth',2,'MarkerSize',10,'DisplayName','Approx. Convex');
semilogx(nlist(1:nGKF),GKF.obj3(1:nGKF)./ref(1:nGKF),'square-','LineWidth',2,'MarkerSize',10,'DisplayName','Greedy');
set( gca, 'FontName','Times New Roman','FontSize',16 );
xlabel('{\it n}')
ylabel('{\it f}_{KF}/{\it f}_{KF}_{-Greedy}')
ylim([0.99 1.01])
grid on
legend
saveas(gcf,'n-fkf.fig')


figure(4)
loglog(nlist(1:nrand),rand.ctime(:),'+-','LineWidth',2,'MarkerSize',10,'DisplayName','Random sel.');
hold on
loglog(nlist(1:nAG),AG.ctime(:),'o-','LineWidth',2,'MarkerSize',10,'DisplayName','Greedy (nondynam.)');
loglog(nlist(1:nGSDP),GSDP.ctime(:),'<-','LineWidth',2,'MarkerSize',10,'DisplayName','SDP (Gform)');
loglog(nlist(1:nWSDPKF),WSDPKF.ctime(:),'>-','LineWidth',2,'MarkerSize',10,'DisplayName','SDP (Wform)');
loglog(nlist(1:nCRconvexKF),CRconvexKF.ctime(:),'^-','LineWidth',2,'MarkerSize',10,'DisplayName','CR Approx. Convex');
loglog(nlist(1:nconvexKF),convexKF.ctime(:),'v-','LineWidth',2,'MarkerSize',10,'DisplayName','Approx. Convex');
loglog(nlist(1:nGKF),GKF.ctime(:),'square-','LineWidth',2,'MarkerSize',10,'DisplayName','Greedy');
set( gca, 'FontName','Times New Roman','FontSize',16 );
xlabel('{\it n}')
ylabel('Time [s]')
ylim([0.9999e-5 1e5])
grid on
legend
saveas(gcf,'n-time.fig')