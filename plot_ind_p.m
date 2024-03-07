clear all 
close all
n=200;
plist=[12, 15, 20, 30, 40, 50, 70, 100];
numave=5;

for i=1:size(plist,2)
  p=plist(i);
  file=['out/obj_p',num2str(p),'n',num2str(n),'m10r110r280numave',num2str(numave),'.mat']
  load(file)
  file=['out/ctime_p',num2str(p),'n',num2str(n),'m10r110r280numave',num2str(numave),'.mat']
  load(file)
  file=['out/obj2_p',num2str(p),'n',num2str(n),'m10r110r280numave',num2str(numave),'.mat']
  load(file)
  file=['out/obj3_p',num2str(p),'n',num2str(n),'m10r110r280numave',num2str(numave),'.mat']
  load(file)

  
  rand.obj(i)=mean(obj.random,2);
  AG.obj(i)=mean(obj.AG,2);
  GSDP.obj(i)=mean(obj.GSDP,2);
  WSDPKF.obj(i)=mean(obj.WSDPKF,2);
  CRconvexKF.obj(i)=mean(obj.CRconvexKF,2);
  convexKF.obj(i)=mean(obj.convexKF,2);
  GKF.obj(i)=mean(obj.GKF,2);

  rand.obj2(i)=mean(obj2.random,2);
  AG.obj2(i)=mean(obj2.AG,2);
  GSDP.obj2(i)=mean(obj2.GSDP,2);
  WSDPKF.obj2(i)=mean(obj2.WSDPKF,2);
  CRconvexKF.obj2(i)=mean(obj2.CRconvexKF,2);
  convexKF.obj2(i)=mean(obj2.convexKF,2);
  GKF.obj2(i)=mean(obj2.GKF,2);

  rand.obj3(i)=mean(obj3.random,2);
  AG.obj3(i)=mean(obj3.AG,2);
  GSDP.obj3(i)=mean(obj3.GSDP,2);
  WSDPKF.obj3(i)=mean(obj3.WSDPKF,2);
  CRconvexKF.obj3(i)=mean(obj3.CRconvexKF,2);
  convexKF.obj3(i)=mean(obj3.convexKF,2);
  GKF.obj3(i)=mean(obj3.GKF,2);

  rand.ctime(i)=mean(ctime.random,2);
  AG.ctime(i)=mean(ctime.AG,2);
  GSDP.ctime(i)=mean(ctime.GSDP,2);
  WSDPKF.ctime(i)=mean(ctime.WSDPKF,2);
  CRconvexKF.ctime(i)=mean(ctime.CRconvexKF,2);
  convexKF.ctime(i)=mean(ctime.convexKF,2);
  GKF.ctime(i)=mean(ctime.GKF,2);
end

ref(:)=GKF.obj(:);
figure(1)
plot(plist(:),rand.obj(:)./ref(:),'LineWidth',2,'DisplayName','Rand');
hold on
plot(plist(:),AG.obj(:)./ref(:),        'LineWidth',2, 'DisplayName','AG');
plot(plist(:),GSDP.obj(:)./ref(:),      'LineWidth',2, 'DisplayName','GSDP');
plot(plist(:),WSDPKF.obj(:)./ref(:),     'LineWidth',2, 'DisplayName','WSDPKF');
plot(plist(:),CRconvexKF.obj(:)./ref(:), 'LineWidth',2,'DisplayName','CRconvexKF');
plot(plist(:),convexKF.obj(:)./ref(:),   'LineWidth',2,'DisplayName','convexKF');
plot(plist(:),GKF.obj(:)./ref(:),        'LineWidth',2,'DisplayName','GKF');
set( gca, 'FontName','Arial','FontSize',16 );
legend

ref(:)=GKF.obj2(:);
figure(2)
plot(plist(:),rand.obj(:)./ref(:),'LineWidth',2,'DisplayName','Rand');
hold on
plot(plist(:),AG.obj2(:)./ref(:),        'LineWidth',2, 'DisplayName','AG');
plot(plist(:),GSDP.obj2(:)./ref(:),      'LineWidth',2, 'DisplayName','GSDP');
plot(plist(:),WSDPKF.obj2(:)./ref(:),     'LineWidth',2, 'DisplayName','WSDPKF');
plot(plist(:),CRconvexKF.obj2(:)./ref(:), 'LineWidth',2,'DisplayName','CRconvexKF');
plot(plist(:),convexKF.obj2(:)./ref(:),   'LineWidth',2,'DisplayName','convexKF');
plot(plist(:),GKF.obj2(:)./ref(:),        'LineWidth',2,'DisplayName','GKF');
set( gca, 'FontName','Arial','FontSize',16 );
legend

ref(:)=GKF.obj3(:);
figure(3)
plot(plist(:),rand.obj(:)./ref(:),'LineWidth',2,'DisplayName','Rand');
hold on
plot(plist(:),AG.obj3(:)./ref(:),        'LineWidth',2, 'DisplayName','AG');
plot(plist(:),GSDP.obj3(:)./ref(:),      'LineWidth',2, 'DisplayName','GSDP');
plot(plist(:),WSDPKF.obj3(:)./ref(:),     'LineWidth',2, 'DisplayName','WSDPKF');
plot(plist(:),CRconvexKF.obj3(:)./ref(:), 'LineWidth',2,'DisplayName','CRconvexKF');
plot(plist(:),convexKF.obj3(:)./ref(:),   'LineWidth',2,'DisplayName','convexKF');
plot(plist(:),GKF.obj3(:)./ref(:),        'LineWidth',2,'DisplayName','GKF');
set( gca, 'FontName','Arial','FontSize',16 );
ylim([0.95 1.05])
legend


figure(4)
semilogy(plist(:),rand.ctime(:),'-k','LineWidth',2,'DisplayName','Rand');
hold on
semilogy(plist(:),AG.ctime(:),        'LineWidth',2,'DisplayName','AG');
semilogy(plist(:),GSDP.ctime(:),      'LineWidth',2,'DisplayName','GSPD');
semilogy(plist(:),WSDPKF.ctime(:),    'LineWidth',2,'DisplayName','WSDPKF');
semilogy(plist(:),CRconvexKF.ctime(:),'LineWidth',2,'DisplayName','CRconvexKF');
semilogy(plist(:),convexKF.ctime(:),  'LineWidth',2,'DisplayName','convexKF');
semilogy(plist(:),GKF.ctime(:),       'LineWidth',2,'DisplayName','GKF');
set( gca, 'FontName','Arial','FontSize',16 );
legend


