clear all, close all
set(0,'DefaultTextFontName','Times',...
'DefaultTextFontSize',18,...
'DefaultAxesFontName','Times',...
'DefaultAxesFontSize',18,...
'DefaultLineLineWidth',1.5,...
'DefaultLineMarkerSize',7.75)
N = [602 2452 5552 8012 9902];

%% banded A and banded B time
BandedHM_time = [1.16714 11.092535 35.417107 44.885422 78.815647];
BandedMe_time = [0.974696 11.821737 34.6623 47.446783 79.522842];
figure()
sz = 65;
c1 = [0, 0.4470, 0.7410];
c2 = [0.8500, 0.3250, 0.0980];
scatter(N,BandedHM_time,sz, c1,'filled');
hold on
scatter(N,BandedMe_time,sz, c2, 'filled');
legend('hm-toolbox time','Our multiplication time')
P1 = polyfit(N,BandedHM_time,2);
P2 = polyfit(N,BandedMe_time,2);
x = linspace(602,9902);
y1 = polyval(P1,x);
plot(x,y1, '--', 'color', c1);
y2 = polyval(P2,x);
plot(x,y2, '--', 'color', c2);
legend('hm-toolbox time','Our multiplication time');
grid on
hold off
xlabel('Matrix size (number of rows)') 
ylabel('Time (s)') 
xlim([602 9902])

%% bytes
BandedHM_time = [2603845 21031877 67603805 110420049 134864645];
BandedMe_time = [2603957 21175277 67111693 108465193 132472285];
figure()
sz = 65;
c1 = [0, 0.4470, 0.7410];
c2 = [0.8500, 0.3250, 0.0980];
scatter(N,BandedHM_time,sz, c1,'filled');
hold on
scatter(N,BandedMe_time,sz, c2, 'filled');
P1 = polyfit(N,BandedHM_time,2);
P2 = polyfit(N,BandedMe_time,2);
x = linspace(602,9902);
y1 = polyval(P1,x);
plot(x,y1, '--', 'color', c1);
y2 = polyval(P2,x);
plot(x,y2, '--', 'color', c2);
legend('hm-toolbox','Our multiplication')

grid on
hold off
xlabel('Matrix size (number of rows)') 
ylabel('bytes') 
xlim([602 9902])
