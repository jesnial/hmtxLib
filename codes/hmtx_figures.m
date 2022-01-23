clear all, close all
set(0,'DefaultTextFontName','Times',...
'DefaultTextFontSize',18,...
'DefaultAxesFontName','Times',...
'DefaultAxesFontSize',18,...
'DefaultLineLineWidth',1.5,...
'DefaultLineMarkerSize',7.75)
N = [1192	1982	602	2972	4162	10922];

%% banded A and banded B time
BandedHM_time = [2.32E-05	2.52E-05	1.05E-05	4.16E-05	4.58E-05	6.87E-05];
BandedMe_time = [7.42E-06	8.66E-06	4.91E-06	4.16E-05	1.62E-05	2.65E-05];
addsub = [2.29E-05	2.53E-05	1.05E-05	4.16E-05	4.58E-05	7.00E-05];
figure()
sz = 65;
c1 = [0, 0.4470, 0.7410];
c2 = [0.8500, 0.3250, 0.0980];
scatter(N,BandedHM_time,sz, c1,'filled');
hold on
scatter(N,BandedMe_time,sz, c2, 'filled');
scatter(N, addsub, sz, 'filled');
legend('Reconstruction/Transposition/Scalar product','Left/Right matrix-vector product', 'Addition/Subtraction')
% P1 = polyfit(N,BandedHM_time,2);
% P2 = polyfit(N,BandedMe_time,2);
% x1 = linspace(512,4096);
% y1 = polyval(P1,x1);
% plot(x1,y1, '--', 'color', c1);
% x2 = linspace(512,4096);
% y2 = polyval(P2,x2);
% plot(x2,y2, '--', 'color', c2);
% legend('hm-toolbox time','Our multiplication time');
grid on
hold off
xlabel('Matrix size (number of rows)') 
ylabel('Relative error') 
xlim([600 11000])

%% banded A and banded B error
N =[1192	1982	602	2972	4162	7142	10922];
BandedHM_err = [1.37E-05	1.56E-05	5.47E-06	2.05E-05	3.24E-05	3.88E-05	4.62E-05];
BandedMe_err = [25.89E-05	5.57E-05	4.08E-05	7.31E-05	8.77E-05	1.22E-04	1.04E-04];
figure()
sz = 65;
c1 = [0, 0.4470, 0.7410];
c2 = [0.8500, 0.3250, 0.0980];
scatter(N,BandedHM_err,sz, c1,'filled');
hold on
scatter(N,BandedMe_err,sz, c2, 'filled');
l=legend('LU factorization','$\mathcal{H}$-matrix product');
set(l, 'interpreter', 'latex')
grid on
hold off
xlabel('Matrix size (number of rows)') 
ylabel('Relative error') 
xlim([600 11000])

%% Actual number of entries/bytes over number of entries/bytes for a MATLAB array of same size
N = [1192	1982	602	2972	4162	7142	10922	13112	15502];
BandedHM_time = [0.712177	0.493956	1.035204	0.340436	0.251806	0.163939	0.118051	9.72E-02	0.083018];
BandedMe_time = [0.396	0.303806	0.589331	0.215484	0.165996	0.112977	0.076787	6.74E-02	0.058327];
figure()
sz = 65;
c1 = [0, 0.4470, 0.7410];
c2 = [0.8500, 0.3250, 0.0980];
scatter(N,BandedHM_time,sz, c1,'filled');
hold on
scatter(N,BandedMe_time,sz, c2, 'filled');
% P1 = polyfit(N,BandedHM_time,2);
% P2 = polyfit(N,BandedMe_time,2);
% x = linspace(602,15502);
% y1 = polyval(P1,x);
% plot(x,y1, '--', 'color', c1);
% x = linspace(602,15502);
% y2 = polyval(P2,x);
% plot(x,y2, '--', 'color', c2);
legend('Memory compression','Entry compression');
grid on
hold off
xlabel('Matrix size (number of rows)') 
ylabel('Compression rate') 
xlim([600 15550])

%% primal functions time
N = [1192	1982	602	2972	4162	7142	10922	13112];
BandedHM_err = [0.034378	0.062712	0.011731	0.065981	0.087611	1.62E-01	0.292209	0.413265];
BandedMe_err = [0.008541	0.030992	0.00366	0.025772	0.028074	0.158557	0.103943	0.159265];
add = [0.096255	0.272046	0.035312	0.277278	0.385542	1.029323	1.534244	1.806641];
LRmv = [0.032083	0.083905	0.011113	0.080969	0.23077	0.309773	0.64879	0.793916];
Trisolv =[0.022338	0.03549	0.012143	0.040784	0.04931	0.37614	0.302624	0.29366];
figure()
sz = 65;
c1 = [0, 0.4470, 0.7410];
c2 = [0.8500, 0.3250, 0.0980];
c3 =[0.9290, 0.6940, 0.1250];
c4 =[0.4940, 0.1840, 0.5560];
c5 = [0.4660, 0.6740, 0.1880];
scatter(N,BandedHM_err,sz, c1,'filled');
hold on
scatter(N,BandedMe_err,sz, c2, 'filled');
scatter(N,add,sz, c3, 'filled');
scatter(N,LRmv,sz, c4, 'filled');
scatter(N,Trisolv,sz, c5, 'filled');
P1 = polyfit(N,BandedHM_err,2);
P2 = polyfit(N,BandedMe_err,2);
P3 = polyfit(N,add,2);
P4 = polyfit(N,LRmv,2);
P5 = polyfit(N,Trisolv,2);
x = linspace(602,13112);
y1 = polyval(P1,x);
plot(x,y1, '--', 'color', c1);
y2 = polyval(P2,x);
plot(x,y2, '--', 'color', c2);
y3 = polyval(P3,x);
plot(x,y3, '--', 'color', c3);
y4 = polyval(P4,x);
plot(x,y4, '--', 'color', c4);
y5 = polyval(P5,x);
plot(x,y5, '--', 'color', c5);
legend('Transposition','Scalar mult.', 'Addition/Subtraction', 'L/R matrix-vector product', 'Matrix-vector solver')
grid on
hold off
xlabel('Matrix size (number of rows)') 
ylabel('Time (s)') 
xlim([602 13112])

%% secondary functions time
N = [1192	1982	602	2972	4162	7142	10922	13112];
BandedHM_err = [6.052756	13.588431	1.58345	21.251853	27.13442	67.869073	1.39E+02	145.806805];
BandedMe_err = [6.299897	14.137629	1.573823	23.427183	30.161221	80.320588	1.52E+02	189.138745];
add = [2.764215	5.423584	0.86886	9.072814	14.450342	3.51E+01	9.28E+01	97.389369];
LRmv = [4.098435	8.115565	1.16487	15.466938	29.243371	5.53E+01	1.20E+02	123.524776];

figure()
sz = 65;
c1 = [0, 0.4470, 0.7410];
c2 = [0.8500, 0.3250, 0.0980];
c3 =[0.9290, 0.6940, 0.1250];
c4 =[0.4940, 0.1840, 0.5560];
c5 = [0.4660, 0.6740, 0.1880];
scatter(N,BandedHM_err,sz, c1,'filled');
hold on
scatter(N,BandedMe_err,sz, c2, 'filled');
scatter(N,add,sz, c3, 'filled');
scatter(N,LRmv,sz, c4, 'filled');

P1 = polyfit(N,BandedHM_err,2);
P2 = polyfit(N,BandedMe_err,2);
P3 = polyfit(N,add,2);
P4 = polyfit(N,LRmv,2);

x = linspace(602,13112);
y1 = polyval(P1,x);
plot(x,y1, '--', 'color', c1);
y2 = polyval(P2,x);
plot(x,y2, '--', 'color', c2);
y3 = polyval(P3,x);
plot(x,y3, '--', 'color', c3);
y4 = polyval(P4,x);
plot(x,y4, '--', 'color', c4);

l2 =legend('$\mathcal{H}$-matrix product','$\mathcal{H}$-matrix inverset', 'LU factorization', 'Triangular system solver')
set(l2, 'interpreter', 'latex')
grid on
hold off
xlabel('Matrix size (number of rows)') 
ylabel('Time (s)') 
xlim([602 13112])