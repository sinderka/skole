





% data(j,1) = k % number of time steps
% data(j,2) = m % number of steps in space
% data(j,3) = n % restartvariable
% data(j,4) = maximum error
% data(j,5) = time
% data(j,6) = prob
% data(j,7) = numf
% data(j,8) = para



s = 0;
if s == 1
    data = zeros(1568,8,10);
    avg = 10;
    
    
    i = 1;
    for para = [1,2,4]
        parpool(para);
        for numf = [1,4,16,32]
            if para <= numf
                for k = [10,50,100,200]
                    for m = [10,50,100,200]
                        for n = [1,4,16,32]
                            if m>n
                                for l = 1:avg
                                    prob = 2;
                                    [ data(i,:,l),~ ] = kpmtester(k,m,n,prob,numf);
                                    data(i,8,l) = para;
                                    i = i+1;
                                    prob = 4;
                                    [ data(i,:,l),~ ] = kpmtester(k,m,n,prob,numf);
                                    data(i,8,l) = para;
                                    i = i-1;
                                end
                                i = i+2;
                            end
                        end
                    end
                end
            end
        end
        delete(gcp);
    end
    
    fprintf('1\n')
    
    for numf = [1,4,16,32]
        for k = [10,50,100,200]
            for m = [10,50,100,200]
                for n = [1,4,16,32]
                    if m>n
                        for l = 1:avg
                            prob = 1;
                            [ data(i,:,l),~ ] = kpmtester(k,m,n,prob,numf);
                            data(i,8,l) = 0;
                            i = i+1;
                            prob = 3;
                            [ data(i,:,l),~ ] = kpmtester(k,m,n,prob,numf);
                            data(i,8,l) = 0;
                            i = i-1;
                        end
                        i = i+2;
                    end
                end
            end
        end
    end
    
    
    fprintf('2\n')
    
    data = sum(data,3 )/avg;
    s = 'data';
    save(s,'data');
end
% Hva er interesangt av resulteter?
% konvergens: for alle metodene. surf plot med k en vei og m en annen vei

load('data','data');



%%%% Må sammenligne kjøretid med restartvariabel og uten KPM

%%%% Må sammenligne tider paralell/ikke paralell med/uten KPM med tider
a1 = intersect(intersect(intersect(intersect(find(data(:,8) == 0), find(data(:,2) == 200)),find(data(:,1) == 200)), find(data(:,6) == 1)), find(data(:,3) == 1) );   % ikkeparalell KPM bortsett fra numf
a2 = intersect(intersect(intersect(intersect(find(data(:,8) == 4), find(data(:,2) == 200)),find(data(:,1) == 200)), find(data(:,6) == 2)), find(data(:,3) == 1) );   % parallell KPM bortsett fra numf
a3 = intersect(intersect(intersect(intersect(find(data(:,8) == 0), find(data(:,2) == 200)),find(data(:,1) == 200)), find(data(:,6) == 3)), find(data(:,3) == 1) );
a4 = intersect(intersect(intersect(intersect(find(data(:,8) == 4), find(data(:,2) == 200)),find(data(:,1) == 200)), find(data(:,6) == 4)), find(data(:,3) == 1) );

figure(1)
%plot(data(a1,7),data(a1,5),'k:h');
%hold on
plot(data(a2,7),data(a2,5),'k:p');
hold on
%loglog(data(a3,7),data(a3,5),'k:s');
loglog(data(a4,7),data(a4,5),'k:d');
ylabel('computation time'); xlabel('numf');
%legend('KPM','KPM, nP = 4','Integration','Integration, nP = 4');
%legend('KPM','Integration');
legend('KPM, nP = 4','Integration, nP = 4');

%%%% Konvergens:
b1 = intersect(intersect(intersect(intersect(find(data(:,8) == 0), find(data(:,2) == 200)),find(data(:,7) == 32)), find(data(:,6) == 1)), find(data(:,3) == 1) );
b2 = intersect(intersect(intersect(intersect(find(data(:,8) == 0), find(data(:,1) == 200)),find(data(:,7) == 32)), find(data(:,6) == 1)), find(data(:,3) == 1) );

figure(2)
loglog(data(b1,1),data(b1,4),'k:h');
hold on
loglog(data(b2,2),data(b2,4),'k:s');
loglog(data(b2,2),1./(data(b2,2)).^2,'k-');
xlabel('m / k'); ylabel('error');
legend('KPM, k changing','KPM, m changing','helpline');


%%%% Sammenligen de forskjellige tidene med forskjllige restartvariabler
c1 = intersect(intersect(intersect(intersect(find(data(:,8) == 0), find(data(:,2) == 200)),find(data(:,7) == 32)), find(data(:,6) == 1)), find(data(:,1) == 200) );
c2 = intersect(intersect(intersect(intersect(find(data(:,8) == 0), find(data(:,2) == 200)),find(data(:,7) == 4)), find(data(:,6) == 1)), find(data(:,1) == 200) );
c3 = intersect(intersect(intersect(intersect(find(data(:,8) == 4), find(data(:,2) == 200)),find(data(:,7) == 32)), find(data(:,6) == 2)), find(data(:,1) == 200) );
c4 = intersect(intersect(intersect(intersect(find(data(:,8) == 4), find(data(:,2) == 200)),find(data(:,7) == 4)), find(data(:,6) == 2)), find(data(:,1) == 200) );

figure(3)
plot(data(c1,3),data(c1,5),'k:h');
hold on
loglog(data(c2,3),data(c2,5),'k:p');
loglog(data(c3,3),data(c3,5),'k:s');
loglog(data(c4,3),data(c4,5),'k:d');
xlabel('n'); ylabel('computation time');
legend('KPM, numf = 32','KPM, numf = 4','KPM, numf = 32, nP = 4', 'KPM, numf = 4, nP = 4');


%%%% Sammenligne KPM med integrator for forskjellige k
d1 = intersect(intersect(intersect(intersect(find(data(:,8) == 0), find(data(:,2) == 200)),find(data(:,7) == 1)), find(data(:,6) == 1)), find(data(:,3) == 1) );
d2 = intersect(intersect(intersect(intersect(find(data(:,8) == 0), find(data(:,2) == 200)),find(data(:,7) == 1)), find(data(:,6) == 1)), find(data(:,3) == 4) );
d3 = intersect(intersect(intersect(intersect(find(data(:,8) == 0), find(data(:,2) == 200)),find(data(:,7) == 1)), find(data(:,6) == 3)), find(data(:,3) == 1) );
figure(5)
loglog(data(d1,1),data(d1,5),'k:h');
hold on
loglog(data(d2,1),data(d2,5),'k:p');
loglog(data(d3,1),data(d3,5),'k:s');
xlabel('k'); ylabel('computation time');
legend('KPM','KPM, n = 4','Integration');
%legend('KPM','Integration');

%%%% Sammenligne KPM med integrator for forskjellige k
e1 = intersect(intersect(intersect(intersect(find(data(:,8) == 0), find(data(:,1) == 200)),find(data(:,7) == 1)), find(data(:,6) == 1)), find(data(:,3) == 1) );
e2 = intersect(intersect(intersect(intersect(find(data(:,8) == 0), find(data(:,1) == 200)),find(data(:,7) == 1)), find(data(:,6) == 1)), find(data(:,3) == 4) );
e3 = intersect(intersect(intersect(intersect(find(data(:,8) == 0), find(data(:,1) == 200)),find(data(:,7) == 1)), find(data(:,6) == 3)), find(data(:,3) == 1) );
figure(6)
loglog(data(e1,2),data(e1,5),'k:h');
hold on
loglog(data(e2,2),data(e2,5),'k:p');
loglog(data(e3,2),data(e3,5),'k:s');
xlabel('m'); ylabel('computation time');
legend('KPM','KPM, n = 4','Integration');
%legend('KPM','Integration');

%%%% Sammenligne KPM med forskjellig antal processorer
f1 = intersect(intersect(intersect(intersect(find(data(:,2) == 200), find(data(:,1) == 200)),find(data(:,7) == 32)), find(data(:,6) == 2)), find(data(:,3) == 1) );
f2 = intersect(intersect(intersect(intersect(find(data(:,2) == 200), find(data(:,1) == 200)),find(data(:,7) == 16)), find(data(:,6) == 2)), find(data(:,3) == 1) );
f3 = intersect(intersect(intersect(intersect(find(data(:,2) == 200), find(data(:,1) == 200)),find(data(:,7) == 4)), find(data(:,6) == 2)), find(data(:,3) == 1) );
figure(7)
loglog(data(f1,8),data(f1,5),'k:h');
hold on
loglog(data(f2,8),data(f2,5),'k:p');
loglog([1,2,4],[4,2,1],'k-');
xlabel('nP'); ylabel('computation time');
legend('KPM, numf = 32','KPM, numf = 16','helpline');

g1 = intersect(intersect(intersect(intersect(intersect(find(data(:,2) == 200), find(data(:,1) == 200)),find(data(:,7) == 32)), find(data(:,6) == 2)), find(data(:,3) == 1)),find(data(:,8) == 1));
g2 = intersect(intersect(intersect(intersect(intersect(find(data(:,2) == 200), find(data(:,1) == 200)),find(data(:,7) == 32)), find(data(:,6) == 2)), find(data(:,3) == 1)),find(data(:,8) == 2));
g3 = intersect(intersect(intersect(intersect(intersect(find(data(:,2) == 200), find(data(:,1) == 200)),find(data(:,7) == 32)), find(data(:,6) == 2)), find(data(:,3) == 1)),find(data(:,8) == 4));
Sp1 = data(g1,5)./data(g2,5);
Sp2 = data(g1,5)/data(g3,5);

eta1 = Sp1/2;
eta2 = Sp2/4;





