%%%pploter
%%% ( k,m,n,prob,para,func)
if 0 % convergence
    %close 1
    j = 1;
    
    p = [5,10,20,50, 100 ,200, 500];
    inndata = zeros(length(p),6); inndatak = zeros(length(p),6);
    utdata  = zeros(length(p),3); utdatak  = zeros(length(p),3);
    for i = p
        [ inndatak(j,:),utdatak(j,:) ] = kpmtester( i,i,1,1,1,1);
        [ inndata(j,:),utdata(j,:) ] = kpmtester( i,i,1,3,1,1);
        j = j + 1;
    end
    figure(1)
    loglog(p,utdatak(:,3),'k:p')
    hold on
    plot(p,utdata(:,3),'k:d')
    plot(p,1./p.^2,'k-')
    legend('KPM','Direct integration','helpline');
    xlabel('m / k'); ylabel('error');
end
if 0 % convergence
    %close 1
    j = 1;
    
    p = [5,10,20,50, 100 ,200, 500];
    inndata = zeros(length(p),6); inndatak = zeros(length(p),6);
    utdata  = zeros(length(p),3); utdatak  = zeros(length(p),3);
    for i = p
        [ inndatak(j,:),utdatak(j,:) ] = kpmtester( i,i,1,1,1,2);
        [ inndata(j,:),utdata(j,:) ] = kpmtester( i,i,1,3,1,2);
        j = j + 1;
    end
    figure(500)
    loglog(p,utdatak(:,3),'k:p')
    hold on
    plot(p,utdata(:,3),'k:d')
    plot(p,1./p.^2,'k-')
    legend('KPM','Direct integration','helpline');
    xlabel('m / k'); ylabel('error');
end

if 0 % ant vs m % bare for funksjon 1
    p = [5,10,20,50, 100 ,200, 500,1000];
    k = 200;
    j = 1;
    inndata = zeros(length(p),6); inndatak = zeros(length(p),6);
    utdata  = zeros(length(p),3); utdatak  = zeros(length(p),3);
    for i = p
        [ inndatak(j,:),utdatak(j,:) ] = kpmtester( k,i,1,1,1,1);
        [ inndataj(j,:),utdataj(j,:) ] = kpmtester( k,i,2,1,1,1);
        [ inndata(j,:),utdata(j,:) ] = kpmtester( k,i,5,1,1,1);
        j = j + 1;
    end
    figure(2)
    plot(p,utdatak(:,1),'k:p')
    hold on
    plot(p,utdataj(:,1),'k:s')
    plot(p,utdata(:,1),'k:d')
    legend('n = 1','n = 2','n = 5');
    xlabel('m'); ylabel('nI');
end

if 0 % ant vs m % bare for funksjon 2
    p = [5,10,20,50, 100 ,200, 500,1000];
    k = 200;
    j = 1;
    inndata = zeros(length(p),6); inndatak = zeros(length(p),6);
    utdata  = zeros(length(p),3); utdatak  = zeros(length(p),3);
    for i = p
        [ inndatak(j,:),utdatak(j,:) ] = kpmtester( k,i,1,1,1,2);
        [ inndataj(j,:),utdataj(j,:) ] = kpmtester( k,i,2,1,1,2);
        [ inndata(j,:),utdata(j,:) ] = kpmtester( k,i,5,1,1,2);
        j = j + 1;
    end
    figure(3)
    plot(p,utdatak(:,1),'k:p')
    hold on
    plot(p,utdataj(:,1),'k:s')
    plot(p,utdata(:,1),'k:d')
    legend('n = 1','n = 2','n = 5');
    xlabel('m'); ylabel('nI');
end


if 0 % time vs m % bare for funksjon 1
    p = [5,10,20,50, 100 ,200, 500,1000];
    avg = 10;
    k = 200;
    j = 1;
    inndata = zeros(length(p),6,avg); utdata  = zeros(length(p),3,avg);
    inndatak = zeros(length(p),6,avg); utdatak  = zeros(length(p),3,avg);
    
    for i = p
        for a = 1:avg
            [ inndata(j,:,a),utdata(j,:,a) ] = kpmtester( i,k,1,1,1,1);
            [ inndatak(j,:,a),utdatak(j,:,a) ] = kpmtester( i,k,1,3,1,1);
        end
        j = j + 1;
        
    end
    utdata = sum(utdata,3)/avg;
    utdatak = sum(utdatak,3)/avg;
    figure(4)
    loglog(p,utdatak(:,2),'k:p')
    hold on
    plot(p,utdata(:,2),'k:s')
    plot(p,p/100,'k-')
    legend('KPM','Direct integration');
    xlabel('m'); ylabel('computation time');
end

if 0 % time vs m % bare for funksjon 2
    p = [5,10,20,50, 100 ,200, 500,1000];
    avg = 10;
    k = 200;
    j = 1;
    inndata = zeros(length(p),6,avg); utdata  = zeros(length(p),3,avg);
    inndatak = zeros(length(p),6,avg); utdatak  = zeros(length(p),3,avg);
    
    for i = p
        for a = 1:avg
            [ inndata(j,:,a),utdata(j,:,a) ] = kpmtester( i,k,1,1,1,2);
            [ inndatak(j,:,a),utdatak(j,:,a) ] = kpmtester( i,k,1,3,1,2);
        end
        j = j + 1;
        
    end
    utdata = sum(utdata,3)/avg;
    utdatak = sum(utdatak,3)/avg;
    figure(5)
    loglog(p,utdatak(:,2),'k:p')
    hold on
    plot(p,utdata(:,2),'k:s')
    plot(p,p/100,'k-')
    legend('KPM','Direct integration');
    xlabel('m'); ylabel('computation time');
end
if 0 % time vs k % bare for funksjon 1
    p = [5,10,20,50, 100 ,200, 500,1000];
    avg = 10;
    k = 200;
    j = 1;
    inndata = zeros(length(p),6,avg); utdata  = zeros(length(p),3,avg);
    inndatak = zeros(length(p),6,avg); utdatak  = zeros(length(p),3,avg);
    
    for i = p
        for a = 1:avg
            [ inndata(j,:,a),utdata(j,:,a) ] = kpmtester( i,k,1,1,1,1);
            [ inndatak(j,:,a),utdatak(j,:,a) ] = kpmtester( i,k,1,3,1,1);
        end
        j = j + 1;
        
    end
    utdata = sum(utdata,3)/avg;
    utdatak = sum(utdatak,3)/avg;
    figure(6)
    loglog(p,utdatak(:,2),'k:p')
    hold on
    plot(p,utdata(:,2),'k:s')
    plot(p,p/100,'k-')
    legend('KPM','Direct integration');
    xlabel('k'); ylabel('computation time');
end
if 0 % time vs k % bare for funksjon 2
    p = [5,10,20,50, 100 ,200, 500,1000];
    avg = 10;
    k = 200;
    j = 1;
    inndata = zeros(length(p),6,avg); utdata  = zeros(length(p),3,avg);
    inndatak = zeros(length(p),6,avg); utdatak  = zeros(length(p),3,avg);
    
    for i = p
        for a = 1:avg
            [ inndata(j,:,a),utdata(j,:,a) ] = kpmtester( i,k,1,1,1,2);
            [ inndatak(j,:,a),utdatak(j,:,a) ] = kpmtester( i,k,1,3,1,2);
        end
        j = j + 1;
        
    end
    utdata = sum(utdata,3)/avg;
    utdatak = sum(utdatak,3)/avg;
    figure(7)
    loglog(p,utdatak(:,2),'k:p')
    hold on
    plot(p,utdata(:,2),'k:s')
    plot(p,p/100,'k-')
    legend('KPM','Direct integration');
    xlabel('k'); ylabel('computation time');
end

if 0 %restartvariabel funksjon 1
    p = [1,2,3,4,5,10,20,30,40,50, 100];
    avg = 10;
    k = 200;
    j = 1;
    inndata = zeros(length(p),6,avg); utdata  = zeros(length(p),3,avg);
    inndatak = zeros(length(p),6,avg); utdatak  = zeros(length(p),3,avg);
    inndataj = zeros(length(p),6,avg); utdataj  = zeros(length(p),3,avg);
    
    for i = p
        for a = 1:avg
            [ inndata(j,:,a),utdata(j,:,a) ] = kpmtester( k,500,i,1,1,1);
            [ inndataj(j,:,a),utdataj(j,:,a) ] = kpmtester( k,100,i,1,1,1);
            [ inndatak(j,:,a),utdatak(j,:,a) ] = kpmtester( k,10,i,1,1,1);
        end
        j = j + 1;
        
    end
    utdata = sum(utdata,3)/avg;
    utdatak = sum(utdatak,3)/avg;
    utdataj = sum(utdataj,3)/avg;
    figure(9)
    loglog(p,utdatak(:,2),'k:h')
    hold on
    loglog(p,utdataj(:,2),'k:s')
    loglog(p,utdata(:,2),'k:p')
    xlabel('n'); ylabel('computation time');
    legend('m = 10','m = 100','m = 500')
    
end
if 0 %restartvariabel
    p = [1,2,3,4,5,10,20,30,40,50, 100 ,200, 500];
    avg = 10;
    k = 200;
    j = 1;
    inndata = zeros(length(p),6,avg); utdata  = zeros(length(p),3,avg);
    inndatak = zeros(length(p),6,avg); utdatak  = zeros(length(p),3,avg);
    inndataj = zeros(length(p),6,avg); utdataj  = zeros(length(p),3,avg);
    
    for i = p
        for a = 1:avg
            [ inndata(j,:,a),utdata(j,:,a) ] = kpmtester( k,500,i,1,1,2);
            [ inndataj(j,:,a),utdataj(j,:,a) ] = kpmtester( k,100,i,1,1,2);
            [ inndatak(j,:,a),utdatak(j,:,a) ] = kpmtester( k,10,i,1,1,2);
        end
        j = j + 1;
        
    end
    utdata = sum(utdata,3)/avg;
    utdatak = sum(utdatak,3)/avg;
    utdataj = sum(utdataj,3)/avg;
    figure(10)
    loglog(p,utdatak(:,2),'k:h')
    hold on
    loglog(p,utdataj(:,2),'k:s')
    loglog(p,utdata(:,2),'k:p')
    xlabel('n'); ylabel('computation time');
    legend('m = 10','m = 100','m = 500')
end

%%%%%%%%%%%%%%%%%burde ha med noe med konvergenskriteriet, men er ikke helt
%%%%%%%%%%%%%%%%%sikker på ha det burde være...%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%må ha med time vs n med forskjellige m og func%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%