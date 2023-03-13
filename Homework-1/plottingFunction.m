function plottingFunction(constraint,font,thickness,t,xk,uk,yk,ref1,ref2)
    s = get(0, 'ScreenSize');
    figure('Position', [10 50 750 650]);
    hold on;
    grid on;
    stairs(0:t,ref1(1:t+1),LineWidth=thickness);
    stairs(0:t,ref2(1:t+1),LineWidth=thickness);
    stairs(0:t,yk',LineWidth=thickness);
    ylim([-2 12]);
    xlabel('Sample time [k]',FontSize=font+4,Interpreter="latex");
    ylabel({'Reference $r(k)$';'and output $y(k)$ [ft/sec]'},FontSize=font+4,Interpreter="latex");
    if(strcmp(constraint,'unconstrained'))
        title({'Reference inputs and outputs for unconstrained MPC'},FontSize=font+4,Interpreter="latex");
    elseif(strcmp(constraint,'Altconstrained'))
        title({'Reference inputs and outputs for the altered constrained MPC'},FontSize=font+4,Interpreter="latex");    
    else
        title({'Reference inputs and outputs for constrained MPC'},FontSize=font+4,Interpreter="latex");
    end
    legend({'Reference airspeed $(v^*)$','Reference climb rate $(h^*)$',...
        'Output airspeed $(v)$','Output climb rate $(h)$'},FontSize=font,Interpreter="latex");

    figure('Position', [10 50 750 650]);
    hold on;
    grid on;
    stairs(0:t,uk(:,1:t+1)',LineWidth=thickness);
    xlabel('Sample time [k]',FontSize=font+4,Interpreter="latex");
    ylabel('Control input $u(k)$ [-]',FontSize=font+4,Interpreter="latex");
    if(strcmp(constraint,'unconstrained'))
        title('Control inputs for unconstrained MPC',FontSize=font+4,Interpreter="latex");
    elseif(strcmp(constraint,'Altconstrained'))
        title('Control inputs for the altered MPC',FontSize=font+4,Interpreter="latex");    
    else
        title('Control inputs for constrained MPC',FontSize=font+4,Interpreter="latex");
    end
    legend({'Elevator $(e)$','Throttle $(\tau)$'},FontSize=font,Interpreter="latex");

    figure('Position', [10 50 750 650]);
    subplot(2,2,1);
    stairs(0:t,xk(1,:),LineWidth=thickness);
    xlabel('Sample time [k]',FontSize=font,Interpreter="latex");
    ylabel('Airspeed $v(k)$ [ft/sec]',FontSize=font,Interpreter="latex");
    grid on;
    subplot(2,2,2);
    stairs(0:t,xk(2,:),LineWidth=thickness);
    xlabel('Sample time [k]',FontSize=font,Interpreter="latex");
    ylabel({'y-axis';'velocity $w(k)$ [ft/sec]'},FontSize=font,Interpreter="latex");
    grid on;
    subplot(2,2,3);
    stairs(0:t,xk(3,:),LineWidth=thickness);
    xlabel('Sample time [k]',FontSize=font,Interpreter="latex");
    ylabel({'Angular velocity';'component $q(k)$ [rad/sec]'},FontSize=font,Interpreter="latex");
    grid on;
    subplot(2,2,4);
    stairs(0:t,xk(4,:),LineWidth=thickness);
    xlabel('Sample time [k]',FontSize=font,Interpreter="latex");
    ylabel({'x-axis angle';'w.r.t. horizontal $\theta(k)$ [rad]'},FontSize=font,Interpreter="latex");
    grid on;
    if(strcmp(constraint,'unconstrained'))
        sgtitle('State vectors for unconstrained MPC',FontSize=font+4,Interpreter="latex");
    elseif(strcmp(constraint,'Altconstrained'))
        sgtitle('State vectors for the altered MPC',FontSize=font+4,Interpreter="latex");    
    else
        sgtitle('State vectors for constrained MPC',FontSize=font+4,Interpreter="latex");
    end
end