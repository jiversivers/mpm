function params = uiGetPhasorParamsv2
    % Default values
    params = struct('f', 80, ...
                    'n', 1, ...
                   'dt', 0.0488, ...
              'tau_ref', 0, ...
                 'calc', false);
    calibration = struct();
    params.w = 2*pi*params.f*params.n*10^6;
    calibration.m_ref = 1/(sqrt(1+params.w*params.tau_ref*10^-9));
    calibration.phi_ref = atan(params.w*params.tau_ref*10^-9);
    
    %% Figure
    % Parent figure
    fig = uifigure;
    fig.Position = [561 529 500 250];
    gl = uigridlayout(fig, [1 2]);
    gl.Padding = [10 10 10 10];
    
    % Input panel
    ip = uipanel(gl);
    ip.Layout.Row = 1;
    ip.Layout.Column = 1;
    ip.BackgroundColor = [0.9 0.9 0.9];
    
    % Calculation panel
    xp = uipanel(gl);
    xp.Layout.Row = 1;
    xp.Layout.Column = 2;
    xp.BackgroundColor = [0.75 0.75 0.75];
    
    %% Headers
    % Prompt
    l = uilabel(ip);
    l.Position = [10 200 215 25];
    l.Text = 'Input Imaging Parameters';
    l.FontWeight = 'bold';
    l.HorizontalAlignment = 'center';
    
    % Calculations
    x = uilabel(xp);
    x.Position = [10 200 215 25];
    x.Text = 'Calculated Imaging Parameters';
    x.FontWeight = 'bold';
    x.HorizontalAlignment = 'center';
    
    %% Okay button
    ok = uibutton(ip);
    ok.Text = 'Load IRF';
    ok.Position = [80 10 125 25];
    ok.ButtonPushedFcn = 'IRF = uigetfile;';
    
    %% Reference lifetime
    rt = uieditfield(ip);
    rt.Position = [105 50 100 25];
    rt.HorizontalAlignment = 'right';
    rt.Value = num2str(params.tau_ref);
    rt.ValueChangedFcn = "params=changeParam(params, str2num(rt.Value), 'tau_ref');";
    
    % Labels
    rl1 = uilabel(ip);
    rl1.Position = [20 55 75 25];
    rl1.HorizontalAlignment = 'right';
    rl1.Text = 'Ref. Lifetime';
    
    rl2 = uilabel(ip);
    rl2.Position = [20 45 75 25];
    rl2.HorizontalAlignment = 'right';
    rl2.Text = '(ns)';
    
    %% Bin width
    bw = uieditfield(ip);
    bw.Position = [105 90 100 25];
    bw.HorizontalAlignment = 'right';
    bw.Value = num2str(params.dt);
    bw.ValueChangedFcn = "params = changeParam(params, str2num(bw.Value), 'dt');";
    
    % Check box
    br = uicheckbox(ip);
    br.Position = [10 95 100 25];
    br.Text = '';
    br.Value = params.calc;
    br.ValueChangedFcn = "params.calc = ~params.calc; switchButton(br, bw);";
    
    % Labels
    bl1 = uilabel(ip);
    bl1.Position = [20 95 75 25];
    bl1.HorizontalAlignment = 'right';
    bl1.Text = 'Bin Width';
    
    bl2 = uilabel(ip);
    bl2.Position = [20 85 75 25];
    bl2.HorizontalAlignment = 'right';
    bl2.Text = '(ns)';
    
    %% Harmonic
    hd = uieditfield(ip);
    hd.Position = [105 130 100 25];
    hd.HorizontalAlignment = 'right';
    hd.Value = num2str(params.n);
    hd.ValueChangedFcn = "params = changeParam(params, str2num(hd.Value), 'n');";
    
    % Label
    hl = uilabel(ip);
    hl.Position = [10 130 85 25];
    hl.HorizontalAlignment = 'right';
    hl.Text = 'Harmonic';
    
    %% Laser Frequency
    fd = uieditfield(ip);
    fd.Position = [105 170 100 25];
    fd.HorizontalAlignment = 'right';
    fd.Value = num2str(params.f);
    fd.ValueChangedFcn = "params = changeParam(params, str2num(fd.Value), 'f');";
    
    % labels
    fl1 = uilabel(ip);
    fl1.Position = [10 175 85 25];
    fl1.Text = 'Frequency';
    fl1.HorizontalAlignment = 'right';
    
    fl2 = uilabel(ip);
    fl2.Position = [10 165 85 25];
    fl2.Text = '(MHz)';
    fl2.HorizontalAlignment = 'right';
    
    %% Calulations %%
    %% Angular Frequency
    wd = uieditfield(xp);
    wd.Position = [105 170 100 25];
    wd.BackgroundColor = [0.7 0.7 0.7];
    wd.HorizontalAlignment = 'right';
    wd.Value = num2str(params.w);
    wd.Editable = "off";
    
    % labels
    xl1 = uilabel(xp);
    xl1.Position = [10 175 85 25];
    xl1.Text = 'Angular';
    xl1.HorizontalAlignment = 'right';
    
    xl2 = uilabel(xp);
    xl2.Position = [10 165 85 25];
    xl2.Text = 'Frequency';
    xl2.HorizontalAlignment = 'right';
    
    %% Ref Modulation
    md = uieditfield(xp);
    md.Position = [105 130 100 25];
    md.BackgroundColor = [0.7 0.7 0.7];
    md.HorizontalAlignment = 'right';
    md.Value = num2str(calibration.m_ref);
    md.Editable = "off";
    
    % Label
    ml = uilabel(xp);
    ml.Position = [10 130 85 25];
    ml.HorizontalAlignment = 'right';
    ml.Text = 'Ref. Modulation';
    
    %% Ref Phase
    md = uieditfield(xp);
    md.Position = [105 90 100 25];
    md.BackgroundColor = [0.7 0.7 0.7];
    md.HorizontalAlignment = 'right';
    md.Value = num2str(calibration.phi_ref);
    md.Editable = "off";
    
    % Label
    ml1 = uilabel(xp);
    ml1.Position = [20 95 75 25];
    ml1.HorizontalAlignment = 'right';
    ml1.Text = 'Ref. Phase';
    
    ml2 = uilabel(xp);
    ml2.Position = [20 85 75 25];
    ml2.HorizontalAlignment = 'right';
    ml2.Text = '(rad)';
    
    %% Calibraiton  %%
    % Label
    % Calculations
    cp = uipanel(xp);
    cp.Position = [10 10 215 70];
    
    % Label
    cl = uilabel(cp);
    cl.Text = 'Calibration';
    cl.FontWeight = 'bold';
    cl.Position = [0 45 215 25];
    cl.HorizontalAlignment = 'center';

%% Controls %%
    function params = changeParam(params, newValue, paramName)
        params.(paramName) = newValue;
        params = updateParams(params);

    end

    function params = updateParams(params)
        params.w = 2*params.n*pi*params.f*10^6;
    end
end
