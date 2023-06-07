classdef rendezvous_game_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                 matlab.ui.Figure
        AddimpulseButton         matlab.ui.control.Button
        ImpulseZEditField        matlab.ui.control.NumericEditField
        ImpulseZEditFieldLabel   matlab.ui.control.Label
        ImpulseYEditField        matlab.ui.control.NumericEditField
        ImpulseYEditFieldLabel   matlab.ui.control.Label
        ImpulseXEditField        matlab.ui.control.NumericEditField
        ImpulseXEditFieldLabel   matlab.ui.control.Label
        PlayButton               matlab.ui.control.Button
        vzEditField              matlab.ui.control.NumericEditField
        vzEditFieldLabel         matlab.ui.control.Label
        vyEditField              matlab.ui.control.NumericEditField
        vyEditFieldLabel         matlab.ui.control.Label
        vxEditField              matlab.ui.control.NumericEditField
        vxEditFieldLabel         matlab.ui.control.Label
        PositionZEditField       matlab.ui.control.NumericEditField
        PositionZEditFieldLabel  matlab.ui.control.Label
        PositionYEditField       matlab.ui.control.NumericEditField
        PositionYEditFieldLabel  matlab.ui.control.Label
        PositionXEditField       matlab.ui.control.NumericEditField
        PositionXEditFieldLabel  matlab.ui.control.Label
        UIAxes                   matlab.ui.control.UIAxes
    end

    
    properties (Access = public)
        X_t_before
        X_t_actual
        hitted
        X_0
        target
    end
   

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function costanti(app)
         
        end

        % Button pushed function: PlayButton
        function PlayButtonPushed(app, event)
            
            h = 500e3; 
            R = 6378.137e3; 
            mu = 3.986004418e14; 
            r = h + R; 
            n = sqrt(mu/r^3); 
            tau =(2*pi)/n; 
            toll = 3; 
            A = [zeros(3), eye(3);[3*n^2 0 0 0 2*n 0;0 0 0 -2*n 0 0;0 0 -n^2 0 0 0]]; 
            phi=@(t) expm(A*t);
            phiX  = @(t,X0) phi(t)*X0;
            
            app.target = [0 0 0 0 0 0];
            plot3(app.UIAxes,app.target(1),app.target(2),app.target(3),'ro','LineWidth',4)
            hold(app.UIAxes,'on')

            x0 = app.PositionXEditField.Value;
            y0 = app.PositionYEditField.Value;
            z0 = app.PositionZEditField.Value;
            vx0 = app.vxEditField.Value;
            vy0 = app.vyEditField.Value;
            vz0 = app.vzEditField.Value;
            X_0=[x0 y0 z0 vx0 vy0 vz0]';
            if app.hitted == 1
                X_0 = app.X_0;
            end

            phiX  = @(t,X0) phi(t)*X0;

            plot3(app.UIAxes,X_0(1)-app.target(1),X_0(2)-app.target(2),X_0(3)-app.target(3),'bo','LineWidth',1.5)
            
            t=linspace(1,tau,200);
            X_t = zeros(length(X_0),length(t));
            X_t_old = zeros(length(X_0),length(t));
            X_t_old_1 = zeros(length(X_0),length(t));
            X_t_old_2 = zeros(length(X_0),length(t));
            X_t_old_3 = zeros(length(X_0),length(t));
            X_t_old_4 = zeros(length(X_0),length(t));
            for i=1:length(t) 
                X_t(:,i)=phiX(t(i),X_0);
            end
            app.X_t_actual = X_t;

            %%%%%%%%%%%%%PARTE ANIMATA%%%%%%%%%%%%%
            i=1;
            control = 0;
            while control < 5

                 if  (abs(X_t(1,i))-app.target(1)+abs(X_t(2,i))-app.target(2)+abs(X_t(3,i))-app.target(3))<toll
                    break
                end

                app.PositionXEditField.Value=X_t(1,i);
                app.PositionYEditField.Value=X_t(2,i);
                app.PositionZEditField.Value=X_t(3,i);

                app.vxEditField.Value=X_t(4,i);
                app.vyEditField.Value=X_t(5,i);
                app.vzEditField.Value=X_t(6,i);

                cla(app.UIAxes);

                plot3(app.UIAxes,app.target(1),app.target(2),app.target(3),'ro','LineWidth',4)
                if app.hitted == 1
                    plot3(app.UIAxes,app.X_t_before(1,:)-app.target(1),app.X_t_before(2,:)-app.target(2),app.X_t_before(3,:)-app.target(3),'b--','LineWidth',0.5)
                end
                hold (app.UIAxes,'on')
                % plot3(app.UIAxes,X_0(1)-app.target(1),X_0(2)-app.target(2),X_0(3)-app.target(3),'bo','LineWidth',1.5)
                plot3(app.UIAxes,X_t(1,i)-app.target(1),X_t(2,i)-app.target(2),X_t(3,i)-app.target(3),'bo','LineWidth',1)
                plot3(app.UIAxes,X_t(1,:)-app.target(1),X_t(2,:)-app.target(2),X_t(3,:)-app.target(3),'g--','LineWidth',0.5)
                plot3(app.UIAxes,X_t_old(1,:)-app.target(1),X_t_old(2,:)-app.target(2),X_t_old(3,:)-app.target(3),'k--','LineWidth',0.5)
                plot3(app.UIAxes,X_t_old_1(1,:)-app.target(1),X_t_old_1(2,:)-app.target(2),X_t_old_1(3,:)-app.target(3),'k--','LineWidth',0.5)
                plot3(app.UIAxes,X_t_old_2(1,:)-app.target(1),X_t_old_2(2,:)-app.target(2),X_t_old_2(3,:)-app.target(3),'k--','LineWidth',0.5)
                plot3(app.UIAxes,X_t_old_3(1,:)-app.target(1),X_t_old_3(2,:)-app.target(2),X_t_old_3(3,:)-app.target(3),'k--','LineWidth',0.5)
                plot3(app.UIAxes,X_t_old_4(1,:)-app.target(1),X_t_old_4(2,:)-app.target(2),X_t_old_4(3,:)-app.target(3),'k--','LineWidth',0.5)
                grid(app.UIAxes,'on')
                xlabel(app.UIAxes,'x')
                ylabel(app.UIAxes,'y')
                zlabel(app.UIAxes,'z')

                pause(0.2)

                if i==length(t)
                    X_0_new = X_t(:,i);
                    X_t_old_4 = X_t_old_3;
                    X_t_old_3 = X_t_old_2;
                    X_t_old_2 = X_t_old_1;
                    X_t_old_1 = X_t_old;
                    X_t_old = X_t;                    
                    X_t = zeros(6,length(t));
                    for j=1:length(t) 
                        X_t(:,j)=phiX(t(j),X_0_new);
                    end
                    i=0;
                    control = control + 1;
                    app.X_t_actual = X_t;
                end                 
                i=i+1;              
            end

        end

        % Value changed function: PositionXEditField
        function PositionXEditFieldValueChanged(app, event)

            
        end

        % Value changed function: PositionYEditField
        function PositionYEditFieldValueChanged(app, event)
           
            
        end

        % Value changed function: PositionZEditField
        function PositionZEditFieldValueChanged(app, event)
            
            
        end

        % Value changed function: vxEditField
        function vxEditFieldValueChanged(app, event)
            
            
        end

        % Value changed function: vyEditField
        function vyEditFieldValueChanged(app, event)
           
            
        end

        % Value changed function: vzEditField
        function vzEditFieldValueChanged(app, event)
           
            
        end

        % Callback function
        function ImpulseXEditFieldValueChanged(app, event)
           
            
        end

        % Value changed function: ImpulseYEditField
        function ImpulseYEditFieldValueChanged(app, event)
            
            
        end

        % Value changed function: ImpulseZEditField
        function ImpulseZEditFieldValueChanged(app, event)
           
            
        end

        % Button pushed function: AddimpulseButton
        function AddimpulseButtonPushed(app, event)
 
            impulsox = app.ImpulseXEditField.Value;
            impulsoy = app.ImpulseYEditField.Value;
            impulsoz = app.ImpulseZEditField.Value;

            x0 = app.PositionXEditField.Value;
            y0 = app.PositionYEditField.Value;
            z0 = app.PositionZEditField.Value;
            vx0 = app.vxEditField.Value;
            vy0 = app.vyEditField.Value;
            vz0 = app.vzEditField.Value;
    
            app.X_0=[x0 y0 z0 vx0+impulsox vy0+impulsoy vz0+impulsoz]';

            for k = 1:size(app.X_t_actual,2)
                if app.X_t_actual(1:3,k) == app.X_0(1:3)
                    imp_t = k;
                end
            end

            app.X_t_before = app.X_t_actual(:,1:imp_t);

            app.hitted = 1;

            app.PlayButtonPushed(@AddimpulseButtonPushed)

        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 650 488];
            app.UIFigure.Name = 'MATLAB App';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'Title')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.XGrid = 'on';
            app.UIAxes.YGrid = 'on';
            app.UIAxes.ZGrid = 'on';
            app.UIAxes.Position = [137 117 477 316];

            % Create PositionXEditFieldLabel
            app.PositionXEditFieldLabel = uilabel(app.UIFigure);
            app.PositionXEditFieldLabel.HorizontalAlignment = 'right';
            app.PositionXEditFieldLabel.Position = [51 441 59 22];
            app.PositionXEditFieldLabel.Text = 'Position X';

            % Create PositionXEditField
            app.PositionXEditField = uieditfield(app.UIFigure, 'numeric');
            app.PositionXEditField.ValueChangedFcn = createCallbackFcn(app, @PositionXEditFieldValueChanged, true);
            app.PositionXEditField.Position = [125 441 53 22];

            % Create PositionYEditFieldLabel
            app.PositionYEditFieldLabel = uilabel(app.UIFigure);
            app.PositionYEditFieldLabel.HorizontalAlignment = 'right';
            app.PositionYEditFieldLabel.Position = [290 441 59 22];
            app.PositionYEditFieldLabel.Text = 'Position Y';

            % Create PositionYEditField
            app.PositionYEditField = uieditfield(app.UIFigure, 'numeric');
            app.PositionYEditField.ValueChangedFcn = createCallbackFcn(app, @PositionYEditFieldValueChanged, true);
            app.PositionYEditField.Position = [364 441 54 22];

            % Create PositionZEditFieldLabel
            app.PositionZEditFieldLabel = uilabel(app.UIFigure);
            app.PositionZEditFieldLabel.HorizontalAlignment = 'right';
            app.PositionZEditFieldLabel.Position = [499 441 58 22];
            app.PositionZEditFieldLabel.Text = 'Position Z';

            % Create PositionZEditField
            app.PositionZEditField = uieditfield(app.UIFigure, 'numeric');
            app.PositionZEditField.ValueChangedFcn = createCallbackFcn(app, @PositionZEditFieldValueChanged, true);
            app.PositionZEditField.Position = [572 441 53 22];

            % Create vxEditFieldLabel
            app.vxEditFieldLabel = uilabel(app.UIFigure);
            app.vxEditFieldLabel.HorizontalAlignment = 'right';
            app.vxEditFieldLabel.Position = [38 344 25 22];
            app.vxEditFieldLabel.Text = 'vx';

            % Create vxEditField
            app.vxEditField = uieditfield(app.UIFigure, 'numeric');
            app.vxEditField.ValueChangedFcn = createCallbackFcn(app, @vxEditFieldValueChanged, true);
            app.vxEditField.Position = [78 344 48 22];

            % Create vyEditFieldLabel
            app.vyEditFieldLabel = uilabel(app.UIFigure);
            app.vyEditFieldLabel.HorizontalAlignment = 'right';
            app.vyEditFieldLabel.Position = [38 301 25 22];
            app.vyEditFieldLabel.Text = 'vy';

            % Create vyEditField
            app.vyEditField = uieditfield(app.UIFigure, 'numeric');
            app.vyEditField.ValueChangedFcn = createCallbackFcn(app, @vyEditFieldValueChanged, true);
            app.vyEditField.Position = [78 301 48 22];

            % Create vzEditFieldLabel
            app.vzEditFieldLabel = uilabel(app.UIFigure);
            app.vzEditFieldLabel.HorizontalAlignment = 'right';
            app.vzEditFieldLabel.Position = [38 256 25 22];
            app.vzEditFieldLabel.Text = 'vz';

            % Create vzEditField
            app.vzEditField = uieditfield(app.UIFigure, 'numeric');
            app.vzEditField.ValueChangedFcn = createCallbackFcn(app, @vzEditFieldValueChanged, true);
            app.vzEditField.Position = [78 256 48 22];

            % Create PlayButton
            app.PlayButton = uibutton(app.UIFigure, 'push');
            app.PlayButton.ButtonPushedFcn = createCallbackFcn(app, @PlayButtonPushed, true);
            app.PlayButton.Position = [38 196 100 23];
            app.PlayButton.Text = 'Play';

            % Create ImpulseXEditFieldLabel
            app.ImpulseXEditFieldLabel = uilabel(app.UIFigure);
            app.ImpulseXEditFieldLabel.HorizontalAlignment = 'right';
            app.ImpulseXEditFieldLabel.Position = [23 79 58 22];
            app.ImpulseXEditFieldLabel.Text = 'Impulse X';

            % Create ImpulseXEditField
            app.ImpulseXEditField = uieditfield(app.UIFigure, 'numeric');
            app.ImpulseXEditField.Position = [96 79 100 22];

            % Create ImpulseYEditFieldLabel
            app.ImpulseYEditFieldLabel = uilabel(app.UIFigure);
            app.ImpulseYEditFieldLabel.HorizontalAlignment = 'right';
            app.ImpulseYEditFieldLabel.Position = [253 79 58 22];
            app.ImpulseYEditFieldLabel.Text = 'Impulse Y';

            % Create ImpulseYEditField
            app.ImpulseYEditField = uieditfield(app.UIFigure, 'numeric');
            app.ImpulseYEditField.ValueChangedFcn = createCallbackFcn(app, @ImpulseYEditFieldValueChanged, true);
            app.ImpulseYEditField.Position = [326 79 100 22];

            % Create ImpulseZEditFieldLabel
            app.ImpulseZEditFieldLabel = uilabel(app.UIFigure);
            app.ImpulseZEditFieldLabel.HorizontalAlignment = 'right';
            app.ImpulseZEditFieldLabel.Position = [452 79 58 22];
            app.ImpulseZEditFieldLabel.Text = 'Impulse Z';

            % Create ImpulseZEditField
            app.ImpulseZEditField = uieditfield(app.UIFigure, 'numeric');
            app.ImpulseZEditField.ValueChangedFcn = createCallbackFcn(app, @ImpulseZEditFieldValueChanged, true);
            app.ImpulseZEditField.Position = [525 79 100 22];

            % Create AddimpulseButton
            app.AddimpulseButton = uibutton(app.UIFigure, 'push');
            app.AddimpulseButton.ButtonPushedFcn = createCallbackFcn(app, @AddimpulseButtonPushed, true);
            app.AddimpulseButton.Position = [270 31 100 23];
            app.AddimpulseButton.Text = 'Add impulse';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = rendezvous_game_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @costanti)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end