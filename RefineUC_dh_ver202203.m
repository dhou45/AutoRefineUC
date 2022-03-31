%% Modified based on B.R. Bickmore Group, http://geol.byu.edu/bickmorewiki/?q=content/refineuc
%% Originally from UnitCell by Tim Holland on Windows

%% need a .csv as input, peak fitting of a few hkl peaks for the in situ datasets.
%% I used LIPRAS for peak fitting, then cleaned the data using excel, then save as .csv 

% read the fitted table, x-axis is Q
    [filename, work_path] = uigetfile('*.csv'); 
    batch_table = readmatrix(fullfile(work_path, filename));
    lamda = 1;
    batch_table_2theta = 2.*asind(batch_table.*lamda./4./pi);
    mkdir('split_csv'); cd('split_csv');

% read the fitted table, x-axis is 2theta
    % batch_table_2theta = readmatrix('/Users/dhou45/Desktop/test_batch.csv');


% initialize the peaks, this case has four peaks with different hkl.
    peak1_index = [0 0 3];
    peak2_index = [1 0 1];
    peak3_index = [1 0 2];
    peak4_index = [1 0 4];

% save each in situ dataset into single .csv file
    for ii = 1 : length (batch_table)
        peak1_index_temp = [peak1_index batch_table_2theta(ii, 1)];
        peak2_index_temp = [peak2_index batch_table_2theta(ii, 2)];
        peak3_index_temp = [peak3_index batch_table_2theta(ii, 3)];
        peak4_index_temp = [peak4_index batch_table_2theta(ii, 4)];
        peak_temp = [peak1_index_temp; peak2_index_temp; peak3_index_temp; peak4_index_temp];
        writematrix(peak_temp, [num2str(ii,  '%02.0f'), '.csv']);
    end

%  read .csv in batch, new folder suggested.
    work_path = uigetdir(pwd,'Select Directory');
    fileList = dir('*.csv');
    LP_a = zeros([length(fileList) 2]);
    LP_c =zeros([length(fileList) 2]);
    LP_vol =zeros([length(fileList) 2]);
    
    for ii = 1 : length (fileList)
        sFullPath = fullfile(work_path, fileList(ii).name);
        RefineUC_dh(sFullPath);
        LP_a(ii,:) = app_UCTable(1, :);
        LP_c(ii,:) = app_UCTable(3, :);
        LP_vol(ii,:) = app_UCTable(7, :);
        % corrected d-space for each hkl.
        % app_cData_final = [cell2table(app_cData) array2table(app_cData_temp)];
    end
    LP = [LP_a, LP_c, LP_vol];
    writematrix(LP, 'RefineUC_results.csv');


%% ##########
%% Change the Wavelength and CrystalSystem for your materials!!!
function RefineUC_dh(sFullPath)
    app.Wavelength = 1;
    app.CrystalSystem = 'Hexagonal';   

    app_lb = []; app_ub = []; app_Aeq = [];
    app_cData = [];
    app_UCTable = [];
    app_cData_temp = [];

    CrystalSystemChanged(app);
    app.lb = app_lb;
    app.ub = app_ub;
    app.Aeq = app_Aeq; 
    app.Beq  =double(zeros(6, 1));

    %Retrieve the data and cull out any that aren't completed.
    mData = getData(app);
    app.cData = app_cData;
 
    %Separate out the different variables into vectors.
    h = mData(:, 1);
    k = mData(:, 2);
    l = mData(:, 3);
    v2Th = mData(:, 4);
    vD = mData(:, 5);
    vQ = CalcQ(app, vD);
    
    %Get the number of peaks.
    iNumPks = length(h);
    %Optimize on Q.
    x0 = [2 2 2 90 90 90];
    vX = fmincon(@fObj1, x0, [], [], app.Aeq, app.Beq, app.lb, app.ub);
    
    %Optimize on 2-theta.
    x0 = vX;
    vX1 = fmincon(@fObj2, x0, [], [], app.Aeq, app.Beq, app.lb, app.ub);

    %%JACKNIFE OPTIMIZATION ON 2-THETA
    mX = zeros(iNumPks-1, 7);
    
    for i = 1:iNumPks-1
        %Perform the deletion and duplication.
        mDataNew = mData;
        mDataNew(i, :) = [];
        
        %Separate out the different variables into vectors.
        h = mDataNew(:, 1);
        k = mDataNew(:, 2);
        l = mDataNew(:, 3);
        v2Th = mDataNew(:, 4);
        vD = mDataNew(:, 5);
        vQ = CalcQ(app, vD);
        
        vX = fmincon(@fObj2, x0, [], [], app.Aeq, app.Beq, app.lb, app.ub);
        mX(i, 1:6) = vX;
        mX(i, 7) = CalcVol(vX(1), vX(2), vX(3), vX(4), vX(5), vX(6));
    end
    
    %Retrieve the data and cull out any that aren't completed.
    mData = getData(app);
    
    %Separate out the different variables into vectors.
    h = mData(:, 1);
    k = mData(:, 2);
    l = mData(:, 3);
    v2Th = mData(:, 4);
    vD = mData(:, 5);
    vQ = CalcQ(app, vD);
    %Get reciprocal lattice parameters.
    [ast, bst, cst, alphast, betast, gammast] = ...
        getRecip(vX1(1), vX1(2), vX1(3), vX1(4), vX1(5), vX1(6));
    

    %Get cell volume.
    dVol = CalcVol(vX1(1), vX1(2), vX1(3), vX1(4), vX1(5), vX1(6));
    
    %Calculate Q
    vQcalc = h.^2 * ast.^2 ...
        + k.^2 * bst.^2 ...
        + l.^2 * cst.^2 ...
        + 2 .* k .* l .* bst .* cst .* cosd(alphast) ...
        + 2 .* l .* h .* cst .* ast .* cosd(betast) ...
        + 2 .* h .* k .* ast .* bst .* cosd(gammast);
    
    %Calculate 2-Theta
    vDcalc = vQcalc .^ (-0.5);
    v2ThCalc = 2 * asind(app.Wavelength./(2*vDcalc));
    
    
    %%Populate the Unit Cell Parameters table.
    mTableData = zeros(7, 2);
    mTableData(:,1) = [vX1'; dVol];
    mTableData(:,2) = 2 * std(mX)';

    % app.UCTable = num2cell(mTableData);
    % app.cData(1:iNumPks,6) = num2cell(vDcalc);

    assignin('base','app_UCTable', mTableData); 
    assignin('base','app_cData_temp', vDcalc); 

        %%FIRST OBJECTIVE FUNCTION FOR OPTIMIZATION ON Q
        function dSSE = fObj1(x)
            %Unpack
            a = x(1);
            b = x(2);
            c = x(3);
            alpha = x(4);
            beta = x(5);
            gamma = x(6);
            %Get reciprocal lattice parameters.
            [ast, bst, cst, alphast, betast, gammast] = ...
                getRecip(a, b, c, alpha, beta, gamma);
            %Calculate Q
            vQcalc = h.^2 * ast.^2 ...
                + k.^2 * bst.^2 ...
                + l.^2 * cst.^2 ...
                + 2 .* k .* l .* bst .* cst .* cosd(alphast) ...
                + 2 .* l .* h .* cst .* ast .* cosd(betast) ...
                + 2 .* h .* k .* ast .* bst .* cosd(gammast);
                
            dSSE = sum((vQ - vQcalc) .^ 2) * 1000;
        end
        
        %%SECOND OBJECTIVE FUNCTION FOR OPTIMIZATION ON 2-THETA
        function dSSE = fObj2(x)
            %Unpack
            a = x(1);
            b = x(2);
            c = x(3);
            alpha = x(4);
            beta = x(5);
            gamma = x(6);
            
            %Get reciprocal lattice parameters.
            [ast, bst, cst, alphast, betast, gammast] = ...
                getRecip(a, b, c, alpha, beta, gamma);
            
            %Calculate Q
            vQcalc = h.^2 * ast.^2 ...
                + k.^2 * bst.^2 ...
                + l.^2 * cst.^2 ...
                + 2 .* k .* l .* bst .* cst .* cosd(alphast) ...
                + 2 .* l .* h .* cst .* ast .* cosd(betast) ...
                + 2 .* h .* k .* ast .* bst .* cosd(gammast);
            
            %Calculate 2-Theta
            vDcalc = vQcalc .^ (-0.5);
            v2ThCalc = 2 * asind(app.Wavelength./(2*vDcalc));
                
            dSSE = sum((v2Th - v2ThCalc) .^ 2) * 1000;
        end
% end



        function mData = getData(app)
            %Retrieve the data and cull out any that aren't completed.
            mData = load(sFullPath);
            mData(:, 1:3) = int8(mData(:, 1:3));
            cTemp = cellstr(num2str(mData, '%i %i %i %f'));
            cData = cell(length(cTemp), 6);           
            for i = 1:length(cTemp)
                cData(i, 1:4) = split(strtrim(cTemp{i}))';
                cData{i, 5} = num2str(CalcD(app, mData(i, 4)));
                % v2Th = mData(i, 4);
                % vD = Wavelength ./ (2*sind(v2Th/2));
                % cData{i, 5} = num2str(vD);
                cData{i, 6} = '0';
            end
            mEmpty = cellfun('isempty', cData);
            vEmpty = logical(sum(mEmpty, 2)); 
            cData(vEmpty,:) = [];
            assignin('base','app_cData', cData);   
            mData = cellfun(@str2num, cData);
             
        end

        %%CALCULATE D
        function vD = CalcD(app, v2Th)
            vD = app.Wavelength ./ (2*sind(v2Th/2));
        end

        %%CALCULATE Q VALUES
        function vQ = CalcQ(app, vD)
            vQ = vD.^(-2);
        end
        
        %%GET METRIC MATRIX
        function mG = getMetMat(a, b, c, alpha, beta, gamma)
            %Construct the metric matrix.
            mG = [a*a a*b*cosd(gamma) a*c*cosd(beta);...
                b*a*cosd(gamma) b*b b*c*cosd(alpha);...
                c*a*cosd(beta) c*b*cosd(alpha) c*c];
        end
        
        %%CALCULATE CELL VOLUME
        function dVol = CalcVol(a, b, c, alpha, beta, gamma)
            mG = getMetMat(a, b, c, alpha, beta, gamma);
            dVol = sqrt(det(mG));
        end
        
        %%CALCULATE RECIPROCAL CELL PARAMETERS
        % function [ast, bst, cst, alphast, betast, gammast] = getRecip(app, a, b, c, alpha, beta, gamma)
        function [ast, bst, cst, alphast, betast, gammast] = getRecip(a, b, c, alpha, beta, gamma)
            dVol = CalcVol(a, b, c, alpha, beta, gamma);
            ast = b * c * sind(alpha)/dVol;
            bst = a * c * sind(beta)/dVol;
            cst = a * b * sind(gamma)/dVol;
            alphast = acosd(cosd(beta) * cosd(gamma) - cosd(alpha) ...
                / (sind(beta) * sind(gamma)));
            betast = acosd(cosd(alpha) * cosd(gamma) - cosd(beta) ...
                / (sind(alpha) * sind(gamma)));
            gammast = acosd(cosd(alpha) * cosd(beta) - cosd(gamma) ...
                / (sind(alpha) * sind(beta)));
        end
            

        % Value changed function: CrystalSystemDropDown
        function CrystalSystemChanged(app)
            % value = app.CrystalSystemDropDown.Value;
            value = app.CrystalSystem;
            %Set the optimization constraints based on the crystal system.
            switch value
                case 'Triclinic'
                    app_lb = [0 0 0 0 0 0];
                    app_ub = [500 500 500 180 180 180];
                    app_Aeq = zeros(6);
                case 'Monoclinic'
                    app_lb = [0 0 0 90 90 90];
                    app_ub = [500 500 500 90 180 90];
                    app_Aeq = zeros(6);
                case 'Orthorhombic'
                    app_lb = [0 0 0 90 90 90];
                    app_ub = [500 500 500 90 90 90];
                    app_Aeq = zeros(6);
                case 'Tetragonal'
                    app_lb = [0 0 0 90 90 90];
                    app_ub = [500 500 500 90 90 90];
                    app_Aeq = zeros(6);
                    app_Aeq(1,:) = [1 -1 0 0 0 0];

                case 'Hexagonal'
                    app_lb = [0 0 0 90 90 120];
                    app_ub = [500 500 500 90 90 120];
                    app_Aeq = zeros(6);
                    app_Aeq(1,:) = [1 -1 0 0 0 0];

                case 'Isometric'
                    app_lb = [0 0 0 90 90 90];
                    app_ub = [500 500 500 90 90 90];
                    app_Aeq = zeros(6);
                    app_Aeq(1:2,:) = [1 -1 0 0 0 0; 1 0 -1 0 0 0];

            end
        % disp(app_lb)
            assignin('base','app_lb', app_lb);
            assignin('base','app_ub', app_ub);
            assignin('base','app_Aeq', app_Aeq);   
        end
end