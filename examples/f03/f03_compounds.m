function [s] = f03_compounds(path)

nC = 28;
nOH = 5;       % hydroxyl
nCH3 = 2;      % terminal 
nCOOH = 2;     % acid
nHCO = 2;      % aldehyde
nCHnCO = 6;    % carbonyl
nCHnO = 2;     % ether
nCHnONO2 = 4;  % nitrate
nCHnOOH = 2;   % hydroperoxide

f = 1;
for j = 5:nC-2 + 1
    for k = 2:nOH + 1
        s(f).nxCH3    = 2;
        s(f).nxCOOH   = 0;
        s(f).nxHCO    = 0;
        
        s(f).nxCH2    = (j-1) - (k-1);
        s(f).nxCH     = (k-1);
        s(f).nxOH     = (k-1);
        s(f).nxCHnCO  = 0;
        s(f).nxCHnO   = 0;
        s(f).nxCHnONO2 = 0;
        s(f).nxCHnOOH = 0;
        s(f).D = 200;
        
        s(f).file = ['CH3_' num2str(s(f).nxCH3) ...
                     '_CH2_' num2str(s(f).nxCH2) ...
                     '_CH_' num2str(s(f).nxCH) ...
                     '_OH_' num2str(s(f).nxOH) ...
                     '_COOH_' num2str(s(f).nxCOOH) ...
                     '_HCO_' num2str(s(f).nxHCO) ...
                     '_CHnCO_' num2str(s(f).nxCHnCO) ...
                     '_CHnO_' num2str(s(f).nxCHnO) ...
                     '_CHnONO2_' num2str(s(f).nxCHnONO2) ...
                     '_CHnOOH_' num2str(s(f).nxCHnOOH)];
   
 
        s(f).groups = zeros(9, 4);
        s(f).groups(1,1) = s(f).nxCH3;
        s(f).groups(1,2) = s(f).nxCH2;
        s(f).groups(1,3) = s(f).nxCH;
        s(f).groups(2,1) = s(f).nxOH;
        s(f).groups(4,2) = s(f).nxCHnCO;
        s(f).groups(5,1) = s(f).nxHCO;
        s(f).groups(6,3) = s(f).nxCHnO;
        s(f).groups(7,1) = s(f).nxCOOH;
        s(f).groups(8,2) = s(f).nxCHnONO2;
        s(f).groups(9,2) = s(f).nxCHnOOH;
                
        of = [path s(f).file];
        dlmwrite(of, s(f).groups,'delimiter','\t','precision', 6);
        f = f+1;

        %% + 2 acids
        if (k-1) == 2 || (k-1) == 3
            s(f).nxCH3    = 0;
            s(f).nxCOOH   = 2;
            s(f).nxHCO    = 0;
            
            s(f).nxCH2    = (j-1) - (k-1);
            s(f).nxCH     = (k-1);
            s(f).nxOH     = (k-1);
            s(f).nxCHnCO  = 0;
            s(f).nxCHnO   = 0;
            s(f).nxCHnONO2 = 0;
            s(f).nxCHnOOH = 0;
            s(f).D = 200;
            
            if  s(f).nxCH2 < 0 
                continue
            end

        
            s(f).file = ['CH3_' num2str(s(f).nxCH3) ...
                         '_CH2_' num2str(s(f).nxCH2) ...
                         '_CH_' num2str(s(f).nxCH) ...
                         '_OH_' num2str(s(f).nxOH) ...
                         '_COOH_' num2str(s(f).nxCOOH) ...
                         '_HCO_' num2str(s(f).nxHCO) ...
                         '_CHnCO_' num2str(s(f).nxCHnCO) ...
                         '_CHnO_' num2str(s(f).nxCHnO) ...
                         '_CHnONO2_' num2str(s(f).nxCHnONO2) ...
                         '_CHnOOH_' num2str(s(f).nxCHnOOH)];
   
 
            s(f).groups = zeros(9, 4);
            s(f).groups(1,1) = s(f).nxCH3;
            s(f).groups(1,2) = s(f).nxCH2;
            s(f).groups(1,3) = s(f).nxCH;
            s(f).groups(2,1) = s(f).nxOH;
            s(f).groups(4,2) = s(f).nxCHnCO;
            s(f).groups(5,1) = s(f).nxHCO;
            s(f).groups(6,3) = s(f).nxCHnO;
            s(f).groups(7,1) = s(f).nxCOOH;
            s(f).groups(8,2) = s(f).nxCHnONO2;
            s(f).groups(9,2) = s(f).nxCHnOOH;
            
            of = [path s(f).file];
            dlmwrite(of, s(f).groups,'delimiter','\t','precision', 6);
            f = f+1;
        end        

        %% + 2 hydroperoxide
        if (k-1) == 2  || (k-1) == 3
            s(f).nxCH3    = 2;
            s(f).nxCOOH   = 0;
            s(f).nxHCO    = 0;
            
            s(f).nxCH2    = (j-1) - (k-1) - 2;
            s(f).nxCH     = (k-1);
            s(f).nxOH     = (k-1);
            s(f).nxCHnCO  = 0;
            s(f).nxCHnO   = 0;
            s(f).nxCHnONO2 = 0;
            s(f).nxCHnOOH = 2;
            s(f).D = 200;

            if  s(f).nxCH2 < 0 
                continue
            end
        
            s(f).file = ['CH3_' num2str(s(f).nxCH3) ...
                         '_CH2_' num2str(s(f).nxCH2) ...
                         '_CH_' num2str(s(f).nxCH) ...
                         '_OH_' num2str(s(f).nxOH) ...
                         '_COOH_' num2str(s(f).nxCOOH) ...
                         '_HCO_' num2str(s(f).nxHCO) ...
                         '_CHnCO_' num2str(s(f).nxCHnCO) ...
                         '_CHnO_' num2str(s(f).nxCHnO) ...
                         '_CHnONO2_' num2str(s(f).nxCHnONO2) ...
                         '_CHnOOH_' num2str(s(f).nxCHnOOH)];
   
 
            s(f).groups = zeros(9, 4);
            s(f).groups(1,1) = s(f).nxCH3;
            s(f).groups(1,2) = s(f).nxCH2;
            s(f).groups(1,3) = s(f).nxCH;
            s(f).groups(2,1) = s(f).nxOH;
            s(f).groups(4,2) = s(f).nxCHnCO;
            s(f).groups(5,1) = s(f).nxHCO;
            s(f).groups(6,3) = s(f).nxCHnO;
            s(f).groups(7,1) = s(f).nxCOOH;
            s(f).groups(8,2) = s(f).nxCHnONO2;
            s(f).groups(9,2) = s(f).nxCHnOOH;
            
            of = [path s(f).file];
            dlmwrite(of, s(f).groups,'delimiter','\t','precision', 6);
            f = f+1;
        end        

        %% + 4 ether
        if (k-1) == 2  || (k-1) == 3
            s(f).nxCH3    = 2;
            s(f).nxCOOH   = 0;
            s(f).nxHCO    = 0;
            
            s(f).nxCH2    = (j-1) - (k-1) - 4;
            s(f).nxCH     = (k-1);
            s(f).nxOH     = (k-1);
            s(f).nxCHnCO  = 0;
            s(f).nxCHnO   = 4;
            s(f).nxCHnONO2 = 0;
            s(f).nxCHnOOH = 0;
            s(f).D = 200;

            if  s(f).nxCH2 < 0 
                continue
            end
        
            s(f).file = ['CH3_' num2str(s(f).nxCH3) ...
                         '_CH2_' num2str(s(f).nxCH2) ...
                         '_CH_' num2str(s(f).nxCH) ...
                         '_OH_' num2str(s(f).nxOH) ...
                         '_COOH_' num2str(s(f).nxCOOH) ...
                         '_HCO_' num2str(s(f).nxHCO) ...
                         '_CHnCO_' num2str(s(f).nxCHnCO) ...
                         '_CHnO_' num2str(s(f).nxCHnO) ...
                         '_CHnONO2_' num2str(s(f).nxCHnONO2) ...
                         '_CHnOOH_' num2str(s(f).nxCHnOOH)];
   
 
            s(f).groups = zeros(9, 4);
            s(f).groups(1,1) = s(f).nxCH3;
            s(f).groups(1,2) = s(f).nxCH2;
            s(f).groups(1,3) = s(f).nxCH;
            s(f).groups(2,1) = s(f).nxOH;
            s(f).groups(4,2) = s(f).nxCHnCO;
            s(f).groups(5,1) = s(f).nxHCO;
            s(f).groups(6,3) = s(f).nxCHnO;
            s(f).groups(7,1) = s(f).nxCOOH;
            s(f).groups(8,2) = s(f).nxCHnONO2;
            s(f).groups(9,2) = s(f).nxCHnOOH;
            
            of = [path s(f).file];
            dlmwrite(of, s(f).groups,'delimiter','\t','precision', 6);
            f = f+1;
        end        

        
        %% + 2 aledyde
        if (k-1) == 2  || (k-1) == 3
            s(f).nxCH3    = 0;
            s(f).nxCOOH   = 0;
            s(f).nxHCO    = 2;
            
            s(f).nxCH2    = (j-1) - (k-1);
            s(f).nxCH     = (k-1);
            s(f).nxOH     = (k-1);
            s(f).nxCHnCO  = 0;
            s(f).nxCHnO   = 0;
            s(f).nxCHnONO2 = 0;
            s(f).nxCHnOOH = 0;
            s(f).D = 200;
        
            s(f).file = ['CH3_' num2str(s(f).nxCH3) ...
                         '_CH2_' num2str(s(f).nxCH2) ...
                         '_CH_' num2str(s(f).nxCH) ...
                         '_OH_' num2str(s(f).nxOH) ...
                         '_COOH_' num2str(s(f).nxCOOH) ...
                         '_HCO_' num2str(s(f).nxHCO) ...
                         '_CHnCO_' num2str(s(f).nxCHnCO) ...
                         '_CHnO_' num2str(s(f).nxCHnO) ...
                         '_CHnONO2_' num2str(s(f).nxCHnONO2) ...
                         '_CHnOOH_' num2str(s(f).nxCHnOOH)];
   
 
            s(f).groups = zeros(9, 4);
            s(f).groups(1,1) = s(f).nxCH3;
            s(f).groups(1,2) = s(f).nxCH2;
            s(f).groups(1,3) = s(f).nxCH;
            s(f).groups(2,1) = s(f).nxOH;
            s(f).groups(4,2) = s(f).nxCHnCO;
            s(f).groups(5,1) = s(f).nxHCO;
            s(f).groups(6,3) = s(f).nxCHnO;
            s(f).groups(7,1) = s(f).nxCOOH;
            s(f).groups(8,2) = s(f).nxCHnONO2;
            s(f).groups(9,2) = s(f).nxCHnOOH;
            
            of = [path s(f).file];
            dlmwrite(of, s(f).groups,'delimiter','\t','precision', 6);
            f = f+1;
        end        

        %% + 3 carbonyls
        if (k-1) == 2  || (k-1) == 3
            s(f).nxCH3    = 2;
            s(f).nxCOOH   = 0;
            s(f).nxHCO    = 0;
            
            s(f).nxCH2    = (j-1) - (k-1) - 6;
            s(f).nxCH     = (k-1);
            s(f).nxOH     = (k-1);
            s(f).nxCHnCO  = 3;
            s(f).nxCHnO   = 0;
            s(f).nxCHnONO2 = 0;
            s(f).nxCHnOOH = 0;
            s(f).D = 200;

            if  s(f).nxCH2 < 0 
                continue
            end
            s(f).file = ['CH3_' num2str(s(f).nxCH3) ...
                         '_CH2_' num2str(s(f).nxCH2) ...
                         '_CH_' num2str(s(f).nxCH) ...
                         '_OH_' num2str(s(f).nxOH) ...
                         '_COOH_' num2str(s(f).nxCOOH) ...
                         '_HCO_' num2str(s(f).nxHCO) ...
                         '_CHnCO_' num2str(s(f).nxCHnCO) ...
                         '_CHnO_' num2str(s(f).nxCHnO) ...
                         '_CHnONO2_' num2str(s(f).nxCHnONO2) ...
                         '_CHnOOH_' num2str(s(f).nxCHnOOH)];
   
 
            s(f).groups = zeros(9, 4);
            s(f).groups(1,1) = s(f).nxCH3;
            s(f).groups(1,2) = s(f).nxCH2;
            s(f).groups(1,3) = s(f).nxCH;
            s(f).groups(2,1) = s(f).nxOH;
            s(f).groups(4,2) = s(f).nxCHnCO;
            s(f).groups(5,1) = s(f).nxHCO;
            s(f).groups(6,3) = s(f).nxCHnO;
            s(f).groups(7,1) = s(f).nxCOOH;
            s(f).groups(8,2) = s(f).nxCHnONO2;
            s(f).groups(9,2) = s(f).nxCHnOOH;
            
            of = [path s(f).file];
            dlmwrite(of, s(f).groups,'delimiter','\t','precision', 6);
            f = f+1;
        end        

        
        %% + 3 nitrate
        if (k-1) == 2  || (k-1) == 3
            s(f).nxCH3    = 2;
            s(f).nxCOOH   = 0;
            s(f).nxHCO    = 0;
            
            s(f).nxCH2    = (j-1) - (k-1) - 5;
            s(f).nxCH     = (k-1);
            s(f).nxOH     = (k-1);
            s(f).nxCHnCO  = 0;
            s(f).nxCHnO   = 0;
            s(f).nxCHnONO2 = 3;
            s(f).nxCHnOOH = 0;
            s(f).D = 200;
        
            if  s(f).nxCH2 < 0 
                continue
            end
            s(f).file = ['CH3_' num2str(s(f).nxCH3) ...
                         '_CH2_' num2str(s(f).nxCH2) ...
                         '_CH_' num2str(s(f).nxCH) ...
                         '_OH_' num2str(s(f).nxOH) ...
                         '_COOH_' num2str(s(f).nxCOOH) ...
                         '_HCO_' num2str(s(f).nxHCO) ...
                         '_CHnCO_' num2str(s(f).nxCHnCO) ...
                         '_CHnO_' num2str(s(f).nxCHnO) ...
                         '_CHnONO2_' num2str(s(f).nxCHnONO2) ...
                         '_CHnOOH_' num2str(s(f).nxCHnOOH)];
   
 
            s(f).groups = zeros(9, 4);
            s(f).groups(1,1) = s(f).nxCH3;
            s(f).groups(1,2) = s(f).nxCH2;
            s(f).groups(1,3) = s(f).nxCH;
            s(f).groups(2,1) = s(f).nxOH;
            s(f).groups(4,2) = s(f).nxCHnCO;
            s(f).groups(5,1) = s(f).nxHCO;
            s(f).groups(6,3) = s(f).nxCHnO;
            s(f).groups(7,1) = s(f).nxCOOH;
            s(f).groups(8,2) = s(f).nxCHnONO2;
            s(f).groups(9,2) = s(f).nxCHnOOH;
            
            of = [path s(f).file];
            dlmwrite(of, s(f).groups,'delimiter','\t','precision', 6);
            f = f+1;
        end        
    
    end
end
