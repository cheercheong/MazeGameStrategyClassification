
% For each field map and corresponding trajectories
for k = 1:8
    
    Trajfile = ['playOnfield' num2str(k) '.txt'];
    Mapfile = ['field' num2str(k) '.txt'];
    
    % Data readin
    InputText = importdata(Trajfile,' ',5);
    InputMap = importdata(Mapfile,' ',5);
    TrajMat = InputText.data;
    MapMat = InputMap.data;

    % Counting StepNum 
    HeaderLoc = find(all(TrajMat == 0,2)); % Find all the starting point
    StepNum = zeros(300,1); % Record the step numbers for every player
    for i = 1:299;
        StepNum(i) = HeaderLoc(i+1) - HeaderLoc(i);
    end
    StepNum(300) = length(TrajMat) - HeaderLoc(300);
    Traj = cell(300,1);     % Record the trajectory for every player
    for i = 1:300
        Traj{i} = TrajMat(HeaderLoc(i):(HeaderLoc(i)+StepNum(i)-1),:);
        % Trajectory of the ith player
    end
    
    % Feature extraction    
    % Dedication: Proportions of steps towards DownRight
    Dedication = zeros(300,1);    
    
    % Diff: Average diffculty value
    Diff = zeros(300,1);
            
    % EasyStep: Proportion of steps towards the easiest direction
    EasyStep = zeros(300,1);
    
    % Laziness: Proportion of turns to easy direction although it is not    
    % the right way
    Laziness = zeros(300,1);
    
    % Memory: Proportion of steps going back
    Memory = zeros(300,1);
    
    % For each player
    for i = 1:300
        
        Down = 1;       % Index of going down
        Right = 1;      % Index of going right
        DirStep = 0;    % # of steps towards downright        
        BackStep = 0;   % # of steps going back
        TurnTime = 0;   % # of turns
        TurnToEasy = 0; % # of turns towards easier but not right direction
        EasyNum = 0;       % # of steps towards the easiest way
        Ngbr = zeros(3,3); % Map of the visible neighboring locations
        Turn = zeros(StepNum(i)-1,2); % Matrix recording each turn
        Turn(1,:) = [Down,Right];     % (1,1)->(3,3)         
        
        % For each step
        for j = 2:StepNum(i)-1
            
            NowRow = Traj{i}(j,1);    % Current location
            NowCol = Traj{i}(j,2);
            NextRow = Traj{i}(j+1,1); % Next location
            NextCol = Traj{i}(j+1,2);
           
            Right = sign(NextCol-NowCol); % Going right if Right=1
            Down = sign(NextRow-NowRow);  % Going down if Down=1
            Turn(j,:) = [Down,Right];
     
            % Counting the steps towards DownRight
            if Right >= 0 && Down >= 0
                DirStep = DirStep+1;
            end
            
            % Finding the difficulty value of current location and the next
            % visible location
            TarRow = NowRow+sign(NextRow-NowRow);
            TarCol = NowCol+sign(NextCol-NowCol);
            NextDiff = MapMat(TarRow,TarCol);
            NowDiff = MapMat(NowRow,NowCol);
            % Calculating the overall difficulty value
            Diff(i) = Diff(i) + NowDiff;
            % Counting steps that towards the easiest direction
            if NowRow>1 && NowRow<512 && NowCol>1 && NowCol<512
                % Map of the neighboring locations
                Ngbr = MapMat((NowRow-1):(NowRow+1),(NowCol-1):(NowCol+1));
                % Exclude the current location
                Ngbr(NowRow,NowCol)=0;
                % Find the easiest locations nearby
                m = min(Ngbr);
                mm = min(m);
                [EasyRow,EasyCol] = find(Ngbr == mm);
                % Count steps towards the easiest location
                ED = ~isempty(find(TarRow==EasyRow,1)); % Easy down
                ER = ~isempty(find(TarCol==EasyCol,1)); % Easy right
                if ED && ER
                    EasyNum = EasyNum + 1;
                end
            end            
            % Count # of turns
            if Down ~= Turn(j-1,1) || Right ~= Turn(j-1,2)
                TurnTime = TurnTime +1;
                % Count stepd towards a easier way, but not downright
                if (NextDiff <= NowDiff) && (Down == -1 || Right == -1)
                    TurnToEasy = TurnToEasy+1;
                end
            end
        end % end j
        
        % Count # of steps going back
        flag = 3;
        while flag <= StepNum(i)-1
            BackTime = 0;
            % Back down
            BD = Turn(flag+BackTime,1) == -1*Turn(flag-1-BackTime,1);
            % Back right
            BR = Turn(flag+BackTime,2) == -1*Turn(flag-1-BackTime,2);
            GoingBack = BD && BR;
            % For avoiding exceeding the range
            BF = flag + BackTime;
            BB = flag - BackTime;
            while (BF <= StepNum(i)-1) &&(BB >= 2) && GoingBack
                 BackStep = BackStep + 1;
                 BackTime = BackTime + 1;
                 BD = Turn(flag+BackTime,1) == -1*Turn(flag-1-BackTime,1);
                 BR = Turn(flag+BackTime,2) == -1*Turn(flag-1-BackTime,2);
                 GoingBack = BD && BR;
                 BF = flag + BackTime;
                 BB = flag - BackTime;
            end
            flag = BF + 1;
        end
        
        % Five features extracted 
        Dedication(i) = DirStep/StepNum(i); 
        Diff(i) = Diff(i)/StepNum(i);
        EasyStep(i) = EasyNum/StepNum(i);
        Laziness(i) = TurnToEasy/TurnTime;
        Memory(i) = BackStep/StepNum(i);

    end % end i
    
    % Write feature data to .txt file
    Datafile = ['data' num2str(k) '.txt'];
    Feature = [Dedication,Diff,EasyStep,Laziness,Memory]; 
    fid = fopen(Datafile,'w');
    fprintf(fid,'%d %d %d\n', 300,4,3);
    for i=1:300
        fprintf(fid,'%.2f %.2f %.2f %.2f\n', Feature(i,:));
    end
    fclose(fid);
end

        

