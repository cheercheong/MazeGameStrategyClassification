 k = 1;
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
    
    % Playfulness: Proportion of steps changing direction
    Playfulness = zeros(300,1);
    
    
    % For each player
    for i = 1:300
        
        Down = 1;       % Index of going down
        Right = 1;      % Index of going right
        DirStep = 0;    % # of steps towards downright        
        BackStep = 0;   % # of steps going back
        TurnTime = 0;   % # of turns
        EasyNum = 0;       % # of steps towards the easiest way
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
            Turn(j,:) = [Down,Right]; % Direction of each step
                
            
            % Counting the steps towards DownRight
            if Right >= 0 && Down >= 0
                DirStep = DirStep+1;
            end
            
            % Calculating the overall difficulty value
            NowDiff = MapMat(NowRow,NowCol);           
            Diff(i) = Diff(i) + NowDiff;
      
            % Count # of turns
            if Down ~= Turn(j-1,1) || Right ~= Turn(j-1,2)
                TurnTime = TurnTime +1;
            end
        end % end j
        
        % Features extracted 
        Dedication(i) = DirStep/StepNum(i);  
        Playfulness(i) = TurnTime/StepNum(i);
        Diff(i) = Diff(i)/StepNum(i);

    end
   EleNum = [180,79,41]; 
   Result = cell(3,1);

    Result{1} = [1	3	4	5	9	10	12	13	15	17	19	20	23	26	27	29	30	34	37	38	39	41	42	43	44	45	47	48	50	51	53	54	56	57	58	59	61	62	63	64	66	67	68	69	70	71	73	74	75	76	77	79	81	82	83	84	85	86	91	92	93	94	95	96	97	98	100	102	107	109	113	114	116	119	122	123	124	125	126	127	129	130	131	132	134	136	137	138	139	141	142	144	145	146	147	149	151	153	154	155	156	162	165	169	170	171	172	173	174	176	179	180	181	182	183	185	186	188	189	192	193	199	200	201	203	204	205	206	209	211	213	214	215	216	217	219	221	227	228	230	231	233	235	236	238	239	241	243	245	250	251	252	253	254	256	257	258	259	260	262	264	265	269	270	272	274	275	277	278	279	280	281	283	284	286	288	291	293	295	296	
    ];
    Result{2} = [2	6	7	8	11	18	21	24	28	31	33	35	49	52	60	78	80	87	90	103	104	106	110	111	112	115	117	118	120	121	128	135	140	148	150	152	158	159	160	161	164	166	167	168	175	177	184	187	190	191	194	195	196	197	198	202	207	210	218	223	224	226	229	244	247	248	249	255	261	263	271	276	282	287	290	292	298	299	300	
    ];
    Result{3} = [14	16	22	25	32	36	40	46	55	65	72	88	89	99	101	105	108	133	143	157	163	178	208	212	220	222	225	232	234	237	240	242	246	266	267	268	273	285	289	294	297	
    ];
    Color = cell(3,1);
    Color{1} = 'r.';
    Color{2} = 'b.';
    Color{3} = 'y.';
    figure(1)
    for i = 1:3;
        for j = 1: EleNum(i)
            plot3(Dedication(Result{i}(j)),Playfulness(Result{i}(j)), Diff(Result{i}(j)),Color{i});
            hold on;
        end        
    end
    xlabel('Dedication');
    ylabel('Playfulness');
    zlabel('Diff');
 
    grid on;
    FigNam = ['Field ' num2str(k) ' Clustering after EM Algorithm'];
    title(FigNam);
    figure(2)
    plot3(Dedication,Playfulness, Diff,'r.');
    xlabel('Dedication');
    ylabel('Playfulness');
    zlabel('Diff');
    grid on;
    FigNam = ['Field ' num2str(k) ' Clustering'];
