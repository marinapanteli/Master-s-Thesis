function C = Coupling

    A = regexp(fileread('parametersTHESIS_after_arabinose.txt'),'\n','split');

    C = {
    find(contains(A,'kfSGCASsg4'))-1,[find(contains(A,'kbSGCASsg4'))-1 find(contains(A,'kfGBSG4'))-1 find(contains(A,'kbGBSG4'))-1 find(contains(A,'kMC4'))-1];
    find(contains(A,'kbSGCASsg4'))-1,[find(contains(A,'kfSGCASsg4'))-1 find(contains(A,'kfGBSG4'))-1 find(contains(A,'kbGBSG4'))-1 find(contains(A,'kMC4'))-1];
    find(contains(A,'kfGBSG4'))-1,[find(contains(A,'kbSGCASsg4'))-1 find(contains(A,'kfSGCASsg4'))-1 find(contains(A,'kbGBSG4'))-1 find(contains(A,'kMC4'))-1];
    find(contains(A,'kbGBSG4'))-1,[find(contains(A,'kbSGCASsg4'))-1 find(contains(A,'kfSGCASsg4'))-1 find(contains(A,'kfGBSG4'))-1 find(contains(A,'kMC4'))-1];
    find(contains(A,'kMC4'))-1,[find(contains(A,'kbSGCASsg4'))-1 find(contains(A,'kfGBSG4'))-1 find(contains(A,'kbGBSG4'))-1 find(contains(A,'kfSGCASsg4'))-1];
    
    find(contains(A,'kfSGCASsg5'))-1,[find(contains(A,'kbSGCASsg5'))-1 find(contains(A,'kfGASG5'))-1 find(contains(A,'kbGASG5'))-1 find(contains(A,'kMC5'))-1];
    find(contains(A,'kbSGCASsg5'))-1,[find(contains(A,'kfSGCASsg5'))-1 find(contains(A,'kfGASG5'))-1 find(contains(A,'kbGASG5'))-1 find(contains(A,'kMC5'))-1];
    find(contains(A,'kfGASG5'))-1,[find(contains(A,'kbSGCASsg5'))-1 find(contains(A,'kfSGCASsg5'))-1 find(contains(A,'kbGASG5'))-1 find(contains(A,'kMC5'))-1];
    find(contains(A,'kbGASG5'))-1,[find(contains(A,'kbSGCASsg5'))-1 find(contains(A,'kfSGCASsg5'))-1 find(contains(A,'kfGASG5'))-1 find(contains(A,'kMC5'))-1];
    find(contains(A,'kMC5'))-1,[find(contains(A,'kbSGCASsg5'))-1 find(contains(A,'kfGASG5'))-1 find(contains(A,'kbGASG5'))-1 find(contains(A,'kfSGCASsg5'))-1];

    find(contains(A,'kfSGCASsg6'))-1,[find(contains(A,'kbSGCASsg6'))-1 find(contains(A,'kfGASG6'))-1 find(contains(A,'kbGASG6'))-1 find(contains(A,'kMB6'))-1];
    find(contains(A,'kbSGCASsg6'))-1,[find(contains(A,'kfSGCASsg6'))-1 find(contains(A,'kfGASG6'))-1 find(contains(A,'kbGASG6'))-1 find(contains(A,'kMB6'))-1];
    find(contains(A,'kfGASG6'))-1,[find(contains(A,'kbSGCASsg6'))-1 find(contains(A,'kfSGCASsg6'))-1 find(contains(A,'kbGASG6'))-1 find(contains(A,'kMB6'))-1];
    find(contains(A,'kbGASG6'))-1,[find(contains(A,'kbSGCASsg6'))-1 find(contains(A,'kfSGCASsg6'))-1 find(contains(A,'kfGASG6'))-1 find(contains(A,'kMB6'))-1];
    find(contains(A,'kMB6'))-1,[find(contains(A,'kbSGCASsg6'))-1 find(contains(A,'kfGASG6'))-1 find(contains(A,'kbGASG6'))-1 find(contains(A,'kfSGCASsg6'))-1]};

end
